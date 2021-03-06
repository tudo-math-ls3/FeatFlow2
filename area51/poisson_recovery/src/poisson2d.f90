!##############################################################################
!# ****************************************************************************
!# <name> poisson2d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple Poisson
!# problem with constant coefficients on a simple domain in 2D.
!# </purpose>
!##############################################################################

module poisson2d

  use bcassembly
  use bilinearformevaluation
  use boundary
  use collection
  use convergencetable
  use cubature
  use derivatives
  use discretebc
  use dofmapping
  use element
  use filtersupport
  use fsystem
  use genoutput
  use linearformevaluation
  use linearsolver
  use linearsystemblock
  use linearsystemscalar
  use matrixfilters
  use multileveloperators
  use multilevelprojection
  use paramlist
  use pprocerror
  use pprocgradients
  use scalarpde
  use sortstrategy
  use sortstrategybase
  use spatialdiscretisation
  use storage
  use triangulation
  use vectorfilters

  use poisson2d_callback

  implicit none

  public :: t_problem
  public :: poisson_initParamTriang
  public :: poisson_initDiscretisation
  public :: poisson_initMatVec
  public :: poisson_initDiscreteBC
  public :: poisson_implementBC
  public :: poisson_solve
  public :: poisson_postprocessing
  public :: poisson_doneMatVec
  public :: poisson_doneBC
  public :: poisson_doneDiscretisation
  public :: poisson_doneParamTriang

  private 

!<types>

!<typeblock description="Type block defining all information about one level">

  type t_problem_lvl
  
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

    ! A scalar matrix that will recieve the prolongation matrix for this level.
    type(t_matrixScalar) :: rmatProl,rmatRest
    
    ! An interlevel projection structure for changing levels
    type(t_interlevelProjectionBlock) :: rprojection

    ! A variable describing the discrete boundary conditions.
    type(t_discreteBC) :: rdiscreteBC
  
    ! Sorting strategy for resorting vectors/matrices.
    type(t_blockSortStrategy) :: rsortStrategy
    
    ! A filter chain to filter the vectors and the matrix during the
    ! solution process.
    type(t_filterChain), dimension(1) :: RfilterChain

    ! Number of filters in the filter chain.
    integer :: nfilters

  end type
  
!</typeblock>

!<typeblock description="Application-specific type block for Poisson problem">

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
  
  subroutine poisson_initParamTriang(rparlist,ilvmin,ilvmax,rproblem)
  
!<description>
  ! This routine initialises the parametrisation and triangulation of the
  ! domain. The corresponding .prm/.tri files are read from disc and
  ! the triangulation is refined as described by the parameter ilv.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist

  ! Minimum refinement level of the mesh; = coarse grid = level 1
  integer, intent(in) :: ilvmin
  
  ! Maximum refinement level
  integer, intent(in) :: ilvmax
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i

  ! Path to the mesh
  character(len=SYS_STRLEN) :: spredir,strifile,sprmfile

  ! Initialise the level in the problem structure
  rproblem%ilvmin = ilvmin
  rproblem%ilvmax = ilvmax
  allocate(rproblem%RlevelInfo(rproblem%ilvmin:rproblem%ilvmax))

  ! Get the path $PREDIR from the environment, where to read .prm/.tri files
  ! from. If that does not exist, write to the directory "./pre".
  if (.not. sys_getenv_string("PREDIR", spredir)) spredir = "./pre"

  ! At first, read in the parametrisation of the boundary and save
  ! it to rboundary.
  call parlst_getvalue_string(rparlist, '', 'prmfile', sprmfile, '')
  if (trim(sprmfile) .ne. '') then
    call boundary_read_prm(rproblem%rboundary, trim(spredir)//'/'//trim(sprmfile))
  end if

  ! Now read in the basic triangulation.
  call parlst_getvalue_string(rparlist, '', 'trifile', strifile)
  if (trim(sprmfile) .ne. '') then
    call tria_readTriFile2D(rproblem%RlevelInfo(rproblem%ilvmin)%rtriangulation, &
        trim(spredir)//'/'//trim(strifile), rproblem%rboundary)
  else
    call tria_readTriFile2D(rproblem%RlevelInfo(rproblem%ilvmin)%rtriangulation, &
        trim(spredir)//'/'//trim(strifile))
  end if

  ! Refine the mesh up to the minimum level
  if (trim(sprmfile) .ne. '') then
    call tria_quickRefine2LevelOrdering(rproblem%ilvmin-1,&
        rproblem%RlevelInfo(rproblem%ilvmin)%rtriangulation,rproblem%rboundary)
  else
    call tria_quickRefine2LevelOrdering(rproblem%ilvmin-1,&
        rproblem%RlevelInfo(rproblem%ilvmin)%rtriangulation)
  end if

  ! Create information about adjacencies and everything one needs from
  ! a triangulation. Afterwards, we have the coarse mesh.
  if (trim(sprmfile) .ne. '') then    
    call tria_initStandardMeshFromRaw(&
        rproblem%RlevelInfo(rproblem%ilvmin)%rtriangulation,rproblem%rboundary)
  else
    call tria_initStandardMeshFromRaw(&
        rproblem%RlevelInfo(rproblem%ilvmin)%rtriangulation)
  end if

  ! Now, refine to level up to nlmax.
  if (trim(sprmfile) .ne. '') then  
    do i=rproblem%ilvmin+1,rproblem%ilvmax
      call tria_refine2LevelOrdering(rproblem%RlevelInfo(i-1)%rtriangulation,&
          rproblem%RlevelInfo(i)%rtriangulation, rproblem%rboundary)
      call tria_initStandardMeshFromRaw(rproblem%RlevelInfo(i)%rtriangulation,&
          rproblem%rboundary)
    end do
  else
    do i=rproblem%ilvmin+1,rproblem%ilvmax
      call tria_refine2LevelOrdering(rproblem%RlevelInfo(i-1)%rtriangulation,&
          rproblem%RlevelInfo(i)%rtriangulation)
      call tria_initStandardMeshFromRaw(rproblem%RlevelInfo(i)%rtriangulation)
    end do
  end if
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine poisson_initDiscretisation(rparlist,rproblem)
  
!<description>
  ! This routine initialises the discretisation structure of the underlying
  ! problem and saves it to the problem structure.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  character(len=SYS_STRLEN) :: sparameter
  integer :: i,ccubaturetype,celementtype
          
  ! Get type of element
  call parlst_getvalue_string(rparlist, '', 'ELEMENTTYPE', sparameter)
  celementtype = elem_igetID(sparameter)
  
  ! Get type of cubature formula
  call parlst_getvalue_string(rparlist, '', 'CUBATURETYPE', sparameter)
  ccubaturetype = cub_igetID(sparameter)
  
  do i=rproblem%ilvmin,rproblem%ilvmax
    
    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we only have one block.
    call spdiscr_initBlockDiscr(rproblem%RlevelInfo(i)%rdiscretisation,1,&
        rproblem%RlevelInfo(i)%rtriangulation, rproblem%rboundary)
    
    ! rproblem%RlevelInfo(i)%rdiscretisation is a list of scalar
    ! discretisation structures for every component of the solution vector.
    ! Initialise the first element of the list to specify the element
    ! for this solution component:
    call spdiscr_initDiscr_simple( &
        rproblem%RlevelInfo(i)%rdiscretisation%RspatialDiscr(1), &
        celementtype, rproblem%RlevelInfo(i)%rtriangulation, rproblem%rboundary)
    
    ! Set up an cubature info structure to tell the code which cubature
    ! formula to use
    ! Create an assembly information structure which tells the code
    ! the cubature formula to use. Standard: Gauss 3x3.
    call spdiscr_createDefCubStructure(&  
        rproblem%RlevelInfo(i)%rdiscretisation%RspatialDiscr(1),&
        rproblem%RlevelInfo(i)%rcubatureInfo,&
        ccubaturetype)
    
  end do
                                   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine poisson_initMatVec(rparlist,rproblem)
  
!<description>
  ! Calculates the system matrix and RHS vector of the linear system
  ! by discretising the problem with the default discretisation structure
  ! in the problem structure.
  ! Sets up a solution vector for the linear system.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  character(len=SYS_STRLEN) :: ssortstrategy
  integer :: i
  
  ! A bilinear and linear form describing the analytic problem to solve
  type(t_bilinearForm) :: rform
  type(t_linearForm) :: rlinform
    
  ! Get sort strategy from parameter list
  call parlst_getvalue_string(rparlist, '', 'SORTSTRATEGY', ssortstrategy, '')
  
  do i=rproblem%ilvmin,rproblem%ilvmax
        
    ! Initialise the block matrix with default values based on
    ! the discretisation.
    call lsysbl_createMatrix(rproblem%RlevelInfo(i)%rdiscretisation,&
        rproblem%RlevelInfo(i)%rmatrix)
    
    ! Save matrix to the collection.
    ! They maybe used later, expecially in nonlinear problems.
    call collct_setvalue_mat(rproblem%rcollection,"LAPLACE",&
        rproblem%RlevelInfo(i)%rmatrix,.true.,i)
    
    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! We create that directly in the block (1,1) of the block matrix
    ! using the discretisation structure of the first block.
    call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,&
        1, 1, LSYSSC_MATRIX9)
    
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
    ! the framework will call the callback routine to get analytical data.
    !
    ! We pass our collection structure as well to this routine,
    ! so the callback routine has access to everything what is
    ! in the collection.
    call bilf_buildMatrixScalar(rform,.true.,&
        rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
        rproblem%RlevelInfo(i)%rcubatureInfo,&
        rcollection=rproblem%rcollection)
    
    ! Apply sort strategy (if any)
    if (trim(ssortstrategy) .ne. '') then
      
      ! Create a sort strategy structure for our discretisation
      call sstrat_initBlockSorting(rproblem%RlevelInfo(i)%rsortStrategy,&
          rproblem%RlevelInfo(i)%rdiscretisation)
      
      ! Calculate the resorting strategy.
      if (trim(ssortstrategy) .eq. 'CM') then
        call sstrat_initCuthillMcKee(&
            rproblem%RlevelInfo(i)%rsortStrategy%p_Rstrategies(1),&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1))
        
      elseif (trim(ssortstrategy) .eq. 'RCM') then
        call sstrat_initRevCuthillMcKee(&
            rproblem%RlevelInfo(i)%rsortStrategy%p_Rstrategies(1),&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1))
        
      elseif (trim(ssortstrategy) .eq. 'XYZCOORD') then
        call sstrat_initXYZsorting(&
            rproblem%RlevelInfo(i)%rsortStrategy%p_Rstrategies(1),&
            rproblem%RlevelInfo(i)%rdiscretisation%RspatialDiscr(1), 0)
        
      elseif (trim(ssortstrategy) .eq. 'ZYXCOORD') then
        call sstrat_initXYZsorting(&
            rproblem%RlevelInfo(i)%rsortStrategy%p_Rstrategies(1),&
            rproblem%RlevelInfo(i)%rdiscretisation%RspatialDiscr(1), 1)
        
      elseif (trim(ssortstrategy) .eq. 'FEAST') then
        call sstrat_initFEASTsorting(&
            rproblem%RlevelInfo(i)%rsortStrategy%p_Rstrategies(1),&
            rproblem%RlevelInfo(i)%rdiscretisation%RspatialDiscr(1))
        
      elseif (trim(ssortstrategy) .eq. 'RANDOM') then
        call sstrat_initRandom(&
            rproblem%RlevelInfo(i)%rsortStrategy%p_Rstrategies(1))
        
      else
        
        call output_line("Invalid sort strategy.", &
            OU_CLASS_ERROR,OU_MODE_STD,"poisson_initMatVec")
        call sys_halt()
        
      end if
      
      ! Attach the sorting strategy to the matrix. The matrix is not yet sorted.
      call lsysbl_setSortStrategy(rproblem%RlevelInfo(i)%rmatrix,&
          rproblem%RlevelInfo(i)%rsortStrategy,&
          rproblem%RlevelInfo(i)%rsortStrategy)
      
      ! Resort the matrix.
      call lsysbl_sortMatrix(rproblem%RlevelInfo(i)%rmatrix,.true.)
    end if
    
  end do
  
  ! Next step: Create a RHS vector and a solution vector and a temporary
  ! vector. All are filled with zero.
  call lsysbl_createVector(&
      rproblem%RlevelInfo(rproblem%ilvmax)%rdiscretisation,rproblem%rrhs,.true.)
  call lsysbl_createVector(&
      rproblem%RlevelInfo(rproblem%ilvmax)%rdiscretisation,rproblem%rvector,.true.)
  
  ! Save the solution/RHS vector to the collection. Might be used
  ! later (e.g. in nonlinear problems)
  call collct_setvalue_vec(rproblem%rcollection,"RHS",rproblem%rrhs,.true.)
  call collct_setvalue_vec(rproblem%rcollection,"SOLUTION",rproblem%rvector,.true.)
  
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
  call linf_buildVectorScalar(&
      rlinform,.true.,rproblem%rrhs%RvectorBlock(1),&
      rproblem%RlevelInfo(rproblem%ilvmax)%rcubatureInfo,&
      coeff_RHS,rproblem%rcollection)
  
  ! Clear the solution vector on the finest level.
  call lsysbl_clearVector(rproblem%rvector)
  
  ! Install the resorting strategy in the RHS- and the solution
  ! vector, but do not resort them yet!
  ! We resort the vectors just before solving.
  if (trim(ssortstrategy) .ne. '') then
    call lsysbl_setSortStrategy(rproblem%rrhs,&
        rproblem%RlevelInfo(rproblem%ilvmax)%rsortStrategy)
    call lsysbl_setSortStrategy(rproblem%rvector,&
        rproblem%RlevelInfo(rproblem%ilvmax)%rsortStrategy)
  end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine poisson_initDiscreteBC(rproblem)
  
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
  integer :: i,iboundComp,iboundSeg
  type(t_boundaryRegion) :: rboundaryRegion
  
  do i=rproblem%ilvmin,rproblem%ilvmax
    
    ! Create a t_discreteBC structure where we store all discretised boundary
    ! conditions.
    call bcasm_initDiscreteBC(rproblem%RlevelInfo(i)%rdiscreteBC)
    
    do iboundComp=1,boundary_igetNBoundComp(rproblem%rboundary)
      
      do iboundSeg=1,boundary_igetNsegments(rproblem%rboundary,iboundComp)
        
        ! We ask the boundary routines to create a "boundary region"
        ! - which is simply a part of the boundary corresponding to
        ! a boundary segment.  A boundary region roughly contains
        ! the type, the min/max parameter value and whether the
        ! endpoints are inside the region or not.
        call boundary_createRegion(rproblem%rboundary,iboundComp,iboundSeg,&
            rboundaryRegion)
        
        ! We use this boundary region and specify that we want to
        ! have Dirichlet boundary there. The following call does the
        ! following:
        ! - Create Dirichlet boundary conditions on the region
        !   rboundaryRegion.  We specify icomponent="1" to indicate
        !   that we set up the Dirichlet BC`s for the first (here:
        !   one and only) component in the solution vector.
        ! - Discretise the boundary condition so that the BC`s can
        !   be applied to matrices and vectors
        ! - Add the calculated discrete BC`s to rdiscreteBC for
        !   later use.
        call bcasm_newDirichletBConRealBD(&
            rproblem%RlevelInfo(i)%rmatrix%p_rblockDiscrTest,1,&
            rboundaryRegion,rproblem%RlevelInfo(i)%rdiscreteBC,&
            getBoundaryValues,rproblem%rcollection)
        
      end do
    end do
  end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine poisson_implementBC (rproblem)
  
!<description>
  ! Implements boundary conditions into the RHS and into a given solution vector.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i
  
  ! Next step is to implement boundary conditions into the RHS,
  ! solution and matrix. This is done using a vector/matrix filter
  ! for discrete boundary conditions.
  call vecfil_discreteBCrhs(rproblem%rrhs,&
      rproblem%RlevelInfo(rproblem%ilvmax)%rdiscreteBC)
  call vecfil_discreteBCsol(rproblem%rvector,&
      rproblem%RlevelInfo(rproblem%ilvmax)%rdiscreteBC)
  
  ! Implement discrete boundary conditions into the matrices on all
  ! levels, too. Call the appropriate matrix filter to modify
  ! all matrices according to the attached discrete boundary conditions.
  do i=rproblem%ilvmin,rproblem%ilvmax
    call matfil_discreteBC(rproblem%RlevelInfo(i)%rmatrix,&
        rproblem%RlevelInfo(i)%rdiscreteBC)
  end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine poisson_solve(rproblem)
  
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

  ! Set up prolongation matrices
  do i=ilvmin+1,ilvmax
    ! Create the matrix structure of the prolongation matrix.
    call mlop_create2LvlMatrixStruct(&
        rproblem%RlevelInfo(i-1)%rdiscretisation%RspatialDiscr(1),&
        rproblem%RlevelInfo(i)%rdiscretisation%RspatialDiscr(1),&
        LSYSSC_MATRIX9, rproblem%RlevelInfo(i)%rmatProl)

    ! Assemble the entries of the prolongation matrix.
    call mlop_build2LvlProlMatrix(&
        rproblem%RlevelInfo(i-1)%rdiscretisation%RspatialDiscr(1),&
        rproblem%RlevelInfo(i)%rdiscretisation%RspatialDiscr(1),&
        .true., rproblem%RlevelInfo(i)%rmatProl,&
        rcubatureInfoCoarse=rproblem%RlevelInfo(i-1)%rcubatureInfo,&
        rcubatureInfoFine=rproblem%RlevelInfo(i)%rcubatureInfo)

    ! Now set up an interlevel projecton structure for this level
    ! based on the Laplace matrix on this level.
    call mlprj_initProjectionMat(rproblem%RlevelInfo(i)%rprojection,&
        rproblem%RlevelInfo(i)%rmatrix)
    call lsyssc_transposeMatrix(rproblem%RlevelInfo(i)%rmatProl,&
        rproblem%RlevelInfo(i)%rmatRest,LSYSSC_TR_ALL)

    ! Initialise the matrix-based projection
    call mlprj_initMatrixProjection(&
        rproblem%RlevelInfo(i)%rprojection%RscalarProjection(1,1),&
        rproblem%RlevelInfo(i)%rmatProl,&
        rmatrixRest=rproblem%RlevelInfo(i)%rmatRest)
  end do

  ! Set up an interlevel projecton structure for the coarse-most level.
  call mlprj_initProjectionMat(rproblem%RlevelInfo(ilvmin)%rprojection,&
      rproblem%RlevelInfo(ilvmin)%rmatrix)

  ! Create a temporary vector we need that for some preparation.
  call lsysbl_createVector(p_rrhs, rvecTmp, .false.)

  ! Now we have to build up the level information for multigrid.
  !
  ! Create a Multigrid-solver. Attach the above filter chain
  ! to the solver, so that the solver automatically filters
  ! the vector during the solution process.
  call linsol_initMultigrid2(p_rsolverNode,ilvmax-ilvmin+1)

  ! Then set up smoothers / coarse grid solver:
  do i=ilvmin,ilvmax

    ! Set up a filter chain for implementing boundary conditions on that level
    call filter_initFilterChain(rproblem%RlevelInfo(i)%RfilterChain,&
        rproblem%RlevelInfo(i)%nfilters)
    call filter_newFilterDiscBCDef(rproblem%RlevelInfo(i)%RfilterChain,&
        rproblem%RlevelInfo(i)%nfilters,rproblem%RlevelInfo(i)%rdiscreteBC)

    ! On the coarsest grid, set up a coarse grid solver and no smoother
    ! On finer grids, set up a smoother but no coarse grid solver.
    nullify(p_rpreconditioner)
    nullify(p_rsmoother)
    nullify(p_rcoarseGridSolver)

    ! Get the level
    call linsol_getMultigrid2Level(p_rsolverNode,1,p_rlevelInfo)

    if (i .eq. ilvmin) then
      ! Set up a BiCGStab solver with ILU preconditioning as coarse grid solver
      ! would be:
      call linsol_initMILUs1x1(p_rpreconditioner,0,0.0_DP)
      call linsol_initBiCGStab(p_rcoarseGridSolver,p_rpreconditioner,&
          rproblem%RlevelInfo(i)%RfilterChain)

      ! Set up UMFPACK coarse grid solver.
      !call linsol_initUMFPACK4(p_rcoarseGridSolver)

    else
      ! Setting up Jacobi smoother for multigrid would be:
      ! call linsol_initJacobi(p_rsmoother)

      ! Setting up VANCA smoother for multigrid would be:
      ! call linsol_initVANKA(p_rsmoother)

      ! Set up an ILU smoother for multigrid with damping parameter 0.7,
      ! 4 smoothing steps:
      call linsol_initMILUs1x1(p_rsmoother,0,0.0_DP)
      call linsol_convertToSmoother(p_rsmoother,4,0.7_DP)
    end if

    ! And add this multi-grid level. We will use the same smoother
    ! for pre- and post-smoothing.
    call linsol_getMultigrid2Level(p_rsolverNode,i-ilvmin+1,p_rlevelInfo)
    p_rlevelInfo%p_rcoarseGridSolver => p_rcoarseGridSolver
    p_rlevelInfo%p_rpresmoother => p_rsmoother
    p_rlevelInfo%p_rpostsmoother => p_rsmoother

    ! Attach the filter chain which imposes boundary conditions on that level.
    p_rlevelInfo%p_RfilterChain => rproblem%RlevelInfo(i)%RfilterChain

    ! Attach our user-defined projection to the level.
    call linsol_initProjMultigrid2Level(p_rlevelInfo,&
        rproblem%RlevelInfo(i)%rprojection)
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
    call lsysbl_duplicateMatrix(rproblem%RlevelInfo(i)%rmatrix,&
        Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
  end do

  call linsol_setMatrices(p_RsolverNode,Rmatrices(ilvmin:ilvmax))

  ! We can release Rmatrices immediately -- as long as we do not
  ! release rproblem%RlevelInfo(i)%rmatrix!
  do i=ilvmin,ilvmax
    call lsysbl_releaseMatrix(Rmatrices(i))
  end do
  deallocate(Rmatrices)

  ! Initialise structure/data of the solver. This allows the
  ! solver to allocate memory / perform some precalculation
  ! to the problem.
  call linsol_initStructure(p_rsolverNode, ierror)

  if (ierror .ne. LINSOL_ERR_NOERROR) then
    call output_line("Matrix structure invalid!",OU_CLASS_ERROR)
    call sys_halt()
  end if

  call linsol_initData(p_rsolverNode, ierror)

  if (ierror .ne. LINSOL_ERR_NOERROR) then
    call output_line("Matrix singular!",OU_CLASS_ERROR)
    call sys_halt()
  end if

  ! Finally solve the system. As we want to solve Ax=b with
  ! b being the real RHS and x being the real solution vector,
  ! we use linsol_solveAdaptively. If b is a defect
  ! RHS and x a defect update to be added to a solution vector,
  ! we would have to use linsol_precondDefect instead.
  call linsol_solveAdaptively(p_rsolverNode,p_rvector,p_rrhs,rvecTmp)

  ! Release solver data and structure
  call linsol_doneData(p_rsolverNode)
  call linsol_doneStructure(p_rsolverNode)

  ! Release the solver node and all subnodes attached to it (if at all):
  call linsol_releaseSolver(p_rsolverNode)

  ! Release the prolongation matrices
  do i=ilvmax,ilvmin+1,-1
    ! Release the projection structure itself
    call mlprj_doneProjection(rproblem%RlevelInfo(i)%rprojection)
    
    ! Release the prolongation matrix
    call lsyssc_releaseMatrix(rproblem%RlevelInfo(i)%rmatProl)
    
    ! Release the restriction matrix
    call lsyssc_releaseMatrix(rproblem%RlevelInfo(i)%rmatRest)
  end do
  
  ! Release the projection structure on the coarse mesh
  call mlprj_doneProjection(rproblem%RlevelInfo(ilvmin)%rprojection)

  ! Release the filter chain
  do i=ilvmin,ilvmax
    call filter_doneFilterChain(rproblem%RlevelInfo(i)%RfilterChain,&
        rproblem%RlevelInfo(i)%nfilters)
  end do

  ! Release the temporary vector
  call lsysbl_releaseVector(rvecTmp) 

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine poisson_postprocessing(rproblem,rtable)
  
!<description>
  ! Writes the solution into a VTK file.
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</input>

!<inputoutput>
  ! Convergence table
  type(t_convergenceTable), intent(inout) :: rtable
!</inputoutput>

!</subroutine>

  ! A temporal block discretisation structure
  type(t_blockDiscretisation) :: rblockDiscr

  ! A temporal vector to store the recovered gradient
  type(t_vectorBlock) :: rvectorBlock

  ! Error of FE function to reference function
  real(DP) :: derror

  ! Add number of cells and vertices
  call ctab_addValue(rtable, "cells",&
      rproblem%rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL)
  call ctab_addValue(rtable, "dofs",&
      dof_igetNDofGlob(rproblem%rvector%RvectorBlock(1)%p_rspatialDiscr))

  ! Calculate the error to the reference function.
  call pperr_scalar(PPERR_L2ERROR,derror,rproblem%rvector%RvectorBlock(1),&
      getReferenceFunction, rcubatureInfo=&
      rproblem%RlevelInfo(rproblem%ilvmax)%rcubatureInfo)
  call ctab_addValue(rtable, "L2-error", derror)

  call pperr_scalar(PPERR_H1ERROR,derror,rproblem%rvector%RvectorBlock(1),&
      getReferenceFunction, rcubatureInfo=&
      rproblem%RlevelInfo(rproblem%ilvmax)%rcubatureInfo)
  call ctab_addValue(rtable, "H1-error", derror)

  ! Recover gradient by superconvergent patch recovery
  call spdiscr_initBlockDiscr(rblockDiscr,2,&
      rproblem%rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation,&
      rproblem%rvector%p_rblockDiscr%p_rboundary)
  call spdiscr_duplicateDiscrSc(rproblem%rvector%RvectorBlock(1)%p_rspatialDiscr,&
      rblockDiscr%RspatialDiscr(1), .true.)
  call spdiscr_duplicateDiscrSc(rproblem%rvector%RvectorBlock(1)%p_rspatialDiscr,&
      rblockDiscr%RspatialDiscr(2), .true.)
  call lsysbl_createVector(rblockDiscr, rvectorBlock, .false.)
  call ppgrd_calcGradient(rproblem%rvector%RvectorBlock(1), rvectorBlock,&
      PPGRD_ZZTECHNIQUE, PPGRD_NODEPATCH)

  ! Calculate the error to the reference DerivX.
  call pperr_scalar(PPERR_L2ERROR,derror,rvectorBlock%RvectorBlock(1),&
      getReferenceDerivX, rcubatureInfo=&
      rproblem%RlevelInfo(rproblem%ilvmax)%rcubatureInfo)
  call ctab_addValue(rtable, "L2-error DerivX", derror)

  call pperr_scalar(PPERR_h1ERROR,derror,rvectorBlock%RvectorBlock(1),&
      getReferenceDerivX, rcubatureInfo=&
      rproblem%RlevelInfo(rproblem%ilvmax)%rcubatureInfo)
  call ctab_addValue(rtable, "H1-error DerivX", derror)

  ! Calculate the error to the reference DerivY.
  call pperr_scalar(PPERR_L2ERROR,derror,rvectorBlock%RvectorBlock(2),&
      getReferenceDerivY, rcubatureInfo=&
      rproblem%RlevelInfo(rproblem%ilvmax)%rcubatureInfo)
  call ctab_addValue(rtable, "L2-error DerivY", derror)

  call pperr_scalar(PPERR_H1ERROR,derror,rvectorBlock%RvectorBlock(2),&
      getReferenceDerivY, rcubatureInfo=&
      rproblem%RlevelInfo(rproblem%ilvmax)%rcubatureInfo)
  call ctab_addValue(rtable, "H1-error DerivY", derror)

  call lsysbl_releaseVector(rvectorBlock)
  call spdiscr_releaseBlockDiscr(rblockDiscr)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine poisson_doneMatVec(rproblem)
  
!<description>
  ! Releases system matrix and vectors.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i

  ! Delete solution/RHS vector
  call lsysbl_releaseVector(rproblem%rvector)
  call lsysbl_releaseVector(rproblem%rrhs)

  ! Delete the variables from the collection.
  call collct_deletevalue(rproblem%rcollection,"RHS")
  call collct_deletevalue(rproblem%rcollection,"SOLUTION")

  ! Release matrices and vectors on all levels
  do i=rproblem%ilvmax,rproblem%ilvmin,-1

    ! Delete the matrix
    call lsysbl_releaseMatrix(rproblem%RlevelInfo(i)%rmatrix)

    ! Delete the variables from the collection.
    call collct_deletevalue(rproblem%rcollection,"LAPLACE",i)

    ! Release the sorting strategy
    call sstrat_doneBlockSorting(rproblem%RlevelInfo(i)%rsortStrategy)

  end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine poisson_doneBC(rproblem)
  
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
      call bcasm_releaseDiscreteBC(rproblem%RlevelInfo(i)%rdiscreteBC)
    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine poisson_doneDiscretisation(rproblem)
  
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
    ! scalar spatial discretisations
    call spdiscr_releaseBlockDiscr(rproblem%RlevelInfo(i)%rdiscretisation)
  end do
    
  end subroutine
    
  ! ***************************************************************************

!<subroutine>

  subroutine poisson_doneParamTriang(rproblem)
  
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
    call tria_done(rproblem%RlevelInfo(i)%rtriangulation)
  end do

  deallocate(rproblem%RlevelInfo)

  ! Finally release the domain.
  call boundary_release(rproblem%rboundary)

  end subroutine

end module poisson2d
