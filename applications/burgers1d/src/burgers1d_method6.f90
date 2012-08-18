!##############################################################################
!# ****************************************************************************
!# <name> burgers1d_method6 </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstation program how to solve a 1D Burgers equation
!# simultaneously in space/time.
!#
!# The 1D burgers equation in space/time is defined by:
!#
!#    $$  u_t  +  u*u_x  -  \nu u_xx  =  0,   u(x,0) = sin(Pi*x)  $$
!#
!# We solve this equation in the space/time domain
!# $(x,t) \in \Omega=[0,1]x[0,1]$. In this example, we use a direct space-time
!# discretisation, i.e. we do not use a separate time-discretisation.
!# Instead, we replace the $t$ variable by $y$-variable of a usual
!# 2D space discretisation, thus resulting in the formula
!#
!#    $$  u_y  +  u*u_x  -  \nu u_xx  =  0,   u(x,0) = sin(Pi*x)  $$
!#
!# For Solving this in the domain $\Omega$, we follow the usual w3ay:
!# The tasks of reading the domain, creating triangulations, discretisation,
!# solving, postprocessing and creanup into different subroutines.
!# The communication between these subroutines is done using an
!# application-specific structure saving problem data.
!#
!# As this problem is nonlinear, we need to invoke a nonlinear solver,
!# which uses a couple of application-spcific callback routines.
!# To provide these routines with necessary data, we build up a collection
!# structure. This is passed through the solver to the callback routines.
!#
!# The preconditioning in this example is done by a Multigrid solver,
!# which performs only a few number of steps, so does not iterate till
!# convergence.
!#
!# </purpose>
!##############################################################################

module burgers1d_method6

  use fsystem
  use storage
  use genoutput
  use linearsolver
  use boundary
  use cubature
  use matrixfilters
  use vectorfilters
  use discretebc
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use scalarpde
  use element
  use bilinearformevaluation
  use linearformevaluation
  use sortstrategybase
  use sortstrategy
  use filtersupport
  use multilevelprojection
  use nonlinearsolver
  use ucd
  
  use matrixio
  
  use collection
    
  use burgers1d_callback
  
  implicit none
  
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

    ! A temporary vector for building the solution when assembling the
    ! matrix on lower levels.
    type(t_vectorBlock) :: rtempVector

    ! A variable describing the discrete boundary conditions.
    type(t_discreteBC) :: rdiscreteBC
  
    ! Sorting strategy for resorting vectors/matrices.
    type(t_blockSortStrategy) :: rsortStrategy

  end type
  
!</typeblock>


!<typeblock description="Application-specific type block for burgers1d problem">

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

  subroutine b1d6_initParamTriang (ilvmin,ilvmax,rproblem)
  
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
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i
  
  ! Path to the mesh
  character(len=SYS_STRLEN) :: spredir

    ! Initialise the level in the problem structure
    rproblem%ilvmin = ilvmin
    rproblem%ilvmax = ilvmax
    
    ! Store the min/max level in the collection to be used in
    ! callback routines
    call collct_setvalue_int(rproblem%rcollection,"NLMIN",ilvmin,.true.)
    call collct_setvalue_int(rproblem%rcollection,"NLMAX",ilvmax,.true.)

    ! Get the path $PREDIR from the environment, where to read .prm/.tri files
    ! from. If that does not exist, write to the directory "./pre".
    if (.not. sys_getenv_string("PREDIR", spredir)) spredir = "./pre"

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call boundary_read_prm(rproblem%rboundary, trim(spredir)//"/QUAD.prm")
        
    ! Now read in the basic triangulation.
    call tria_readTriFile2D (rproblem%RlevelInfo(rproblem%ilvmin)%rtriangulation, &
        trim(spredir)//"/QUAD.tri", rproblem%rboundary)
    
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

  subroutine b1d6_initDiscretisation (rproblem)
  
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
    integer :: i
  
    ! An object for saving the domain:
    type(t_boundary), pointer :: p_rboundary
    
    ! An object for saving the triangulation on the domain
    type(t_triangulation), pointer :: p_rtriangulation

    ! An object for the block discretisation on one level
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
    do i=rproblem%ilvmin,rproblem%ilvmax
      ! Ask the problem structure to give us the boundary and triangulation.
      ! We need it for the discretisation.
      p_rboundary => rproblem%rboundary
      p_rtriangulation => rproblem%RlevelInfo(i)%rtriangulation
      
      ! Now we can start to initialise the discretisation. At first, set up
      ! a block discretisation structure that specifies the blocks in the
      ! solution vector. In this simple problem, we only have one block.
      allocate(p_rdiscretisation)
      call spdiscr_initBlockDiscr (rproblem%RlevelInfo(i)%rdiscretisation,1,&
          p_rtriangulation, p_rboundary)

      ! Get the discretisation
      p_rdiscretisation => rproblem%RlevelInfo(i)%rdiscretisation

      ! p_rdiscretisation%Rdiscretisations is a list of scalar
      ! discretisation structures for every component of the solution vector.
      ! Initialise the first element of the list to specify the element
      ! and cubature rule for this solution component:
      call spdiscr_initDiscr_simple ( &
                  p_rdiscretisation%RspatialDiscr(1), &
                  EL_Q1, p_rtriangulation, p_rboundary)

      ! Set up an cubature info structure to tell the code which cubature
      ! formula to use
      ! Create an assembly information structure which tells the code
      ! the cubature formula to use. Standard: Gauss 3x3.
      call spdiscr_createDefCubStructure(&  
          p_rdiscretisation%RspatialDiscr(1),rproblem%RlevelInfo(i)%rcubatureInfo,&
          CUB_GEN_AUTO_G2)

      call collct_setvalue_cubinfo(rproblem%rcollection,"CUBINFO",&
          rproblem%RlevelInfo(i)%rcubatureInfo,.true.,i)

    end do
                                   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine b1d6_initMatVec (rproblem)
  
!<description>
  ! Calculates the RHS-ector, set up the solution vetor.
  ! Set up the structure of the system matrix/matrices of the linear
  ! system.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i
    
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_linearForm) :: rlinform
    
    ! A pointer to the system matrix and the RHS/solution vectors.
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_vectorBlock), pointer :: p_rrhs,p_rvector,p_rtempVector

    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
  
    ! A cubature information structure
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo

    do i=rproblem%ilvmin,rproblem%ilvmax
      ! Ask the problem structure to give us the discretisation structure
      p_rdiscretisation => rproblem%RlevelInfo(i)%rdiscretisation
      
      ! Get the matrix structure; we want to build a template matrix
      ! on the level, which receives the entries later.
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      
      ! Get the cubature information structure
      p_rcubatureInfo => rproblem%RlevelInfo(i)%rcubatureInfo
            
      ! Initialise the block matrix with default values based on
      ! the discretisation.
      call lsysbl_createMatBlockByDiscr (p_rdiscretisation,p_rmatrix)
      
      ! Save matrix to the collection.
      ! They maybe used later, expecially in nonlinear problems.
      call collct_setvalue_mat(rproblem%rcollection,"SYSTEMMAT",p_rmatrix,.true.,i)

      ! Now using the discretisation, we can start to generate
      ! the structure of the system matrix which is to solve.
      ! We create that directly in the block (1,1) of the block matrix
      ! using the discretisation structure of the first block.
      call bilf_createMatrixStructure (&
          p_rdiscretisation%RspatialDiscr(1),LSYSSC_MATRIX9,&
          p_rmatrix%RmatrixBlock(1,1))

      ! Allocate memory for the matrix, do not calculate the entries.
      ! Remember hat we have a nonlinear matrix, which entries must be build
      ! in evey step of the nonlinear iteration!
      ! We fill the matrix with 1. This is necessary, as the UMFPACK solver
      ! needs nonzero matrix entries for the symbolic factorisation!
      call lsyssc_allocEmptyMatrix(p_rmatrix%RmatrixBlock(1,1),LSYSSC_SETM_ONE)
      
      ! Create a sort strategy structure for our discretisation
      call sstrat_initBlockSorting (rproblem%RlevelInfo(1)%rsortStrategy,p_rdiscretisation)

      ! Calculate the resorting strategy.
      call sstrat_initCuthillMcKee (rproblem%RlevelInfo(1)%rsortStrategy%p_Rstrategies(1),&
          p_rmatrix%RmatrixBlock(1,1))

      ! Attach the sorting strategy to the matrix. Thematrix is not yet sorted.
      call lsysbl_setSortStrategy (p_rmatrix,&
          rproblem%RlevelInfo(1)%rsortStrategy,&
          rproblem%RlevelInfo(1)%rsortStrategy)

      ! Now on all levels except for the maximum one, create a temporary
      ! vector on that level, based on the matrix template.
      ! It is used for building the matrices on lower levels.
      if (i .lt. rproblem%ilvmax) then
        p_rtempVector => rproblem%RlevelInfo(i)%rtempVector
        call lsysbl_createVectorBlock (p_rdiscretisation,p_rtempVector,.false.)
        
        ! Add the temp vector to the collection on level i
        ! for use in the callback routine
        call collct_setvalue_vec(rproblem%rcollection,"RTEMPVEC",p_rtempVector,&
                                .true.,i)
      end if
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
    !
    ! Note that the vector is unsorted after calling this routine!
    call linf_buildVectorScalar (&
              rlinform,.true.,p_rrhs%RvectorBlock(1),p_rcubatureInfo,&
              coeff_RHS,rproblem%rcollection)
                                
    ! Clear the solution vector on the finest level.
    call lsysbl_clearVector(rproblem%rvector)
    
    ! Install the resorting strategy in the RHS- and the solution
    ! vector, but do not resort them yet!
    ! We resort the vectors just before solving.
    call lsysbl_setSortStrategy (p_rrhs,rproblem%RlevelInfo(1)%rsortStrategy)
    call lsysbl_setSortStrategy (p_rvector,rproblem%RlevelInfo(1)%rsortStrategy)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine b1d6_initDiscreteBC (rproblem)
  
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
    type(t_boundaryRegion) :: rboundaryRegion

    ! A pointer to the system matrix and the RHS vector as well as
    ! the discretisation
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_vectorBlock), pointer :: p_rrhs,p_rvector
    type(t_blockDiscretisation), pointer :: p_rdiscretisation

    ! Pointer to structure for saving discrete BC`s:
    type(t_discreteBC), pointer :: p_rdiscreteBC
      
    ! A pointer to the domain
    type(t_boundary), pointer :: p_rboundary

    ! Get the domain from the problem structure
    p_rboundary => rproblem%rboundary

    do i=rproblem%ilvmin,rproblem%ilvmax
    
      ! Get our matrix from the problem structure.
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      
      ! From the matrix or the RHS we have access to the discretisation and the
      ! analytic boundary conditions.
      p_rdiscretisation => p_rmatrix%p_rblockDiscrTest
      
      ! Now we have the raw problem. What is missing is the definition of the boundary
      ! conditions.
      ! For implementing boundary conditions, we use a `filter technique with
      ! discretised boundary conditions`. This means, we first have to calculate
      ! a discrete version of the analytic BC, which we can implement into the
      ! solution/RHS vectors using the corresponding filter.
      !
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
      call boundary_createRegion(p_rboundary,1,1,rboundaryRegion)
      
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
          getBoundaryValues,rproblem%rcollection)
                               
      ! Now to the edge 2 of boundary component 1 the domain.
      call boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
         rboundaryRegion,rproblem%RlevelInfo(i)%rdiscreteBC,&
         getBoundaryValues,rproblem%rcollection)
                               
      ! Edge 3 of boundary component 1.
      ! Ege 3 must be set up as Neumann boundary, which is realised as
      ! simple "do-$nothing"-boundary conditions. So we do not do anything with edge 3
      ! CALL boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
      ! CALL bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
      !     rboundaryRegion,rproblem%RlevelInfo(i)%rdiscreteBC,&
      !     getBoundaryValues,rproblem%rcollection)
      
      ! Edge 4 of boundary component 1. That is it.
      call boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
          rboundaryRegion,rproblem%RlevelInfo(i)%rdiscreteBC,&
          getBoundaryValues,rproblem%rcollection)
                               
      ! Hang the pointer into the vectors and the matrix - more precisely,
      ! to the first block matrix and the first subvector. That way, these
      ! boundary conditions are always connected to that matrix and that
      ! vector.
      p_rdiscreteBC => rproblem%RlevelInfo(i)%rdiscreteBC
      
      p_rmatrix%p_rdiscreteBC => p_rdiscreteBC
      
    end do

    ! On the finest level, attach the discrete BC also
    ! to the solution and RHS vector. They need it to be compatible
    ! to the matrix on the finest level.
    p_rdiscreteBC => rproblem%RlevelInfo(rproblem%ilvmax)%rdiscreteBC
    
    p_rrhs    => rproblem%rrhs
    p_rvector => rproblem%rvector
    
    call lsysbl_assignDiscreteBC(p_rrhs,p_rdiscreteBC)
    call lsysbl_assignDiscreteBC(p_rvector,p_rdiscreteBC)
                
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine b1d6_implementBC (rproblem)
  
!<description>
  ! Implements boundary conditions into the RHS and into a given solution vector.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: ilvmax
  
    ! A pointer to the RHS/solution vector
    type(t_vectorBlock), pointer :: p_rrhs,p_rvector
    
    ! Get our the right hand side and solution from the problem structure
    ! on the finest level
    ilvmax = rproblem%ilvmax
    p_rrhs    => rproblem%rrhs
    p_rvector => rproblem%rvector
    
    ! Next step is to implement boundary conditions into the RHS and
    ! solution. This is done using a vector filter for discrete boundary
    ! conditions.
    ! The discrete boundary conditions are already attached to the
    ! vectors. Call the appropriate vector filter that modifies the vectors
    ! according to the boundary conditions.
    call vecfil_discreteBCrhs (p_rrhs)
    call vecfil_discreteBCsol (p_rvector)

  end subroutine

  ! ***************************************************************************
    subroutine b1d6_getDefect (ite,rx,rb,rd,p_rcollection)
  
    use linearsystemblock
    use collection
    
  !<description>
    ! FOR NONLINEAR ITERATION:
    ! Defect vector calculation callback routine. Based on the current iteration
    ! vector rx and the right hand side vector rb, this routine has to compute the
    ! defect vector rd. The routine accepts a pointer to a collection structure
    ! p_rcollection, which allows the routine to access information from the
    ! main application (e.g. system matrices).
  !</description>

  !<input>
    ! Number of current iteration. 0=build initial defect
    integer, intent(in)                           :: ite

    ! Current iteration vector
    type(t_vectorBlock), intent(in),target        :: rx

    ! Right hand side vector of the equation.
    type(t_vectorBlock), intent(in), target       :: rb
  !</input>
               
  !<inputoutput>
    ! Pointer to collection structure of the application. Points to NULL()
    ! if there is none.
    type(t_collection), pointer                   :: p_rcollection

    ! Defect vector b-A(x)x. This must be filled by the callback routine
    ! with data.
    type(t_vectorBlock), intent(inout), target    :: rd
  !</inputoutput>

      ! local variables
      type(t_bilinearForm) :: rform
      integer :: ilvmax
      type(t_matrixBlock), pointer :: p_rmatrix

      ! A cubature information structure
      type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo

      ! Get maximum level from the collection
      ilvmax = collct_getvalue_int (p_rcollection,"NLMAX")

      ! Get the system matrix on the maximum level
      p_rmatrix => collct_getvalue_mat (p_rcollection,"SYSTEMMAT",ilvmax)
      
      ! Get the cubature information structure
      p_rcubatureInfo => collct_getvalue_cubinfo (p_rcollection,"CUBINFO",ilvmax)

      ! Put a reference to rx into the collection. This way, we inform the callback
      ! routine of the matrix assembly about the solution vector to use
      ! fot the nonlinear term.
      call collct_setvalue_vec(p_rcollection,"RX",rx,.true.)
      
      ! Build the entries with the discretisation routine.
      
      ! For assembling of the entries,
      ! we need a bilinear form, which first has to be set up manually.
      ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
      ! scalar system matrix in 2D.
      
      rform%itermCount = 4
      rform%Idescriptors(1,1) = DER_DERIV_Y   ! u_t
      rform%Idescriptors(2,1) = DER_FUNC
      
      rform%Idescriptors(1,2) = DER_DERIV_X   ! u_x
      rform%Idescriptors(2,2) = DER_FUNC
      
      rform%Idescriptors(1,3) = DER_DERIV_X   ! -u_xx  -> u_x phi_x
      rform%Idescriptors(2,3) = DER_DERIV_X

      ! The 4th last term u_yy is actually not needed. By setting the coefficient
      ! in front of thisterm to 0 (see matrix assembly callback routine), the term
      ! can be switched off. Nevertheless, we add it here for having the
      ! possibility to use it (by setting the cofficient to a value not equal
      ! to 0), which serves as stabilisation for the problem!

      rform%Idescriptors(1,4) = DER_DERIV_Y   ! -u_xx  -> u_x phi_x
      rform%Idescriptors(2,4) = DER_DERIV_Y

      ! In the standard case, we have constant coefficients.
      ! Theoretically, there are some coefficients constant - but for
      ! simplicity, we define them all in the callback routine of the matrix
      ! assembly.
      rform%ballCoeffConstant = .false.
      rform%BconstantCoeff = .false.

      ! Now we can build the matrix entries.
      ! We specify the callback function coeff_burgers for the coefficients.
      ! As long as we use constant coefficients, this routine is not used.
      ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
      ! the framework will call the callback routine to get analytical data.
      !
      ! We pass our collection structure as well to this routine,
      ! so the callback routine has access to everything what is
      ! in the collection.
      call bilf_buildMatrixScalar (rform,.true.,p_rmatrix%RmatrixBlock(1,1),&
          p_rcubatureInfo,coeff_burgers,p_rcollection)

      ! Remove the vector from the collection - not necessary anymore.
      call collct_deletevalue (p_rcollection,"RX")
      
      ! Implement discrete boundary conditions into the matrix.
      ! Call the appropriate matrix filter to modify the system matrix
      ! according to the attached discrete boundary conditions.
      call matfil_discreteBC (p_rmatrix)
      
      ! Build the defect: d=b-Ax
      call lsysbl_copyVector (rb,rd)
      call lsysbl_blockMatVec (p_rmatrix, rx, rd, -1.0_DP, 1.0_DP)
    
      ! Apply the defect-vector filter for discrete boundary conditions
      ! to modify the defect vector according to the (discrete) boundary
      ! conditions.
      call vecfil_discreteBCdef (rd)

      ! That is it
  
    end subroutine
    
  ! ***************************************************************************

    subroutine b1d6_precondDefect (ite,rd,rx,rb,domega,bsuccess,p_rcollection)
  
    use linearsystemblock
    use collection
    
  !<description>
    ! FOR NONLINEAR ITERATION:
    ! Defect vector calculation callback routine. Based on the current iteration
    ! vector rx and the right hand side vector rb, this routine has to compute the
    ! defect vector rd. The routine accepts a pointer to a collection structure
    ! p_rcollection, which allows the routine to access information from the
    ! main application (e.g. system matrices).
  !</description>

  !<inputoutput>
    ! Number of current iteration.
    integer, intent(in)                           :: ite

    ! Defect vector b-A(x)x. This must be replaced by J^{-1} rd by a preconditioner.
    type(t_vectorBlock), intent(inout), target    :: rd

    ! Pointer to collection structure of the application. Points to NULL()
    ! if there is none.
    type(t_collection), pointer                   :: p_rcollection
    
    ! Damping parameter. Is set to rsolverNode%domega (usually = 1.0_DP)
    ! on the first call to the callback routine.
    ! The callback routine can modify this parameter according to any suitable
    ! algorithm to calculate an "optimal damping" parameter. The nonlinear loop
    ! will then use this for adding rd to the solution vector:
    ! $$ x_{n+1} = x_n + domega*rd $$
    ! domega will stay at this value until it is changed again.
    real(DP), intent(inout)                       :: domega
  
    ! If the preconditioning was a success. Is normally automatically set to
    ! TRUE. If there is an error in the preconditioner, this flag can be
    ! set to FALSE. In this case, the nonlinear solver breaks down with
    ! the error flag set to "preconditioner broke down".
    logical, intent(inout)                        :: bsuccess
  !</inputoutput>
  
  !<input>
    ! Current iteration vector
    type(t_vectorBlock), intent(in), target       :: rx

    ! Current right hand side of the nonlinear system
    type(t_vectorBlock), intent(in), target       :: rb
  !</input>
  
    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), pointer :: p_rmatrix
    integer :: ierror,ilvmin,ilvmax,i
    type(t_linsolNode), pointer :: p_rsolverNode
    type(t_vectorBlock), pointer :: p_rvectorFine,p_rvectorCoarse
    type(t_interlevelProjectionBlock), pointer :: p_rprojection
    type(t_vectorScalar), pointer :: p_rvectorTemp
    type(t_bilinearForm) :: rform
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo

      ! Our "parent" (the caller of the nonlinear solver) has prepared
      ! a preconditioner node for us (a linear solver with symbolically
      ! factorised matrices). Get this from the collection.
      p_rsolverNode => collct_getvalue_linsol(p_rcollection,"LINSOLVER")

      ! The matrices are all attached to the solver node, and
      ! symbolic factorisation was already performed. The only thing that
      ! is missing is the numerical factorisation... no wonder,
      ! the matrix entries do not yet exist!
      !
      ! So the task is now: Calculate the matrix entries and perform
      ! a numerical factorisation / initData in the solver!
      !
      ! The matrix on the maximum level is already prepared.
      ! This was done for calculating the nonlinear defect. The matrices
      ! on the other levels are still missing.
      !
      ! Get maximum/minimum level from the collection
      ilvmin = collct_getvalue_int (p_rcollection,"NLMIN")
      ilvmax = collct_getvalue_int (p_rcollection,"NLMAX")
      
      ! Get the interlevel projection structure and the temporary vector
      ! from the collection.
      ! Our "parent" prepared there how to interpolate the solution on the
      ! fine grid to coarser grids.
      p_rprojection => collct_getvalue_ilvp(p_rcollection,"ILVPROJECTION")
      p_rvectorTemp => collct_getvalue_vecsca(p_rcollection,"RTEMPSCALAR")
      
      ! Prepare the matrix assembly on level < NLMAX.
      !
      ! For assembling of the entries,
      ! we need a bilinear form, which first has to be set up manually.
      ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
      ! scalar system matrix in 2D.
      
      rform%itermCount = 4
      rform%Idescriptors(1,1) = DER_DERIV_Y   ! u_t
      rform%Idescriptors(2,1) = DER_FUNC
      
      rform%Idescriptors(1,2) = DER_DERIV_X   ! u_x
      rform%Idescriptors(2,2) = DER_FUNC
      
      rform%Idescriptors(1,3) = DER_DERIV_X   ! -u_xx  -> u_x phi_x
      rform%Idescriptors(2,3) = DER_DERIV_X

      ! The 4th last term u_yy is actually not needed. By setting the coefficient
      ! in front of thisterm to 0 (see matrix assembly callback routine), the term
      ! can be switched off. Nevertheless, we add it here for having the
      ! possibility to use it (by setting the cofficient to a value not equal
      ! to 0), which serves as stabilisation for the problem!

      rform%Idescriptors(1,4) = DER_DERIV_Y   ! -u_xx  -> u_x phi_x
      rform%Idescriptors(2,4) = DER_DERIV_Y

      ! In the standard case, we have constant coefficients.
      ! Theoretically, there are some coefficients constant - but for
      ! simplicity, we define them all in the callback routine of the matrix
      ! assembly.
      rform%ballCoeffConstant = .false.
      rform%BconstantCoeff = .false.
      
      ! Loop through all the levels to set up the matrices.
      ! Note that we loop from the last but one level to the minimum
      ! level, because we must interpolate the solution vector
      ! from the finest to the coarsest for getting the nonlinarity
      ! on all the levels.
      do i=ilvmax-1,ilvmin,-1

        ! Get the cubature information structure
        p_rcubatureInfo => collct_getvalue_cubinfo (p_rcollection,"CUBINFO",i)
      
        ! Get the destination matrix on that level
        p_rmatrix => collct_getvalue_mat (p_rcollection,"SYSTEMMAT",i)
        
        ! Get the temporary vector on level i. Will receive the solution
        ! vector on that level.
        p_rvectorCoarse => collct_getvalue_vec (p_rcollection,"RTEMPVEC",i)
        
        ! Get the solution vector on level i+1. This is either the temporary
        ! vector on that level, or the solution vector on the maximum level.
        if (i .lt. ilvmax-1) then
          p_rvectorFine => collct_getvalue_vec (p_rcollection,"RTEMPVEC",i+1)
        else
          p_rvectorFine => rx
        end if
        
        ! Interpolate the solution from the finer grid to the coarser grid.
        ! The interpolation is configured in the interlevel projection
        ! structure we got from the collection.
        call mlprj_performInterpolation (p_rprojection,p_rvectorCoarse, &
                                         p_rvectorFine,p_rvectorTemp)
        
        ! Now we have the solution vector on the current level.
        ! Next step is to build the matrix entries with the discretisation routine.
        !
        ! Get the system matrix on the current level.
        p_rmatrix => collct_getvalue_mat (p_rcollection,"SYSTEMMAT",i)
        
        ! Put a reference to the solution vector on the current level into
        ! the collection. This way, we inform the callback
        ! routine of the matrix assembly about the solution vector to use
        ! fot the nonlinear term.
        call collct_setvalue_vec(p_rcollection,"RX",p_rvectorCoarse,.true.)
        
        ! We specify the callback function coeff_Laplace for the coefficients.
        ! As long as we use constant coefficients, this routine is not used.
        ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
        ! the framework will call the callback routine to get analytical data.
        !
        ! We pass our collection structure as well to this routine,
        ! so the callback routine has access to everything what is
        ! in the collection.
        call bilf_buildMatrixScalar (rform,.true.,p_rmatrix%RmatrixBlock(1,1),&
            p_rcubatureInfo,coeff_burgers,p_rcollection)
        
        ! Remove RX from the collection, not needed there anymore.
        call collct_deletevalue(p_rcollection,"RX")
        
        ! Implement discrete boundary conditions into the matrix.
        ! Call the appropriate matrix filter to modify the system matrix
        ! according to the attached discrete boundary conditions.
        call matfil_discreteBC (p_rmatrix)
        
        ! Sort the matrix according to the attached sorting strategy -
        ! if there is a sorting strategy attached at all
        ! (this is prepared by the application).
        ! The sorting is important,
        ! - to make the matrices compatible to the vector rd (which is also
        !   resorted later)
        ! - to make the matrices consistent to those we put into the
        !   linear solver; we put sorted matrices into the linear solver!
        call lsyssc_sortMatrix (p_rmatrix%RmatrixBlock(1,1),.true.,.true.)
        
      end do

      ! Sort the matrix on the maximum level - do not forget this!
      p_rmatrix => collct_getvalue_mat (p_rcollection,"SYSTEMMAT",ilvmax)
      call lsyssc_sortMatrix (p_rmatrix%RmatrixBlock(1,1),.true.,.true.)

      ! Ok, system matrices on all levels are assembled now.
      ! Now we turn to invokle the linear solver for preconditioning...
      !
      ! Resort the RHS and solution vector according to the resorting
      ! strategy given in the matrix.
      call lsysbl_sortVector (rd,.true.,p_rvectorTemp)
      
      ! Initialise data of the solver. This in fact performs a numeric
      ! factorisation of the matrices in UMFPACK-like solvers.
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
      call linsol_precondDefect (p_rsolverNode,rd)

      ! Release the numeric factorisation of the matrix.
      ! We do not release the symbolic factorisation, as we can use them
      ! for the next iteration.
      call linsol_doneData (p_rsolverNode)

      ! Unsort the vectors again in case they were resorted before calling
      ! the solver.
      ! We use the first subvector of rvecTmp as temporary data; it is
      ! large enough, as we only have one block.
      call lsysbl_sortVector (rd,.false.,p_rvectorTemp)

      ! Unsort the structure of all matrices without unsorting the entries.
      ! This of course means throwing away all matrices, but we do not need
      ! them anymore - they are reassembled in the next sweep.
      do i=ilvmin,ilvmax
        p_rmatrix => collct_getvalue_mat (p_rcollection,"SYSTEMMAT",i)
        call lsyssc_sortMatrix (p_rmatrix%RmatrixBlock(1,1),.false.,.false.)
      end do

    end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine b1d6_solve (rproblem)
  
!<description>
  ! Solves the given problem by applying a nonlinear solver.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_vectorBlock), target :: rtempBlock
    type(t_vectorScalar), target :: rtempVectorSc
    type(t_matrixBlock), pointer :: p_rmatrix
    
    type(t_nlsolNode) :: rnlSol
    integer :: ierror,ilvmin,ilvmax,i
    type(t_linsolNode), pointer :: p_rsolverNode,p_rpreconditioner,p_rsmoother
    type(t_linsolNode), pointer :: p_rcoarseGridSolver

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(rproblem%ilvmax) :: Rmatrices

    ! A filter chain to filter the vectors and the matrix during the
    ! solution process.
    type(t_filterChain), dimension(1), target :: RfilterChain

    ! One level of multigrid
    type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo
    
    ! An interlevel projection structure for changing levels
    type(t_interlevelProjectionBlock) :: rprojection
    
    integer :: imaxmem
    
    ! Min/Max level?
    ilvmin = rproblem%ilvmin
    ilvmax = rproblem%ilvmax
    
    ! Now we have to build up the level information for multigrid.
    !
    ! At first, initialise a standard interlevel projection structure. We
    ! can use the same structure for all levels. Therefore it is enough
    ! to initialise one structure using the RHS vector on the finest
    ! level to specify the shape of the PDE-discretisation.
    call mlprj_initProjectionVec (rprojection,rproblem%rrhs)
    
    ! For solving the problem, we need to invoke a nonlinear solver.
    ! This nonlinear solver
    !  - builds in every step a linearisation of the system matrix
    !  - calls a linear solver for preconditioning
    ! We can save some time in this situation, if we prepare the
    ! linear solver in-advance. This means:
    !  - allocate memory in advance and
    !  - perform a symbolic factorisation of the matrix in advance.
    ! The point is that the entries of the matrix change in every
    ! iteration - but not the matrix structure! So this is something
    ! we can prepare once and use it through the whole solution process!
    !
    ! At first, set up the linear solver as usual.
    !
    ! During the linear solver, the boundary conditions must
    ! frequently be imposed to the vectors. This is done using
    ! a filter chain. As the linear solver does not work with
    ! the actual solution vectors but with defect vectors instead,
    ! a filter for implementing the real boundary conditions
    ! would be wrong.
    ! Therefore, create a filter chain with one filter only,
    ! which implements Dirichlet-conditions into a defect vector.
    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

    ! Create a Multigrid-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    call linsol_initMultigrid2 (p_rsolverNode,ilvmax-ilvmin+1,RfilterChain)
    
    ! Set the output level of the solver for some output
    p_rsolverNode%ioutputLevel = 1

    ! Configure MG to gain only one digit. We can do this, as MG is only
    ! used as preconditioner inside of another solver. Furthermore, perform
    ! at least two, at most 10 steps.
    p_rsolverNode%depsRel = 1E-1_DP
    p_rsolverNode%depsAbs = 0.0_DP
    p_rsolverNode%nminIterations = 2
    p_rsolverNode%nmaxIterations = 10
    
    ! Add the interlevel projection structure to the collection; we can
    ! use it later for setting up nonlinear matrices.
    call collct_setvalue_ilvp(rproblem%rcollection,"ILVPROJECTION",&
                              rprojection,.true.)
    
    ! Then set up smoothers / coarse grid solver:
    imaxmem = 0
    do i=ilvmin,ilvmax
      
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

        ! Set up an ILU smoother for multigrid with damping parameter 0.7,
        ! 2 smoothing steps (pre- and postsmoothing). As the problem is
        ! badly conditioned, we need even ILU(4) to get this problem solved,
        ! otherwise the smoother is diverging!
        ! note that if the trapezoidal rule is used for setting up the matrix,
        ! one would even need ILU(6)!!!
        call linsol_initMILUs1x1 (p_rsmoother,4,0.0_DP)
        call linsol_convertToSmoother (p_rsmoother,2,0.7_DP)
      end if
    
      ! Add the level.
      ! And add this multi-grid level. We will use the same smoother
      ! for pre- and post-smoothing.
      call linsol_getMultigrid2Level (p_rsolverNode,i-ilvmin+1,p_rlevelInfo)
      p_rlevelInfo%p_rcoarseGridSolver => p_rcoarseGridSolver
      p_rlevelInfo%p_rpresmoother => p_rsmoother
      p_rlevelInfo%p_rpostsmoother => p_rsmoother

      ! How much memory is necessary for performing the level change?
      ! We ourself must build nonlinear matrices on multiple levels and have
      ! to interpolate the solution vector from finer level to coarser ones.
      ! We need temporary memory for this purpose...
      if (i .gt. ilvmin) then
        ! Pass the system metrices on the coarse/fine grid to
        ! mlprj_getTempMemoryMat to specify the discretisation structures
        ! of all equations in the PDE there.
        imaxmem = max(imaxmem,mlprj_getTempMemoryMat (rprojection,&
                              rproblem%RlevelInfo(i-1)%rmatrix,&
                              rproblem%RlevelInfo(i)%rmatrix))
      end if

    end do
    
    ! Set up a scalar temporary vector that we need for building up nonlinear
    ! matrices. It must be at least as large as MAXMEM and NEQ(finest level),
    ! as we use it for resorting vectors as well.
    call lsyssc_createVector (rtempVectorSc,max(imaxmem,rproblem%rrhs%NEQ),.false.)
    call collct_setvalue_vecsca(rproblem%rcollection,"RTEMPSCALAR",rtempVectorSc,.true.)
    
    ! Before attaching the matrices to the solver and the initialisation of
    ! the problem structure, sort the matrix structure on all levels
    ! according top the associated permutation. Do not sort the
    ! entries - there are none!
    do i=ilvmin,ilvmax
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      call lsyssc_sortMatrix (p_rmatrix%RmatrixBlock(1,1),.true.,.false.)
    end do
    
    ! Attach the system matrices to the solver.
    !
    ! We copy our matrices to a big matrix array and transfer that
    ! to the setMatrices routines. This intitialises then the matrices
    ! on all levels according to that array.
    Rmatrices(ilvmin:ilvmax) = rproblem%RlevelInfo(ilvmin:ilvmax)%rmatrix
    call linsol_setMatrices(p_RsolverNode,Rmatrices(ilvmin:ilvmax))
    
    ! Initialise structure of the solver. This allows the
    ! solver to allocate memory / perform some precalculation
    ! to the problem.
    ! In fact, solvers like UMFPACK use this for a symbolic factorisation
    ! of the matrix.
    call linsol_initStructure (p_rsolverNode, ierror)
    
    if (ierror .ne. LINSOL_ERR_NOERROR) then
      call output_line("Matrix structure invalid!",OU_CLASS_ERROR)
      call sys_halt()
    end if

    ! Unsort the matrix structure again. The matrices stay in unsorted form
    ! until the entries are assembled.
    ! Remark: This makes the matrices inconsistent to those attached to the
    !  linear solver! So before invoking the linear solver, the matrices
    !  must be sorted to make their entries consistent!!!
    do i=ilvmin,ilvmax
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      call lsyssc_sortMatrix (p_rmatrix%RmatrixBlock(1,1),.false.,.false.)
    end do
    
    ! Put the prepared solver node to the collection for later use.
    call collct_setvalue_linsol(rproblem%rcollection,"LINSOLVER",p_rsolverNode,.true.)
    
    ! Create a temporary vector we need for the nonliner iteration.
    call lsysbl_createVecBlockIndirect (rproblem%rrhs, rtempBlock, .false.)

    ! The nonlinear solver structure rnlSol is initialised by the default
    ! initialisation with all necessary information to solve the problem.
    ! We call the nonlinear solver directly. For preconditioning
    ! and defect calculation, we use our own callback routine.
    rnlSol%ioutputLevel = 2
    call nlsol_performSolve(rnlSol,rproblem%rvector,rproblem%rrhs,rtempBlock,&
                            b1d6_getDefect,b1d6_precondDefect,&
                            rcollection=rproblem%rcollection)

    ! Release the temporary vector(s)
    call lsysbl_releaseVector (rtempBlock)
    call lsyssc_releaseVector (rtempVectorSc)
    
    ! Remove the solver node from the collection - not needed anymore there
    call collct_deletevalue(rproblem%rcollection,"LINSOLVER")
    
    ! Remove the temporary vector from the collection
    call collct_deletevalue(rproblem%rcollection,"RTEMPSCALAR")
    
    ! Remove the interlevel projection structure
    call collct_deletevalue(rproblem%rcollection,"ILVPROJECTION")
    
    ! Clean up the linear solver, release all memory, remove the solver node
    ! from memory.
    call linsol_releaseSolver (p_rsolverNode)
    
    ! Release the multilevel projection structure.
    call mlprj_doneProjection (rprojection)
    
    call output_lbrk()
    call output_line ("Nonlinear solver statistics")
    call output_line ("---------------------------")
    call output_line ("Initial defect: "//trim(sys_sdEL(rnlSol%DinitialDefect(1),15)))
    call output_line ("Final defect:  "//trim(sys_sdEL(rnlSol%DfinalDefect(1),15)))
    call output_line ("#Iterations:   "//trim(sys_siL(rnlSol%iiterations,10)))

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine b1d6_postprocessing (rproblem)
  
!<description>
  ! Writes the solution into a GMV file.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  
    ! We need some more variables for postprocessing
    real(DP), dimension(:), pointer :: p_Ddata
    character(len=SYS_STRLEN) :: sucddir
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport

    ! A cubature information structure
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo

    ! A pointer to the solution vector and to the triangulation.
    type(t_vectorBlock), pointer :: p_rvector
    type(t_triangulation), pointer :: p_rtriangulation

    ! Get the solution vector from the problem structure.
    p_rvector => rproblem%rvector
    
    ! From the attached discretisation, get the underlying triangulation
    p_rtriangulation => &
      p_rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation
    
    ! Get the cubature information structure
    p_rcubatureInfo => rproblem%RlevelInfo(rproblem%ilvmax)%rcubatureInfo
    
    ! p_rvector now contains our solution. We can now
    ! start the postprocessing.
    !
    ! Get the path for writing postprocessing files from the environment variable
    ! $UCDDIR. If that does not exist, write to the directory "./gmv".
    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = "./gmv"

    ! Start UCD export to VTK file:
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,p_rtriangulation,trim(sucddir)//"/u6.vtk")
    
    call lsyssc_getbase_double (p_rvector%RvectorBlock(1),p_Ddata)
    call ucd_addVariableVertexBased (rexport,"sol",UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that is it.
    call ucd_write (rexport)
    call ucd_release (rexport)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine b1d6_doneMatVec (rproblem)
  
!<description>
  ! Releases system matrix and vectors.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    integer :: i

    ! Release matrices and vectors on all levels
    do i=rproblem%ilvmax,rproblem%ilvmin,-1
      ! Delete the matrix
      call lsysbl_releaseMatrix (rproblem%RlevelInfo(i)%rmatrix)

      ! Delete the variables from the collection.
      call collct_deletevalue (rproblem%rcollection,"SYSTEMMAT",i)
      
      ! Remove the temp vector that was used for interpolating the solution
      ! from higher to lower levels in the nonlinear iteration.
      if (i .lt. rproblem%ilvmax) then
        call lsysbl_releaseVector(rproblem%RlevelInfo(i)%rtempVector)
        call collct_deletevalue(rproblem%rcollection,"RTEMPVEC",i)
      end if
      
    end do

    ! Delete solution/RHS vector
    call lsysbl_releaseVector (rproblem%rvector)
    call lsysbl_releaseVector (rproblem%rrhs)

    ! Release the sorting strategy
    call sstrat_doneBlockSorting (rproblem%RlevelInfo(1)%rsortStrategy)

    ! Delete the variables from the collection.
    call collct_deletevalue (rproblem%rcollection,"RHS")
    call collct_deletevalue (rproblem%rcollection,"SOLUTION")

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine b1d6_doneBC (rproblem)
  
!<description>
  ! Releases discrete boundary conditions from the heap.
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

  subroutine b1d6_doneDiscretisation (rproblem)
  
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
      call collct_deletevalue (rproblem%rcollection,"CUBINFO",i)

      ! Delete the block discretisation together with the associated
      ! scalar spatial discretisations.
      call spdiscr_releaseBlockDiscr(rproblem%RlevelInfo(i)%rdiscretisation)

    end do
    
  end subroutine
    
  ! ***************************************************************************

!<subroutine>

  subroutine b1d6_doneParamTriang (rproblem)
  
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
    
    ! Finally release the domain.
    call boundary_release (rproblem%rboundary)
    
    call collct_deleteValue(rproblem%rcollection,"NLMAX")
    call collct_deleteValue(rproblem%rcollection,"NLMIN")

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine burgers1d6
  
!<description>
  ! This is a "separated" burgers1d solver for solving a Burgers-1D
  ! problem. The different tasks of the problem are separated into
  ! subroutines. The problem uses a problem-specific structure for the
  ! communication: All subroutines add their generated information to the
  ! structure, so that the other subroutines can work with them.
  ! For the communication to callback routines of black-box subroutines
  ! (matrix-assembly, preconditioning), a collection is used.
  !
  ! The following tasks are performed by the subroutines:
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

    ! NLMIN receives the minimal level where to discretise for supporting
    ! the solution process.
    ! NLMAX receives the level where we want to solve.
    integer :: NLMIN,NLMAX
    
    ! A problem structure for our problem
    type(t_problem), pointer :: p_rproblem
    
    integer :: i
    
    ! Ok, let us start.
    !
    ! We want to solve our Laplace problem on level...

    NLMIN = 3
    NLMAX = 7
    
    ! Allocate the problem structure -- it is rather large
    allocate(p_rproblem)
    allocate(p_rproblem%RlevelInfo(1:NLMAX))

    ! Initialise the collection
    call collct_init (p_rproblem%rcollection)
    do i=1,NLMAX
      call collct_addlevel (p_rproblem%rcollection)
    end do

    ! So now the different steps - one after the other.
    
    ! Initialisation.
    call b1d6_initParamTriang (NLMIN,NLMAX,p_rproblem)
    call b1d6_initDiscretisation (p_rproblem)
    call b1d6_initMatVec (p_rproblem)
    call b1d6_initDiscreteBC (p_rproblem)
    
    ! Implementation of boundary conditions
    call b1d6_implementBC (p_rproblem)
    
    ! Solve the problem
    call b1d6_solve (p_rproblem)
    
    ! Postprocessing
    call b1d6_postprocessing (p_rproblem)
    
    ! Cleanup
    call b1d6_doneMatVec (p_rproblem)
    call b1d6_doneBC (p_rproblem)
    call b1d6_doneDiscretisation (p_rproblem)
    call b1d6_doneParamTriang (p_rproblem)
    
    ! Print some statistical data about the collection - anything forgotten?
    print *
    print *,"Remaining collection statistics:"
    print *,"--------------------------------"
    print *
    call collct_printStatistics (p_rproblem%rcollection)
    
    ! Finally release the collection and the problem structure.
    call collct_done (p_rproblem%rcollection)
    
    deallocate(p_rproblem%RlevelInfo)
    deallocate(p_rproblem)
    
  end subroutine

end module
