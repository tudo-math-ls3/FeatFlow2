!##############################################################################
!# ****************************************************************************
!# <name> burgers1d_method5 </name>
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
!# For Solving this in the domain $\Omega$, we follow the usual way:
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
!# The preconditioning in this example is done by a UMFPACK4 Gauss
!# elimination.
!#
!# </purpose>
!##############################################################################

module burgers1d_method5

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
  use sortstrategy
  use nonlinearsolver
  use ucd
  
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

    ! A variable describing the discrete boundary conditions.
    type(t_discreteBC) :: rdiscreteBC
  
  end type
  
!</typeblock>


!<typeblock description="Application-specific type block for burgers1d problem">

  type t_problem
  
    ! Maximum refinement level = level where the system is solved
    integer :: ilvmax

    ! An object for saving the domain:
    type(t_boundary) :: rboundary

    ! A solution vector and a RHS vector on the finest level.
    type(t_vectorBlock) :: rvector,rrhs

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode

    ! An array of t_problem_lvl structures, each corresponding
    ! to one level of the discretisation. There is currently
    ! only one level supported, identified by LV!
    type(t_problem_lvl) :: rlevelInfo
    
    ! A collection structure with problem-dependent data
    type(t_collection) :: rcollection
    
  end type

!</typeblock>

!</types>
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine b1d5_initParamTriang (ilvmax,rproblem)
  
!<description>
  ! This routine initialises the parametrisation and triangulation of the
  ! domain. The corresponding .prm/.tri files are read from disc and
  ! the triangulation is refined as described by the parameter ilv.
!</description>

!<input>
  ! Maximum refinement level
  integer, intent(in) :: ilvmax
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

!</subroutine>

    ! Path to the mesh
    character(len=SYS_STRLEN) :: spredir

    ! Initialise the level in the problem structure
    rproblem%ilvmax = ilvmax
    
    ! Store the min/max level in the collection to be used in
    ! callback routines
    call collct_setvalue_int(rproblem%rcollection,"NLMAX",ilvmax,.true.)

    ! Get the path $PREDIR from the environment, where to read .prm/.tri files
    ! from. If that does not exist, write to the directory "./pre".
    if (.not. sys_getenv_string("PREDIR", spredir)) spredir = "./pre"

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call boundary_read_prm(rproblem%rboundary, trim(spredir)//"/QUAD.prm")
        
    ! Now read in the basic triangulation.
    call tria_readTriFile2D (rproblem%rlevelInfo%rtriangulation, &
        trim(spredir)//"/QUAD.tri", rproblem%rboundary)
    
    ! Refine the mesh up to the maximum level
    call tria_quickRefine2LevelOrdering(rproblem%ilvmax-1,&
        rproblem%rlevelInfo%rtriangulation,rproblem%rboundary)
    
    ! Create information about adjacencies and everything one needs from
    ! a triangulation. Afterwards, we have the coarse mesh.
    call tria_initStandardMeshFromRaw (&
        rproblem%rlevelInfo%rtriangulation,rproblem%rboundary)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine b1d5_initDiscretisation (rproblem)
  
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

    ! An object for the spatial discretisation
    type(t_blockDiscretisation), pointer :: p_rdiscretisation

    i=rproblem%ilvmax

    ! Ask the problem structure to give us the boundary and triangulation.
    ! We need it for the discretisation.
    p_rboundary => rproblem%rboundary
    p_rtriangulation => rproblem%rlevelInfo%rtriangulation
    
    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we only have one block.
    allocate(p_rdiscretisation)
    call spdiscr_initBlockDiscr (rproblem%rlevelInfo%rdiscretisation,1,&
                                 p_rtriangulation, p_rboundary)

    ! Get the discretisation
    p_rdiscretisation => rproblem%rlevelInfo%rdiscretisation

    ! p_rdiscretisation%Rdiscretisations is a list of scalar
    ! discretisation structures for every component of the solution vector.
    ! Initialise the first element of the list to specify the element
    ! for this solution component:
    call spdiscr_initDiscr_simple ( &
                 p_rdiscretisation%RspatialDiscr(1), &
                 EL_Q1,p_rtriangulation, p_rboundary)

    ! Set up an cubature info structure to tell the code which cubature
    ! formula to use
    ! Create an assembly information structure which tells the code
    ! the cubature formula to use. Standard: Gauss 3x3.
    call spdiscr_createDefCubStructure(&  
        p_rdiscretisation%RspatialDiscr(1),rproblem%rlevelInfo%rcubatureInfo,&
        CUB_GEN_AUTO_G3)
        
    call collct_setvalue_cubinfo(rproblem%rcollection,"CUBINFO",&
        rproblem%rlevelInfo%rcubatureInfo,.true.,i)
                                   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine b1d5_initMatVec (rproblem)
  
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
    type(t_vectorBlock), pointer :: p_rrhs,p_rvector

    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
  
    ! A cubature information structure
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo

    ! Arrays for the Cuthill McKee renumbering strategy
    integer, dimension(1) :: H_Iresort
    integer, dimension(:), pointer :: p_Iresort

    i=rproblem%ilvmax

    ! Ask the problem structure to give us the discretisation structure
    p_rdiscretisation => rproblem%rlevelInfo%rdiscretisation
    
    p_rmatrix => rproblem%rlevelInfo%rmatrix
    
    ! Get the cubature information structure
    p_rcubatureInfo => rproblem%rlevelInfo%rcubatureInfo
           
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
    ! in every step of the nonlinear iteration!
    ! We fill the matrix with 1. This is necessary, as the UMFPACK solver
    ! needs nonzero matrix entries for the symbolic factorisation!
    call lsyssc_allocEmptyMatrix(p_rmatrix%RmatrixBlock(1,1),LSYSSC_SETM_ONE)
    
    ! Allocate an array for holding the resorting strategy.
    call storage_new ("b1d5_initMatVec", "Iresort", &
          p_rmatrix%RmatrixBlock(1,1)%NEQ*2, ST_INT, h_Iresort(1), ST_NEWBLOCK_ZERO)
    call storage_getbase_int(h_Iresort(1),p_Iresort)
    
    ! Calculate the resorting strategy.
    call sstrat_calcCuthillMcKee (p_rmatrix%RmatrixBlock(1,1),p_Iresort)
    
    ! Save the handle of the resorting strategy to the collection.
    call collct_setvalue_int(rproblem%rcollection,"LAPLACE-CM",h_Iresort(1),.true.,i)
    
    ! We do not resort the matrix yet - this is done later when the entries
    ! are assembled.
      
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
    ! Note that the vector is unsorted when this call finishes!
    call linf_buildVectorScalar (&
        rlinform,.true.,p_rrhs%RvectorBlock(1),p_rcubatureInfo,&
        coeff_RHS,rproblem%rcollection)
                                
    ! Clear the solution vector on the finest level.
    call lsysbl_clearVector(rproblem%rvector)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine b1d5_initDiscreteBC (rproblem)
  
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
    integer :: ilvmax
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
    
    ilvmax=rproblem%ilvmax
    
    ! Get our matrix from the problem structure.
    p_rmatrix => rproblem%rlevelInfo%rmatrix
    
    ! From the matrix or the RHS we have access to the discretisation and the
    ! analytic boundary conditions.
    p_rdiscretisation => p_rmatrix%p_rblockDiscrTest
    
    ! Get the domain from the problem structure
    p_rboundary => rproblem%rboundary

    ! Now we have the raw problem. What is missing is the definition of the boundary
    ! conditions.
    ! For implementing boundary conditions, we use a `filter technique with
    ! discretised boundary conditions`. This means, we first have to calculate
    ! a discrete version of the analytic BC, which we can implement into the
    ! solution/RHS vectors using the corresponding filter.
    !
    ! Create a t_discreteBC structure where we store all discretised boundary
    ! conditions.
    call bcasm_initDiscreteBC(rproblem%rlevelInfo%rdiscreteBC)
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
       rboundaryRegion,rproblem%rlevelInfo%rdiscreteBC,&
       getBoundaryValues,rproblem%rcollection)
                             
    ! Now to the edge 2 of boundary component 1 the domain.
    call boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    call bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
       rboundaryRegion,rproblem%rlevelInfo%rdiscreteBC,&
       getBoundaryValues,rproblem%rcollection)
                             
    ! Edge 3 of boundary component 1.
    ! Edge 3 must be set up as Neumann boundary, which is realised as
    ! simple "do-$nothing"-boundary conditions. So we do not do anything with edge 3!
    ! CALL boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
    ! CALL bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
    !    rboundaryRegion,rproblem%rlevelInfo%rdiscreteBC,&
    !    getBoundaryValues,rproblem%rcollection)
    
    ! Edge 4 of boundary component 1. That is it.
    call boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
       rboundaryRegion,rproblem%rlevelInfo%rdiscreteBC,&
       getBoundaryValues,rproblem%rcollection)
                             
    ! Hang the pointer into the vectors and the matrix. That way, these
    ! boundary conditions are always connected to that matrix and that
    ! vector.
    p_rdiscreteBC => rproblem%rlevelInfo%rdiscreteBC
    
    p_rmatrix%p_rdiscreteBC => p_rdiscreteBC
      
    ! On the finest level, attach the discrete BC also
    ! to the solution and RHS vector. They need it to be compatible
    ! to the matrix on the finest level.
    p_rdiscreteBC => rproblem%rlevelInfo%rdiscreteBC
    
    p_rrhs    => rproblem%rrhs
    p_rvector => rproblem%rvector
    
    call lsysbl_assignDiscreteBC(p_rrhs,p_rdiscreteBC)
    call lsysbl_assignDiscreteBC(p_rvector,p_rdiscreteBC)
                
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine b1d5_implementBC (rproblem)
  
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
    subroutine b1d5_getDefect (ite,rx,rb,rd,p_rcollection)
  
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

    subroutine b1d5_precondDefect (ite,rd,rx,rb,domega,bsuccess,p_rcollection)
  
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
    integer :: ierror
    type(t_linsolNode), pointer :: p_rsolverNode

      ! A cubature information structure
      type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo

      ! Our "parent" (the caller of the nonlinear solver) has prepared
      ! a preconditioner node for us (a linear solver with symbolically
      ! factorised matrices). Get this from the collection.
      
      p_rsolverNode => collct_getvalue_linsol(p_rcollection,"LINSOLVER")

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

    end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine b1d5_solve (rproblem)
  
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
    
    type(t_nlsolNode) :: rnlSol
    integer :: ierror
    type(t_linsolNode), pointer :: p_rsolverNode
    type(t_matrixBlock), dimension(1) :: Rmatrices
    type(t_matrixBlock), pointer :: p_rmatrix
    
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
    ! At first, set up the linear solver as usual:
    call linsol_initUMFPACK4 (p_rsolverNode)

    ! Get the system matrix on the finest level...
    p_rmatrix => rproblem%rlevelInfo%rmatrix

    ! And associate it to the solver
    Rmatrices = (/p_rmatrix/)
    call linsol_setMatrices(p_rsolverNode,Rmatrices)

    ! Initialise structure of the solver. This allows the
    ! solver to allocate memory / perform some precalculation
    ! to the problem.
    ! In fact, solvers like UMFPACK use this for a symbolic factorisation
    ! of the matrix.
    call linsol_initStructure (p_rsolverNode, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) stop
  
    ! Put the prepared solver node to the collection for later use.
    call collct_setvalue_linsol(rproblem%rcollection,"LINSOLVER",p_rsolverNode,.true.)
    
    ! Create a temporary vector we need for the nonlinera iteration.
    call lsysbl_createVecBlockIndirect (rproblem%rrhs, rtempBlock, .false.)

    ! The nonlinear solver structure rnlSol is initialised by the default
    ! initialisation with all necessary information to solve the problem.
    ! We call the nonlinear solver directly. For preconditioning
    ! and defect calculation, we use our own callback routine.
    !rnlSol%domega = 0.25_DP
    rnlSol%ioutputLevel = 2
    call nlsol_performSolve(rnlSol,rproblem%rvector,rproblem%rrhs,rtempBlock,&
                            b1d5_getDefect,b1d5_precondDefect,&
                            rcollection=rproblem%rcollection)

    ! Release the temporary vector
    call lsysbl_releaseVector (rtempBlock)
    
    ! Remove the solver node from the collection - not needed anymore there
    call collct_deletevalue(rproblem%rcollection,"LINSOLVER")
    
    ! Clean up the linear solver, release all memory, remove the solver node
    ! from memory.
    call linsol_releaseSolver (p_rsolverNode)
    
    call output_lbrk()
    call output_line ("Nonlinear solver statistics")
    call output_line ("---------------------------")
    call output_line ("Initial defect: "//trim(sys_sdEL(rnlSol%DinitialDefect(1),15)))
    call output_line ("Final defect:  "//trim(sys_sdEL(rnlSol%DfinalDefect(1),15)))
    call output_line ("#Iterations:   "//trim(sys_siL(rnlSol%iiterations,10)))

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine b1d5_postprocessing (rproblem)
  
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
    
    ! Output block for UCD output to VTK file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir

    ! A pointer to the solution vector and to the triangulation.
    type(t_vectorBlock), pointer :: p_rvector
    type(t_triangulation), pointer :: p_rtriangulation
    
    ! Cubature information structure
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo

    ! Get the solution vector from the problem structure.
    p_rvector => rproblem%rvector
    
    ! From the attached discretisation, get the underlying triangulation
    p_rtriangulation => &
      p_rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation
    
    ! Get the cubature information structure
    p_rcubatureInfo => rproblem%rlevelInfo%rcubatureInfo
    
    ! p_rvector now contains our solution. We can now
    ! start the postprocessing.
    !
    ! Get the path for writing postprocessing files from the environment variable
    ! $UCDDIR. If that does not exist, write to the directory "./gmv".
    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = "./gmv"

    ! Start UCD export to VTK file:
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,p_rtriangulation,trim(sucddir)//"/u5.vtk")
    
    call lsyssc_getbase_double (p_rvector%RvectorBlock(1),p_Ddata)
    call ucd_addVariableVertexBased (rexport,"sol",UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that is it.
    call ucd_write (rexport)
    call ucd_release (rexport)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine b1d5_doneMatVec (rproblem)
  
!<description>
  ! Releases system matrix and vectors.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    integer :: ihandle,ilvmax

    ! Release matrix and vectors on all levels
    ilvmax=rproblem%ilvmax

    ! Delete the matrix
    call lsysbl_releaseMatrix (rproblem%rlevelInfo%rmatrix)

    ! Delete the variables from the collection.
    call collct_deletevalue (rproblem%rcollection,"SYSTEMMAT",ilvmax)
    
    ! Release the permutation for sorting matrix/vectors
    ihandle = collct_getvalue_int (rproblem%rcollection,"LAPLACE-CM",ilvmax)
    call storage_free (ihandle)
    call collct_deletevalue (rproblem%rcollection,"LAPLACE-CM",ilvmax)

    ! Delete solution/RHS vector
    call lsysbl_releaseVector (rproblem%rvector)
    call lsysbl_releaseVector (rproblem%rrhs)

    ! Delete the variables from the collection.
    call collct_deletevalue (rproblem%rcollection,"RHS")
    call collct_deletevalue (rproblem%rcollection,"SOLUTION")

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine b1d5_doneBC (rproblem)
  
!<description>
  ! Releases discrete boundary conditions from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: ilvmax

    ilvmax=rproblem%ilvmax
      
    ! Release our discrete version of the boundary conditions
    call bcasm_releaseDiscreteBC (rproblem%rlevelInfo%rdiscreteBC)
    
  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine b1d5_doneDiscretisation (rproblem)
  
!<description>
  ! Releases the discretisation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: ilvmax

    ilvmax=rproblem%ilvmax

    ! Release the cubature info structure.
    call spdiscr_releaseCubStructure(rproblem%rlevelInfo%rcubatureInfo)
    call collct_deletevalue (rproblem%rcollection,"CUBINFO",ilvmax)

    ! Delete the block discretisation together with the associated
    ! scalar spatial discretisations....
    call spdiscr_releaseBlockDiscr(rproblem%rlevelInfo%rdiscretisation)

  end subroutine
    
  ! ***************************************************************************

!<subroutine>

  subroutine b1d5_doneParamTriang (rproblem)
  
!<description>
  ! Releases the triangulation and parametrisation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    ! Release the triangulation
    call tria_done (rproblem%rlevelInfo%rtriangulation)
    
    ! Finally release the domain.
    call boundary_release (rproblem%rboundary)
    
    call collct_deleteValue(rproblem%rcollection,"NLMAX")

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine burgers1d5
  
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

    ! NLMAX receives the level where we want to solve
    integer :: NLMAX
    
    ! A problem structure for our problem
    type(t_problem), pointer :: p_rproblem
    
    integer :: i
    
    ! Ok, let us start.
    !
    ! We want to solve our Laplace problem on level...
    NLMAX = 7
    
    ! Allocate the problem structure -- it is rather large
    allocate(p_rproblem)

    ! Initialise the collection
    call collct_init (p_rproblem%rcollection)
    do i=1,NLMAX
      call collct_addlevel (p_rproblem%rcollection)
    end do

    ! So now the different steps - one after the other.
    
    ! Initialisation.
    call b1d5_initParamTriang (NLMAX,p_rproblem)
    call b1d5_initDiscretisation (p_rproblem)
    call b1d5_initMatVec (p_rproblem)
    call b1d5_initDiscreteBC (p_rproblem)
    
    ! Implementation of boundary conditions
    call b1d5_implementBC (p_rproblem)
    
    ! Solve the problem
    call b1d5_solve (p_rproblem)
    
    ! Postprocessing
    call b1d5_postprocessing (p_rproblem)
    
    ! Cleanup
    call b1d5_doneMatVec (p_rproblem)
    call b1d5_doneBC (p_rproblem)
    call b1d5_doneDiscretisation (p_rproblem)
    call b1d5_doneParamTriang (p_rproblem)

    ! Print some statistical data about the collection - anything forgotten?
    print *
    print *,"Remaining collection statistics:"
    print *,"--------------------------------"
    print *
    call collct_printStatistics (p_rproblem%rcollection)
    
    ! Finally release the collection and the problem structure.
    call collct_done (p_rproblem%rcollection)

    deallocate(p_rproblem)
    
  end subroutine

end module
