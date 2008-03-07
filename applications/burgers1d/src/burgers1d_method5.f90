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
!# discretisation, i.e. we don't use a separate time-discretisation.
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

MODULE burgers1d_method5

  USE fsystem
  USE storage
  USE linearsolver
  USE boundary
  USE bilinearformevaluation
  USE linearformevaluation
  USE cubature
  USE matrixfilters
  USE vectorfilters
  USE bcassembly
  USE triangulation
  USE spatialdiscretisation
  USE sortstrategy
  USE nonlinearsolver
  USE ucd
  
  USE collection
    
  USE burgers1d_callback
  
  IMPLICIT NONE
  
!<types>

!<typeblock description="Type block defining all information about one level">

  TYPE t_problem_lvl
  
    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation (trial/test functions,...)
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
    
    ! A system matrix for that specific level. The matrix will receive the 
    ! discrete Laplace operator.
    TYPE(t_matrixBlock) :: rmatrix

    ! A variable describing the discrete boundary conditions.    
    TYPE(t_discreteBC) :: rdiscreteBC
  
  END TYPE
  
!</typeblock>


!<typeblock description="Application-specific type block for burgers1d problem">

  TYPE t_problem
  
    ! Maximum refinement level = level where the system is solved
    INTEGER :: ilvmax

    ! An object for saving the domain:
    TYPE(t_boundary), POINTER :: p_rboundary

    ! A solution vector and a RHS vector on the finest level. 
    TYPE(t_vectorBlock) :: rvector,rrhs

    ! A solver node that accepts parameters for the linear solver    
    TYPE(t_linsolNode), POINTER :: p_rsolverNode

    ! An array of t_problem_lvl structures, each corresponding
    ! to one level of the discretisation. There is currently
    ! only one level supported, identified by LV!
    TYPE(t_problem_lvl) :: rlevelInfo
    
    ! A collection structure with problem-dependent data
    TYPE(t_collection) :: rcollection
    
  END TYPE

!</typeblock>

!</types>
  
CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE b1d5_initParamTriang (ilvmax,rproblem)
  
!<description>
  ! This routine initialises the parametrisation and triangulation of the
  ! domain. The corresponding .prm/.tri files are read from disc and
  ! the triangulation is refined as described by the parameter ilv.
!</description>

!<input>
  ! Maximum refinement level
  INTEGER, INTENT(IN) :: ilvmax
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

    ! Initialise the level in the problem structure
    rproblem%ilvmax = ilvmax
    
    ! Store the min/max level in the collection to be used in 
    ! callback routines
    CALL collct_setvalue_int(rproblem%rcollection,'NLMAX',ilvmax,.TRUE.)

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    ! Set p_rboundary to NULL() to create a new structure.
    NULLIFY(rproblem%p_rboundary)
    CALL boundary_read_prm(rproblem%p_rboundary, './pre/QUAD.prm')
        
    ! Now read in the basic triangulation.
    CALL tria_readTriFile2D (rproblem%rlevelInfo%rtriangulation, &
        './pre/QUAD.tri', rproblem%p_rboundary)
    
    ! Refine the mesh up to the maximum level
    CALL tria_quickRefine2LevelOrdering(rproblem%ilvmax-1,&
        rproblem%rlevelInfo%rtriangulation,rproblem%p_rboundary)
    
    ! Create information about adjacencies and everything one needs from
    ! a triangulation. Afterwards, we have the coarse mesh.
    CALL tria_initStandardMeshFromRaw (&
        rproblem%rlevelInfo%rtriangulation,rproblem%p_rboundary)
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE b1d5_initDiscretisation (rproblem)
  
!<description>
  ! This routine initialises the discretisation structure of the underlying
  ! problem and saves it to the problem structure.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: I
  
    ! An object for saving the domain:
    TYPE(t_boundary), POINTER :: p_rboundary
    
    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation), POINTER :: p_rtriangulation
    
    ! An object for the spatial discretisation
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation

    i=rproblem%ilvmax
      
    ! Ask the problem structure to give us the boundary and triangulation.
    ! We need it for the discretisation.
    p_rboundary => rproblem%p_rboundary
    p_rtriangulation => rproblem%rlevelInfo%rtriangulation
    
    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we only have one block.
    ALLOCATE(p_rdiscretisation)
    CALL spdiscr_initBlockDiscr2D (p_rdiscretisation,1,&
                                   p_rtriangulation, p_rboundary)
                                   
    ! Save the discretisation structure to our local LevelInfo structure
    ! for later use.
    rproblem%rlevelInfo%p_rdiscretisation => p_rdiscretisation

    ! p_rdiscretisation%Rdiscretisations is a list of scalar 
    ! discretisation structures for every component of the solution vector.
    ! Initialise the first element of the list to specify the element
    ! and cubature rule for this solution component:
    CALL spdiscr_initDiscr_simple ( &
                 p_rdiscretisation%RspatialDiscretisation(1), &
                 EL_E011,CUB_G2X2, &
                 p_rtriangulation, p_rboundary)
                                   
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE b1d5_initMatVec (rproblem)
  
!<description>
  ! Calculates the RHS-ector, set up the solution vetor.
  ! Set up the structure of the system matrix/matrices of the linear
  ! system.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: i
    
    ! A bilinear and linear form describing the analytic problem to solve
    TYPE(t_linearForm) :: rlinform
    
    ! A pointer to the system matrix and the RHS/solution vectors.
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    TYPE(t_vectorBlock), POINTER :: p_rrhs,p_rvector

    ! A pointer to the discretisation structure with the data.
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
  
    ! Arrays for the Cuthill McKee renumbering strategy
    INTEGER, DIMENSION(1) :: H_Iresort 
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Iresort
    
    i=rproblem%ilvmax
      
    ! Ask the problem structure to give us the discretisation structure
    p_rdiscretisation => rproblem%rlevelInfo%p_rdiscretisation
    
    p_rmatrix => rproblem%rlevelInfo%rmatrix
    
    ! Initialise the block matrix with default values based on
    ! the discretisation.
    CALL lsysbl_createMatBlockByDiscr (p_rdiscretisation,p_rmatrix)    

    ! Save matrix to the collection.
    ! They maybe used later, expecially in nonlinear problems.
    CALL collct_setvalue_mat(rproblem%rcollection,'SYSTEMMAT',p_rmatrix,.TRUE.,i)

    ! Now using the discretisation, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! We create that directly in the block (1,1) of the block matrix
    ! using the discretisation structure of the first block.
    CALL bilf_createMatrixStructure (&
              p_rdiscretisation%RspatialDiscretisation(1),LSYSSC_MATRIX9,&
              p_rmatrix%RmatrixBlock(1,1))
                                    
    ! Update the structural information of the block matrix, as we manually
    ! changed one of the submatrices:
    CALL lsysbl_updateMatStrucInfo (p_rmatrix)
    
    ! Allocate memory for the matrix, don't calculate the entries.
    ! Remember hat we have a nonlinear matrix, which entries must be build
    ! in every step of the nonlinear iteration!
    ! We fill the matrix with 1. This is necessary, as the UMFPACK solver
    ! needs nonzero matrix entries for the symbolic factorisation!
    CALL lsyssc_allocEmptyMatrix(p_rmatrix%RmatrixBlock(1,1),LSYSSC_SETM_ONE)
    
    ! Allocate an array for holding the resorting strategy.
    CALL storage_new ('b1d5_initMatVec', 'Iresort', &
          p_rmatrix%RmatrixBlock(1,1)%NEQ*2, ST_INT, h_Iresort(1), ST_NEWBLOCK_ZERO)
    CALL storage_getbase_int(h_Iresort(1),p_Iresort)
    
    ! Calculate the resorting strategy.
    CALL sstrat_calcCuthillMcKee (p_rmatrix%RmatrixBlock(1,1),p_Iresort)
    
    ! Save the handle of the resorting strategy to the collection.
    CALL collct_setvalue_int(rproblem%rcollection,'LAPLACE-CM',h_Iresort(1),.TRUE.,i)
    
    ! We don't resort the matrix yet - this is done later when the entries
    ! are assembled.
      
    ! (Only) on the finest level, we need to calculate a RHS vector
    ! and to allocate a solution vector.
    
    p_rrhs    => rproblem%rrhs   
    p_rvector => rproblem%rvector

    ! Although we could manually create the solution/RHS vector,
    ! the easiest way to set up the vector structure is
    ! to create it by using our matrix as template:
    CALL lsysbl_createVecBlockIndMat (p_rmatrix,p_rrhs, .FALSE.)
    CALL lsysbl_createVecBlockIndMat (p_rmatrix,p_rvector, .FALSE.)

    ! Save the solution/RHS vector to the collection. Might be used
    ! later (e.g. in nonlinear problems)
    CALL collct_setvalue_vec(rproblem%rcollection,'RHS',p_rrhs,.TRUE.)
    CALL collct_setvalue_vec(rproblem%rcollection,'SOLUTION',p_rvector,.TRUE.)
    
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
    CALL linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscretisation(1),rlinform,.TRUE.,&
              p_rrhs%RvectorBlock(1),coeff_RHS,&
              rproblem%rcollection)
                                
    ! Clear the solution vector on the finest level.
    CALL lsysbl_clearVector(rproblem%rvector)
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE b1d5_initDiscreteBC (rproblem)
  
!<description>
  ! This calculates the discrete version of the boundary conditions and
  ! assigns it to the system matrix and RHS vector.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: ilvmax
    TYPE(t_boundaryRegion) :: rboundaryRegion

    ! A pointer to the system matrix and the RHS vector as well as 
    ! the discretisation
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    TYPE(t_vectorBlock), POINTER :: p_rrhs,p_rvector
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation

    ! Pointer to structure for saving discrete BC's:
    TYPE(t_discreteBC), POINTER :: p_rdiscreteBC
      
    ! A pointer to the domain
    TYPE(t_boundary), POINTER :: p_rboundary
    
    ilvmax=rproblem%ilvmax
    
    ! Get our matrix from the problem structure.
    p_rmatrix => rproblem%rlevelInfo%rmatrix
    
    ! From the matrix or the RHS we have access to the discretisation and the
    ! analytic boundary conditions.
    p_rdiscretisation => p_rmatrix%p_rblockDiscretisation
    
    ! Get the domain from the problem structure
    p_rboundary => rproblem%p_rboundary

    ! Now we have the raw problem. What is missing is the definition of the boudary
    ! conditions.
    ! For implementing boundary conditions, we use a 'filter technique with
    ! discretised boundary conditions'. This means, we first have to calculate
    ! a discrete version of the analytic BC, which we can implement into the
    ! solution/RHS vectors using the corresponding filter.
    !
    ! Create a t_discreteBC structure where we store all discretised boundary
    ! conditions.
    CALL bcasm_initDiscreteBC(rproblem%rlevelInfo%rdiscreteBC)
    !
    ! We 'know' already (from the problem definition) that we have four boundary
    ! segments in the domain. Each of these, we want to use for enforcing
    ! some kind of boundary condition.
    !
    ! We ask the bondary routines to create a 'boundary region' - which is
    ! simply a part of the boundary corresponding to a boundary segment.
    ! A boundary region roughly contains the type, the min/max parameter value
    ! and whether the endpoints are inside the region or not.
    CALL boundary_createRegion(p_rboundary,1,1,rboundaryRegion)
    
    ! We use this boundary region and specify that we want to have Dirichlet
    ! boundary there. The following call does the following:
    ! - Create Dirichlet boundary conditions on the region rboundaryRegion.
    !   We specify icomponent='1' to indicate that we set up the
    !   Dirichlet BC's for the first (here: one and only) component in the 
    !   solution vector.
    ! - Discretise the boundary condition so that the BC's can be applied
    !   to matrices and vectors
    ! - Add the calculated discrete BC's to rdiscreteBC for later use.
    CALL bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
       rboundaryRegion,rproblem%rlevelInfo%rdiscreteBC,&
       getBoundaryValues,rproblem%rcollection)
                             
    ! Now to the edge 2 of boundary component 1 the domain. We use the
    ! same two routines to add the boundary condition to p_rboundaryConditions.
    CALL boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    CALL bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
       rboundaryRegion,rproblem%rlevelInfo%rdiscreteBC,&
       getBoundaryValues,rproblem%rcollection)
                             
    ! Edge 3 of boundary component 1.
    ! Edge 3 must be set up as Neumann boundary, which is realised as
    ! simple 'do-$nothing'-boundary conditions. So we don't do anything with edge 3!
    ! CALL boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
    ! CALL bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
    !    rboundaryRegion,rproblem%rlevelInfo%rdiscreteBC,&
    !    getBoundaryValues,rproblem%rcollection)
    
    ! Edge 4 of boundary component 1. That's it.
    CALL boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
    CALL bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
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
    p_rrhs    => rproblem%rrhs   
    p_rvector => rproblem%rvector
    
    p_rrhs%p_rdiscreteBC => p_rdiscreteBC
    p_rvector%p_rdiscreteBC => p_rdiscreteBC
                
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE b1d5_implementBC (rproblem)
  
!<description>
  ! Implements boundary conditions into the RHS and into a given solution vector.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: ilvmax
  
    ! A pointer to the RHS/solution vector 
    TYPE(t_vectorBlock), POINTER :: p_rrhs,p_rvector
    
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
    CALL vecfil_discreteBCrhs (p_rrhs)
    CALL vecfil_discreteBCsol (p_rvector)

  END SUBROUTINE

  ! ***************************************************************************
    SUBROUTINE b1d5_getDefect (ite,rx,rb,rd,p_rcollection)
  
    USE linearsystemblock
    USE collection
    
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
    INTEGER, INTENT(IN)                           :: ite

    ! Current iteration vector
    TYPE(t_vectorBlock), INTENT(IN),TARGET        :: rx

    ! Right hand side vector of the equation.
    TYPE(t_vectorBlock), INTENT(IN)               :: rb
  !</input>
               
  !<inputoutput>
    ! Pointer to collection structure of the application. Points to NULL()
    ! if there is none.
    TYPE(t_collection), POINTER                   :: p_rcollection

    ! Defect vector b-A(x)x. This must be filled by the callback routine
    ! with data.
    TYPE(t_vectorBlock), INTENT(INOUT)            :: rd
  !</inputoutput>

      ! local variables
      TYPE(t_bilinearForm) :: rform
      INTEGER :: ilvmax
      TYPE(t_matrixBlock), POINTER :: p_rmatrix

      ! Get maximum level from the collection
      ilvmax = collct_getvalue_int (p_rcollection,'NLMAX')

      ! Get the system matrix on the maximum level
      p_rmatrix => collct_getvalue_mat (p_rcollection,'SYSTEMMAT',ilvmax)
      
      ! Put a reference to rx into the collection. This way, we inform the callback
      ! routine of the matrix assembly about the solution vector to use
      ! fot the nonlinear term.
      CALL collct_setvalue_vec(p_rcollection,'RX',rx,.TRUE.)
      
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
      rform%ballCoeffConstant = .FALSE.
      rform%BconstantCoeff = .FALSE.

      ! Now we can build the matrix entries.
      ! We specify the callback function coeff_burgers for the coefficients.
      ! As long as we use constant coefficients, this routine is not used.
      ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
      ! the framework will call the callback routine to get analytical data.
      !
      ! We pass our collection structure as well to this routine, 
      ! so the callback routine has access to everything what is
      ! in the collection.
      CALL bilf_buildMatrixScalar (rform,.TRUE.,&
                                   p_rmatrix%RmatrixBlock(1,1),coeff_burgers,&
                                   p_rcollection)

      ! Remove the vector from the collection - not necessary anymore.
      CALL collct_deletevalue (p_rcollection,'RX')
      
      ! Implement discrete boundary conditions into the matrix. 
      ! Call the appropriate matrix filter to modify the system matrix
      ! according to the attached discrete boundary conditions.
      CALL matfil_discreteBC (p_rmatrix)
      
      ! Build the defect: d=b-Ax
      CALL lsysbl_copyVector (rb,rd)
      CALL lsysbl_blockMatVec (p_rmatrix, rx, rd, -1.0_DP, 1.0_DP)
    
      ! Apply the defect-vector filter for discrete boundary conditions
      ! to modify the defect vector according to the (discrete) boundary
      ! conditions.
      CALL vecfil_discreteBCdef (rd)
      
      ! That's it
  
    END SUBROUTINE
    
  ! ***************************************************************************

    SUBROUTINE b1d5_precondDefect (ite,rd,rx,rb,domega,bsuccess,p_rcollection)
  
    USE linearsystemblock
    USE collection
    
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
    INTEGER, INTENT(IN)                           :: ite

    ! Defect vector b-A(x)x. This must be replaced by J^{-1} rd by a preconditioner.
    TYPE(t_vectorBlock), INTENT(INOUT)            :: rd

    ! Pointer to collection structure of the application. Points to NULL()
    ! if there is none.
    TYPE(t_collection), POINTER                   :: p_rcollection
    
    ! Damping parameter. Is set to rsolverNode%domega (usually = 1.0_DP)
    ! on the first call to the callback routine.
    ! The callback routine can modify this parameter according to any suitable
    ! algorithm to calculate an 'optimal damping' parameter. The nonlinear loop
    ! will then use this for adding rd to the solution vector:
    ! $$ x_{n+1} = x_n + domega*rd $$
    ! domega will stay at this value until it's changed again.
    REAL(DP), INTENT(INOUT)                       :: domega

    ! If the preconditioning was a success. Is normally automatically set to
    ! TRUE. If there is an error in the preconditioner, this flag can be
    ! set to FALSE. In this case, the nonlinear solver breaks down with
    ! the error flag set to 'preconditioner broke down'.
    LOGICAL, INTENT(INOUT)                        :: bsuccess
  !</inputoutput>
  
  !<input>
    ! Current iteration vector
    TYPE(t_vectorBlock), INTENT(IN)               :: rx

    ! Current right hand side of the nonlinear system
    TYPE(t_vectorBlock), INTENT(IN), TARGET       :: rb
  !</input>
  
    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    INTEGER :: ierror
    TYPE(t_linsolNode), POINTER :: p_rsolverNode
  
      ! Our 'parent' (the caller of the nonlinear solver) has prepared
      ! a preconditioner node for us (a linear solver with symbolically
      ! factorised matrices). Get this from the collection.
      
      p_rsolverNode => collct_getvalue_linsol(p_rcollection,'LINSOLVER')

      ! Initialise data of the solver. This in fact performs a numeric
      ! factorisation of the matrices in UMFPACK-like solvers.
      CALL linsol_initData (p_rsolverNode, ierror)
      IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
      
      ! Finally solve the system. As we want to solve Ax=b with
      ! b being the real RHS and x being the real solution vector,
      ! we use linsol_solveAdaptively. If b is a defect
      ! RHS and x a defect update to be added to a solution vector,
      ! we would have to use linsol_precondDefect instead.
      CALL linsol_precondDefect (p_rsolverNode,rd)

      ! Release the numeric factorisation of the matrix.
      ! We don't release the symbolic factorisation, as we can use them
      ! for the next iteration.
      CALL linsol_doneData (p_rsolverNode)

    END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE b1d5_solve (rproblem)
  
!<description>
  ! Solves the given problem by applying a nonlinear solver.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    TYPE(t_vectorBlock), TARGET :: rtempBlock
    
    TYPE(t_nlsolNode) :: rnlSol
    INTEGER :: ierror
    TYPE(t_linsolNode), POINTER :: p_rsolverNode
    TYPE(t_matrixBlock), DIMENSION(1) :: Rmatrices
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    
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
    CALL linsol_initUMFPACK4 (p_rsolverNode)

    ! Get the system matrix on the finest level...
    p_rmatrix => rproblem%rlevelInfo%rmatrix

    ! And associate it to the solver
    Rmatrices = (/p_rmatrix/)
    CALL linsol_setMatrices(p_rsolverNode,Rmatrices)

    ! Initialise structure of the solver. This allows the
    ! solver to allocate memory / perform some precalculation
    ! to the problem.
    ! In fact, solvers like UMFPACK use this for a symbolic factorisation
    ! of the matrix.
    CALL linsol_initStructure (p_rsolverNode, ierror)
    IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
  
    ! Put the prepared solver node to the collection for later use.
    CALL collct_setvalue_linsol(rproblem%rcollection,'LINSOLVER',p_rsolverNode,.TRUE.)
    
    ! Create a temporary vector we need for the nonlinera iteration.
    CALL lsysbl_createVecBlockIndirect (rproblem%rrhs, rtempBlock, .FALSE.)

    ! The nonlinear solver structure rnlSol is initialised by the default
    ! initialisation with all necessary information to solve the problem.
    ! We call the nonlinear solver directly. For preconditioning
    ! and defect calculation, we use our own callback routine.
    !rnlSol%domega = 0.25_DP
    rnlSol%ioutputLevel = 2
    CALL nlsol_performSolve(rnlSol,rproblem%rvector,rproblem%rrhs,rtempBlock,&
                            b1d5_getDefect,b1d5_precondDefect,&
                            rcollection=rproblem%rcollection)

    ! Release the temporary vector
    CALL lsysbl_releaseVector (rtempBlock)
    
    ! Remove the solver node from the collection - not needed anymore there
    CALL collct_deletevalue(rproblem%rcollection,'LINSOLVER')
    
    ! Clean up the linear solver, release all memory, remove the solver node
    ! from memory.
    CALL linsol_releaseSolver (p_rsolverNode)
    
    CALL output_lbrk()
    CALL output_line ('Nonlinear solver statistics')
    CALL output_line ('---------------------------')
    CALL output_line ('Intial defect: '//TRIM(sys_sdEL(rnlSol%DinitialDefect(1),15)))
    CALL output_line ('Final defect:  '//TRIM(sys_sdEL(rnlSol%DfinalDefect(1),15)))
    CALL output_line ('#Iterations:   '//TRIM(sys_siL(rnlSol%iiterations,10)))

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE b1d5_postprocessing (rproblem)
  
!<description>
  ! Writes the solution into a GMV file.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  
    ! We need some more variables for postprocessing
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    
    ! Output block for UCD output to GMV file
    TYPE(t_ucdExport) :: rexport

    ! A pointer to the solution vector and to the triangulation.
    TYPE(t_vectorBlock), POINTER :: p_rvector
    TYPE(t_triangulation), POINTER :: p_rtriangulation

    ! Get the solution vector from the problem structure.
    p_rvector => rproblem%rvector
    
    ! From the attached discretisation, get the underlying triangulation
    p_rtriangulation => &
      p_rvector%RvectorBlock(1)%p_rspatialDiscretisation%p_rtriangulation
    
    ! p_rvector now contains our solution. We can now
    ! start the postprocessing. 
    ! Start UCD export to GMV file:
    CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,'gmv/u5.gmv')
    
    CALL lsyssc_getbase_double (p_rvector%RvectorBlock(1),p_Ddata)
    CALL ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that's it.
    CALL ucd_write (rexport)
    CALL ucd_release (rexport)
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE b1d5_doneMatVec (rproblem)
  
!<description>
  ! Releases system matrix and vectors.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

    INTEGER :: ihandle,ilvmax

    ! Release matrix and vectors on all levels
    ilvmax=rproblem%ilvmax
      
    ! Delete the matrix
    CALL lsysbl_releaseMatrix (rproblem%rlevelInfo%rmatrix)

    ! Delete the variables from the collection.
    CALL collct_deletevalue (rproblem%rcollection,'SYSTEMMAT',ilvmax)
    
    ! Release the permutation for sorting matrix/vectors
    ihandle = collct_getvalue_int (rproblem%rcollection,'LAPLACE-CM',ilvmax)
    CALL storage_free (ihandle)
    CALL collct_deletevalue (rproblem%rcollection,'LAPLACE-CM',ilvmax)

    ! Delete solution/RHS vector
    CALL lsysbl_releaseVector (rproblem%rvector)
    CALL lsysbl_releaseVector (rproblem%rrhs)

    ! Delete the variables from the collection.
    CALL collct_deletevalue (rproblem%rcollection,'RHS')
    CALL collct_deletevalue (rproblem%rcollection,'SOLUTION')

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE b1d5_doneBC (rproblem)
  
!<description>
  ! Releases discrete boundary conditions from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: ilvmax

    ilvmax=rproblem%ilvmax
      
    ! Release our discrete version of the boundary conditions
    CALL bcasm_releaseDiscreteBC (rproblem%rlevelInfo%rdiscreteBC)
    
  END SUBROUTINE


  ! ***************************************************************************

!<subroutine>

  SUBROUTINE b1d5_doneDiscretisation (rproblem)
  
!<description>
  ! Releases the discretisation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: ilvmax

    ilvmax=rproblem%ilvmax
      
    ! Delete the block discretisation together with the associated
    ! scalar spatial discretisations....
    CALL spdiscr_releaseBlockDiscr(rproblem%rlevelInfo%p_rdiscretisation)

    ! and remove the allocated block discretisation structure from the heap.
    DEALLOCATE(rproblem%rlevelInfo%p_rdiscretisation)
    
  END SUBROUTINE
    
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE b1d5_doneParamTriang (rproblem)
  
!<description>
  ! Releases the triangulation and parametrisation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

    ! Release the triangulation
    CALL tria_done (rproblem%rlevelInfo%rtriangulation)
    
    ! Finally release the domain.
    CALL boundary_release (rproblem%p_rboundary)
    
    CALL collct_deleteValue(rproblem%rcollection,'NLMAX')

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE burgers1d5
  
!<description>
  ! This is a 'separated' burgers1d solver for solving a Burgers-1D
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
    INTEGER :: NLMAX
    
    ! A problem structure for our problem
    TYPE(t_problem), POINTER :: p_rproblem
    
    INTEGER :: i
    
    ! Ok, let's start. 
    !
    ! We want to solve our Laplace problem on level...
    NLMAX = 7
    
    ! Allocate the problem structure -- it's rather large
    ALLOCATE(p_rproblem)
    
    ! Initialise the collection
    CALL collct_init (p_rproblem%rcollection)
    DO i=1,NLMAX
      CALL collct_addlevel_all (p_rproblem%rcollection)
    END DO

    ! So now the different steps - one after the other.
    
    ! Initialisation.
    CALL b1d5_initParamTriang (NLMAX,p_rproblem)
    CALL b1d5_initDiscretisation (p_rproblem)    
    CALL b1d5_initMatVec (p_rproblem)    
    CALL b1d5_initDiscreteBC (p_rproblem)
    
    ! Implementation of boundary conditions
    CALL b1d5_implementBC (p_rproblem)
    
    ! Solve the problem
    CALL b1d5_solve (p_rproblem)
    
    ! Postprocessing
    CALL b1d5_postprocessing (p_rproblem)
    
    ! Cleanup
    CALL b1d5_doneMatVec (p_rproblem)
    CALL b1d5_doneBC (p_rproblem)
    CALL b1d5_doneDiscretisation (p_rproblem)
    CALL b1d5_doneParamTriang (p_rproblem)
    
    ! Print some statistical data about the collection - anything forgotten?
    PRINT *
    PRINT *,'Remaining collection statistics:'
    PRINT *,'--------------------------------'
    PRINT *
    CALL collct_printStatistics (p_rproblem%rcollection)
    
    ! Finally release the collection and the problem structure.
    CALL collct_done (p_rproblem%rcollection)
    
    DEALLOCATE(p_rproblem)
    
  END SUBROUTINE

END MODULE
