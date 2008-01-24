!##############################################################################
!# ****************************************************************************
!# <name> poisson2d_method2 </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple Poisson
!# problem with constant coefficients on a simple domain.
!#
!# The routine splits up the tasks of reading the domain, creating 
!# triangulations, discretisation, solving, postprocessing and cleanup into
!# different subroutines. The communication between these subroutines
!# is done using an application-specific structure saving problem data
!# as well as a collection structure for the communication with callback
!# routines.
!# </purpose>
!##############################################################################

MODULE poisson2d_method2

  USE fsystem
  USE genoutput
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
  USE ucd
  USE pprocerror
  
  USE collection
    
  USE poisson2d_callback
  
  IMPLICIT NONE
  
!<types>

!<typeblock description="Type block defining all information about one level">

  TYPE t_problem_lvl
  
    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation (trial/test functions,...)
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
    
    ! A matrix and vector. The vector accepts the RHS of the problem
    ! in scalar form. The matrix will receive the discrete Laplace operator.
    TYPE(t_matrixBlock) :: rmatrix
    TYPE(t_vectorBlock) :: rvector,rrhs

    ! A variable describing the discrete boundary conditions.    
    TYPE(t_discreteBC), POINTER :: p_rdiscreteBC
  
  END TYPE
  
!</typeblock>


!<typeblock description="Application-specific type block for poisson problem">

  TYPE t_problem
  
    ! NLMAX receives the level where we want to solve
    INTEGER :: NLMAX

    ! An object for saving the domain:
    TYPE(t_boundary), POINTER :: p_rboundary

    ! A variable describing the analytic boundary conditions.    
    TYPE(t_boundaryConditions), POINTER :: p_rboundaryConditions

    ! A solver node that accepts parameters for the linear solver    
    TYPE(t_linsolNode), POINTER :: p_rsolverNode

    ! An array of t_problem_lvl structures, each corresponding
    ! to one level of the discretisation. There is currently
    ! only one level supported, identified by NLMAX!
    TYPE(t_problem_lvl), DIMENSION(1) :: RlevelInfo
    
    ! A collection object that saves structural data and some 
    ! problem-dependent information which is e.g. passed to 
    ! callback routines.
    TYPE(t_collection) :: rcollection
    
  END TYPE

!</typeblock>

!</types>
  
CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE pm3_initParamTriang (ilv,rproblem)
  
!<description>
  ! This routine initialises the parametrisation and triangulation of the
  ! domain. The corresponding .prm/.tri files are read from disc and
  ! the triangulation is refined as described by the parameter ilv.
!</description>

!<input>
  ! The level up to where we refine the coarse mesh.
  INTEGER, INTENT(IN) :: ilv
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
!</inputoutput>

  ! local variables
  
    ! Initialise the level in the problem structure
    rproblem%NLMAX = ilv

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    ! Set p_rboundary to NULL() to create a new structure.
    NULLIFY(rproblem%p_rboundary)
    CALL boundary_read_prm(rproblem%p_rboundary, './pre/QUAD.prm')
        
    ! Now read in the basic triangulation.
    CALL tria_readTriFile2D (rproblem%RlevelInfo(1)%rtriangulation, &
        './pre/QUAD.tri', rproblem%p_rboundary)
    
    ! Refine it.
    CALL tria_quickRefine2LevelOrdering (rproblem%NLMAX-1, &
        rproblem%RlevelInfo(1)%rtriangulation,rproblem%p_rboundary)
    
    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    CALL tria_initStandardMeshFromRaw (rproblem%RlevelInfo(1)%rtriangulation, &
        rproblem%p_rboundary)
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE pm3_initDiscretisation (rproblem)
  
!<description>
  ! This routine initialises the discretisation structure of the underlying
  ! problem and saves it to the problem structure.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

  ! local variables
  
    ! An object for saving the domain:
    TYPE(t_boundary), POINTER :: p_rboundary
    
    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation), POINTER :: p_rtriangulation
    
    ! An object for the spatial discretisation
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
    
    ! Ask the problem structure to give us the boundary and triangulation.
    ! We need it for the discretisation.
    p_rboundary => rproblem%p_rboundary
    p_rtriangulation => rproblem%RlevelInfo(1)%rtriangulation
    
    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we only have one block.
    ALLOCATE(p_rdiscretisation)
    CALL spdiscr_initBlockDiscr2D (p_rdiscretisation,1,&
                                   p_rtriangulation, p_rboundary)
                                   
    ! SAve the discretisation structure to our local LevelInfo structure
    ! for later use.
    rproblem%RlevelInfo(1)%p_rdiscretisation => p_rdiscretisation

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

  SUBROUTINE pm3_initMatVec (rproblem)
  
!<description>
  ! Calculates the system matrix and RHS vector of the linear system
  ! by discretising the problem with the default discretisation structure
  ! in the problem structure.
  ! Sets up a solution vector for the linear system.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

  ! local variables
  
    ! A bilinear and linear form describing the analytic problem to solve
    TYPE(t_bilinearForm) :: rform
    TYPE(t_linearForm) :: rlinform
    
    ! A pointer to the system matrix and the RHS/solution vectors.
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    TYPE(t_vectorBlock), POINTER :: p_rrhs,p_rvector

    ! A pointer to the discretisation structure with the data.
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
  
    ! Ask the problem structure to give us the discretisation structure
    p_rdiscretisation => rproblem%RlevelInfo(1)%p_rdiscretisation
    
    p_rmatrix => rproblem%RlevelInfo(1)%rmatrix
    p_rrhs    => rproblem%RlevelInfo(1)%rrhs   
    p_rvector => rproblem%RlevelInfo(1)%rvector
    
    ! Initialise the block matrix with default values based on
    ! the discretisation.
    CALL lsysbl_createMatBlockByDiscr (p_rdiscretisation,p_rmatrix)    
    
    ! Save matrix and vectors to the collection.
    ! They maybe used later, expecially in nonlinear problems.
    CALL collct_setvalue_vec(rproblem%rcollection,'RHS',p_rrhs,.TRUE.)
    CALL collct_setvalue_vec(rproblem%rcollection,'SOLUTION',p_rvector,.TRUE.)
    CALL collct_setvalue_mat(rproblem%rcollection,'LAPLACE',p_rmatrix,.TRUE.)

    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! We create that directly in the block (1,1) of the block matrix
    ! using the discretisation structure of the first block.
    CALL bilf_createMatrixStructure (&
              p_rdiscretisation%RspatialDiscretisation(1),LSYSSC_MATRIX9,&
              p_rmatrix%RmatrixBlock(1,1))
    
    ! Update the structural information of the block matrix, as we manually
    ! changed one of the submatrices:
    CALL lsysbl_updateMatStrucInfo (p_rmatrix)
    
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
    rform%ballCoeffConstant = .TRUE.
    rform%BconstantCoeff = .TRUE.
    rform%Dcoefficients(1)  = 1.0 
    rform%Dcoefficients(2)  = 1.0 

    ! Now we can build the matrix entries.
    ! We specify the callback function coeff_Laplace for the coefficients.
    ! As long as we use constant coefficients, this routine is not used.
    ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
    ! the framework will call the callback routine to get analytical data.
    !
    ! We pass our collection structure as well to this routine, 
    ! so the callback routine has access to everything what is
    ! in the collection.
    CALL bilf_buildMatrixScalar (rform,.TRUE.,&
                                 p_rmatrix%RmatrixBlock(1,1),coeff_Laplace_2D,&
                                 rproblem%rcollection)
    
    ! Now we want to build up the right hand side. At first we need a block
    ! vector of the right structure. Although we could manually create
    ! that vector, the easiest way to set up the vector structure is
    ! to create it by using our matrix as template:
    CALL lsysbl_createVecBlockIndMat (p_rmatrix,p_rrhs, .FALSE.)
    
    ! The vector structure is done but the entries are missing. 
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
    CALL linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscretisation(1),rlinform,.TRUE.,&
              p_rrhs%RvectorBlock(1),coeff_RHS_2D,&
              rproblem%rcollection)
    
    ! Now we have block vectors for the RHS and the matrix. What we
    ! need additionally is a block vector for the solution. 
    ! Create them using the RHS as template.
    ! Fill the solution vector with 0:
    CALL lsysbl_createVecBlockIndirect (p_rrhs, p_rvector, .TRUE.)
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE pm3_initAnalyticBC (rproblem)
  
!<description>
  ! This initialises the analytic bonudary conditions of the problem
  ! and saves them to the problem structure.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

  ! local variables

    ! A set of variables describing the analytic boundary conditions.    
    TYPE(t_boundaryRegion) :: rboundaryRegion
    TYPE(t_bcRegion), POINTER :: p_rbcRegion
    
    ! A pointer to the discretisation structure with the data.
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
    
    ! A pointer to the domain
    TYPE(t_boundary), POINTER :: p_rboundary
  
    ! Ask the problem structure to give us the discretisation structure and
    p_rdiscretisation => rproblem%RlevelInfo(1)%p_rdiscretisation
    
    ! Get the domain from the discretisation
    p_rboundary => p_rdiscretisation%p_rboundary
    
    ! For implementing boundary conditions, we use a 'filter technique with
    ! discretised boundary conditions'. This means, we first have to calculate
    ! a discrete version of the analytic BC, which we can implement into the
    ! solution/RHS vectors using the corresponding filter.
    !
    ! At first, we need the analytic description of the boundary conditions.
    ! Initialise a structure for boundary conditions, which accepts this,
    ! on the heap.
    !
    ! Set p_rboundaryConditions to NULL() to create a new structure on the heap.
    NULLIFY (rproblem%p_rboundaryConditions)
    CALL bcond_initBC (rproblem%p_rboundaryConditions,&
                       p_rdiscretisation%p_rboundary)
    
    ! We 'know' already (from the problem definition) that we have four boundary
    ! segments in the domain. Each of these, we want to use for inforcing
    ! some kind of boundary condition.
    !
    ! We ask the bondary routines to create a 'boundary region' - which is
    ! simply a part of the boundary corresponding to a boundary segment.
    ! A boundary region roughly contains the type, the min/max parameter value
    ! and whether the endpoints are inside the region or not.
    CALL boundary_createRegion(p_rboundary,1,1,rboundaryRegion)
    
    ! We use this boundary region and specify that we want to have Dirichlet
    ! boundary there. The following routine adds a new 'boundary condition region'
    ! for the first segment to the boundary condition structure.
    ! The region will be set up as 'Dirichlet boundary'.
    ! We specify icomponent='1' to indicate that we set up the
    ! Dirichlet BC's for the first (here: one and only) component in the solution
    ! vector.
    ! The routine also returns the created object in p_rbcRegion so that we can
    ! modify it - but accept it as it is, so we can ignore that.
    CALL bcond_newDirichletBConRealBD (rproblem%p_rboundaryConditions,1,&
                                       rboundaryRegion,p_rbcRegion)
                             
    ! Now to the edge 2 of boundary component 1 the domain. We use the
    ! same two routines to add the boundary condition to p_rboundaryConditions.
    CALL boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
    CALL bcond_newDirichletBConRealBD (rproblem%p_rboundaryConditions,1,&
                                       rboundaryRegion,p_rbcRegion)
                             
    ! Edge 3 of boundary component 1.
    CALL boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
    CALL bcond_newDirichletBConRealBD (rproblem%p_rboundaryConditions,1,&
                                       rboundaryRegion,p_rbcRegion)
    
    ! Edge 4 of boundary component 1. That's it.
    CALL boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
    CALL bcond_newDirichletBConRealBD (rproblem%p_rboundaryConditions,1,&
                                       rboundaryRegion,p_rbcRegion)
                             
    ! The boundary conditions are set up, but still the discretisation
    ! does not know about it. So inform the discretisation which
    ! analytic boundary conditions to use:
    p_rdiscretisation%p_rboundaryConditions => rproblem%p_rboundaryConditions
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE pm3_initDiscreteBC (rproblem)
  
!<description>
  ! This calculates the discrete version of the boundary conditions and
  ! assigns it to the system matrix and RHS vector.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

  ! local variables

    ! A pointer to the system matrix and the RHS vector as well as 
    ! the discretisation
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    TYPE(t_vectorBlock), POINTER :: p_rrhs,p_rvector
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
    
    ! Pointer to structure for saving discrete BC's:
    TYPE(t_discreteBC), POINTER :: p_rdiscreteBC

    ! Get our matrix and right hand side from the problem structure.
    p_rrhs    => rproblem%RlevelInfo(1)%rrhs   
    p_rvector => rproblem%RlevelInfo(1)%rvector
    p_rmatrix => rproblem%RlevelInfo(1)%rmatrix
    
    ! From the matrix or the RHS we have access to the discretisation and the
    ! analytic boundary conditions.
    p_rdiscretisation => p_rmatrix%p_rblockDiscretisation
    
    ! For the discrete problem, we need a discrete version of the above
    ! boundary conditions. So we have to discretise them.
    ! The following routine gives back p_rdiscreteBC, a pointer to a
    ! discrete version of the boundary conditions. Remark that
    ! the pointer has to be nullified before calling the routine,
    ! otherwise, the routine tries to update the boundary conditions
    ! in p_rdiscreteBC!
    ! getBoundaryValues is a callback routine that specifies the
    ! values on the boundary. We pass our collection structure as well
    ! to this routine, so the callback routine has access to everything what is
    ! in the collection.
    NULLIFY(rproblem%RlevelInfo(1)%p_rdiscreteBC)
    CALL bcasm_discretiseBC (p_rdiscretisation,rproblem%RlevelInfo(1)%p_rdiscreteBC, &
                             .FALSE.,getBoundaryValues_2D,rproblem%rcollection)
                             
    ! Hang the pointer into the vectors and the matrix. That way, these
    ! boundary conditions are always connected to that matrix and that
    ! vector.
    p_rdiscreteBC => rproblem%RlevelInfo(1)%p_rdiscreteBC
    
    p_rmatrix%p_rdiscreteBC => p_rdiscreteBC
    p_rrhs%p_rdiscreteBC => p_rdiscreteBC
    p_rvector%p_rdiscreteBC => p_rdiscreteBC
                
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE pm3_implementBC (rproblem)
  
!<description>
  ! Implements boundary conditions into the RHS and into a given solution vector.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

  ! local variables
  
    ! A pointer to the system matrix and the RHS vector as well as 
    ! the discretisation
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    TYPE(t_vectorBlock), POINTER :: p_rrhs,p_rvector

    ! Get our matrix and right hand side from the problem structure.
    p_rrhs    => rproblem%RlevelInfo(1)%rrhs   
    p_rvector => rproblem%RlevelInfo(1)%rvector
    p_rmatrix => rproblem%RlevelInfo(1)%rmatrix
    
    ! Next step is to implement boundary conditions into the RHS,
    ! solution and matrix. This is done using a vector/matrix filter
    ! for discrete boundary conditions.
    ! The discrete boundary conditions are already attached to the
    ! vectors/matrix. Call the appropriate vector/matrix filter that
    ! modifies the vectors/matrix according to the boundary conditions.
    CALL vecfil_discreteBCrhs (p_rrhs)
    CALL vecfil_discreteBCsol (p_rvector)
    CALL matfil_discreteBC (p_rmatrix)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE pm3_solve (rproblem)
  
!<description>
  ! Solves the given problem by applying a linear solver.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

  ! local variables

    ! Error indicator during initialisation of the solver
    INTEGER :: ierror    
  
    ! A filter chain to filter the vectors and the matrix during the
    ! solution process.
    TYPE(t_filterChain), DIMENSION(1), TARGET :: RfilterChain
    TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain

    ! A pointer to the system matrix and the RHS vector as well as 
    ! the discretisation
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    TYPE(t_vectorBlock), POINTER :: p_rrhs,p_rvector
    TYPE(t_vectorBlock), TARGET :: rtempBlock

    ! A solver node that accepts parameters for the linear solver    
    TYPE(t_linsolNode), POINTER :: p_rsolverNode,p_rpreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    TYPE(t_matrixBlock), DIMENSION(1) :: Rmatrices

    ! Get our matrix and right hand side from the problem structure.
    p_rrhs    => rproblem%RlevelInfo(1)%rrhs   
    p_rvector => rproblem%RlevelInfo(1)%rvector
    p_rmatrix => rproblem%RlevelInfo(1)%rmatrix
    
    ! Create a temporary vector for the solver - it needs that.
    CALL lsysbl_createVecBlockIndirect (p_rrhs, rtempBlock, .FALSE.)
    
    ! During the linear solver, the boundary conditions must
    ! frequently be imposed to the vectors. This is done using
    ! a filter chain. As the linear solver does not work with 
    ! the actual solution vectors but with defect vectors instead,
    ! a filter for implementing the real boundary conditions 
    ! would be wrong.
    ! Therefore, create a filter chain with one filter only,
    ! which implements Dirichlet-conditions into a defect vector.
    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

    ! Create a BiCGStab-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    p_RfilterChain => RfilterChain
    NULLIFY(p_rpreconditioner)
    CALL linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,p_RfilterChain)
    
    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2

    ! Attach the system matrix to the solver.
    ! First create an array with the matrix data (on all levels, but we
    ! only have one level here), then call the initialisation 
    ! routine to attach all these matrices.
    ! Remark: Don't make a call like
    !    CALL linsol_setMatrices(p_RsolverNode,(/p_rmatrix/))
    ! This doesn't work on all compilers, since the compiler would have
    ! to create a temp array on the stack - which does not always work!
    Rmatrices = (/p_rmatrix/)
    CALL linsol_setMatrices(p_RsolverNode,Rmatrices)
    
    ! Initialise structure/data of the solver. This allows the
    ! solver to allocate memory / perform some precalculation
    ! to the problem.
    CALL linsol_initStructure (p_rsolverNode,ierror)
    IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
    CALL linsol_initData (p_rsolverNode,ierror)
    IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
    
    ! Finally solve the system. As we want to solve Ax=b with
    ! b being the real RHS and x being the real solution vector,
    ! we use linsol_solveAdaptively. If b is a defect
    ! RHS and x a defect update to be added to a solution vector,
    ! we would have to use linsol_precondDefect instead.
    CALL linsol_solveAdaptively (p_rsolverNode,p_rvector,p_rrhs,rtempBlock)
    
    ! Release solver data and structure
    CALL linsol_doneData (p_rsolverNode)
    CALL linsol_doneStructure (p_rsolverNode)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    CALL linsol_releaseSolver (p_rsolverNode)
    
    ! Release the temporary vector
    CALL lsysbl_releaseVector (rtempBlock)
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE pm3_postprocessing (rproblem)
  
!<description>
  ! Writes the solution into a GMV file.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

  ! local variables
  
    ! We need some more variables for postprocessing.
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    
    ! Output block for UCD output to GMV file
    TYPE(t_ucdExport) :: rexport

    ! A pointer to the solution vector and to the triangulation.
    TYPE(t_vectorBlock), POINTER :: p_rvector
    TYPE(t_triangulation), POINTER :: p_rtriangulation

    ! Error of FE function to reference function
    REAL(DP) :: derror

    ! Get the solution vector from the problem structure.
    p_rvector => rproblem%RlevelInfo(1)%rvector
    
    ! From the attached discretisation, get the underlying triangulation
    p_rtriangulation => &
      p_rvector%RvectorBlock(1)%p_rspatialDiscretisation%p_rtriangulation
    
    ! p_rvector now contains our solution. We can now
    ! start the postprocessing. 
    ! Start UCD export to GMV file:
    CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,&
                       'gmv/u2d_2.gmv')
    
    CALL lsyssc_getbase_double (p_rvector%RvectorBlock(1),p_Ddata)
    CALL ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that's it.
    CALL ucd_write (rexport)
    CALL ucd_release (rexport)
    
    ! Calculate the error to the reference function.
    CALL pperr_scalar (p_rvector%RvectorBlock(1),PPERR_L2ERROR,derror,&
                       getReferenceFunction_2D)
    CALL output_line ('L2-error: ' // sys_sdEL(derror,10) )

    CALL pperr_scalar (p_rvector%RvectorBlock(1),PPERR_H1ERROR,derror,&
                       getReferenceFunction_2D)
    CALL output_line ('H1-error: ' // sys_sdEL(derror,10) )
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE pm3_doneMatVec (rproblem)
  
!<description>
  ! Releases system matrix and vectors.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

    ! Release matrix and vectors
    CALL lsysbl_releaseVector (rproblem%RlevelInfo(1)%rvector)
    CALL lsysbl_releaseVector (rproblem%RlevelInfo(1)%rrhs)
    CALL lsysbl_releaseMatrix (rproblem%RlevelInfo(1)%rmatrix)

    ! Delete the variables from the collection.
    CALL collct_deletevalue (rproblem%rcollection,'RHS')
    CALL collct_deletevalue (rproblem%rcollection,'SOLUTION')
    CALL collct_deletevalue (rproblem%rcollection,'LAPLACE')

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE pm3_doneBC (rproblem)
  
!<description>
  ! Releases discrete and analytic boundary conditions from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

    ! Release our discrete version of the boundary conditions
    CALL bcasm_releaseDiscreteBC (rproblem%RlevelInfo(1)%p_rdiscreteBC)

    ! ...and also the corresponding analytic description.
    CALL bcond_doneBC (rproblem%p_rboundaryConditions)
    
  END SUBROUTINE


  ! ***************************************************************************

!<subroutine>

  SUBROUTINE pm3_doneDiscretisation (rproblem)
  
!<description>
  ! Releases the discretisation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

    ! Delete the block discretisation together with the associated
    ! scalar spatial discretisations...
    CALL spdiscr_releaseBlockDiscr(rproblem%RlevelInfo(1)%p_rdiscretisation)
    
    ! and remove the allocated block discretisation structure from the heap.
    DEALLOCATE(rproblem%RlevelInfo(1)%p_rdiscretisation)
    
  END SUBROUTINE
    
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE pm3_doneParamTriang (rproblem)
  
!<description>
  ! Releases the triangulation and parametrisation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

    ! Release the triangulation
    CALL tria_done (rproblem%RlevelInfo(1)%rtriangulation)
    
    ! Finally release the domain.
    CALL boundary_release (rproblem%p_rboundary)
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE poisson2d_2
  
!<description>
  ! This is a 'separated' poisson solver for solving a Poisson
  ! problem. The different tasks of the problem are separated into
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
  ! 7.) Write solution to GMV file
  ! 8.) Release all variables, finish
!</description>

!</subroutine>

    ! NLMAX receives the level where we want to solve.
    INTEGER :: NLMAX
    
    ! A problem structure for our problem
    TYPE(t_problem), TARGET :: rproblem
    
    ! Ok, let's start. 
    ! We want to solve our Poisson problem on level...

    NLMAX = 7
    
    ! Initialise the collection.
    CALL collct_init (rproblem%rcollection)

    ! So now the different steps - one after the other.
    !
    ! Initialisation
    CALL pm3_initParamTriang (NLMAX,rproblem)
    CALL pm3_initDiscretisation (rproblem)    
    CALL pm3_initMatVec (rproblem)
    CALL pm3_initAnalyticBC (rproblem)   
    CALL pm3_initDiscreteBC (rproblem)
    
    ! Implementation of boundary conditions
    CALL pm3_implementBC (rproblem)
    
    ! Solve the problem
    CALL pm3_solve (rproblem)
    
    ! Postprocessing
    CALL pm3_postprocessing (rproblem)
    
    ! Cleanup
    CALL pm3_doneMatVec (rproblem)
    CALL pm3_doneBC (rproblem)
    CALL pm3_doneDiscretisation (rproblem)
    CALL pm3_doneParamTriang (rproblem)

    ! Print some statistical data about the collection - anything forgotten?
    CALL output_lbrk ()
    CALL output_line ('Remaining collection statistics:')
    CALL output_line ('--------------------------------')
    CALL output_lbrk ()
    CALL collct_printStatistics (rproblem%rcollection)
    
    ! Finally release the collection.
    CALL collct_done (rproblem%rcollection)
    
  END SUBROUTINE

END MODULE
