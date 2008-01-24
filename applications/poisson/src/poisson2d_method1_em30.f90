!##############################################################################
!# ****************************************************************************
!# <name> poisson_method7 </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple Poisson
!# problem with constant coefficients on a simple domain.
!#
!# The module is the same as poisson_method1 but with different element
!# implementation that is, EL_EM30. The different focuses on writing the
!# solution vector into the way that GMV can read it.
!# </purpose>
!##############################################################################

MODULE poisson2d_method1_em30

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
  USE matrixio

  USE poisson2d_callback
  USE spdiscprojection

  IMPLICIT NONE

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE poisson2d_1_em30

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
  ! 7.) Convert the solution vector into GMV's style
  ! 8.) Write solution to GMV file
  ! 9.) Release all variables, finish
!</description>

!</subroutine>

    ! Definitions of variables.
    !
    ! We need a couple of variables for this problem. Let's see...
    !
    ! An object for saving the domain:
    TYPE(t_boundary), POINTER :: p_rboundary

    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    TYPE(t_blockDiscretisation) :: rdiscretisation

    ! A bilinear and linear form describing the analytic problem to solve
    TYPE(t_bilinearForm) :: rform
    TYPE(t_linearForm) :: rlinform

    ! A scalar matrix and vector. The vector accepts the RHS of the problem
    ! in scalar form.
    TYPE(t_matrixScalar) :: rmatrix
    TYPE(t_vectorScalar) :: rrhs

    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    TYPE(t_matrixBlock) :: rmatrixBlock
    TYPE(t_vectorBlock)  :: rvectorBlock,rrhsBlock,rtempBlock

    ! A set of variables describing the analytic and discrete boundary
    ! conditions.
    TYPE(t_boundaryConditions), POINTER :: p_rboundaryConditions
    TYPE(t_boundaryRegion) :: rboundaryRegion
    TYPE(t_bcRegion), POINTER :: p_rbcRegion
    TYPE(t_discreteBC), POINTER :: p_rdiscreteBC

    ! A solver node that accepts parameters for the linear solver
    TYPE(t_linsolNode), POINTER :: p_rsolverNode,p_rpreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    TYPE(t_matrixBlock), DIMENSION(1) :: Rmatrices

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    TYPE(t_filterChain), DIMENSION(1), TARGET :: RfilterChain
    TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain

    ! NLMAX receives the level where we want to solve.
    INTEGER :: NLMAX

    ! Error indicator during initialisation of the solver
    INTEGER :: ierror

    ! We need some more variables for pre/postprocessing.
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    
    ! Output block for UCD output to GMV file
    TYPE(t_ucdExport) :: rexport

    ! Error of FE function to reference function
    REAL(DP) :: derror

    ! Declarations for projecting a vector to the Q1 space for GMV export
    TYPE(t_vectorBlock) :: rprjVector
    TYPE(t_blockDiscretisation) :: rprjDiscretisation
    TYPE(t_discreteBC), POINTER ::p_rdiscreteBC_Q1

    ! Ok, let's start.
    !
    ! We want to solve our Poisson problem on level...
    NLMAX = 7

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    ! Set p_rboundary to NULL to create a new structure on the heap.
    NULLIFY(p_rboundary)
    CALL boundary_read_prm(p_rboundary, './pre/QUAD.prm')

    ! Now read in the basic triangulation.
    CALL tria_readTriFile2D (rtriangulation, './pre/QUAD.tri', p_rboundary)
    
    ! Refine it.
    CALL tria_quickRefine2LevelOrdering (NLMAX-1,rtriangulation,p_rboundary)
    
    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    CALL tria_initStandardMeshFromRaw (rtriangulation,p_rboundary)
    
    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we only have one block.
    CALL spdiscr_initBlockDiscr (rdiscretisation,1,&
                                   rtriangulation, p_rboundary)

    ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
    ! structures for every component of the solution vector.
    ! Initialise the first element of the list to specify the element
    ! and cubature rule for this solution component:
    CALL spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscretisation(1), &
                                   EL_EM30,CUB_G2X2,rtriangulation, p_rboundary)

    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! We create a scalar matrix, based on the discretisation structure
    ! for our one and only solution component.
    CALL bilf_createMatrixStructure (rdiscretisation%RspatialDiscretisation(1),&
                                     LSYSSC_MATRIX9,rmatrix)

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
    ! the framework will call the callback routine to get analytical
    ! data.
    CALL bilf_buildMatrixScalar (rform,.TRUE.,rmatrix,coeff_Laplace_2D)

    ! The same has to be done for the right hand side of the problem.
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC

    ! ... and then discretise the RHS to get a discrete version of it.
    ! Again we simply create a scalar vector based on the one and only
    ! discretisation structure.
    ! This scalar vector will later be used as the one and only first
    ! component in a block vector.
    CALL linf_buildVectorScalar (rdiscretisation%RspatialDiscretisation(1),&
                                 rlinform,.TRUE.,rrhs,coeff_RHS_2D)

    ! The linear solver only works for block matrices/vectors - but above,
    ! we created scalar ones. So the next step is to make a 1x1 block
    ! system from the matrices/vectors above which the linear solver
    ! understands.
    CALL lsysbl_createMatFromScalar (rmatrix,rmatrixBlock,rdiscretisation)
    CALL lsysbl_createVecFromScalar (rrhs,rrhsBlock,rdiscretisation)

    ! Now we have the raw problem. What is missing is the definition of the boudary
    ! conditions.
    ! For implementing boundary conditions, we use a 'filter technique with
    ! discretised boundary conditions'. This means, we first have to calculate
    ! a discrete version of the analytic BC, which we can implement into the
    ! solution/RHS vectors using the corresponding filter.
    !
    ! At first, we need the analytic description of the boundary conditions.
    ! Initialise a structure for boundary conditions which accepts this:
    !
    ! Set p_rboundaryConditions to create a new structure on the heap.
    NULLIFY (p_rboundaryConditions)
    CALL bcond_initBC (p_rboundaryConditions,p_rboundary)

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
    CALL bcond_newDirichletBConRealBD (p_rboundaryConditions,1,&
                                       rboundaryRegion,p_rbcRegion)

    ! Now to the edge 2 of boundary component 1 the domain. We use the
    ! same two routines to add the boundary condition to p_rboundaryConditions.
    CALL boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
    CALL bcond_newDirichletBConRealBD (p_rboundaryConditions,1,&
                                       rboundaryRegion,p_rbcRegion)

    ! Edge 3 of boundary component 1.
    CALL boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
    CALL bcond_newDirichletBConRealBD (p_rboundaryConditions,1,&
                                       rboundaryRegion,p_rbcRegion)

    ! Edge 4 of boundary component 1. That's it.
    CALL boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
    CALL bcond_newDirichletBConRealBD (p_rboundaryConditions,1,&
                                       rboundaryRegion,p_rbcRegion)

    ! The boundary conditions are set up, but still the block discretisation
    ! does not know about it. So inform the discretisation which
    ! analytic boundary conditions to use:
    rdiscretisation%p_rboundaryConditions => p_rboundaryConditions

    ! For the discrete problem, we need a discrete version of the above
    ! boundary conditions. So we have to discretise them.
    ! The following routine gives back p_rdiscreteBC, a pointer to a
    ! discrete version of the boundary conditions. Remark that
    ! the pointer has to be nullified before calling the routine,
    ! otherwise, the routine tries to update the boundary conditions
    ! in p_rdiscreteBC!
    NULLIFY(p_rdiscreteBC)
    CALL bcasm_discretiseBC (rdiscretisation,p_rdiscreteBC,.FALSE., &
                             getBoundaryValues_2D)

    ! Hang the pointer into the vector and matrix. That way, these
    ! boundary conditions are always connected to that matrix and that
    ! vector.
    rmatrixBlock%p_rdiscreteBC => p_rdiscreteBC
    rrhsBlock%p_rdiscreteBC => p_rdiscreteBC

    ! Now we have block vectors for the RHS and the matrix. What we
    ! need additionally is a block vector for the solution and
    ! temporary data. Create them using the RHS as template.
    ! Fill the solution vector with 0:
    CALL lsysbl_createVecBlockIndirect (rrhsBlock, rvectorBlock, .TRUE.)
    CALL lsysbl_createVecBlockIndirect (rrhsBlock, rtempBlock, .FALSE.)

    ! Next step is to implement boundary conditions into the RHS,
    ! solution and matrix. This is done using a vector/matrix filter
    ! for discrete boundary conditions.
    ! The discrete boundary conditions are already attached to the
    ! vectors/matrix. Call the appropriate vector/matrix filter that
    ! modifies the vectors/matrix according to the boundary conditions.
    CALL vecfil_discreteBCrhs (rrhsBlock)
    CALL vecfil_discreteBCsol (rvectorBlock)
    CALL matfil_discreteBC (rmatrixBlock)

    ! During the linear solver, the boundary conditions are also
    ! frequently imposed to the vectors. But as the linear solver
    ! does not work with the actual solution vectors but with
    ! defect vectors instead.
    ! So, set up a filter chain that filters the defect vector
    ! during the solution process to implement discrete boundary conditions.
    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

    ! Create a BiCGStab-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    p_RfilterChain => RfilterChain
    NULLIFY(p_rpreconditioner)
    !CALL linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,p_RfilterChain)
    CALL linsol_initUMFPACK4 (p_rsolverNode)

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
    Rmatrices = (/rmatrixBlock/)
    CALL linsol_setMatrices(p_RsolverNode,Rmatrices)

    ! Initialise structure/data of the solver. This allows the
    ! solver to allocate memory / perform some precalculation
    ! to the problem.
    CALL linsol_initStructure (p_rsolverNode, ierror)
    IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
    CALL linsol_initData (p_rsolverNode, ierror)
    IF (ierror .NE. LINSOL_ERR_NOERROR) STOP

    ! Finally solve the system. As we want to solve Ax=b with
    ! b being the real RHS and x being the real solution vector,
    ! we use linsol_solveAdaptively. If b is a defect
    ! RHS and x a defect update to be added to a solution vector,
    ! we would have to use linsol_precondDefect instead.
    CALL linsol_solveAdaptively (p_rsolverNode,rvectorBlock,rrhsBlock,rtempBlock)

    ! -------------------------------------------------------------------------
    ! Solution vector projection for projecting a solution vector of
    ! an arbitrary element to Q1 such that it can be written into a GMV file
    ! -------------------------------------------------------------------------

    ! Now, Our vector block is off the way GMV's style. We need to convert the vector
    ! to Q1 as that is the format, GMV's understands. So the task is to 
    ! create a Q1 solution from rvectorBlock.
    !
    ! Step 1: Create a discretisation structure for Q1, based on our
    ! previous discretisation structure:
    CALL spdiscr_duplicateBlockDiscr(rvectorBlock%p_rblockDiscretisation,rprjDiscretisation)
    CALL spdiscr_deriveSimpleDiscrSc (&
                 rvectorBlock%p_rblockDiscretisation%RspatialDiscretisation(1), &
                 EL_Q1, CUB_G2X2, rprjDiscretisation%RspatialDiscretisation(1))

    ! Step 2: Setup a new solution vector based on this discretisation,
    ! allocate memory.
    CALL lsysbl_createVecBlockByDiscr (rprjDiscretisation,rprjVector,.FALSE.)
 
    ! Step 3: Take the original solution vector and convert it according to the
    ! new discretisation:
    CALL spdp_projectSolution(rvectorBlock,rprjVector)

    ! Step 4: Discretise the boundary condition according to the Q1
    ! discretisation:
    NULLIFY(p_rdiscreteBC_Q1)
    CALL bcasm_discretiseBC (rprjDiscretisation,p_rdiscreteBC_Q1, &
                            .FALSE.,getBoundaryValues_2D)
                            
    ! Connect the vector to the BC's
    rprjVector%p_rdiscreteBC => p_rdiscreteBC_Q1

    ! Step 5: Set up a boundary condition filter for Dirichtley boundary conditions
    ! and pass the vetor through it. This finally implement the Dirichtley 
    ! boundary conditions into the output vector
    RfilterChain(1)%ifilterType = FILTER_DISCBCSOLREAL
    CALL filter_applyFilterChainVec (rprjVector, RfilterChain)

    ! Now we have a Q1 solution in rprjVector.
    !
    ! Step 6: Write the GMV file.
    ! That's it, rvectorBlock now contains our solution in GMV's style. 
    ! We can now start the postprocessing. 
      
    ! Call the GMV library to write out a GMV file for our solution.
    ! Start UCD export to GMV file:
    CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,'gmv/u7.gmv')
    
    CALL lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
    CALL ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that's it.
    CALL ucd_write (rexport)
    CALL ucd_release (rexport)

    ! Step 7: Release our discrete version of the boundary conditions of Q1.
    CALL bcasm_releaseDiscreteBC(p_rdiscreteBC_Q1)
    
    ! Release the block matrix/vectors rprjVector
    CALL lsysbl_releaseVector(rprjVector)
    
    ! Release the discretisation structure.
    CALL spdiscr_releaseBlockDiscr (rprjDiscretisation)
    
    ! -------------------------------------------------------------------------
    ! Projection and GMV export finished.
    ! -------------------------------------------------------------------------

    ! Calculate the error to the reference function.
    CALL pperr_scalar (rvectorBlock%RvectorBlock(1),PPERR_L2ERROR,derror,&
                       getReferenceFunction_2D)
    CALL output_line ('L2-error: ' // sys_sdEL(derror,10) )

    CALL pperr_scalar (rvectorBlock%RvectorBlock(1),PPERR_H1ERROR,derror,&
                       getReferenceFunction_2D)
    CALL output_line ('H1-error: ' // sys_sdEL(derror,10) )
    
    ! We are finished - but not completely!
    ! Now, clean up so that all the memory is available again.
    !
    ! Release solver data and structure
    CALL linsol_doneData (p_rsolverNode)
    CALL linsol_doneStructure (p_rsolverNode)

    ! Release the solver node and all subnodes attached to it (if at all):
    CALL linsol_releaseSolver (p_rsolverNode)

    ! Release the block matrix/vectors
    CALL lsysbl_releaseVector (rtempBlock)
    CALL lsysbl_releaseVector (rvectorBlock)
    CALL lsysbl_releaseVector (rrhsBlock)
    CALL lsysbl_releaseMatrix (rmatrixBlock)

    ! Release the scalar matrix/rhs vector which were used to create
    ! the block matrices/vectors. These must exist as long as the
    ! block matrices/vectors exist, as the block matrices/vectors are
    ! only 'copies' of the scalar ones, sharing the same handles!
    CALL lsyssc_releaseVector (rrhs)
    CALL lsyssc_releaseMatrix (rmatrix)

    ! Release our discrete version of the boundary conditions
    CALL bcasm_releaseDiscreteBC (p_rdiscreteBC)

    ! ...and also the corresponding analytic description.
    CALL bcond_doneBC (p_rboundaryConditions)

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    CALL spdiscr_releaseBlockDiscr(rdiscretisation)

    ! Release the triangulation. 
    CALL tria_done (rtriangulation)

    ! Finally release the domain, that's it.
    CALL boundary_release (p_rboundary)

  END SUBROUTINE

END MODULE
