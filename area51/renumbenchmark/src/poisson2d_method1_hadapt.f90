!##############################################################################
!# ****************************************************************************
!# <name> poisson2d_method1_hadapt </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple Poisson
!# problem with constant coefficients on a simple domain.
!#
!# The triangulation is adaptively refined by the h-adaptivity module,
!# mixing triangle and quadrilateral elements. As solver, the Gauss elimination
!# UMFPACK solver is used.
!# </purpose>
!##############################################################################

MODULE poisson2d_method1_hadapt

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
  USE hadaptaux
  USE hadaptivity
    
  USE poisson2d_callback
  
  IMPLICIT NONE

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE poisson2d_1_hadapt
  
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
  ! 7.) Write solution to GMV file
  ! 8.) Release all variables, finish
!</description>

!</subroutine>

    ! Definitions of variables.
    !
    ! We need a couple of variables for this problem. Let's see...
    !
    ! An object for saving the domain:
    TYPE(t_boundary) :: rboundary
    
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
    TYPE(t_vectorBlock) :: rvectorBlock,rrhsBlock,rtempBlock

    ! A set of variables describing the analytic and discrete boundary
    ! conditions.    
    TYPE(t_boundaryRegion) :: rboundaryRegion
    TYPE(t_discreteBC), TARGET :: rdiscreteBC

    ! A solver node that accepts parameters for the linear solver    
    TYPE(t_linsolNode), POINTER :: p_rsolverNode,p_rpreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    TYPE(t_matrixBlock), DIMENSION(1) :: Rmatrices

    ! Adaptivity structure
    TYPE(t_hadapt) :: rhadapt
    TYPE(t_vectorScalar) :: rindicator

    ! NLMIN defines the pre-refinement level of the coarse mesh.
    INTEGER :: NLMIN
    
    ! NLMAXhRefinement defines how much refinement levels are used by
    ! the adaptive refinement
    INTEGER :: NLMAXhRefinement
    
    ! Maximum number of h-adaptivity cycles
    INTEGER :: MAXhRefinementSteps

    ! Error indicator during initialisation of the solver
    INTEGER :: ierror    
    
    ! Error of FE function to reference function
    REAL(DP) :: derror
    
    ! Output block for UCD output to GMV file
    TYPE(t_ucdExport) :: rexport
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata

    ! Ok, let's start. 
    !
    ! At first we define...
    ! 1.) the minimum refinement level that is used to get the coarse mesh. 
    NLMIN = 2
    
    ! 2.) the maximum refinement level of the h-adaptive refinement routines.
    ! This is relative to NLMIN: Starting from level NLMIN, the h-adaptivity
    ! routines refine the mesh up to NLMAXhRefinement times.
    NLMAXhRefinement = 6

    ! 3.) the maximum number of refinement sweeps.
    MAXhRefinementSteps = 8
    
    ! +------------------------------------------------------------------------
    ! | BOUNDARY AND TRIANGULATION
    ! +------------------------------------------------------------------------

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    ! Set rboundary to NULL to create a new structure on the heap.
    CALL boundary_read_prm(rboundary, './pre/QUAD.prm')
        
    ! Now read in the basic triangulation.
    CALL tria_readTriFile2D (rtriangulation, './pre/QUAD.tri', rboundary)
    
    ! Refine it to get the coarse mesh.
    CALL tria_quickRefine2LevelOrdering (NLMIN-1,rtriangulation,rboundary)
    
    ! Create information about adjacencies and everything one needs from
    ! a triangulation.
    CALL tria_initStandardMeshFromRaw (rtriangulation,rboundary)

    ! +------------------------------------------------------------------------
    ! | SETUP OF THE H-ADAPTIVITY
    ! +------------------------------------------------------------------------
       
    rhadapt%iSpec               = 2**1
    rhadapt%nsubdividemax       = NLMAXhRefinement
    rhadapt%iadaptationStrategy = HADAPT_REDGREEN
    
    ! The following tolerance values are the bounds for the monitor function
    ! gethadaptMonitorFunction below. If this monitor function returns a value
    ! > drefinementTolerance to an element, that element is refined.
    ! A value < dcoarseningTolerance on the other hand results in coarsening.
    rhadapt%drefinementTolerance = 0.5
    rhadapt%dcoarseningTolerance = 0.1
    CALL hadapt_initFromTriangulation(rhadapt,rtriangulation)

    ! +------------------------------------------------------------------------
    ! | SETUP OF THE DISCRETISATION
    ! +------------------------------------------------------------------------

    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we only have one block.
    CALL spdiscr_initBlockDiscr2D (rdiscretisation,1,&
                                   rtriangulation, rboundary)

    ! Repeat the procedure until the maximum number of refinement
    ! steps has been reached. This will be checked below.
    DO 
      
      ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
      ! structures for every component of the solution vector.
      ! Initialise the first element of the list to specify the element
      ! and cubature rule for this solution component:
      CALL spdiscr_initDiscr_triquad (rdiscretisation%RspatialDiscr(1), &
          EL_P1,EL_Q1,CUB_G3_T,CUB_G2X2,rtriangulation, &
          rboundary)
      
      ! Now as the discretisation is set up, we can start to generate
      ! the structure of the system matrix which is to solve.
      ! We create a scalar matrix, based on the discretisation structure
      ! for our one and only solution component.
      CALL bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
          LSYSSC_MATRIX9,rmatrix)
      
      ! +------------------------------------------------------------------------
      ! | DISCRETISATION OF MATRICES/VECTORS
      ! +------------------------------------------------------------------------
      
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
      CALL linf_buildVectorScalar (rdiscretisation%RspatialDiscr(1),&
          rlinform,.TRUE.,rrhs,coeff_RHS_2D)
      
      ! The linear solver only works for block matrices/vectors - but above,
      ! we created scalar ones. So the next step is to make a 1x1 block
      ! system from the matrices/vectors above which the linear solver
      ! understands.
      CALL lsysbl_createMatFromScalar (rmatrix,rmatrixBlock,rdiscretisation)
      CALL lsysbl_createVecFromScalar (rrhs,rrhsBlock,rdiscretisation)
      
      ! +------------------------------------------------------------------------
      ! | DISCRETISATION AND IMPLEMENTATION OF BOUNDARY CONDITIONS
      ! +------------------------------------------------------------------------
      
      ! Now we have the raw problem. What is missing is the definition of the boudary
      ! conditions.
      ! For implementing boundary conditions, we use a 'filter technique with
      ! discretised boundary conditions'. This means, we first have to calculate
      ! a discrete version of the analytic BC, which we can implement into the
      ! solution/RHS vectors using the corresponding filter.
      !
      ! Create a t_discreteBC structure where we store all discretised boundary
      ! conditions.
      CALL bcasm_initDiscreteBC(rdiscreteBC)
      !
      ! We 'know' already (from the problem definition) that we have four boundary
      ! segments in the domain. Each of these, we want to use for enforcing
      ! some kind of boundary condition.
      !
      ! We ask the bondary routines to create a 'boundary region' - which is
      ! simply a part of the boundary corresponding to a boundary segment.
      ! A boundary region roughly contains the type, the min/max parameter value
      ! and whether the endpoints are inside the region or not.
      CALL boundary_createRegion(rboundary,1,1,rboundaryRegion)
      
      ! We use this boundary region and specify that we want to have Dirichlet
      ! boundary there. The following call does the following:
      ! - Create Dirichlet boundary conditions on the region rboundaryRegion.
      !   We specify icomponent='1' to indicate that we set up the
      !   Dirichlet BC's for the first (here: one and only) component in the 
      !   solution vector.
      ! - Discretise the boundary condition so that the BC's can be applied
      !   to matrices and vectors
      ! - Add the calculated discrete BC's to rdiscreteBC for later use.
      CALL bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues_2D)
                               
      ! Now to the edge 2 of boundary component 1 the domain. 
      CALL boundary_createRegion(rboundary,1,2,rboundaryRegion)
      CALL bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues_2D)
                               
      ! Edge 3 of boundary component 1.
      CALL boundary_createRegion(rboundary,1,3,rboundaryRegion)
      CALL bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues_2D)
      
      ! Edge 4 of boundary component 1. That's it.
      CALL boundary_createRegion(rboundary,1,4,rboundaryRegion)
      CALL bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues_2D)
      
      ! Hang the pointer into the vector and matrix. That way, these
      ! boundary conditions are always connected to that matrix and that
      ! vector.
      rmatrixBlock%p_rdiscreteBC => rdiscreteBC
      rrhsBlock%p_rdiscreteBC => rdiscreteBC
      
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
      
      ! +------------------------------------------------------------------------
      ! | INVOKE THE SOLVER
      ! +------------------------------------------------------------------------
      
      NULLIFY(p_rpreconditioner)
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
      
      ! Do we have to perform one step of h-adaptivity?
      IF (rhadapt%nRefinementSteps .GE. MAXhRefinementSteps) EXIT
      
      ! +----------------------------------------------------------------------
      ! | COMPUTE INDICATOR FOR H-ADAPTIVITY
      ! +----------------------------------------------------------------------     

      ! Perform a posteriori error estimation
      CALL lsyssc_createVector(rindicator,rtriangulation%NEL,.TRUE.)
      CALL gethadaptMonitorFunction_2D(rtriangulation,rvectorBlock%RvectorBlock(1),&
          -1,3,rindicator)
      
      ! Output error
      CALL lsyssc_getbase_double(rindicator,p_Ddata)
      CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,'gmv/error8.'//&
          TRIM(sys_siL(rhadapt%nRefinementSteps,3))//'.gmv')
      CALL ucd_addVariableElementBased (rexport,'error',UCD_VAR_STANDARD, p_Ddata)
      CALL ucd_write (rexport)
      CALL ucd_release (rexport)

      ! Perform one step h-adaptivity
      CALL hadapt_refreshAdaptation(rhadapt,rtriangulation)
      CALL hadapt_performAdaptation(rhadapt,rindicator)
      
      ! Release the indicator vector again
      CALL lsyssc_releaseVector(rindicator)
      
      ! Generate raw mesh from adaptivity structure
      CALL hadapt_generateRawMesh(rhadapt,rtriangulation)
      
      ! Create information about adjacencies and everything one needs from
      ! a triangulation.
      CALL tria_initStandardMeshFromRaw (rtriangulation,rboundary)
      
      ! +----------------------------------------------------------------------
      ! | TEMPORAL CLEANUP
      ! +----------------------------------------------------------------------

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
      CALL bcasm_releaseDiscreteBC (rdiscreteBC)
      
    END DO
    
    ! +------------------------------------------------------------------------
    ! | POSTPROCESSING
    ! +------------------------------------------------------------------------
    
    ! That's it, rvectorBlock now contains our solution. We can now
    ! start the postprocessing. 
    ! Start UCD export to GMV file:
    CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       'gmv/u2d_1_hadapt.gmv')
    
    CALL lsyssc_getbase_double (rvectorBlock%RvectorBlock(1),p_Ddata)
    CALL ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
      
    ! Write the file to disc, that's it.
    CALL ucd_write (rexport)
    CALL ucd_release (rexport)
    
    ! Calculate the error to the reference function.
    CALL pperr_scalar (rvectorBlock%RvectorBlock(1),PPERR_L2ERROR,derror,&
                       getReferenceFunction_2D)
    CALL output_line ('L2-error: ' // sys_sdEL(derror,10) )

    CALL pperr_scalar (rvectorBlock%RvectorBlock(1),PPERR_H1ERROR,derror,&
                       getReferenceFunction_2D)
    CALL output_line ('H1-error: ' // sys_sdEL(derror,10) )
    
    ! +------------------------------------------------------------------------
    ! | CLEANUP
    ! +------------------------------------------------------------------------

    ! We are finished - but not completely!
    ! Now, clean up so that all the memory is available again.
    !
    ! Release adaptivity structure
    CALL hadapt_releaseAdaptation(rhadapt)

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
    CALL bcasm_releaseDiscreteBC (rdiscreteBC)

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    CALL spdiscr_releaseBlockDiscr(rdiscretisation)
    
    ! Release the triangulation. 
    CALL tria_done (rtriangulation)
    
    ! Finally release the domain, that's it.
    CALL boundary_release (rboundary)
    
  END SUBROUTINE

END MODULE

