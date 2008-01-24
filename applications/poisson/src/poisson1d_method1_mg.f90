!##############################################################################
!# ****************************************************************************
!# <name> poisson1d_method1_mg </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple 1D Poisson
!# problem with constant coefficients on a simple domain.
!#
!# This module is based on poisson1d_method0_simple, but using a multi-grid 
!# solver.
!# </purpose>
!##############################################################################

MODULE poisson1d_method1_mg

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
  USE genoutput
  USE matrixio
  USE vectorio
    
  USE poisson1d_callback
  
  IMPLICIT NONE
  
!<types>

!<typeblock description="Type block defining all information about one level">

  TYPE t_level
  
    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation (structure of the
    ! solution, trial/test functions,...)
    TYPE(t_blockDiscretisation) :: rdiscretisation
    
    ! A system matrix for that specific level. The matrix will receive the 
    ! discrete Laplace operator.
    TYPE(t_matrixBlock) :: rmatrix

    ! A variable describing the discrete boundary conditions.    
    TYPE(t_discreteBC), POINTER :: p_rdiscreteBC
  
  END TYPE
  
!</typeblock>

!</types>

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE poisson1d_1_mg
  
!<description>
  ! This is an all-in-one poisson solver for directly solving a Poisson
  ! problem without making use of special features like collections
  ! and so on. The routine performs the following tasks:
  !
  ! 1.) Create triangulation
  ! 2.) Set up RHS
  ! 3.) Set up matrix
  ! 4.) Create solver structure
  ! 5.) Solve the problem
  ! 6.) Write solution to GMV file
  ! 7.) Release all variables, finish
!</description>

!</subroutine>

    ! Definitions of variables.
    !
    ! We need a couple of variables for this problem. Let's see...
    !
    ! An array of problem levels for the multigrid solver
    TYPE(t_level), DIMENSION(:), ALLOCATABLE :: Rlevels
    
    ! A bilinear and linear form describing the analytic problem to solve
    TYPE(t_bilinearForm) :: rform
    TYPE(t_linearForm) :: rlinform
    
    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    TYPE(t_vectorBlock) :: rvectorBlock,rrhsBlock,rtempBlock

    ! A solver node that accepts parameters for the linear solver    
    TYPE(t_linsolNode), POINTER :: p_rsolverNode,p_rcoarseGridSolver,p_rsmoother

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    TYPE(t_matrixBlock), DIMENSION(:),ALLOCATABLE :: Rmatrices

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    TYPE(t_filterChain), DIMENSION(1), TARGET :: RfilterChain
    TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
    
    ! An interlevel projection structure for changing levels
    TYPE(t_interlevelProjectionBlock) :: rprojection

    ! One level of multigrid
    TYPE(t_linsolMGLevelInfo), POINTER :: p_rlevelInfo

    ! Error indicator during initialisation of the solver
    INTEGER :: ierror    
    
    ! Error of FE function to reference function
    REAL(DP) :: derror

    ! Output block for UCD output to GMV file
    TYPE(t_ucdExport) :: rexport
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata

    ! Some temporary variables
    INTEGER :: i
    
    ! The number of sub-intervals for the discretisation
    ! This defines the coarse mesh on level 1!
    INTEGER :: nintervals = 16
    
    ! The number of levels. The coarse mesh on level 1 with nintervals intervals is
    ! refined nlevel-1 times to get the fine mesh
    INTEGER :: nlevels = 4
    
    ! Ok, let's start. 
    !
    ! Allocate memory for all levels
    ALLOCATE(Rlevels(nlevels))
    
    ! At first, create the basic coarse grid triangulation.
    ! Our domain is [0, 1], divided into nintervals sub-intervals.
    CALL tria_createRawTria1D(Rlevels(1)%rtriangulation, 0.0_DP, 1.0_DP, nintervals)
    
    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    CALL tria_initStandardMeshFromRaw (Rlevels(1)%rtriangulation)

    ! Now refine the grid for the fine levels.
    DO i = 2, nlevels
    
      ! Refine the grid using the 2-Level-Ordering algorithm
      CALL tria_refine2LevelOrdering(&
          Rlevels(i-1)%rtriangulation, Rlevels(i)%rtriangulation)
      
      ! Create a standard mesh
      CALL tria_initStandardMeshFromRaw(Rlevels(i)%rtriangulation)
    
    END DO
    
    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we only have one block.
    ! Do this for all levels
    DO i = 1, nlevels
      CALL spdiscr_initBlockDiscr (&
          Rlevels(i)%rdiscretisation, 1, Rlevels(i)%rtriangulation)
    END DO
    
    ! In the next step, we will define the element type and the cubature
    ! formula that is to be used. For this 1D poisson-example we currently
    ! have 2 possible element types: linear and quadratic ones.
    ! For linear elements the trapezoidal formula satisfies all our
    ! needs, for quadratic elements we should choose a 3-point Gauss-formula.
    !
    ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
    ! structures for every component of the solution vector.
    ! Initialise the first element of the list to specify the element
    ! and cubature rule for this solution component:
    DO i = 1, nlevels
      CALL spdiscr_initDiscr_simple (&
          Rlevels(i)%rdiscretisation%RspatialDiscretisation(1), &
      ! Setting up a linear element and trapezoidal rule would be...
          EL_P1_1D,CUB_TRZ_1D,Rlevels(i)%rtriangulation)
      ! Setting up a quadratic element and 3-point Gauss rule would be...
          !EL_P2_1D,CUB_G3_1D,Rlevels(i)%rtriangulation)
    END DO

    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! We create a scalar matrix, based on the discretisation structure
    ! for our one and only solution component.
    DO i = 1, nlevels
      
      ! Initialise the block matrix with default values based on
      ! the discretisation.
      CALL lsysbl_createMatBlockByDiscr (&
          Rlevels(i)%rdiscretisation,Rlevels(i)%rmatrix)    

      ! Now as the discretisation is set up, we can start to generate
      ! the structure of the system matrix which is to solve.
      ! We create that directly in the block (1,1) of the block matrix
      ! using the discretisation structure of the first block.
      CALL bilf_createMatrixStructure ( &
           Rlevels(i)%rdiscretisation%RspatialDiscretisation(1),&
           LSYSSC_MATRIX9,Rlevels(i)%rmatrix%RmatrixBlock(1,1))

      ! Update the structural information of the block matrix, as we manually
      ! changed one of the submatrices:
      CALL lsysbl_updateMatStrucInfo (Rlevels(i)%rmatrix)
      
      ! And now to the entries of the matrix. For assembling of the entries,
      ! we need a bilinear form, which first has to be set up manually.
      ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
      ! scalar system matrix in 1D.
      rform%itermCount = 1
      rform%Idescriptors(1,1) = DER_DERIV_X
      rform%Idescriptors(2,1) = DER_DERIV_X

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
      CALL bilf_buildMatrixScalar (rform,.TRUE.,&
           Rlevels(i)%rmatrix%RmatrixBlock(1,1),coeff_Laplace_1D)

    END DO
    ! Keep in mind that now p_rdisc and p_rmatrix are pointers to the
    ! discretisation and matrix structures on the finest level!
    
    ! Although we could manually create the solution/RHS vector,
    ! the easiest way to set up the vector structure is
    ! to create it by using our matrix as template:
    CALL lsysbl_createVecBlockIndMat (Rlevels(nlevels)%rmatrix,rrhsBlock, .FALSE.)

    ! The vector structure is ready but the entries are missing. 
    ! So the next thing is to calculate the content of that vector.
    !
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    
    ! ... and then discretise the RHS to get a discrete version of it.
    CALL linf_buildVectorScalar (&
        Rlevels(nlevels)%rdiscretisation%RspatialDiscretisation(1),&
        rlinform,.TRUE.,rrhsBlock%RvectorBlock(1),coeff_RHS_1D)
    
    ! Now we have the raw problem. What is missing is the definition of the boudary
    ! conditions.
    ! For implementing boundary conditions, we use a 'filter technique with
    ! discretised boundary conditions'. This means, we first have to calculate
    ! a discrete version of the analytic BC, which we can implement into the
    ! solution/RHS vectors using the corresponding filter.
    !
    ! In contrast to the 2D poisson examples, we will directly set the
    ! dirichlet boundary conditions by hand instead of discretising an analytic
    ! boundary condition function using a boundary structure.
    !
    ! Set p_rdiscreteBC to NULL -- bcasm_initDirichletBC_1D will allocate it.
    DO i = 1, nlevels
    
      NULLIFY(Rlevels(i)%p_rdiscreteBC)

      CALL bcasm_initDirichletBC_1D(&
          Rlevels(i)%rdiscretisation, Rlevels(i)%p_rdiscreteBC,0.0_DP, 0.0_DP)

      ! Hang the pointer into the matrix. That way, these
      ! boundary conditions are always connected to that matrix.
      Rlevels(i)%rmatrix%p_rdiscreteBC => Rlevels(i)%p_rdiscreteBC
      
      ! Also implement the boundary conditions into the matrix.
      CALL matfil_discreteBC (Rlevels(i)%rmatrix)
  
    END DO
    
    ! Our one and only right-hand-side also needs to know the boundary conditions.
    rrhsBlock%p_rdiscreteBC => Rlevels(nlevels)%p_rdiscreteBC
                             
    ! Now we have block vectors for the RHS and the matrix. What we
    ! need additionally is a block vector for the solution and
    ! temporary data. Create them using the RHS as template.
    ! Fill the solution vector with 0:
    CALL lsysbl_createVecBlockIndirect (rrhsBlock, rvectorBlock, .TRUE.)
    CALL lsysbl_createVecBlockIndirect (rrhsBlock, rtempBlock, .FALSE.)
    
    ! Next step is to implement boundary conditions into the RHS and
    ! solution vectors. This is done using a vector filter
    ! for discrete boundary conditions.
    ! The discrete boundary conditions are already attached to the
    ! vectors. Call the appropriate vector filter that
    ! modifies the vectors according to the boundary conditions.
    CALL vecfil_discreteBCrhs (rrhsBlock)
    CALL vecfil_discreteBCsol (rvectorBlock)

    ! During the linear solver, the boundary conditions are also
    ! frequently imposed to the vectors. But as the linear solver
    ! does not work with the actual solution vectors but with
    ! defect vectors instead.
    ! So, set up a filter chain that filters the defect vector
    ! during the solution process to implement discrete boundary conditions.
    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

    ! Now we have to build up the level information for multigrid.
    !
    ! At first, initialise a standard interlevel projection structure. We
    ! can use the same structure for all levels.
    CALL mlprj_initProjectionMat (rprojection,Rlevels(nlevels)%rmatrix)

    ! Create a Multigrid-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    p_RfilterChain => RfilterChain
    CALL linsol_initMultigrid (p_rsolverNode,p_RfilterChain)
    
    ! Set up a coarse grid solver.
    CALL linsol_initUMFPACK4 (p_rcoarsegridSolver)
    
    ! Add the coarse grid level.
    CALL linsol_addMultiGridLevel(p_rlevelInfo,p_rsolverNode,rprojection,&
                                  NULL(), NULL(), p_rcoarseGridSolver)

    ! Now set up the other levels...
    DO i = 2, nlevels
    
      ! Create a Jacobi smoother
      CALL linsol_initJacobi(p_rsmoother)
      
      ! We will use 4 smoothing steps with damping parameter 0.7
      CALL linsol_convertToSmoother(p_rsmoother, 4, 0.7_DP)
      
      ! And add this multi-grid level. We will use the same smoother
      ! for pre- and post-smoothing.
      CALL linsol_addMultiGridLevel(p_rlevelInfo,p_rsolverNode,rprojection,&
                                    p_rsmoother, p_rsmoother, NULL())
      
    END DO

    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2
    
    ! Attach the system matrices to the solver.
    !
    ! We copy our matrices to a big matrix array and transfer that
    ! to the setMatrices routines. This intitialises then the matrices
    ! on all levels according to that array. Note that this does not
    ! allocate new memory, we create only 'links' to existing matrices
    ! into Rmatrices(:)!
    ALLOCATE(Rmatrices(nlevels))
    DO i=1,nlevels
      CALL lsysbl_duplicateMatrix (Rlevels(i)%rmatrix,&
          Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    END DO
    
    CALL linsol_setMatrices(p_RsolverNode,Rmatrices(1:nlevels))

    ! We can release Rmatrices immediately -- as long as we don't
    ! release rproblem%RlevelInfo(i)%rmatrix!
    DO i=1,nlevels
      CALL lsysbl_releaseMatrix (Rmatrices(i))
    END DO
    DEALLOCATE(Rmatrices)
    
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
    
    ! That's it, rvectorBlock now contains our solution. We can now
    ! start the postprocessing. 
    CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,&
        Rlevels(nlevels)%rtriangulation,'gmv/u0_mg.gmv')
    
    CALL lsyssc_getbase_double (rvectorBlock%RvectorBlock(1),p_Ddata)
    CALL ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that's it.
    CALL ucd_write (rexport)
    CALL ucd_release (rexport)
    
    ! Calculate the error to the reference function.
    CALL pperr_scalar (rvectorBlock%RvectorBlock(1),PPERR_L2ERROR,derror,&
                       getReferenceFunction_1D)
    CALL output_line ('L2-error: ' // sys_sdEL(derror,10) )
    CALL pperr_scalar (rvectorBlock%RvectorBlock(1),PPERR_H1ERROR,derror,&
                       getReferenceFunction_1D)
    CALL output_line ('H1-error: ' // sys_sdEL(derror,10) )

    ! We are finished - but not completely!
    ! Now, clean up so that all the memory is available again.
    !
    ! Release solver data and structure
    CALL linsol_doneData (p_rsolverNode)
    CALL linsol_doneStructure (p_rsolverNode)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    CALL linsol_releaseSolver (p_rsolverNode)
    
    ! Release the multilevel projection structure.
    CALL mlprj_doneProjection (rprojection)
    
    ! Release the block matrices/vectors
    CALL lsysbl_releaseVector (rtempBlock)
    CALL lsysbl_releaseVector (rvectorBlock)
    CALL lsysbl_releaseVector (rrhsBlock)
    DO i = nlevels, 1, -1
      CALL lsysbl_releaseMatrix (Rlevels(i)%rmatrix)
    END DO

    ! Release our discrete version of the boundary conditions
    DO i = nlevels, 1, -1
      CALL bcasm_releaseDiscreteBC (Rlevels(i)%p_rdiscreteBC)
    END DO

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    DO i = nlevels, 1, -1
      CALL spdiscr_releaseBlockDiscr(Rlevels(i)%rdiscretisation)
    END DO
    
    ! Release the triangulation. 
    DO i = nlevels, 1, -1
      CALL tria_done (Rlevels(i)%rtriangulation)
    END DO

    ! Release level information
    DEALLOCATE(Rlevels)
    
    ! That's it!
    
  END SUBROUTINE

END MODULE
