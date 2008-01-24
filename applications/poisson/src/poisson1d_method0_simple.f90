!##############################################################################
!# ****************************************************************************
!# <name> poisson1d_method0_simple </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple 1D Poisson
!# problem with constant coefficients on a simple domain.
!# </purpose>
!##############################################################################

MODULE poisson1d_method0_simple

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

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE poisson1d_0_simple
  
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
    
    ! Error indicator during initialisation of the solver
    INTEGER :: ierror    
    
    ! Error of FE function to reference function
    REAL(DP) :: derror

    ! Output block for UCD output to GMV file
    TYPE(t_ucdExport) :: rexport
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata

    ! The number of sub-intervals for the discretisation
    INTEGER :: nintervals = 16
    
    ! Ok, let's start. 
    ! At first, create the basic triangulation.
    ! Our domain is [0, 1], divided into nintervals sub-intervals.
    CALL tria_createRawTria1D(rtriangulation, 0.0_DP, 1.0_DP, nintervals)
    
    ! As the tria_createRawTria1D routine always generates a grid
    ! with sub-intervals of equal length, we can optionally disturb
    ! the mesh. This will result in an unsymmetric matrix.
    !CALL meshmod_disturbMesh(rtriangulation, 0.2_DP)

    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    CALL tria_initStandardMeshFromRaw (rtriangulation)
    
    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we only have one block.
    CALL spdiscr_initBlockDiscr (rdiscretisation, 1, rtriangulation)
    
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
    CALL spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscretisation(1), &
    ! Setting up a linear element and trapezoidal rule would be...
                                   EL_P1_1D,CUB_TRZ_1D,rtriangulation)
    ! Setting up a quadratic element and 3-point Gauss rule would be...
                                   !EL_P2_1D,CUB_G3_1D,rtriangulation)
    ! Setting up a cubic spline element and 4-point Gauss rule would be...
                                   !EL_S31_1D,CUB_G4_1D,rtriangulation)
                                   
                 
    ! We will set the evaluation cubature formula to 3-point Gauss.
    ! If we don't do this, then the L2-error, which is calculated in the
    ! post-processing phase would be beyond machine exactness...
    rdiscretisation%RspatialDiscretisation(1)%RelementDistribution(1)%ccubTypeEval = &
      CUB_G6_1D

    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! We create a scalar matrix, based on the discretisation structure
    ! for our one and only solution component.
    CALL bilf_createMatrixStructure (rdiscretisation%RspatialDiscretisation(1),&
                                     LSYSSC_MATRIX9,rmatrix)
    
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
    CALL bilf_buildMatrixScalar (rform,.TRUE.,rmatrix,coeff_Laplace_1D)
    
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
                                 rlinform,.TRUE.,rrhs,coeff_RHS_1D)
    
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
    ! In contrast to the 2D poisson examples, we will directly set the
    ! dirichlet boundary conditions by hand instead of discretising an analytic
    ! boundary condition function using a boundary structure.
    !
    ! Set p_rdiscreteBC to NULL -- bcasm_initDirichletBC_1D will allocate it.
    NULLIFY(p_rdiscreteBC)
    CALL bcasm_initDirichletBC_1D(rdiscretisation, p_rdiscreteBC, 0.0_DP, 0.0_DP)
                             
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

    ! Attach the above filter chain to the solver, so that the solver
    ! automatically filters the vector during the solution process.
    p_RfilterChain => RfilterChain
    
    ! We now have the option to create a preconditioner for the solver.
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Preconditioner Remark:
    ! ----------------------
    ! Please keep in mind that if the 1D grid for the discretisation was
    ! created using the tria_createRawTria1D routine (which is the default
    ! case in this example) and the grid was NOT refined afterwards (which is
    ! also the default case in this example), the resulting matrix will be
    ! tri-diagonal for linear elements.
    ! In this case, a LU-decomposition would produce no fill-in - therefore the
    ! ILU(0) preconditioner computes a complete LU-decomposition of the system
    ! matrix. So do not wonder if you set ILU(0) as a preconditioner and the
    ! solver always converges after 1 iteration... ^_^
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    NULLIFY(p_rpreconditioner)
    ! Setting up a Jacobi-Preconditioner would be..
    !CALL linsol_initJacobi (p_rpreconditioner)
    ! Setting up a Jin-Wei-Tam-Preconditioner would be...
    !CALL linsol_initJinWeiTam (p_rpreconditioner)
    ! Setting up a SOR[1.2]-Preconditioner would be...
    !CALL linsol_initSOR (p_rpreconditioner, 1.2_DP)
    ! Setting up a SSOR[1.2]-Preconditioner would be...
    !CALL linsol_initSSOR (p_rpreconditioner, 1.2_DP)
    ! Setting up a ILU(0)-Preconditioner would be...
    !CALL linsol_initMILUs1x1(p_rpreconditioner, 0, 0.0_DP)
    
    ! We now need to create a solver for the linear system.
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Solver Remark:
    ! --------------
    ! Please keep in mind that the CG solver needs a symmetric preconditioner
    ! and will (most probably) not work with unsymmetric preconditioners as
    ! SOR or (M)ILU(s)
    ! Also remember that the CG solver might diverge if the grid was disturbed
    ! using the 'meshmod_disturbMesh' routine.
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Setting up a Defect-Correction-Solver would be...
    !CALL linsol_initDefCorr (p_rsolverNode,p_rpreconditioner,p_RfilterChain)
    ! Setting up a PCG-Solver would be...
    !CALL linsol_initCG (p_rsolverNode,p_rpreconditioner,p_RfilterChain)
    ! Setting up a BiCGStab-Solver would be...
    CALL linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,p_RfilterChain)
    ! Setting up a GMRES(17)-Solver would be...
    !CALL linsol_initGMRES (p_rsolverNode,17,p_rpreconditioner,p_RfilterChain)
    
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
    
    ! That's it, rvectorBlock now contains our solution. We can now
    ! start the postprocessing. 
    CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,'gmv/u0.gmv')
    
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

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    CALL spdiscr_releaseBlockDiscr(rdiscretisation)
    
    ! Release the triangulation. 
    CALL tria_done (rtriangulation)
    
  END SUBROUTINE

END MODULE
