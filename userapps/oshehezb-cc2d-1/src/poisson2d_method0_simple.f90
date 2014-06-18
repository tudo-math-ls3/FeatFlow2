!##############################################################################
!# ****************************************************************************
!# <name> poisson2d_method0_simple </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple Poisson
!# problem with constant coefficients on a simple domain.
!# </purpose>
!##############################################################################

MODULE poisson2d_method0_simple

  USE fsystem
  USE genoutput
  USE genoutput
  USE storage
  USE linearsolver
  USE boundary
  USE linearformevaluation
  USE cubature
  USE matrixfilters
  USE vectorfilters
  use discretebc
  USE bcassembly
  USE triangulation
  USE linearsystemblock
  USE spatialdiscretisation
  USE ucd
  USE pprocerror
  use scalarpde
  USE bilinearformevaluation
  USE trilinearformevaluation
  use filtersupport
  use linearsystemscalar
  use linearsystemblock
    
  USE poisson2d_callback
  
  IMPLICIT NONE


CONTAINS

 
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE poisson2d_0_simple_core(rvector,rvelocity,dtstep,dtime)
  
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
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
    TYPE(t_boundary), POINTER :: p_rboundary
    TYPE(t_triangulation), POINTER :: p_rtriangulation
      
    ! A bilinear,trilinear and linear form describing the analytic problem to solve
    TYPE(t_bilinearForm) :: rform
    TYPE(t_trilinearForm) :: rform1
    TYPE(t_linearForm) :: rlinform

    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    TYPE(t_matrixBlock) :: rmatrix
    
    TYPE(t_vectorBlock) :: rtempBlock
    TYPE(t_vectorBlock),INTENT(INOUT), TARGET :: rvector
    TYPE(t_vectorBlock), INTENT(IN):: rvelocity
    TYPE(t_vectorBlock) :: rrhs
    
    ! A set of variables describing the discrete boundary conditions.
    TYPE(t_boundaryRegion) :: rboundaryRegion
    TYPE(t_discreteBC), TARGET :: rdiscreteBC
    TYPE(t_discreteBC), TARGET :: p_rdiscreteBC
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
    
    ! Error of FE function to reference function
    REAL(DP) :: derror
    
    
    ! Output block for UCD output to GMV file
    TYPE(t_ucdExport) :: rexport
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_Ddata1
    
    !Time-depend values
    REAL(DP), INTENT(IN):: dtime,dtstep
    
    ! Ok, let's start.
    !
    ! We want to solve our Poisson problem on level...
    NLMAX = 7
   
    CALL lsyssc_getbase_double (rvector%RvectorBlock(1),p_Ddata)
    p_rdiscretisation => rvector%p_rblockDiscr
    
    p_rboundary => rvector%p_rblockDiscr%p_rboundary
    
    p_rtriangulation => rvector%p_rblockDiscr%p_rtriangulation
    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! We create a scalar matrix, based on the discretisation structure
    ! for our one and only solution component.
    CALL lsysbl_createMatBlockByDiscr (p_rdiscretisation,rmatrix)
    CALL bilf_createMatrixStructure (p_rdiscretisation%RspatialDiscr(1),&
                                     LSYSSC_MATRIX9,rmatrix%RmatrixBlock(1,1))
    
        
    ! The same has to be done for the right hand side of the problem.
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    
    ! ... and then discretise the RHS to get a discrete version of it.
    ! Again we simply create a scalar vector based on the one and only
    ! discretisation structure.
    ! This scalar vector will later be used as the one and only first
    ! component in a block vector.
    CALL lsysbl_createVecBlockByDiscr (p_rdiscretisation,rrhs,.TRUE.)
    CALL linf_buildVectorScalar (p_rdiscretisation%RspatialDiscr(1),&
                                rlinform,.TRUE.,rrhs%RvectorBlock(1),coeff_RHS_2D)
                                    
    rform%itermCount = 1
    rform%Idescriptors(1,1) = DER_FUNC
    rform%Idescriptors(2,1) = DER_FUNC
   
    rform%ballCoeffConstant = .TRUE.
    rform%BconstantCoeff = .TRUE.
    rform%Dcoefficients(1)  = 1.0/dtstep
    
    CALL bilf_buildMatrixScalar (rform,.TRUE.,rmatrix%RmatrixBlock(1,1))
    CALL lsysbl_blockMatVec(rmatrix,rvector,rrhs,1.0_DP,1.0_DP)
    
    
    !Set up the correspponding trilinear form
    rform1%itermCount = 1
    rform1%Idescriptors(1,1) = DER_FUNC
    rform1%Idescriptors(2,1) = DER_DERIV_X
    rform1%Idescriptors(3,1) = DER_FUNC
    
    
    ! In the standard case, we have constant coefficients:
    rform1%ballCoeffConstant = .TRUE.
    rform1%BconstantCoeff = .TRUE.
    rform1%Dcoefficients(1)  = 1.0
    
    
    CALL trilf_buildMatrixScalar (rform1,.FALSE.,rmatrix%RmatrixBlock(1,1),rvelocity%RvectorBlock(1))
    
    rform1%Idescriptors(1,1) = DER_FUNC
    rform1%Idescriptors(2,1) = DER_DERIV_Y
    
    CALL trilf_buildMatrixScalar (rform1,.FALSE.,rmatrix%RmatrixBlock(1,1),rvelocity%RvectorBlock(2))
    
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
    rform%Dcoefficients(1)  = 0.001
    rform%Dcoefficients(2)  = 0.001

    ! Now we can build the matrix entries.
    ! We specify the callback function coeff_Laplace for the coefficients.
    ! As long as we use constant coefficients, this routine is not used.
    ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
    ! the framework will call the callback routine to get analytical
    ! data.
    CALL bilf_buildMatrixScalar (rform,.FALSE.,rmatrix%RmatrixBlock(1,1),coeff_Laplace_2D)
   
    ! Now we have the raw problem. What is missing is the definition of the boundary
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
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
                             
    ! Now to the edge 2 of boundary component 1 the domain.
    !CALL boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
    CALL bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
                             
    CALL boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
    CALL bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
    
    ! Edge 4 of boundary component 1. That's it.
    CALL boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
    CALL bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
    ! Edge 5 of boundary component 2. That's it.
    CALL boundary_createRegion(p_rboundary,2,1,rboundaryRegion)
    CALL bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
    ! Hang the pointer into the vector and matrix. That way, these
    ! boundary conditions are always connected to that matrix and that
    ! vector.
    rmatrix%p_rdiscreteBC => rdiscreteBC
    rrhs%p_rdiscreteBC => rdiscreteBC
    rvector%p_rdiscreteBC => rdiscreteBC
    ! Now we have block vectors for the RHS and the matrix. What we
    ! need additionally is a block vector for the solution and
    ! temporary data. Create them using the RHS as template.
    ! Fill the solution vector with 0:
    !CALL lsysbl_createVecBlockIndirect (rrhs, rvector, .TRUE.)
    CALL lsysbl_createVecBlockIndirect (rrhs, rtempBlock, .FALSE.)
    
    ! Next step is to implement boundary conditions into the RHS,
    ! solution and matrix. This is done using a vector/matrix filter
    ! for discrete boundary conditions.
    ! The discrete boundary conditions are already attached to the
    ! vectors/matrix. Call the appropriate vector/matrix filter that
    ! modifies the vectors/matrix according to the boundary conditions.
   
    CALL vecfil_discreteBCrhs (rrhs)
    CALL vecfil_discreteBCsol (rvector)
    CALL matfil_discreteBC (rmatrix)
   
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
    !CALL linsol_initSSOR (p_rpreconditioner,1.0_DP)
    !CALL linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,p_RfilterChain)
    CALL linsol_initUmfpack4 (p_rsolverNode)
    
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
    Rmatrices = (/rmatrix/)
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
    CALL linsol_solveAdaptively (p_rsolverNode,rvector,rrhs,rtempBlock)
    
    ! That's it, rvector now contains our solution. We can now
    ! start the postprocessing.
    ! Start UCD export to GMV file:
    CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,&
                    'gmv/u2d_0_simple.gmv')
    
    CALL lsyssc_getbase_double (rvector%RvectorBlock(1),p_Ddata)
    CALL ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that's it.
    CALL ucd_write (rexport)
    CALL ucd_release (rexport)
    
    ! Calculate the error to the reference function.
    CALL pperr_scalar (rvector%RvectorBlock(1),PPERR_L2ERROR,derror,&
                    getReferenceFunction_2D)
    CALL output_line ('L2-error: ' // sys_sdEL(derror,10) )

    CALL pperr_scalar (rvector%RvectorBlock(1),PPERR_H1ERROR,derror,&
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
    CALL lsysbl_releaseVector (rrhs)
    CALL lsysbl_releaseMatrix (rmatrix)

    ! Release the scalar matrix/rhs vector which were used to create
    ! the block matrices/vectors. These must exist as long as the
    ! block matrices/vectors exist, as the block matrices/vectors are
    ! only 'copies' of the scalar ones, sharing the same handles!
    !CALL lsysbl_releaseVector (rrhs)
    !CALL lsysbl_releaseMatrix (rmatrix)
  
    ! Release our discrete version of the boundary conditions
    CALL bcasm_releaseDiscreteBC (rdiscreteBC)

    CALL lsyssc_getbase_double (rvector%RvectorBlock(1),p_Ddata)
  END SUBROUTINE

END MODULE
