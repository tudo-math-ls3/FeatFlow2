!##############################################################################
!# ****************************************************************************
!# <name> stokes2d_method0_simple </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstation program how to solve a Stokes
!# problem on a simple domain.
!#
!# The routine uses the BiCGStab solver with a simple-VANCA preconditioner
!# for 2D saddle point problems, Jacobi-Type.
!# </purpose>
!##############################################################################

MODULE stokes2d_method0_simple

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
  USE coarsegridcorrection
  USE spdiscprojection
  USE ucd
  
  USE stokes2d_callback
  
  IMPLICIT NONE

CONTAINS
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE stokes2d_0_simple
  
!<description>
  ! This is an all-in-one stokes solver for directly solving a stokes
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
    TYPE(t_boundary), POINTER :: p_rboundary

    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    TYPE(t_blockDiscretisation) :: rdiscretisation,rprjDiscretisation
    
    ! A bilinear and linear form describing the analytic problem to solve
    TYPE(t_bilinearForm) :: rform
    TYPE(t_linearForm) :: rlinform

    ! A scalar matrix and vector. The vector accepts the RHS of the problem
    ! in scalar form.
    TYPE(t_matrixScalar) :: rmatrixB1, rmatrixB2

    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    TYPE(t_matrixBlock) :: rmatrix
    TYPE(t_vectorBlock) :: rvector,rrhs,rtempBlock,rprjVector
    
    ! A set of variables describing the analytic and discrete boundary
    ! conditions.    
    TYPE(t_boundaryConditions), POINTER :: p_rboundaryConditions
    TYPE(t_boundaryRegion) :: rboundaryRegion
    TYPE(t_bcRegion), POINTER :: p_rbcRegion
    TYPE(t_discreteBC), POINTER :: p_rdiscreteBC, p_rprjDiscreteBC

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
    
    ! Viscosity parameter nu = 1/Re
    REAL(DP) :: dnu
    
    ! Error indicator during initialisation of the solver
    INTEGER :: ierror    
    
    ! Output block for UCD output to GMV file
    TYPE(t_ucdExport) :: rexport
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_Ddata2

    ! Ok, let's start. 
    !
    ! We want to solve our Poisson problem on level...
    ! As we do not use a multigrid solver here, we'll set the level
    ! to 5 instead of 7.
    NLMAX = 5
    
    ! Viscosity parameter:
    dnu = 1.0_DP

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    ! Set p_rboundary to NULL() to create a new structure.
    NULLIFY(p_rboundary)
    CALL boundary_read_prm(p_rboundary, './pre/QUAD.prm')
        
    ! Now read in the basic triangulation.
    CALL tria_readTriFile2D (rtriangulation, './pre/QUAD.tri', p_rboundary)
    
    ! Refine the mesh up to the minimum level
    CALL tria_quickRefine2LevelOrdering(NLMAX-1,rtriangulation,p_rboundary)
    
    ! Create information about adjacencies and everything one needs from
    ! a triangulation. Afterwards, we have the coarse mesh.
    CALL tria_initStandardMeshFromRaw (rtriangulation,p_rboundary)

    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies 3 blocks in the
    ! solution vector.
    CALL spdiscr_initBlockDiscr2D (rdiscretisation,3,&
                                   rtriangulation, p_rboundary)

    ! rdiscretisation%RspatialDiscretisation is a list of scalar 
    ! discretisation structures for every component of the solution vector.
    ! We have a solution vector with three components:
    !  Component 1 = X-velocity
    !  Component 2 = Y-velocity
    !  Component 3 = Pressure
    ! For simplicity, we set up one discretisation structure for the 
    ! velocity...
    CALL spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscretisation(1),&
                EL_EM30, CUB_G2X2, rtriangulation, p_rboundary)
                
    ! ...and copy this structure also to the discretisation structure
    ! of the 2nd component (Y-velocity). This needs no additional memory, 
    ! as both structures will share the same dynamic information afterwards.
    CALL spdiscr_duplicateDiscrSc(rdiscretisation%RspatialDiscretisation(1),&
        rdiscretisation%RspatialDiscretisation(2))

    ! For the pressure (3rd component), we set up a separate discretisation 
    ! structure, as this uses different finite elements for trial and test
    ! functions.
    CALL spdiscr_initDiscr_combined (rdiscretisation%RspatialDiscretisation(3),&
                EL_Q0,EL_EM30,CUB_G2X2,rtriangulation,p_rboundary)

    ! Initialise the block matrix with default values based on
    ! the discretisation.
    CALL lsysbl_createMatBlockByDiscr (rdiscretisation,rmatrix)    
    
    ! Inform the matrix that we build a saddle-point problem.
    ! Normally, imatrixSpec has the value LSYSBS_MSPEC_GENERAL,
    ! but probably some solvers can use the special structure later.
    rmatrix%imatrixSpec = LSYSBS_MSPEC_SADDLEPOINT
    
    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! We create that directly in the block (1,1) of the block matrix
    ! using the discretisation structure of the first block.
    !
    ! Create the matrix structure of the X-velocity.
    CALL bilf_createMatrixStructure (rdiscretisation%RspatialDiscretisation(1),&
                                     LSYSSC_MATRIX9, rmatrix%RmatrixBlock(1,1))

    ! In the Stokes problem, the matrix for the Y-velocity is identical to
    ! the matrix for the X-velocity; both are Laplace-matrices!
    ! Therefore, we can simply make a copy of the matrix for the X-velocity.
    ! This we do later after the entries are created.
    !
    ! In the global system, there are two coupling matrices B1 and B2.
    ! Both have the same structure.
    !
    !    ( A         B1 )
    !    (      A    B2 )
    !    ( B1^T B2^T    )
    !
    ! Create the matrices structure of the pressure using the 3rd
    ! spatial discretisation structure in p_rdiscretisation%RspatialDiscretisation.
    CALL bilf_createMatrixStructure (rdiscretisation%RspatialDiscretisation(3),&
                                     LSYSSC_MATRIX9, rmatrixB1)
              
    ! Duplicate the B1 matrix structure to the B2 matrix, so use
    ! lsyssc_duplicateMatrix to create B2. Share the matrix 
    ! structure between B1 and B2 (B1 is the parent and B2 the child). 
    ! Don't create a content array yet, it will be created by 
    ! the assembly routines later.
    CALL lsyssc_duplicateMatrix (rmatrixB1, rmatrixB2, LSYSSC_DUP_COPY,&
                                 LSYSSC_DUP_REMOVE)
                                     
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
    rform%Dcoefficients(1)  = dnu
    rform%Dcoefficients(2)  = dnu

    ! Now we can build the matrix entries.
    ! We specify the callback function coeff_Laplace for the coefficients.
    ! As long as we use constant coefficients, this routine is not used.
    ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
    ! the framework will call the callback routine to get analytical data.
    !
    ! We pass our collection structure as well to this routine, 
    ! so the callback routine has access to everything what is
    ! in the collection.
    !
    ! Build the X-velocity matrix:
    CALL bilf_buildMatrixScalar (rform,.TRUE.,rmatrix%RmatrixBlock(1,1),&
                                 coeff_Stokes_2D)
    
    ! Duplicate the matrix to the Y-velocity matrix, share structure and
    ! content between them (as the matrices are the same).
    CALL lsyssc_duplicateMatrix (rmatrix%RmatrixBlock(1,1),&
             rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    
    ! Manually change the discretisation structure of the Y-velocity 
    ! matrix to the Y-discretisation structure.
    ! Ok, we use the same discretisation structure for both, X- and Y-velocity,
    ! so this is not really necessary - we do this for sure...
    rmatrix%RmatrixBlock(2,2)%p_rspatialDiscretisation => &
      rdiscretisation%RspatialDiscretisation(2)
                                
    ! Build the first pressure matrix B1.
    ! Again first set up the bilinear form, then call the matrix assembly.
    rform%itermCount = 1
    rform%Idescriptors(1,1) = DER_FUNC
    rform%Idescriptors(2,1) = DER_DERIV_X

    ! In the standard case, we have constant coefficients:
    rform%ballCoeffConstant = .TRUE.
    rform%BconstantCoeff = .TRUE.
    rform%Dcoefficients(1)  = -1.0_DP
    
    CALL bilf_buildMatrixScalar (rform,.TRUE.,rmatrixB1,coeff_Pressure_2D)

    ! Build the second pressure matrix B2.
    ! Again first set up the bilinear form, then call the matrix assembly.
    rform%itermCount = 1
    rform%Idescriptors(1,1) = DER_FUNC
    rform%Idescriptors(2,1) = DER_DERIV_Y

    ! In the standard case, we have constant coefficients:
    rform%ballCoeffConstant = .TRUE.
    rform%BconstantCoeff = .TRUE.
    rform%Dcoefficients(1)  = -1.0_DP
    
    CALL bilf_buildMatrixScalar (rform,.TRUE.,rmatrixB2,coeff_Pressure_2D)
                                
    ! The B1/B2 matrices exist up to now only in our local problem structure.
    ! Put a copy of them into the block matrix.
    !
    ! Note that we share the structure of B1/B2 with those B1/B2 of the
    ! block matrix, while we create copies of the entries. The reason is
    ! that these matrices are modified for boundary conditions later.
    CALL lsyssc_duplicateMatrix (rmatrixB1, rmatrix%RmatrixBlock(1,3),&
                                 LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

    CALL lsyssc_duplicateMatrix (rmatrixB2, rmatrix%RmatrixBlock(2,3),&
                                 LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
    
    ! Furthermore, put B1^T and B2^T to the block matrix.
    CALL lsyssc_transposeMatrix (rmatrixB1, rmatrix%RmatrixBlock(3,1),&
                                 LSYSSC_TR_VIRTUAL)

    CALL lsyssc_transposeMatrix (rmatrixB2, rmatrix%RmatrixBlock(3,2),&
                                 LSYSSC_TR_VIRTUAL)

    ! Update the structural information of the block matrix, as we manually
    ! changed the submatrices:
    CALL lsysbl_updateMatStrucInfo (rmatrix)

    ! Although we could manually create the solution/RHS vector,
    ! the easiest way to set up the vector structure is
    ! to create it by using our matrix as template:
    CALL lsysbl_createVecBlockIndMat (rmatrix,rrhs, .FALSE.)
    CALL lsysbl_createVecBlockIndMat (rmatrix,rvector, .FALSE.)

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
    ! Note that the vector is unsorted after calling this routine!
    CALL linf_buildVectorScalar (rdiscretisation%RspatialDiscretisation(1),&
                  rlinform,.TRUE.,rrhs%RvectorBlock(1),coeff_RHS_2D)
                                
    ! The third subvector must be zero - as it represents the RHS of
    ! the equation "div(u) = 0".
    CALL lsyssc_clearVector(rrhs%RvectorBlock(3))
                                
    ! Clear the solution vector on the finest level.
    CALL lsysbl_clearVector(rvector)

    ! For implementing boundary conditions, we use a 'filter technique with
    ! discretised boundary conditions'. This means, we first have to calculate
    ! a discrete version of the analytic BC, which we can implement into the
    ! solution/RHS vectors using the corresponding filter.
    !
    ! At first, we need the analytic description of the boundary conditions.
    ! Initialise a structure for boundary conditions, which accepts this,
    ! on the heap.
    !
    ! We first set up the boundary conditions for the X-velocity, then those
    ! of the Y-velocity.
    !
    ! Set p_rboundaryConditions to NULL() to create a new structure on the heap.
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
    
    ! The endpoint of this segment should also be Dirichlet. We set this by
    ! changing the region properties in rboundaryRegion.
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    
    ! We use this boundary region and specify that we want to have Dirichlet
    ! boundary there. The following routine adds a new 'boundary condition region'
    ! for the first segment to the boundary condition structure.
    ! The region will be set up as 'Dirichlet boundary'.
    ! We specify icomponent='1' to indicate that we set up the
    ! Dirichlet BC's for the firstcomponent in the solution vector,
    ! the X-velocity.
    ! The routine also returns the created object in p_rbcRegion so that we can
    ! modify it - but accept it as it is, so we can ignore that.
    CALL bcond_newDirichletBConRealBD (p_rboundaryConditions,1,&
                                       rboundaryRegion,p_rbcRegion)
                              
    ! Now to the edge 2 of boundary component 1 the domain. We use the
    ! same two routines to add the boundary condition to p_rboundaryConditions.
    !
    ! Edge 2 is Neumann boudary
    !CALL boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
    !CALL bcond_newDirichletBConRealBD (p_rboundaryConditions,1,&
    !                                   rboundaryRegion,p_rbcRegion)
                              
    ! Edge 3 of boundary component 1.
    CALL boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
    CALL bcond_newDirichletBConRealBD (p_rboundaryConditions,1,&
                                       rboundaryRegion,p_rbcRegion)
    
    ! Edge 4 of boundary component 1. That's it.
    CALL boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
    CALL bcond_newDirichletBConRealBD (p_rboundaryConditions,1,&
                                       rboundaryRegion,p_rbcRegion)
      
    ! Now continue with defining the boundary conditions of the Y-velocity:
    !
    ! Define edge 1.
    CALL boundary_createRegion(p_rboundary,1,1,rboundaryRegion)
    
    ! As we define the Y-velocity, we now set icomponent=2 in the following call.
    CALL bcond_newDirichletBConRealBD (p_rboundaryConditions,2,&
                                       rboundaryRegion,p_rbcRegion)
     
    ! Define edge 2 - Neumann boundary                         
    !CALL boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
    !CALL bcond_newDirichletBConRealBD (p_rboundaryConditions,2,&
    !                                   rboundaryRegion,p_rbcRegion)
                              
    ! Define Edge 3
    CALL boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
    CALL bcond_newDirichletBConRealBD (p_rboundaryConditions,2,&
                                       rboundaryRegion,p_rbcRegion)
    
    ! Define Edge 4. That's it.
    CALL boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
    CALL bcond_newDirichletBConRealBD (p_rboundaryConditions,2,&
                                       rboundaryRegion,p_rbcRegion)

    ! The boundary conditions are set up, but still the block discretisation
    ! does not know about it. So inform the discretisation which
    ! analytic boundary conditions to use:
    rdiscretisation%p_rboundaryConditions => p_rboundaryConditions

    ! The pressure does not need boundary conditions.

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
    rmatrix%p_rdiscreteBC => p_rdiscreteBC
    rrhs%p_rdiscreteBC => p_rdiscreteBC
    rvector%p_rdiscreteBC => p_rdiscreteBC
    
    ! Next step is to implement boundary conditions into the RHS,
    ! solution and matrix. This is done using a vector/matrix filter
    ! for discrete boundary conditions.
    ! The discrete boundary conditions are already attached to the
    ! vectors/matrix. Call the appropriate vector/matrix filter that
    ! modifies the vectors/matrix according to the boundary conditions.
    CALL vecfil_discreteBCrhs (rrhs)
    CALL vecfil_discreteBCsol (rvector)
    CALL matfil_discreteBC (rmatrix)

    ! During the linear solver, the boundary conditions must
    ! frequently be imposed to the vectors. This is done using
    ! a filter chain. As the linear solver does not work with 
    ! the actual solution vectors but with defect vectors instead,
    ! a filter for implementing the real boundary conditions 
    ! would be wrong.
    ! Therefore, create a filter chain with one filter only,
    ! which implements Dirichlet-conditions into a defect vector.
    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

    ! Create a BiCGStab-solver with VANCA preconditioner.
    ! Attach the above filter chain to the solver, so that the solver
    ! automatically filters the vector during the solution process.
    p_RfilterChain => RfilterChain
    NULLIFY(p_rpreconditioner)
    CALL linsol_initVANCA (p_rpreconditioner,1.0_DP,LINSOL_VANCA_2DNAVST)
    CALL linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,p_RfilterChain)

    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2

    ! We will allow the solver to perform 200 iterations
    p_rsolverNode%nmaxIterations = 200

    ! Attach the system matrix to the solver.
    ! First create an array with the matrix data (on all levels, but we
    ! only have one level here), then call the initialisation 
    ! routine to attach all these matrices.
    ! Remark: Don't make a call like
    !    CALL linsol_setMatrices(p_RsolverNode,(/p_rmatrix/))
    ! This doesn't work on all compilers, since the compiler would have
    ! to create a temp array on the stack - which does not always work!
    Rmatrices = (/rmatrix/)
    CALL linsol_setMatrices(p_rsolverNode,Rmatrices)
    
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

    ! The solution vector is probably not in the way GMV likes it!
    ! GMV for example does not understand Q1~ vectors!
    ! Therefore, we first have to convert the vector to a form that
    ! GMV understands.
    ! GMV understands only Q1 solutions! So the task is now to create
    ! a Q1 solution from p_rvector and write that out.
    !
    ! For this purpose, first create a 'derived' simple discretisation
    ! structure based on Q1 by copying the main guiding block discretisation
    ! structure and modifying the discretisation structures of the
    ! two velocity subvectors:
    
    CALL spdiscr_duplicateBlockDiscr (rdiscretisation,rprjDiscretisation)
    
    CALL spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscretisation(1), &
                 EL_Q1, CUB_G2X2, rprjDiscretisation%RspatialDiscretisation(1))

    CALL spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscretisation(2), &
                 EL_Q1, CUB_G2X2, rprjDiscretisation%RspatialDiscretisation(2))
                 
    ! The pressure discretisation substructure stays the old.
    !
    ! Now set up a new solution vector based on this discretisation,
    ! allocate memory.
    CALL lsysbl_createVecBlockByDiscr (rprjDiscretisation,rprjVector,.FALSE.)
    
    ! Then take our original solution vector and convert it according to the
    ! new discretisation:
    CALL spdp_projectSolution (rvector,rprjVector)
    
    ! Discretise the boundary conditions according to the Q1/Q1/Q0 
    ! discretisation:
    NULLIFY(p_rprjDiscreteBC)
    CALL bcasm_discretiseBC (rprjDiscretisation,p_rprjDiscreteBC,.FALSE.,&
                             getBoundaryValues_2D)
                            
    ! Connect the vector to the BC's
    rprjVector%p_rdiscreteBC => p_rprjDiscreteBC
    
    ! Send the vector to the boundary-condition implementation filter.
    ! This modifies the vector according to the attached discrete boundary
    ! conditions.
    CALL vecfil_discreteBCsol (rprjVector)
    
    ! Now we have a Q1/Q1/Q0 solution in rprjVector.
    ! We can now start the postprocessing. 
    ! Start UCD export to GMV file:
    CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
        'gmv/u2d_0_simple.gmv')

    ! Write velocity field
    CALL lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
    CALL lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
    
    ! In case we use the VTK exporter, which supports vector output, we will
    ! pass the X- and Y-velocity at once to the ucd module.
    CALL ucd_addVarVertBasedVec(rexport,'velocity',p_Ddata,p_Ddata2)

    ! If we use the GMV exporter, we might replace the line above by the
    ! following two lines:
    !CALL ucd_addVariableVertexBased (rexport,'X-vel',UCD_VAR_XVELOCITY, p_Ddata)
    !CALL ucd_addVariableVertexBased (rexport,'Y-vel',UCD_VAR_YVELOCITY, p_Ddata2)
        
    ! Write pressure
    CALL lsyssc_getbase_double (rprjVector%RvectorBlock(3),p_Ddata)
    CALL ucd_addVariableElementBased (rexport,'pressure',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that's it.
    CALL ucd_write (rexport)
    CALL ucd_release (rexport)

    ! We are finished - but not completely!
    ! Now, clean up so that all the memory is available again.
    !
    ! Release solver data and structure
    CALL linsol_doneData (p_rsolverNode)
    CALL linsol_doneStructure (p_rsolverNode)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    CALL linsol_releaseSolver (p_rsolverNode)
    
    ! Release the block matrix/vectors
    CALL lsysbl_releaseVector (rprjVector)
    CALL lsysbl_releaseVector (rtempBlock)
    CALL lsysbl_releaseVector (rvector)
    CALL lsysbl_releaseVector (rrhs)
    CALL lsysbl_releaseMatrix (rmatrix)
    
    ! Release B1 and B2 matrix
    CALL lsyssc_releaseMatrix (rmatrixB2)
    CALL lsyssc_releaseMatrix (rmatrixB1)
    
    ! Release our discrete version of the boundary conditions
    CALL bcasm_releaseDiscreteBC (p_rprjDiscreteBC)
    CALL bcasm_releaseDiscreteBC (p_rdiscreteBC)

    ! ...and also the corresponding analytic description.
    CALL bcond_doneBC (p_rboundaryConditions)
    
    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    CALL spdiscr_releaseBlockDiscr(rprjDiscretisation)
    CALL spdiscr_releaseBlockDiscr(rdiscretisation)
    
    ! Release the triangulation. 
    CALL tria_done (rtriangulation)
    
    ! Finally release the domain, that's it.
    CALL boundary_release (p_rboundary)

  END SUBROUTINE

END MODULE