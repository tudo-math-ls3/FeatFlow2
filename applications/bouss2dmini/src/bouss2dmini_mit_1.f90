!##############################################################################
!# ****************************************************************************
!# <name> bouss2dmini_mit_1 </name>
!# ****************************************************************************
!#
!# <purpose>
!# Todo
!# </purpose>
!##############################################################################

MODULE bouss2dmini_mit_1

  USE fsystem
  USE storage
  USE linearsolver
  USE boundary
  USE trilinearformevaluation
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
  USE stdoperators
  USE convection
  USE statistics
  
  USE bouss2dmini_callback
  
  IMPLICIT NONE

!<types>

!<typeblock description="Type block defining all information about one level">

  TYPE t_level
  
    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation
    TYPE(t_blockDiscretisation) :: rdiscretisation
    
    ! A system matrix for that specific level. The matrix will receive the 
    ! discrete Laplace operator.
    TYPE(t_matrixBlock) :: rmatrix
    
    ! The solution vector for this level.
    TYPE(t_vectorBlock) :: rvecSol
    
    ! A-matrix for that specific level.
    TYPE(t_matrixScalar) :: rmatrixA

    ! B1-matrix for that specific level. 
    TYPE(t_matrixScalar) :: rmatrixB1

    ! B2-matrix for that specific level. 
    TYPE(t_matrixScalar) :: rmatrixB2
    
    ! M1- / M2-matrix for that specific level.
    TYPE(t_matrixScalar) :: rmatrixM
    
    ! Laplace-Matrix of the secondary system for that specific level.
    TYPE(t_matrixScalar) :: rmatrixN

    ! A variable describing the discrete boundary conditions.    
    TYPE(t_discreteBC) :: rdiscreteBC
  
  END TYPE
  
!</typeblock>

!</types>

CONTAINS
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE b2dm_mit_1
  
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
    ! An array of problem levels for the multigrid solver
    TYPE(t_level), DIMENSION(:), TARGET, ALLOCATABLE :: Rlevels

    ! An object for saving the domain:
    TYPE(t_boundary) :: rboundary

    ! A trilinear form describing the convection.
    TYPE(t_trilinearform) :: rtriform1, rtriform2
    
    ! A streamline diffusion object for the convection.
    TYPE(t_convStreamlineDiffusion) :: rsd

    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    TYPE(t_vectorBlock) :: rrhs,rtmp,rdef
    
    ! A set of variables describing the analytic and discrete boundary
    ! conditions.    
    TYPE(t_boundaryRegion) :: rboundaryRegion

    ! A solver node that accepts parameters for the linear solver    
    TYPE(t_linsolNode), POINTER :: p_rsolver,p_rpreconditioner,&
                                   p_rcoarseGridSolver,p_rsmoother

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    TYPE(t_matrixBlock), DIMENSION(:), ALLOCATABLE :: Rmatrices

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    TYPE(t_filterChain), DIMENSION(1), TARGET :: RfilterChain
    TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
    
    ! An interlevel projection structure for changing levels
    TYPE(t_interlevelProjectionBlock) :: rprojection

    ! One level of multigrid
    TYPE(t_linsolMGLevelInfo), POINTER :: p_rlevelInfo
    
    ! NLMIN receives the level of the coarse grid.
    INTEGER :: NLMIN

    ! NLMAX receives the level where we want to solve.
    INTEGER :: NLMAX
    
    ! Viscosity parameter nu = 1/Re
    REAL(DP) :: dnu
    
    ! Thermal conductivity eta
    REAL(DP) :: deta
    
    ! Gravity direction (dgrav1, dgrav2)
    REAL(DP) :: dgrav1, dgrav2
    
    ! Non-linear damping parameters
    REAL(DP) :: dnlDamping
    
    ! A counter variable for the non-linear loop
    INTEGER :: NLiter
    
    ! Non-linear residuals
    REAL(DP) :: dnlRes

    ! Error indicator during initialisation of the solver
    INTEGER :: ierror    
    
    ! Output block for UCD output to GMV file
    TYPE(t_ucdExport) :: rexport
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_Ddata2

    ! A counter variable
    INTEGER :: i
    LOGICAL :: bError
    
    ! A timer for time measuring
    TYPE(t_timer) :: rtimer

    ! Ok, let's start. 
    !
    ! We want to solve our Stokes problem on level...
    NLMIN = 3
    NLMAX = 6
    
    ! Viscosity parameter:
    dnu = 0.04_DP
    
    ! Thermal conductivity parameter:
    deta = 0.01_DP
    
    ! Gravity direction:
    dgrav1 =  0.0_DP
    dgrav2 = -1.0_DP
    
    ! Set the damping parameter for the non-linear loop.
    dnlDamping = 1.0_DP

    ! Allocate memory for all levels
    ALLOCATE(Rlevels(NLMIN:NLMAX))

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    CALL boundary_read_prm(rboundary, './pre/mit1.prm')
        
    ! Now read in the basic triangulation.
    CALL tria_readTriFile2D (Rlevels(NLMIN)%rtriangulation, &
                             './pre/mit1.tri', rboundary)
    
    ! Refine the mesh up to the minimum level
    CALL tria_quickRefine2LevelOrdering (NLMIN-1,&
        Rlevels(NLMIN)%rtriangulation,rboundary)
    
    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    CALL tria_initStandardMeshFromRaw (Rlevels(NLMIN)%rtriangulation,&
        rboundary)
    
    ! Now refine the grid for the fine levels.
    DO i = NLMIN+1, NLMAX

      ! Refine the grid using the 2-Level-Ordering algorithm
      CALL tria_refine2LevelOrdering(Rlevels(i-1)%rtriangulation,&
          Rlevels(i)%rtriangulation,rboundary)
      
      ! Create a standard mesh
      CALL tria_initStandardMeshFromRaw(Rlevels(i)%rtriangulation,&
        rboundary)
    
    END DO

    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies 4 blocks in the
    ! solution vector.
    DO i = NLMIN, NLMAX
      CALL spdiscr_initBlockDiscr2D (Rlevels(i)%rdiscretisation, 4, &
                                     Rlevels(i)%rtriangulation, rboundary)
    END DO

    ! rdiscretisation%RspatialDiscr is a list of scalar 
    ! discretisation structures for every component of the solution vector.
    ! We have a solution vector with three components:
    !  Component 1 = X-velocity
    !  Component 2 = Y-velocity
    !  Component 3 = Pressure
    !  Component 4 = Temperature
    DO i = NLMIN, NLMAX
      ! For simplicity, we set up one discretisation structure for the 
      ! velocity...
      CALL spdiscr_initDiscr_simple (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
          EL_EM30, CUB_G2X2, Rlevels(i)%rtriangulation, rboundary)
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory, 
      ! as both structures will share the same dynamic information afterwards.
      CALL spdiscr_duplicateDiscrSc (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
          Rlevels(i)%rdiscretisation%RspatialDiscr(2))

      ! For the pressure (3rd component), we set up a separate discretisation 
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      CALL spdiscr_deriveSimpleDiscrSc (Rlevels(i)%rdiscretisation%RspatialDiscr(1), &
          EL_Q0, CUB_G2X2, Rlevels(i)%rdiscretisation%RspatialDiscr(3))
    
      CALL spdiscr_deriveSimpleDiscrSc (Rlevels(i)%rdiscretisation%RspatialDiscr(1), &
          EL_Q1, CUB_G2X2, Rlevels(i)%rdiscretisation%RspatialDiscr(4))

    END DO

    DO i = NLMIN, NLMAX
    
      ! Initialise the block matrix with default values based on
      ! the discretisation.
      CALL lsysbl_createMatBlockByDiscr (Rlevels(i)%rdiscretisation,&
                                         Rlevels(i)%rmatrix)    
      
      ! Inform the matrix that we build a saddle-point problem.
      ! Normally, imatrixSpec has the value LSYSBS_MSPEC_GENERAL,
      ! but probably some solvers can use the special structure later.
      !Rlevels(i)%rmatrix%imatrixSpec = LSYSBS_MSPEC_SADDLEPOINT
      
      ! Now as the discretisation is set up, we can start to generate
      ! the structure of the system matrix which is to solve.
      ! We create that directly in the block (1,1) of the block matrix
      ! using the discretisation structure of the first block.
      !
      ! Create the matrix structure of the X-velocity.
      CALL bilf_createMatrixStructure (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
          LSYSSC_MATRIX9, Rlevels(i)%rmatrixA)

      ! Set up the matrix structure for the B1/B2-matrices
      CALL bilf_createMatrixStructure (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(3),&
          LSYSSC_MATRIX9, Rlevels(i)%rmatrixB1,&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1))

      CALL lsyssc_duplicateMatrix (Rlevels(i)%rmatrixB1, Rlevels(i)%rmatrixB2,&
                                   LSYSSC_DUP_SHARE, LSYSSC_DUP_REMOVE)
      
      ! Set up the matrix structure for the M1/M2 matrices
      CALL bilf_createMatrixStructure (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(4),&
          LSYSSC_MATRIX9, Rlevels(i)%rmatrixM,&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1))

      ! Create the matrix structure for the N-matrix
      CALL bilf_createMatrixStructure (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(4),&
          LSYSSC_MATRIX9, Rlevels(i)%rmatrixN)

      ! Build the A-Matrix.
      CALL stdop_assembleLaplaceMatrix(Rlevels(i)%rmatrixA, .TRUE., dnu)

      ! Build the B1-matrix entries.
      CALL stdop_assembleSimpleMatrix(Rlevels(i)%rmatrixB1, DER_FUNC, DER_DERIV_X,&
                                      -1.0_DP, .TRUE.)

      ! Build the B2-matrix entries.
      CALL stdop_assembleSimpleMatrix(Rlevels(i)%rmatrixB2, DER_FUNC, DER_DERIV_Y,&
                                      -1.0_DP, .TRUE.)
 
      ! Build the M1-matrix entries.
      CALL stdop_assembleSimpleMatrix(Rlevels(i)%rmatrixM, DER_FUNC, DER_FUNC,&
                                      1.0_DP, .TRUE.)
                                      
      ! Build the N-Matrix.
      CALL stdop_assembleLaplaceMatrix(Rlevels(i)%rmatrixN, .TRUE., deta)

      ! Put the A-matrices into the block matrix
      CALL lsyssc_duplicateMatrix (Rlevels(i)%rmatrixA, &
          Rlevels(i)%rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      CALL lsyssc_duplicateMatrix (Rlevels(i)%rmatrixA, &
          Rlevels(i)%rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

      ! Put the B1/B2-matrices into the block matrix
      CALL lsyssc_duplicateMatrix (Rlevels(i)%rmatrixB1, &
          Rlevels(i)%rmatrix%RmatrixBlock(1,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      CALL lsyssc_duplicateMatrix (Rlevels(i)%rmatrixB2, &
          Rlevels(i)%rmatrix%RmatrixBlock(2,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      
      ! Furthermore, put B1^T and B2^T to the block matrix.
      CALL lsyssc_transposeMatrix (Rlevels(i)%rmatrixB1, &
          Rlevels(i)%rmatrix%RmatrixBlock(3,1),LSYSSC_TR_VIRTUAL)

      CALL lsyssc_transposeMatrix (Rlevels(i)%rmatrixB2, &
          Rlevels(i)%rmatrix%RmatrixBlock(3,2),LSYSSC_TR_VIRTUAL)
      
      ! At this point, we have a little problem:
      ! We have assembled the M-matrix, and now we need to "copy" it into our
      ! global block-matrix. But now our block-matrix submatrices M1 and M2
      ! are defined as:
      ! M1 := dgrav1 * M
      ! M2 := dgrav2 * M
      !
      ! We could simply copy the data from M to M1 and M2 and afterwards scale
      ! the M1 and M2 matrices:
      !
      !  ! Put the M1/M2-matrices into the block matrix
      !  CALL lsyssc_duplicateMatrix (Rlevels(i)%rmatrixM, &
      !      Rlevels(i)%rmatrix%RmatrixBlock(1,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      !  CALL lsyssc_duplicateMatrix (Rlevels(i)%rmatrixM, &
      !      Rlevels(i)%rmatrix%RmatrixBlock(2,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      !
      !  ! Scale the M1/M2 matrices
      !  CALL lsyssc_scaleMatrix(Rlevels(i)%rmatrix%RmatrixBlock(1,4),dgrav1)
      !  CALL lsyssc_scaleMatrix(Rlevels(i)%rmatrix%RmatrixBlock(2,4),dgrav2)
      !
      ! But fortunately, the Boussinesq-Vanka supports scaled matrices, so we
      ! will give the M1 and M2 matrices a reference to the data array of M and
      ! simply adjust the scaling factors of M1 and M2:

      ! Put the M1/M2-matrices into the block matrix
      CALL lsyssc_duplicateMatrix (Rlevels(i)%rmatrixM, &
          Rlevels(i)%rmatrix%RmatrixBlock(1,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      CALL lsyssc_duplicateMatrix (Rlevels(i)%rmatrixM, &
          Rlevels(i)%rmatrix%RmatrixBlock(2,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      
      ! Change the scaling factors of the M1/M2 matrices.
      Rlevels(i)%rmatrix%RmatrixBlock(1,4)%dscaleFactor = dgrav1
      Rlevels(i)%rmatrix%RmatrixBlock(2,4)%dscaleFactor = dgrav2
      
      ! Copy the N-matrix into the block matrix
      CALL lsyssc_duplicateMatrix (Rlevels(i)%rmatrixN, &
          Rlevels(i)%rmatrix%RmatrixBlock(4,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      
      ! Update the structural information of the block matrix, as we manually
      ! changed the submatrices:
      CALL lsysbl_updateMatStrucInfo (Rlevels(i)%rmatrix)
      
      ! Create a solution vector on this level
      CALL lsysbl_createVecBlockIndMat (Rlevels(i)%rmatrix,Rlevels(i)%rvecSol,.FALSE.)
    
    END DO

    ! Although we could manually create the solution/RHS vector,
    ! the easiest way to set up the vector structure is
    ! to create it by using our matrix as template:
    CALL lsysbl_createVecBlockIndMat (Rlevels(NLMAX)%rmatrix,rrhs,.FALSE.)
    CALL lsysbl_createVecBlockIndMat (Rlevels(NLMAX)%rmatrix,rdef,.FALSE.)
    CALL lsysbl_createVecBlockIndMat (Rlevels(NLMAX)%rmatrix,rtmp,.FALSE.)

    ! Clear both RHS vectors and the solution vectors on the finest level.
    CALL lsysbl_clearVector(rrhs)
    CALL lsysbl_clearVector(Rlevels(NLMAX)%rvecSol)

    DO i = NLMIN, NLMAX
    
      ! Now implement the boundary conditions.
      CALL bcasm_initDiscreteBC(Rlevels(i)%rdiscreteBC)

      ! Get the bottom edge - Dirichlet(0) BCs for X- and Y-velocities, but
      ! Neumann for the temperature.
      CALL boundary_createRegion(rboundary,1,1,rboundaryRegion)
      
      ! X-/Y-velocity
      CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
          rboundaryRegion, Rlevels(i)%rdiscreteBC, getBoundaryValues)
      CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,2,&
          rboundaryRegion, Rlevels(i)%rdiscreteBC, getBoundaryValues)
      
      ! Get the left edge - Dirichlet(0) BC for the X- and Y-velocities and
      ! Dirichlet(-1/2) for the temperature.
      CALL boundary_createRegion(rboundary,1,2,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
      
      ! X-/Y-velocity
      CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
          rboundaryRegion, Rlevels(i)%rdiscreteBC, getBoundaryValues)
      CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,2,&
          rboundaryRegion, Rlevels(i)%rdiscreteBC, getBoundaryValues)
                               
      ! temperature
      CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,4,&
          rboundaryRegion, Rlevels(i)%rdiscreteBC, getBoundaryValues)

      ! Get the top edge - Dirichlet(0) BCs for X- and Y-velocities, but
      ! Neumann for the temperature.
      CALL boundary_createRegion(rboundary,1,3,rboundaryRegion)

      ! X-/Y-velocity
      CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
          rboundaryRegion, Rlevels(i)%rdiscreteBC, getBoundaryValues)
      CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,2,&
          rboundaryRegion, Rlevels(i)%rdiscreteBC, getBoundaryValues)
      
      ! Get the right edge - Dirichlet(0) BC for the X- and Y-velocities and
      ! Dirichlet(1/2) for the temperature.
      CALL boundary_createRegion(rboundary,1,4,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND

      ! X-/Y-velocity
      CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
          rboundaryRegion, Rlevels(i)%rdiscreteBC, getBoundaryValues)
      CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,2,&
          rboundaryRegion, Rlevels(i)%rdiscreteBC, getBoundaryValues)

      ! temperature
      CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,4,&
          rboundaryRegion, Rlevels(i)%rdiscreteBC, getBoundaryValues)

      ! Hang in the discrete BCs into the system matrices and solution vectors.
      Rlevels(i)%rmatrix%p_rdiscreteBC => Rlevels(i)%rdiscreteBC
      Rlevels(i)%rvecSol%p_rdiscreteBC => Rlevels(i)%rdiscreteBC
      
    END DO

    ! Also implement the discrete boundary conditions on the finest level
    ! onto our right-hand-side and solution vectors.
    rrhs%p_rdiscreteBC => Rlevels(NLMAX)%rdiscreteBC
    CALL vecfil_discreteBCrhs (rrhs)
    CALL vecfil_discreteBCsol (Rlevels(NLMAX)%rvecSol)

    ! During the linear solver, the boundary conditions must
    ! frequently be imposed to the vectors. This is done using
    ! a filter chain. As the linear solver does not work with 
    ! the actual solution vectors but with defect vectors instead,
    ! a filter for implementing the real boundary conditions 
    ! would be wrong.
    ! Therefore, create a filter chain with one filter only,
    ! which implements Dirichlet-conditions into a defect vector.
    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

    ! Now we have to build up the level information for multigrid.
    !
    ! At first, initialise a standard interlevel projection structure. We
    ! can use the same structure for all levels.
    CALL mlprj_initProjectionMat (rprojection,Rlevels(NLMAX)%rmatrix)

    ! Create a Multigrid-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    p_RfilterChain => RfilterChain
    CALL linsol_initMultigrid (p_rsolver,p_RfilterChain)
    
    ! Set up a BiCGStab solver with VANKA preconditioning as coarse grid solver:
    !CALL linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_BOUSS2D_DIAG)
    !CALL linsol_initBiCGStab (p_rcoarseGridSolver,p_rpreconditioner,p_RfilterChain)
    CALL linsol_initUMFPACK4(p_rcoarseGridSolver)
    
    ! Add the coarse grid level.
    CALL linsol_addMultiGridLevel(p_rlevelInfo,p_rsolver,rprojection,&
                                  NULL(), NULL(), p_rcoarseGridSolver)

    ! Now set up the other levels...
    DO i = NLMIN+1, NLMAX
    
      ! Set up the VANKA smoother.
      CALL linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_BOUSS2D_DIAG)
      
      ! We will use 4 smoothing steps with damping parameter 1.0
      CALL linsol_convertToSmoother(p_rsmoother, 4, 1.0_DP)
      
      ! And add this multi-grid level. We will use the same smoother
      ! for pre- and post-smoothing.
      CALL linsol_addMultiGridLevel(p_rlevelInfo,p_rsolver,rprojection,&
                                    p_rsmoother, p_rsmoother, NULL())
      
    END DO
    
    ! Set the output level of the solver to 2 for some output
    p_rsolver%ioutputLevel = -1
    p_rsolver%iresNorm = 0
    p_rsolver%depsRel = 1E-5_DP
    p_rsolver%depsAbs = 0.0_DP !1E-9_DP
    p_rsolver%istoppingCriterion = LINSOL_STOP_ONEOF
    p_rsolver%nminIterations = 5
    p_rsolver%nmaxIterations = 5

    ! Attach the system matrix to the solver.
    !
    ! We copy our matrices to a big matrix array and transfer that
    ! to the setMatrices routines. This intitialises then the matrices
    ! on all levels according to that array. Note that this does not
    ! allocate new memory, we create only 'links' to existing matrices
    ! into Rmatrices(:)!
    ALLOCATE(Rmatrices(NLMIN:NLMAX))
    DO i = NLMIN, NLMAX
      CALL lsysbl_duplicateMatrix (Rlevels(i)%rmatrix,&
          Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    END DO
    
    CALL linsol_setMatrices(p_Rsolver,Rmatrices(NLMIN:NLMAX))

    ! We can release Rmatrices immediately -- as long as we don't
    ! release Rlevels(i)%rmatrix!
    DO i=NLMIN,NLMAX
      CALL lsysbl_releaseMatrix (Rmatrices(i))
    END DO
    DEALLOCATE(Rmatrices)
    
    ! Initialise only the structure of the solvers - the data is initialised
    ! in the non-linear loop below.
    CALL linsol_initStructure (p_rsolver, ierror)
    IF (ierror .NE. LINSOL_ERR_NOERROR) STOP

    ! Now prepare the trilinear forms for the convection in the secondary system.
    ! u1 * d_x * T
    rtriform1%itermCount = 1
    rtriform1%ballCoeffConstant = .TRUE.
    rtriform1%BconstantCoeff = .TRUE.
    rtriform1%Dcoefficients = 1.0_DP
    rtriform1%Idescriptors(1,1) = DER_FUNC
    rtriform1%Idescriptors(2,1) = DER_DERIV_X
    rtriform1%Idescriptors(3,1) = DER_FUNC
    ! u2 * d_y * T
    rtriform2%itermCount = 1
    rtriform2%ballCoeffConstant = .TRUE.
    rtriform2%BconstantCoeff = .TRUE.
    rtriform2%Dcoefficients = 1.0_DP
    rtriform2%Idescriptors(1,1) = DER_FUNC
    rtriform2%Idescriptors(2,1) = DER_DERIV_Y
    rtriform2%Idescriptors(3,1) = DER_FUNC
    
    ! And prepare the streamline diffusion object for the convection in the
    ! primary system.
    rsd%dnu = dnu

    ! Start the timer
    CALL stat_clearTimer(rtimer)
    CALL stat_startTimer(rtimer)

    ! Now comes the fun part - the non-linear loop.
    DO NLiter = 0, 50
    
      ! The very first thing is we need to update the non-linear parts of
      ! the system matrices, but to do this, we need the solution vectors
      ! on all levels - so let's restrict the solution vectors.
      DO i = NLMAX, NLMIN+1, -1
        
        ! Restrict the solution vector
        CALL mlprj_performInterpolation(rprojection, Rlevels(i-1)%rvecSol,&
                         Rlevels(i)%rvecSol, rtmp%RvectorBlock(1))

        
        ! And filter the restricted vector.
        CALL vecfil_discreteBCsol (Rlevels(i-1)%rvecSol)

      END DO
      
      ! Now we can update the matrices - so go for it.
      DO i = NLMIN, NLMAX
      
        ! Copy A-matrix into the X-velocity block of the system.
        CALL lsyssc_duplicateMatrix(Rlevels(i)%rmatrixA,&
            Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
            LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPY)
        
        ! Assemble the convection of the system using streamline diffusion.
        CALL conv_streamlineDiffusion2d(Rlevels(i)%rvecSol, &
            Rlevels(i)%rvecSol, 1.0_DP, 0.0_DP, rsd, CONV_MODMATRIX, &
            Rlevels(i)%rmatrix%RmatrixBlock(1,1))
        
        ! The Y-velocity block is automatically updated, as it is just a shared
        ! copy of the X-velocity block.
        
        ! Copy N-matrix into the system block matrix.
        CALL lsyssc_duplicateMatrix(Rlevels(i)%rmatrixN,&
            Rlevels(i)%rmatrix%RmatrixBlock(4,4),&
            LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPY)

        ! Assemble the convection of the secondary system.
        ! Assemble u_1 * d_x T
        CALL trilf_buildMatrixScalar(rtriform1,.FALSE.,&
             Rlevels(i)%rmatrix%RmatrixBlock(4,4), &
             Rlevels(i)%rvecSol%RvectorBlock(1))
        ! Assemble u_2 * d_y T
        CALL trilf_buildMatrixScalar(rtriform2,.FALSE.,&
             Rlevels(i)%rmatrix%RmatrixBlock(4,4), &
             Rlevels(i)%rvecSol%RvectorBlock(2))
        
        ! And filter the system matrix.
        CALL matfil_discreteBC (Rlevels(i)%rmatrix)

      END DO
      
      ! Calculate the non-linear defect.
      CALL lsysbl_copyVector(rrhs, rdef)
      CALL lsysbl_blockMatVec(Rlevels(NLMAX)%rmatrix, &
          Rlevels(NLMAX)%rvecSol, rdef, -1.0_DP, 1.0_DP)
      
      ! Filter the defect vector.
      CALL vecfil_discreteBCdef(rdef)
      
      ! Calculate the residual.
      dnlRes = lsysbl_vectorNorm(rdef, LINALG_NORMEUCLID)
      
      ! Print the residual
      !CALL output_separator(OU_SEP_MINUS)
      CALL output_line('NL-Iteration: ' // TRIM(sys_si(NLiter,2)) // &
                       '   |RES| = ' // TRIM(sys_sdEP(dnlRes,20,12)))
      !CALL output_separator(OU_SEP_MINUS)
      
      ! Is this precise enough?
      IF(dnlRes .LE. 1E-8_DP) EXIT
      
      ! Initialise the data of the solver.
      CALL linsol_initData(p_rsolver, ierror)
      IF(ierror .NE. LINSOL_ERR_NOERROR) STOP
      
      ! And call the solver
      CALL linsol_precondDefect(p_rsolver, rdef)
      bError = (p_rsolver%iresult .NE. 0)
      
      ! Release the solver data.
      CALL linsol_doneData(p_rsolver)
      
      ! Did the solver break down?
      IF(bError) THEN

        ! Print an error
        CALL output_separator(OU_SEP_STAR)
        CALL output_line('NL-Iteration: ERROR: linear solver broke down')
        CALL output_separator(OU_SEP_STAR)
        
        ! Exit the non-linear loop
        EXIT
      
      END IF
      
      ! Now update the solution vector.
      CALL lsysbl_vectorLinearComb(rdef, Rlevels(NLMAX)%rvecSol, dnlDamping, 1.0_DP)

      ! That's it for this non-linear iteration...
          
    END DO ! NLiter
    
    ! Stop the timer and print the elapsed time
    CALL stat_stopTimer(rtimer)
    CALL output_line('Total time for non-linear loop: ' // stat_sgetTime(rtimer))

    ! We can now start the postprocessing. 
    ! Start UCD export to GMV file:
    CALL ucd_startVTK (rexport,UCD_FLAG_STANDARD,&
        Rlevels(NLMAX)%rtriangulation,'gmv/u_mit_1.vtk')
    
    ! Project X- and Y-velocity to the vertices of the mesh.
    NULLIFY(p_Ddata)
    NULLIFY(p_Ddata2)
    CALL spdp_projectToVertices(Rlevels(NLMAX)%rvecSol%RvectorBlock(1), p_Ddata)
    CALL spdp_projectToVertices(Rlevels(NLMAX)%rvecSol%RvectorBlock(2), p_Ddata2)

    ! In case we use the VTK exporter, which supports vector output, we will
    ! pass the X- and Y-velocity at once to the ucd module.
    CALL ucd_addVarVertBasedVec(rexport,'velocity',p_Ddata,p_Ddata2)
    
    ! Release the temporary memory.
    DEALLOCATE(p_Ddata)
    DEALLOCATE(p_Ddata2)
        
    ! Write pressure
    CALL lsyssc_getbase_double (Rlevels(NLMAX)%rvecSol%RvectorBlock(3),p_Ddata)
    CALL ucd_addVariableElementBased (rexport,'pressure',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write temperature
    CALL lsyssc_getbase_double (Rlevels(NLMAX)%rvecSol%RvectorBlock(4),p_Ddata)
    CALL ucd_addVariableVertexBased (rexport,'temperature',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that's it.
    CALL ucd_write (rexport)
    CALL ucd_release (rexport)

    ! We are finished - but not completely!
    ! Now, clean up so that all the memory is available again.
    !
    ! Release solver data and structure
    CALL linsol_doneStructure (p_rsolver)
    
    ! Release the interlevel projection structure
    CALL mlprj_doneProjection (rprojection)

    ! Release the solver node and all subnodes attached to it (if at all):
    CALL linsol_releaseSolver (p_rsolver)
    
    ! Release the block matrix/vectors
    CALL lsysbl_releaseVector (rtmp)
    CALL lsysbl_releaseVector (rdef)
    CALL lsysbl_releaseVector (rrhs)
    DO i = NLMAX, NLMIN, -1
      CALL lsysbl_releaseVector (Rlevels(i)%rvecSol)
      CALL lsysbl_releaseMatrix (Rlevels(i)%rmatrix)
    END DO
    
    ! Release B1 and B2 matrix
    DO i = NLMAX, NLMIN, -1
      CALL lsyssc_releaseMatrix (Rlevels(i)%rmatrixM)
      CALL lsyssc_releaseMatrix (Rlevels(i)%rmatrixN)
      CALL lsyssc_releaseMatrix (Rlevels(i)%rmatrixB2)
      CALL lsyssc_releaseMatrix (Rlevels(i)%rmatrixB1)
      CALL lsyssc_releaseMatrix (Rlevels(i)%rmatrixA)
    END DO
    
    ! Release our discrete version of the boundary conditions
    DO i = NLMAX, NLMIN, -1
      CALL bcasm_releaseDiscreteBC (Rlevels(i)%rdiscreteBC)
    END DO

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    DO i = NLMAX, NLMIN, -1
      CALL spdiscr_releaseBlockDiscr(Rlevels(i)%rdiscretisation)
    END DO
    
    ! Release the triangulation. 
    DO i = NLMAX, NLMIN, -1
      CALL tria_done (Rlevels(i)%rtriangulation)
    END DO
    
    DEALLOCATE(Rlevels)
    
    ! Finally release the domain, that's it.
    CALL boundary_release (rboundary)

  END SUBROUTINE

END MODULE