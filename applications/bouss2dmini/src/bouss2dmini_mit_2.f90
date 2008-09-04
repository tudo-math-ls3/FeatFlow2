!##############################################################################
!# ****************************************************************************
!# <name> bouss2dmini_mit_2 </name>
!# ****************************************************************************
!#
!# <purpose>
!# Todo
!#
!# System:
!#
!# / A(u) 0    B1    M1   \   / u_1 \
!# | 0    A(u) B2    M2   | * | u_2 | = 0
!# | B1^T B2^T 0     0    |   |  p  |
!# \ 0    0    0     N(u) /   \  T  /
!#
!# u = 0 on the complete boundary
!# T =  1/2 on the left boundary region
!# T = -1/2 on the right boundary region
!# </purpose>
!##############################################################################

MODULE bouss2dmini_mit_2

  USE fsystem
  USE storage
  USE linearsolver
  USE boundary
  USE trilinearformevaluation
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

    ! A block discretisation for the primary system (Q1~/Q1~/Q0)
    TYPE(t_blockDiscretisation) :: rdiscrPrimary
    
    ! A block discretisation for the secondary system (Q1)
    TYPE(t_blockDiscretisation) :: rdiscrSecondary
    
    ! A primary system matrix which recieves the Navier-Stokes system
    ! for that specific level.
    TYPE(t_matrixBlock) :: rmatrixPrimary
    
    ! A secondary system matrix which recieves the Convection-Diffusion system
    ! of the temperature for that specific level.
    TYPE(t_matrixBlock) :: rmatrixSecondary
    
    ! Laplace-matrix of the primary system for that specific level.
    TYPE(t_matrixScalar) :: rmatrixA

    ! B1-matrix for that specific level. 
    TYPE(t_matrixScalar) :: rmatrixB1

    ! B2-matrix for that specific level. 
    TYPE(t_matrixScalar) :: rmatrixB2
    
    ! M1- / M2-matrix for that specific level.
    TYPE(t_matrixScalar) :: rmatrixM
    
    ! Laplace-Matrix of the secondary system for that specific level.
    TYPE(t_matrixScalar) :: rmatrixN
    
    ! Primary solution vector for that specific level.
    TYPE(t_vectorBlock) :: rvecSolPrimary
    
    ! Secondary solution vector for that specific level.
    TYPE(t_vectorBlock) :: rvecSolSecondary

    ! A variable describing the discrete boundary conditions for the
    ! primary system.
    TYPE(t_discreteBC) :: rdiscreteBCprimary

    ! A variable describing the discrete boundary conditions for the
    ! secondary system.
    TYPE(t_discreteBC) :: rdiscreteBCsecondary
  
  END TYPE
  
!</typeblock>

!</types>

CONTAINS
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE b2dm_mit_2
  
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
    TYPE(t_vectorBlock) :: rrhsPrimary,rrhsSecondary,&
        rvecTmpPrimary,rvecTmpSecondary,rvecDefPrimary,rvecDefSecondary
    
    ! A set of variables describing the analytic and discrete boundary
    ! conditions.    
    TYPE(t_boundaryRegion) :: rboundaryRegion

    ! A solver node that accepts parameters for the linear solver    
    TYPE(t_linsolNode), POINTER :: p_rsolverPrimary,p_rsolverSecondary,&
        p_rcoarseGridSolver,p_rsmoother,p_rpreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    TYPE(t_matrixBlock), DIMENSION(:), ALLOCATABLE :: Rmatrices

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    TYPE(t_filterChain), DIMENSION(1), TARGET :: RfilterChain
    TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
    
    ! An interlevel projection structure for changing levels
    TYPE(t_interlevelProjectionBlock) :: rprojPrimary, rprojSecondary

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
    REAL(DP) :: dnlDampPrimary, dnlDampSecondary
    
    ! Error indicator during initialisation of the solver
    INTEGER :: ierror
    
    ! A counter variable for the non-linear loop
    INTEGER :: NLiter
    
    ! Non-linear residuals
    REAL(DP) :: dnlResPrimary, dnlResSecondary, dnlResTotal !, dnlResScale
    
    ! Output block for UCD output to GMV file
    TYPE(t_ucdExport) :: rexport
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_Ddata2

    ! A counter variable
    INTEGER :: i
    LOGICAL :: bError
    
    ! A timer object for time measuring
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
    
    ! Set the damping parameters for the systems.
    dnlDampPrimary = 1.0_DP
    dnlDampSecondary = 1.0_DP

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

    ! Now we can start to initialise the discretisation.
    DO i = NLMIN, NLMAX
    
      ! Our primary discretisation has 3 components (Q1~/Q1~/Q0)
      CALL spdiscr_initBlockDiscr2D (Rlevels(i)%rdiscrPrimary, 3, &
                                     Rlevels(i)%rtriangulation, rboundary)
      
      ! For simplicity, we set up one discretisation structure for the 
      ! velocity...
      CALL spdiscr_initDiscr_simple (Rlevels(i)%rdiscrPrimary%RspatialDiscr(1),&
          EL_EM30, CUB_G2X2, Rlevels(i)%rtriangulation, rboundary)
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory, 
      ! as both structures will share the same dynamic information afterwards.
      CALL spdiscr_duplicateDiscrSc (&
          Rlevels(i)%rdiscrPrimary%RspatialDiscr(1),&
          Rlevels(i)%rdiscrPrimary%RspatialDiscr(2))

      ! For the pressure (3rd component), we set up a separate discretisation 
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      CALL spdiscr_deriveSimpleDiscrSc (Rlevels(i)%rdiscrPrimary%RspatialDiscr(1), &
          EL_Q0, CUB_G2X2, Rlevels(i)%rdiscrPrimary%RspatialDiscr(3))

      ! Our secondary discretisation has 1 component (Q1)
      CALL spdiscr_initBlockDiscr2D (Rlevels(i)%rdiscrSecondary, 1, &
                                     Rlevels(i)%rtriangulation, rboundary)

      CALL spdiscr_initDiscr_simple (Rlevels(i)%rdiscrSecondary%RspatialDiscr(1),&
          EL_Q1, CUB_G2X2, Rlevels(i)%rtriangulation, rboundary)
    
    END DO

    DO i = NLMIN, NLMAX
    
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! Initialise primary system for this level
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    
      ! Initialise a block matrix based on the primary discretisation.
      CALL lsysbl_createMatBlockByDiscr (Rlevels(i)%rdiscrPrimary,&
                                         Rlevels(i)%rmatrixPrimary)    
      
      ! Inform the matrix that we build a saddle-point problem.
      Rlevels(i)%rmatrixPrimary%imatrixSpec = LSYSBS_MSPEC_SADDLEPOINT
      
      ! Create the matrix structure of the A-matrices of the primary system.
      CALL bilf_createMatrixStructure (Rlevels(i)%rdiscrPrimary%RspatialDiscr(1),&
          LSYSSC_MATRIX9, Rlevels(i)%rmatrixA)

      ! Create the matrix structure of the B-matrices
      CALL bilf_createMatrixStructure (Rlevels(i)%rdiscrPrimary%RspatialDiscr(3),&
          LSYSSC_MATRIX9, Rlevels(i)%rmatrixB1,&
          Rlevels(i)%rdiscrPrimary%RspatialDiscr(1))

      ! Since B2 has the same structure as B1, we'll simply get a copy of it.
      CALL lsyssc_duplicateMatrix (Rlevels(i)%rmatrixB1, Rlevels(i)%rmatrixB2,&
                                   LSYSSC_DUP_SHARE, LSYSSC_DUP_REMOVE)
      
      ! Build the Laplace-Matrix of the primary system.
      CALL stdop_assembleLaplaceMatrix(Rlevels(i)%rmatrixA, .TRUE., dnu)
      
      ! Build the B1-matrix entries
      CALL stdop_assembleSimpleMatrix(Rlevels(i)%rmatrixB1, DER_FUNC, DER_DERIV_X,&
                                      -1.0_DP, .TRUE.)

      ! Build the B2-matrix entries
      CALL stdop_assembleSimpleMatrix(Rlevels(i)%rmatrixB2, DER_FUNC, DER_DERIV_Y,&
                                      -1.0_DP, .TRUE.)
      
      ! Now let's set up our primary system block matrix.
      
      ! Copy the A-matrix into our primary system block matrix.
      CALL lsyssc_duplicateMatrix (Rlevels(i)%rmatrixA, &
          Rlevels(i)%rmatrixPrimary%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      CALL lsyssc_duplicateMatrix (Rlevels(i)%rmatrixPrimary%RmatrixBlock(1,1), &
          Rlevels(i)%rmatrixPrimary%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      
      ! Copy the B1-/B2-matrices into our primary system block matrix.
      CALL lsyssc_duplicateMatrix (Rlevels(i)%rmatrixB1, &
          Rlevels(i)%rmatrixPrimary%RmatrixBlock(1,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      CALL lsyssc_duplicateMatrix (Rlevels(i)%rmatrixB2, &
          Rlevels(i)%rmatrixPrimary%RmatrixBlock(2,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      
      ! Furthermore, put B1^T and B2^T to the block matrix.
      CALL lsyssc_transposeMatrix (Rlevels(i)%rmatrixB1, &
          Rlevels(i)%rmatrixPrimary%RmatrixBlock(3,1),LSYSSC_TR_VIRTUAL)

      CALL lsyssc_transposeMatrix (Rlevels(i)%rmatrixB2, &
          Rlevels(i)%rmatrixPrimary%RmatrixBlock(3,2),LSYSSC_TR_VIRTUAL)

      ! Update the structural information of the block matrix.
      CALL lsysbl_updateMatStrucInfo (Rlevels(i)%rmatrixPrimary)

      ! Finally, create a solution vector based on the matrix.
      CALL lsysbl_createVecBlockIndMat (Rlevels(i)%rmatrixPrimary,&
                                        Rlevels(i)%rvecSolPrimary,.FALSE.)

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! Initialise secondary system for this level
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      ! Initialise a block matrix based on the secondary discretisation.
      CALL lsysbl_createMatBlockByDiscr (Rlevels(i)%rdiscrSecondary,&
                                         Rlevels(i)%rmatrixSecondary)

      ! Create the matrix structure of the secondary system.
      CALL bilf_createMatrixStructure (Rlevels(i)%rdiscrSecondary%RspatialDiscr(1),&
          LSYSSC_MATRIX9, Rlevels(i)%rmatrixN)

      ! Build the Laplace-Matrix of the secondary system.
      CALL stdop_assembleLaplaceMatrix(Rlevels(i)%rmatrixN, .TRUE., deta)

      ! Copy the N-matrix into our secondary system block matrix.
      CALL lsyssc_duplicateMatrix (Rlevels(i)%rmatrixN, &
          Rlevels(i)%rmatrixSecondary%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

      ! Update the structural information of the block matrix.
      CALL lsysbl_updateMatStrucInfo (Rlevels(i)%rmatrixSecondary)

      ! Finally, create a solution vector based on the matrix.
      CALL lsysbl_createVecBlockIndMat (Rlevels(i)%rmatrixSecondary,&
                                        Rlevels(i)%rvecSolSecondary,.FALSE.)

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! Initialise coupling matrices for this level
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      ! Create the matrix structure of the M-matrices
      CALL bilf_createMatrixStructure (Rlevels(i)%rdiscrSecondary%RspatialDiscr(1),&
          LSYSSC_MATRIX9, Rlevels(i)%rmatrixM,&
          Rlevels(i)%rdiscrPrimary%RspatialDiscr(1))

      ! Build a mass matrix.
      CALL stdop_assembleSimpleMatrix(Rlevels(i)%rmatrixM, DER_FUNC, DER_FUNC,&
                                      1.0_DP, .TRUE.)
      
      ! That's it for now.
    
    END DO

    ! Create the RHS vectors on the finest level.
    CALL lsysbl_createVecBlockIndMat (Rlevels(NLMAX)%rmatrixPrimary,&
                                      rrhsPrimary, .FALSE.)
    CALL lsysbl_createVecBlockIndMat (Rlevels(NLMAX)%rmatrixPrimary,&
                                      rvecTmpPrimary, .FALSE.)
    CALL lsysbl_createVecBlockIndMat (Rlevels(NLMAX)%rmatrixPrimary,&
                                      rvecDefPrimary, .FALSE.)
    CALL lsysbl_createVecBlockIndMat (Rlevels(NLMAX)%rmatrixSecondary,&
                                      rrhsSecondary, .FALSE.)
    CALL lsysbl_createVecBlockIndMat (Rlevels(NLMAX)%rmatrixSecondary,&
                                      rvecTmpSecondary, .FALSE.)
    CALL lsysbl_createVecBlockIndMat (Rlevels(NLMAX)%rmatrixSecondary,&
                                      rvecDefSecondary, .FALSE.)

    ! Clear both RHS vectors and the solution vectors on the finest level.
    CALL lsysbl_clearVector(rrhsPrimary)
    CALL lsysbl_clearVector(rrhsSecondary)
    CALL lsysbl_clearVector(Rlevels(NLMAX)%rvecSolPrimary)
    CALL lsysbl_clearVector(Rlevels(NLMAX)%rvecSolSecondary)

    DO i = NLMIN, NLMAX
    
      ! Now implement the boundary conditions.
      CALL bcasm_initDiscreteBC(Rlevels(i)%rdiscreteBCprimary)
      CALL bcasm_initDiscreteBC(Rlevels(i)%rdiscreteBCsecondary)

      ! Get the bottom edge - Dirichlet(0) BCs for X- and Y-velocities, but
      ! Neumann for the temperature.
      CALL boundary_createRegion(rboundary,1,1,rboundaryRegion)
      
      ! X-/Y-velocity
      CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscrPrimary,1,&
          rboundaryRegion, Rlevels(i)%rdiscreteBCprimary, getBoundaryValuesPrimary)
      CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscrPrimary,2,&
          rboundaryRegion, Rlevels(i)%rdiscreteBCprimary, getBoundaryValuesPrimary)
      
      ! Get the left edge - Dirichlet(0) BC for the X- and Y-velocities and
      ! Dirichlet(-1/2) for the temperature.
      CALL boundary_createRegion(rboundary,1,2,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
      
      ! X-/Y-velocity
      CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscrPrimary,1,&
          rboundaryRegion, Rlevels(i)%rdiscreteBCprimary, getBoundaryValuesPrimary)
      CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscrPrimary,2,&
          rboundaryRegion, Rlevels(i)%rdiscreteBCprimary, getBoundaryValuesPrimary)
                               
      ! temperature
      CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscrSecondary,1,&
          rboundaryRegion, Rlevels(i)%rdiscreteBCsecondary, getBoundaryValuesSecondary)

      ! Get the top edge - Dirichlet(0) BCs for X- and Y-velocities, but
      ! Neumann for the temperature.
      CALL boundary_createRegion(rboundary,1,3,rboundaryRegion)

      ! X-/Y-velocity
      CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscrPrimary,1,&
          rboundaryRegion, Rlevels(i)%rdiscreteBCprimary, getBoundaryValuesPrimary)
      CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscrPrimary,2,&
          rboundaryRegion, Rlevels(i)%rdiscreteBCprimary, getBoundaryValuesPrimary)
      
      ! Get the right edge - Dirichlet(0) BC for the X- and Y-velocities and
      ! Dirichlet(1/2) for the temperature.
      CALL boundary_createRegion(rboundary,1,4,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND

      ! X-/Y-velocity
      CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscrPrimary,1,&
          rboundaryRegion, Rlevels(i)%rdiscreteBCprimary, getBoundaryValuesPrimary)
      CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscrPrimary,2,&
          rboundaryRegion, Rlevels(i)%rdiscreteBCprimary, getBoundaryValuesPrimary)

      ! temperature
      CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscrSecondary,1,&
          rboundaryRegion, Rlevels(i)%rdiscreteBCsecondary, getBoundaryValuesSecondary)

      ! Hang in the discrete BCs into the system matrices and solution vectors.
      Rlevels(i)%rmatrixPrimary%p_rdiscreteBC => Rlevels(i)%rdiscreteBCprimary
      Rlevels(i)%rmatrixSecondary%p_rdiscreteBC => Rlevels(i)%rdiscreteBCsecondary
      Rlevels(i)%rvecSolPrimary%p_rdiscreteBC => Rlevels(i)%rdiscreteBCprimary
      Rlevels(i)%rvecSolSecondary%p_rdiscreteBC => Rlevels(i)%rdiscreteBCsecondary
      
    END DO

    ! Hang in the BCs into the RHS vectors.
    rrhsPrimary%p_rdiscreteBC => Rlevels(NLMAX)%rdiscreteBCprimary
    rvecDefPrimary%p_rdiscreteBC => Rlevels(NLMAX)%rdiscreteBCprimary
    rrhsSecondary%p_rdiscreteBC => Rlevels(NLMAX)%rdiscreteBCsecondary
    rvecDefSecondary%p_rdiscreteBC => Rlevels(NLMAX)%rdiscreteBCsecondary
    
    ! And filter the RHS and solution vectors on the finest level.
    CALL vecfil_discreteBCrhs (rrhsPrimary)
    CALL vecfil_discreteBCrhs (rrhsSecondary)
    CALL vecfil_discreteBCsol (Rlevels(NLMAX)%rvecSolPrimary)
    CALL vecfil_discreteBCsol (Rlevels(NLMAX)%rvecSolSecondary)

    ! Initialise a filter chain.
    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

    ! Set up two inter-level-projection structures - one for the primary and
    ! one for the secondary system.
    CALL mlprj_initProjectionMat (rprojPrimary,Rlevels(NLMAX)%rmatrixPrimary)
    CALL mlprj_initProjectionMat (rprojSecondary,Rlevels(NLMAX)%rmatrixSecondary)

    ! Create two multi-grid solvers.
    p_RfilterChain => RfilterChain
    CALL linsol_initMultigrid (p_rsolverPrimary, p_RfilterChain)
    CALL linsol_initMultigrid (p_rsolverSecondary, p_RfilterChain)
    
    ! Set up a BiCGStab solver with VANKA preconditioner as coarse grid solver
    ! for the primary system.
    !CALL linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DNAVST)
    !CALL linsol_initBiCGStab (p_rcoarseGridSolver,p_rpreconditioner,p_RfilterChain)
    CALL linsol_initUMFPACK4(p_rcoarseGridSolver)
    CALL linsol_addMultiGridLevel(p_rlevelInfo,p_rsolverPrimary,rprojPrimary,&
                                  NULL(), NULL(), p_rcoarseGridSolver)

    ! Set up a BiCGStab solver with SSOR preconditioner as coarse grid solver
    ! for the secondary system.
    !CALL linsol_initSSOR (p_rpreconditioner)
    !CALL linsol_initBiCGStab (p_rcoarseGridSolver,p_rpreconditioner,p_RfilterChain)
    CALL linsol_initUMFPACK4(p_rcoarseGridSolver)
    CALL linsol_addMultiGridLevel(p_rlevelInfo,p_rsolverSecondary,rprojSecondary,&
                                  NULL(), NULL(), p_rcoarseGridSolver)

    ! Now set up the other levels...
    DO i = NLMIN+1, NLMAX
    
      ! Set up a VANKA smoother for the primary system.
      CALL linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DNAVST)
      CALL linsol_convertToSmoother(p_rsmoother, 4, 1.0_DP)
      CALL linsol_addMultiGridLevel(p_rlevelInfo,p_rsolverPrimary,rprojPrimary,&
                                    p_rsmoother, p_rsmoother, NULL())
      
      ! Set up a SSOR smoother for the secondary system.
      CALL linsol_initSSOR (p_rsmoother)
      CALL linsol_convertToSmoother(p_rsmoother, 4, 1.0_DP)
      CALL linsol_addMultiGridLevel(p_rlevelInfo,p_rsolverSecondary,rprojSecondary,&
                                    p_rsmoother, p_rsmoother, NULL())

    END DO
    
    ! Set the output level of the solvers to 2 for some output
    p_rsolverPrimary%ioutputLevel = -1
    p_rsolverSecondary%ioutputLevel = -1
    
    ! And set the correct tolerance parameters
    p_rsolverPrimary%iresNorm = 0
    p_rsolverPrimary%depsRel = 1E-5_DP
    p_rsolverPrimary%depsAbs = 0.0_DP !1E-9_DP
    p_rsolverPrimary%istoppingCriterion = LINSOL_STOP_ONEOF
    p_rsolverPrimary%nminIterations = 3
    p_rsolverPrimary%nmaxIterations = 5

    p_rsolverSecondary%iresNorm = 0
    p_rsolverSecondary%depsRel = 1E-5_DP
    p_rsolverSecondary%depsAbs = 0.0_DP !1E-9_DP
    p_rsolverSecondary%istoppingCriterion = LINSOL_STOP_ONEOF
    p_rsolverSecondary%nminIterations = 3
    p_rsolverSecondary%nmaxIterations = 5

    ! Attach the system matrices to the solvers.
    ALLOCATE(Rmatrices(NLMIN:NLMAX))
    ! Primary System
    DO i = NLMIN, NLMAX
      CALL lsysbl_duplicateMatrix (Rlevels(i)%rmatrixPrimary,&
          Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    END DO
    CALL linsol_setMatrices(p_RsolverPrimary,Rmatrices(NLMIN:NLMAX))
    DO i=NLMIN,NLMAX
      CALL lsysbl_releaseMatrix (Rmatrices(i))
    END DO

    ! Secondary System
    DO i = NLMIN, NLMAX
      CALL lsysbl_duplicateMatrix (Rlevels(i)%rmatrixSecondary,&
          Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    END DO
    CALL linsol_setMatrices(p_RsolverSecondary,Rmatrices(NLMIN:NLMAX))
    DO i=NLMIN,NLMAX
      CALL lsysbl_releaseMatrix (Rmatrices(i))
    END DO

    DEALLOCATE(Rmatrices)
    
    ! Initialise only the structure of the solvers - the data is initialised
    ! in the non-linear loop below.
    CALL linsol_initStructure (p_rsolverPrimary, ierror)
    IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
    CALL linsol_initStructure (p_rsolverSecondary, ierror)
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
        
        ! Restrict primary solution vector.
        CALL mlprj_performInterpolation(rprojPrimary, Rlevels(i-1)%rvecSolPrimary,&
                         Rlevels(i)%rvecSolPrimary, rvecTmpPrimary%RvectorBlock(1))

        ! Restrict secondary solution vector.
        CALL mlprj_performInterpolation(rprojSecondary, Rlevels(i-1)%rvecSolSecondary,&
                         Rlevels(i)%rvecSolSecondary, rvecTmpSecondary%RvectorBlock(1))
        
        ! And filter the restricted vectors.
        CALL vecfil_discreteBCsol (Rlevels(i-1)%rvecSolPrimary)
        CALL vecfil_discreteBCsol (Rlevels(i-1)%rvecSolSecondary)

      END DO
      
      ! Now we can update the matrices - so go for it.
      DO i = NLMIN, NLMAX
      
        ! Copy A-matrix into the X-velocity block of the primary system.
        CALL lsyssc_duplicateMatrix(Rlevels(i)%rmatrixA,&
            Rlevels(i)%rmatrixPrimary%RmatrixBlock(1,1),&
            LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPY)
        
        ! Assemble the convection of the primary system using streamline diffusion.
        CALL conv_streamlineDiffusion2d(Rlevels(i)%rvecSolPrimary, &
            Rlevels(i)%rvecSolPrimary, 1.0_DP, 0.0_DP, rsd, CONV_MODMATRIX, &
            Rlevels(i)%rmatrixPrimary%RmatrixBlock(1,1))
        
        ! And filter the primary system matrix.
        CALL matfil_discreteBC (Rlevels(i)%rmatrixPrimary)
        
        ! The Y-velocity block is automatically updated, as it is just a shared
        ! copy of the X-velocity block.
        
        ! Copy N-matrix into the secondary system block matrix.
        CALL lsyssc_duplicateMatrix(Rlevels(i)%rmatrixN,&
            Rlevels(i)%rmatrixSecondary%RmatrixBlock(1,1),&
            LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPY)

        ! Assemble the convection of the secondary system.
        ! Assemble u_1 * d_x T
        CALL trilf_buildMatrixScalar(rtriform1,.FALSE.,&
             Rlevels(i)%rmatrixSecondary%RmatrixBlock(1,1), &
             Rlevels(i)%rvecSolPrimary%RvectorBlock(1))
        ! Assemble u_2 * d_y T
        CALL trilf_buildMatrixScalar(rtriform2,.FALSE.,&
             Rlevels(i)%rmatrixSecondary%RmatrixBlock(1,1), &
             Rlevels(i)%rvecSolPrimary%RvectorBlock(2))
        
        ! And filter the secondary system matrix.
        CALL matfil_discreteBC (Rlevels(i)%rmatrixSecondary)

      END DO
      
      ! Okay, the system matrices are up-to-date. Now we need to calculate
      ! the non-linear defect. But as we do not have a global system matrix,
      ! we need to calculate the defect of both the primary and secondary
      ! systems.
      ! So let's begin with the secondary system.
      CALL lsysbl_copyVector(rrhsSecondary, rvecDefSecondary)
      CALL lsysbl_blockMatVec(Rlevels(NLMAX)%rmatrixSecondary, &
          Rlevels(NLMAX)%rvecSolSecondary, rvecDefSecondary, -1.0_DP, 1.0_DP)
      
      ! Filter the defect vector.
      CALL vecfil_discreteBCdef(rvecDefSecondary)
      
      ! Calculate the secondary residual.
      dnlResSecondary = lsysbl_vectorNorm(rvecDefSecondary, LINALG_NORMEUCLID)
      
      ! Calculate the defect of the primary system.
      CALL lsysbl_copyVector(rrhsPrimary, rvecDefPrimary)
      CALL lsysbl_blockMatVec(Rlevels(NLMAX)%rmatrixPrimary, &
          Rlevels(NLMAX)%rvecSolPrimary, rvecDefPrimary, -1.0_DP, 1.0_DP)
      
      ! Don't forget the M-matrices which couple our primary system to
      ! the secondary system.
      CALL lsyssc_scalarMatVec(Rlevels(NLMAX)%rmatrixM,&
          Rlevels(NLMAX)%rvecSolSecondary%RvectorBlock(1),&
          rvecDefPrimary%RvectorBlock(1), -dgrav1, 1.0_DP)
      CALL lsyssc_scalarMatVec(Rlevels(NLMAX)%rmatrixM,&
          Rlevels(NLMAX)%rvecSolSecondary%RvectorBlock(1),&
          rvecDefPrimary%RvectorBlock(2), -dgrav2, 1.0_DP)
      
      ! Filter the defect vector
      CALL vecfil_discreteBCdef(rvecDefPrimary)
      
      ! Calculate the primary residual.
      dnlResPrimary = lsysbl_vectorNorm(rvecDefPrimary, LINALG_NORMEUCLID)
      
      ! And calculate the total residual.
      dnlResTotal = SQRT(dnlResPrimary**2 + dnlResSecondary**2)
      
      ! Print the residual
      !CALL output_separator(OU_SEP_MINUS)
      CALL output_line('NL-Iteration: ' // TRIM(sys_si(NLiter,2)) // &
                       '   |RES| = ' // TRIM(sys_sdEP(dnlResTotal,20,12)))
      !CALL output_separator(OU_SEP_MINUS)
      
      ! Is this precise enough?
      IF(dnlResTotal .LE. 1E-8_DP) EXIT
      
      ! Okay, we need to perform another non-linear iteration.
      
      !IF(dnlResSecondary .GT. 1E-10_DP) THEN
      IF(.TRUE.) THEN
      
        ! So let's first solve the secondary system to get a new temperature.

        ! Initialise the data of the secondary solver.
        CALL linsol_initData(p_rsolverSecondary, ierror)
        IF(ierror .NE. LINSOL_ERR_NOERROR) STOP
        
        ! And call the secondary solver
        !CALL output_line('Solving secondary system...')
        !CALL output_line('---------------------------')
        CALL linsol_precondDefect(p_rsolverSecondary, rvecDefSecondary)
        bError = (p_rsolverSecondary%iresult .NE. 0)
        
        ! Release the secondary solver data.
        CALL linsol_doneData(p_rsolverSecondary)
        
        ! Did the secondary solver break down?
        IF(bError) THEN

          ! Print an error
          CALL output_separator(OU_SEP_STAR)
          CALL output_line('NL-Iteration: ERROR: secondary solver broke down')
          CALL output_separator(OU_SEP_STAR)
          
          ! Exit the non-linear loop
          EXIT
        
        END IF
        
        ! Now update the secondary solution vector - which is the temperature.
        CALL lsysbl_vectorLinearComb(rvecDefSecondary, &
            Rlevels(NLMAX)%rvecSolSecondary, dnlDampSecondary, 1.0_DP)
        
        ! Now recalculate the defect of the primary system as we have updated
        ! the temperature.
        CALL lsysbl_copyVector(rrhsPrimary, rvecDefPrimary)
        CALL lsysbl_blockMatVec(Rlevels(NLMAX)%rmatrixPrimary, &
            Rlevels(NLMAX)%rvecSolPrimary, rvecDefPrimary, -1.0_DP, 1.0_DP)
        CALL lsyssc_scalarMatVec(Rlevels(NLMAX)%rmatrixM,&
            Rlevels(NLMAX)%rvecSolSecondary%RvectorBlock(1),&
            rvecDefPrimary%RvectorBlock(1), -dgrav1, 1.0_DP)
        CALL lsyssc_scalarMatVec(Rlevels(NLMAX)%rmatrixM,&
            Rlevels(NLMAX)%rvecSolSecondary%RvectorBlock(1),&
            rvecDefPrimary%RvectorBlock(2), -dgrav2, 1.0_DP)
        
        ! Filter the defect vector
        CALL vecfil_discreteBCdef(rvecDefPrimary)

        ! Calculate the primary residual.
        dnlResPrimary = lsysbl_vectorNorm(rvecDefPrimary, LINALG_NORMEUCLID)
      
      END IF
      
      !IF (dnlResPrimary .GT. 1E-10_DP) THEN
      IF(.TRUE.) THEN

        ! Initialise the data of the primary solver.
        CALL linsol_initData(p_rsolverPrimary, ierror)
        IF(ierror .NE. LINSOL_ERR_NOERROR) STOP
        
        ! And call the primary solver
        !CALL output_line('Solving primary system...')
        !CALL output_line('---------------------------')
        CALL linsol_precondDefect(p_rsolverPrimary, rvecDefPrimary)
        bError = (p_rsolverPrimary%iresult .NE. 0)
        
        ! Release the primary solver data.
        CALL linsol_doneData(p_rsolverPrimary)
        
        ! Did the primary solver break down?
        IF(bError) THEN

          ! Print an error
          CALL output_separator(OU_SEP_STAR)
          CALL output_line('NL-Iteration: ERROR: Primary solver broke down')
          CALL output_separator(OU_SEP_STAR)
          
          ! Exit the non-linear loop
          EXIT
        
        END IF
        
        ! Now update the primary solution vector.
        CALL lsysbl_vectorLinearComb(rvecDefPrimary, &
            Rlevels(NLMAX)%rvecSolPrimary, dnlDampPrimary, 1.0_DP)
      
      END IF
      
      ! That's it for this non-linear iteration...
          
    END DO ! NLiter
    
    ! Stop the timer and print the elapsed time
    CALL stat_stopTimer(rtimer)
    CALL output_line('Total time for non-linear loop: ' // stat_sgetTime(rtimer))

    ! We can now start the postprocessing. 
    ! Start UCD export to GMV file:
    CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,&
        Rlevels(NLMAX)%rtriangulation,'gmv/u_mit_2.gmv')
    
    ! Project X- and Y-velocity to the vertices of the mesh.
    NULLIFY(p_Ddata)
    NULLIFY(p_Ddata2)
    CALL spdp_projectToVertices(Rlevels(NLMAX)%rvecSolPrimary%RvectorBlock(1), p_Ddata)
    CALL spdp_projectToVertices(Rlevels(NLMAX)%rvecSolPrimary%RvectorBlock(2), p_Ddata2)

    ! In case we use the VTK exporter, which supports vector output, we will
    ! pass the X- and Y-velocity at once to the ucd module.
    CALL ucd_addVarVertBasedVec(rexport,'velocity',p_Ddata,p_Ddata2)
    
    ! Release the temporary memory.
    DEALLOCATE(p_Ddata)
    DEALLOCATE(p_Ddata2)
        
    ! Write pressure
    CALL lsyssc_getbase_double (Rlevels(NLMAX)%rvecSolPrimary%RvectorBlock(3),p_Ddata)
    CALL ucd_addVariableElementBased (rexport,'pressure',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write temperature
    CALL lsyssc_getbase_double (Rlevels(NLMAX)%rvecSolSecondary%RvectorBlock(1),p_Ddata)
    CALL ucd_addVariableVertexBased (rexport,'temperature',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that's it.
    CALL ucd_write (rexport)
    CALL ucd_release (rexport)

    ! We are finished - but not completely!
    ! Now, clean up so that all the memory is available again.
    !
    ! Release solver data and structure
   ! CALL linsol_doneData (p_rsolverNode)
    CALL linsol_doneStructure (p_rsolverSecondary)
    CALL linsol_doneStructure (p_rsolverPrimary)
    
    ! Release the interlevel projection structure
    CALL mlprj_doneProjection (rprojSecondary)
    CALL mlprj_doneProjection (rprojPrimary)

    ! Release the solver node and all subnodes attached to it (if at all):
    CALL linsol_releaseSolver (p_rsolverSecondary)
    CALL linsol_releaseSolver (p_rsolverPrimary)
    
    ! Release the block matrix/vectors
    CALL lsysbl_releaseVector (rvecDefPrimary)
    CALL lsysbl_releaseVector (rvecDefSecondary)
    CALL lsysbl_releaseVector (rvecTmpPrimary)
    CALL lsysbl_releaseVector (rvecTmpSecondary)
    CALL lsysbl_releaseVector (rrhsPrimary)
    CALL lsysbl_releaseVector (rrhsSecondary)
    DO i = NLMAX, NLMIN, -1
      CALL lsysbl_releaseVector (Rlevels(i)%rvecSolSecondary)
      CALL lsysbl_releaseVector (Rlevels(i)%rvecSolPrimary)
      CALL lsysbl_releaseMatrix (Rlevels(i)%rmatrixSecondary)
      CALL lsysbl_releaseMatrix (Rlevels(i)%rmatrixPrimary)
    END DO
    
    ! Release scalar matrices
    DO i = NLMAX, NLMIN, -1
      CALL lsyssc_releaseMatrix (Rlevels(i)%rmatrixM)
      CALL lsyssc_releaseMatrix (Rlevels(i)%rmatrixN)
      CALL lsyssc_releaseMatrix (Rlevels(i)%rmatrixB2)
      CALL lsyssc_releaseMatrix (Rlevels(i)%rmatrixB1)
      CALL lsyssc_releaseMatrix (Rlevels(i)%rmatrixA)
    END DO
    
    ! Release our discrete version of the boundary conditions
    DO i = NLMAX, NLMIN, -1
      CALL bcasm_releaseDiscreteBC (Rlevels(i)%rdiscreteBCsecondary)
      CALL bcasm_releaseDiscreteBC (Rlevels(i)%rdiscreteBCprimary)
    END DO

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    DO i = NLMAX, NLMIN, -1
      CALL spdiscr_releaseBlockDiscr(Rlevels(i)%rdiscrSecondary)
      CALL spdiscr_releaseBlockDiscr(Rlevels(i)%rdiscrPrimary)
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