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

module bouss2dmini_mit_2

  use fsystem
  use storage
  use genoutput
  use boundary
  use cubature
  use derivatives
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use element
  use spatialdiscretisation
  use linearalgebra
  use linearsystemscalar
  use linearsystemblock
  use spdiscprojection
  use scalarpde
  use trilinearformevaluation
  use bilinearformevaluation
  use linearformevaluation
  use stdoperators
  use discretebc
  use filtersupport
  use coarsegridcorrection
  use multilevelprojection
  use linearsolver
  use ucd
  use convection
  use statistics
  
  use bouss2dmini_callback
  
  implicit none

!<types>

!<typeblock description="Type block defining all information about one level">

  type t_level
  
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! A block discretisation for the primary system (Q1~/Q1~/Q0)
    type(t_blockDiscretisation) :: rdiscrPrimary
    
    ! A block discretisation for the secondary system (Q1)
    type(t_blockDiscretisation) :: rdiscrSecondary
    
    ! A primary system matrix which recieves the Navier-Stokes system
    ! for that specific level.
    type(t_matrixBlock) :: rmatrixPrimary
    
    ! A secondary system matrix which recieves the Convection-Diffusion system
    ! of the temperature for that specific level.
    type(t_matrixBlock) :: rmatrixSecondary
    
    ! Laplace-matrix of the primary system for that specific level.
    type(t_matrixScalar) :: rmatrixA

    ! B1-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixB1

    ! B2-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixB2
    
    ! M1- / M2-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixM
    
    ! Laplace-Matrix of the secondary system for that specific level.
    type(t_matrixScalar) :: rmatrixN
    
    ! Primary solution vector for that specific level.
    type(t_vectorBlock) :: rvecSolPrimary
    
    ! Secondary solution vector for that specific level.
    type(t_vectorBlock) :: rvecSolSecondary

    ! A variable describing the discrete boundary conditions for the
    ! primary system.
    type(t_discreteBC) :: rdiscreteBCprimary

    ! A variable describing the discrete boundary conditions for the
    ! secondary system.
    type(t_discreteBC) :: rdiscreteBCsecondary
  
  end type
  
!</typeblock>

!</types>

contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine b2dm_mit_2
  
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
    type(t_level), dimension(:), pointer :: Rlevels

    ! An object for saving the domain:
    type(t_boundary) :: rboundary
    
    ! A trilinear form describing the convection.
    type(t_trilinearform) :: rtriform1, rtriform2
    
    ! A streamline diffusion object for the convection.
    type(t_convStreamlineDiffusion) :: rsd

    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_vectorBlock) :: rrhsPrimary,rrhsSecondary,&
        rvecTmpPrimary,rvecTmpSecondary,rvecDefPrimary,rvecDefSecondary
    
    ! A set of variables describing the analytic and discrete boundary
    ! conditions.
    type(t_boundaryRegion) :: rboundaryRegion

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverPrimary,p_rsolverSecondary,&
        p_rcoarseGridSolver,p_rsmoother,p_rpreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(:), pointer :: Rmatrices

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    type(t_filterChain), dimension(1), target :: RfilterChain
    type(t_filterChain), dimension(:), pointer :: p_RfilterChain
    
    ! An interlevel projection structure for changing levels
    type(t_interlevelProjectionBlock) :: rprojPrimary, rprojSecondary

    ! One level of multigrid
    type(t_linsolMGLevelInfo), pointer :: p_rlevelInfo
    
    ! NLMIN receives the level of the coarse grid.
    integer :: NLMIN

    ! NLMAX receives the level where we want to solve.
    integer :: NLMAX
    
    ! Viscosity parameter nu = 1/Re
    real(DP) :: dnu
    
    ! Thermal conductivity eta
    real(DP) :: deta
    
    ! Gravity direction (dgrav1, dgrav2)
    real(DP) :: dgrav1, dgrav2
    
    ! Non-linear damping parameters
    real(DP) :: dnlDampPrimary, dnlDampSecondary
    
    ! Error indicator during initialisation of the solver
    integer :: ierror
    
    ! A counter variable for the non-linear loop
    integer :: NLiter
    
    ! Non-linear residuals
    real(DP) :: dnlResPrimary, dnlResSecondary, dnlResTotal !, dnlResScale
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    real(DP), dimension(:), pointer :: p_Ddata,p_Ddata2

    ! A counter variable
    integer :: i
    logical :: bError
    
    ! A timer object for time measuring
    type(t_timer) :: rtimer

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
    allocate(Rlevels(NLMIN:NLMAX))

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call boundary_read_prm(rboundary, './pre/mit1.prm')
        
    ! Now read in the basic triangulation.
    call tria_readTriFile2D (Rlevels(NLMIN)%rtriangulation, &
                             './pre/mit1.tri', rboundary)
    
    ! Refine the mesh up to the minimum level
    call tria_quickRefine2LevelOrdering (NLMIN-1,&
        Rlevels(NLMIN)%rtriangulation,rboundary)
    
    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    call tria_initStandardMeshFromRaw (Rlevels(NLMIN)%rtriangulation,&
        rboundary)
    
    ! Now refine the grid for the fine levels.
    do i = NLMIN+1, NLMAX

      ! Refine the grid using the 2-Level-Ordering algorithm
      call tria_refine2LevelOrdering(Rlevels(i-1)%rtriangulation,&
          Rlevels(i)%rtriangulation,rboundary)
      
      ! Create a standard mesh
      call tria_initStandardMeshFromRaw(Rlevels(i)%rtriangulation,&
        rboundary)
    
    end do

    ! Now we can start to initialise the discretisation.
    do i = NLMIN, NLMAX
    
      ! Our primary discretisation has 3 components (Q1~/Q1~/Q0)
      call spdiscr_initBlockDiscr (Rlevels(i)%rdiscrPrimary, 3, &
                                   Rlevels(i)%rtriangulation, rboundary)
      
      ! For simplicity, we set up one discretisation structure for the
      ! velocity...
      call spdiscr_initDiscr_simple (Rlevels(i)%rdiscrPrimary%RspatialDiscr(1),&
          EL_EM30, CUB_G2X2, Rlevels(i)%rtriangulation, rboundary)
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory,
      ! as both structures will share the same dynamic information afterwards.
      call spdiscr_duplicateDiscrSc (&
          Rlevels(i)%rdiscrPrimary%RspatialDiscr(1),&
          Rlevels(i)%rdiscrPrimary%RspatialDiscr(2))

      ! For the pressure (3rd component), we set up a separate discretisation
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      call spdiscr_deriveSimpleDiscrSc (Rlevels(i)%rdiscrPrimary%RspatialDiscr(1), &
          EL_Q0, CUB_G2X2, Rlevels(i)%rdiscrPrimary%RspatialDiscr(3))

      ! Our secondary discretisation has 1 component (Q1)
      call spdiscr_initBlockDiscr (Rlevels(i)%rdiscrSecondary, 1, &
                                   Rlevels(i)%rtriangulation, rboundary)

      call spdiscr_initDiscr_simple (Rlevels(i)%rdiscrSecondary%RspatialDiscr(1),&
          EL_Q1, CUB_G2X2, Rlevels(i)%rtriangulation, rboundary)
    
    end do

    do i = NLMIN, NLMAX
    
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! Initialise primary system for this level
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    
      ! Initialise a block matrix based on the primary discretisation.
      call lsysbl_createMatBlockByDiscr (Rlevels(i)%rdiscrPrimary,&
                                         Rlevels(i)%rmatrixPrimary)
      
      ! Inform the matrix that we build a saddle-point problem.
      Rlevels(i)%rmatrixPrimary%imatrixSpec = LSYSBS_MSPEC_SADDLEPOINT
      
      ! Create the matrix structure of the A-matrices of the primary system.
      call bilf_createMatrixStructure (Rlevels(i)%rdiscrPrimary%RspatialDiscr(1),&
          LSYSSC_MATRIX9, Rlevels(i)%rmatrixA)

      ! Create the matrix structure of the B-matrices
      call bilf_createMatrixStructure (Rlevels(i)%rdiscrPrimary%RspatialDiscr(3),&
          LSYSSC_MATRIX9, Rlevels(i)%rmatrixB1,&
          Rlevels(i)%rdiscrPrimary%RspatialDiscr(1))

      ! Since B2 has the same structure as B1, we'll simply get a copy of it.
      call lsyssc_duplicateMatrix (Rlevels(i)%rmatrixB1, Rlevels(i)%rmatrixB2,&
                                   LSYSSC_DUP_SHARE, LSYSSC_DUP_REMOVE)
      
      ! Build the Laplace-Matrix of the primary system.
      call stdop_assembleLaplaceMatrix(Rlevels(i)%rmatrixA, .true., dnu)
      
      ! Build the B1-matrix entries
      call stdop_assembleSimpleMatrix(Rlevels(i)%rmatrixB1, DER_FUNC, DER_DERIV_X,&
                                      -1.0_DP, .true.)

      ! Build the B2-matrix entries
      call stdop_assembleSimpleMatrix(Rlevels(i)%rmatrixB2, DER_FUNC, DER_DERIV_Y,&
                                      -1.0_DP, .true.)
      
      ! Now let's set up our primary system block matrix.
      
      ! Copy the A-matrix into our primary system block matrix.
      call lsyssc_duplicateMatrix (Rlevels(i)%rmatrixA, &
          Rlevels(i)%rmatrixPrimary%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      call lsyssc_duplicateMatrix (Rlevels(i)%rmatrixPrimary%RmatrixBlock(1,1), &
          Rlevels(i)%rmatrixPrimary%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      
      ! Copy the B1-/B2-matrices into our primary system block matrix.
      call lsyssc_duplicateMatrix (Rlevels(i)%rmatrixB1, &
          Rlevels(i)%rmatrixPrimary%RmatrixBlock(1,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      call lsyssc_duplicateMatrix (Rlevels(i)%rmatrixB2, &
          Rlevels(i)%rmatrixPrimary%RmatrixBlock(2,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      
      ! Furthermore, put B1^T and B2^T to the block matrix.
      call lsyssc_transposeMatrix (Rlevels(i)%rmatrixB1, &
          Rlevels(i)%rmatrixPrimary%RmatrixBlock(3,1),LSYSSC_TR_VIRTUAL)

      call lsyssc_transposeMatrix (Rlevels(i)%rmatrixB2, &
          Rlevels(i)%rmatrixPrimary%RmatrixBlock(3,2),LSYSSC_TR_VIRTUAL)

      ! Update the structural information of the block matrix.
      call lsysbl_updateMatStrucInfo (Rlevels(i)%rmatrixPrimary)

      ! Finally, create a solution vector based on the matrix.
      call lsysbl_createVecBlockIndMat (Rlevels(i)%rmatrixPrimary,&
                                        Rlevels(i)%rvecSolPrimary,.false.)

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! Initialise secondary system for this level
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      ! Initialise a block matrix based on the secondary discretisation.
      call lsysbl_createMatBlockByDiscr (Rlevels(i)%rdiscrSecondary,&
                                         Rlevels(i)%rmatrixSecondary)

      ! Create the matrix structure of the secondary system.
      call bilf_createMatrixStructure (Rlevels(i)%rdiscrSecondary%RspatialDiscr(1),&
          LSYSSC_MATRIX9, Rlevels(i)%rmatrixN)

      ! Build the Laplace-Matrix of the secondary system.
      call stdop_assembleLaplaceMatrix(Rlevels(i)%rmatrixN, .true., deta)

      ! Copy the N-matrix into our secondary system block matrix.
      call lsyssc_duplicateMatrix (Rlevels(i)%rmatrixN, &
          Rlevels(i)%rmatrixSecondary%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

      ! Update the structural information of the block matrix.
      call lsysbl_updateMatStrucInfo (Rlevels(i)%rmatrixSecondary)

      ! Finally, create a solution vector based on the matrix.
      call lsysbl_createVecBlockIndMat (Rlevels(i)%rmatrixSecondary,&
                                        Rlevels(i)%rvecSolSecondary,.false.)

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! Initialise coupling matrices for this level
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      ! Create the matrix structure of the M-matrices
      call bilf_createMatrixStructure (Rlevels(i)%rdiscrSecondary%RspatialDiscr(1),&
          LSYSSC_MATRIX9, Rlevels(i)%rmatrixM,&
          Rlevels(i)%rdiscrPrimary%RspatialDiscr(1))

      ! Build a mass matrix.
      call stdop_assembleSimpleMatrix(Rlevels(i)%rmatrixM, DER_FUNC, DER_FUNC,&
                                      1.0_DP, .true.)
      
      ! That's it for now.
    
    end do

    ! Create the RHS vectors on the finest level.
    call lsysbl_createVecBlockIndMat (Rlevels(NLMAX)%rmatrixPrimary,&
                                      rrhsPrimary, .false.)
    call lsysbl_createVecBlockIndMat (Rlevels(NLMAX)%rmatrixPrimary,&
                                      rvecTmpPrimary, .false.)
    call lsysbl_createVecBlockIndMat (Rlevels(NLMAX)%rmatrixPrimary,&
                                      rvecDefPrimary, .false.)
    call lsysbl_createVecBlockIndMat (Rlevels(NLMAX)%rmatrixSecondary,&
                                      rrhsSecondary, .false.)
    call lsysbl_createVecBlockIndMat (Rlevels(NLMAX)%rmatrixSecondary,&
                                      rvecTmpSecondary, .false.)
    call lsysbl_createVecBlockIndMat (Rlevels(NLMAX)%rmatrixSecondary,&
                                      rvecDefSecondary, .false.)

    ! Clear both RHS vectors and the solution vectors on the finest level.
    call lsysbl_clearVector(rrhsPrimary)
    call lsysbl_clearVector(rrhsSecondary)
    call lsysbl_clearVector(Rlevels(NLMAX)%rvecSolPrimary)
    call lsysbl_clearVector(Rlevels(NLMAX)%rvecSolSecondary)

    do i = NLMIN, NLMAX
    
      ! Now implement the boundary conditions.
      call bcasm_initDiscreteBC(Rlevels(i)%rdiscreteBCprimary)
      call bcasm_initDiscreteBC(Rlevels(i)%rdiscreteBCsecondary)

      ! Get the bottom edge - Dirichlet(0) BCs for X- and Y-velocities, but
      ! Neumann for the temperature.
      call boundary_createRegion(rboundary,1,1,rboundaryRegion)
      
      ! X-/Y-velocity
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscrPrimary,1,&
          rboundaryRegion, Rlevels(i)%rdiscreteBCprimary, getBoundaryValuesPrimary)
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscrPrimary,2,&
          rboundaryRegion, Rlevels(i)%rdiscreteBCprimary, getBoundaryValuesPrimary)
      
      ! Get the left edge - Dirichlet(0) BC for the X- and Y-velocities and
      ! Dirichlet(-1/2) for the temperature.
      call boundary_createRegion(rboundary,1,2,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
      
      ! X-/Y-velocity
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscrPrimary,1,&
          rboundaryRegion, Rlevels(i)%rdiscreteBCprimary, getBoundaryValuesPrimary)
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscrPrimary,2,&
          rboundaryRegion, Rlevels(i)%rdiscreteBCprimary, getBoundaryValuesPrimary)
                               
      ! temperature
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscrSecondary,1,&
          rboundaryRegion, Rlevels(i)%rdiscreteBCsecondary, getBoundaryValuesSecondary)

      ! Get the top edge - Dirichlet(0) BCs for X- and Y-velocities, but
      ! Neumann for the temperature.
      call boundary_createRegion(rboundary,1,3,rboundaryRegion)

      ! X-/Y-velocity
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscrPrimary,1,&
          rboundaryRegion, Rlevels(i)%rdiscreteBCprimary, getBoundaryValuesPrimary)
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscrPrimary,2,&
          rboundaryRegion, Rlevels(i)%rdiscreteBCprimary, getBoundaryValuesPrimary)
      
      ! Get the right edge - Dirichlet(0) BC for the X- and Y-velocities and
      ! Dirichlet(1/2) for the temperature.
      call boundary_createRegion(rboundary,1,4,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND

      ! X-/Y-velocity
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscrPrimary,1,&
          rboundaryRegion, Rlevels(i)%rdiscreteBCprimary, getBoundaryValuesPrimary)
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscrPrimary,2,&
          rboundaryRegion, Rlevels(i)%rdiscreteBCprimary, getBoundaryValuesPrimary)

      ! temperature
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscrSecondary,1,&
          rboundaryRegion, Rlevels(i)%rdiscreteBCsecondary, getBoundaryValuesSecondary)

      ! Hang in the discrete BCs into the system matrices and solution vectors.
      Rlevels(i)%rmatrixPrimary%p_rdiscreteBC => Rlevels(i)%rdiscreteBCprimary
      Rlevels(i)%rmatrixSecondary%p_rdiscreteBC => Rlevels(i)%rdiscreteBCsecondary
      Rlevels(i)%rvecSolPrimary%p_rdiscreteBC => Rlevels(i)%rdiscreteBCprimary
      Rlevels(i)%rvecSolSecondary%p_rdiscreteBC => Rlevels(i)%rdiscreteBCsecondary
      
    end do

    ! Hang in the BCs into the RHS vectors.
    rrhsPrimary%p_rdiscreteBC => Rlevels(NLMAX)%rdiscreteBCprimary
    rvecDefPrimary%p_rdiscreteBC => Rlevels(NLMAX)%rdiscreteBCprimary
    rrhsSecondary%p_rdiscreteBC => Rlevels(NLMAX)%rdiscreteBCsecondary
    rvecDefSecondary%p_rdiscreteBC => Rlevels(NLMAX)%rdiscreteBCsecondary
    
    ! And filter the RHS and solution vectors on the finest level.
    call vecfil_discreteBCrhs (rrhsPrimary)
    call vecfil_discreteBCrhs (rrhsSecondary)
    call vecfil_discreteBCsol (Rlevels(NLMAX)%rvecSolPrimary)
    call vecfil_discreteBCsol (Rlevels(NLMAX)%rvecSolSecondary)

    ! Initialise a filter chain.
    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

    ! Set up two inter-level-projection structures - one for the primary and
    ! one for the secondary system.
    call mlprj_initProjectionMat (rprojPrimary,Rlevels(NLMAX)%rmatrixPrimary)
    call mlprj_initProjectionMat (rprojSecondary,Rlevels(NLMAX)%rmatrixSecondary)

    ! Create two multi-grid solvers.
    p_RfilterChain => RfilterChain
    call linsol_initMultigrid (p_rsolverPrimary, p_RfilterChain)
    call linsol_initMultigrid (p_rsolverSecondary, p_RfilterChain)
    
    ! Set up a BiCGStab solver with VANKA preconditioner as coarse grid solver
    ! for the primary system.
    !CALL linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DNAVST)
    !CALL linsol_initBiCGStab (p_rcoarseGridSolver,p_rpreconditioner,p_RfilterChain)
    call linsol_initUMFPACK4(p_rcoarseGridSolver)
    call linsol_addMultiGridLevel(p_rlevelInfo,p_rsolverPrimary,rprojPrimary,&
                                  null(), null(), p_rcoarseGridSolver)

    ! Set up a BiCGStab solver with SSOR preconditioner as coarse grid solver
    ! for the secondary system.
    !CALL linsol_initSSOR (p_rpreconditioner)
    !CALL linsol_initBiCGStab (p_rcoarseGridSolver,p_rpreconditioner,p_RfilterChain)
    call linsol_initUMFPACK4(p_rcoarseGridSolver)
    call linsol_addMultiGridLevel(p_rlevelInfo,p_rsolverSecondary,rprojSecondary,&
                                  null(), null(), p_rcoarseGridSolver)

    ! Now set up the other levels...
    do i = NLMIN+1, NLMAX
    
      ! Set up a VANKA smoother for the primary system.
      call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DNAVST)
      call linsol_convertToSmoother(p_rsmoother, 4, 1.0_DP)
      call linsol_addMultiGridLevel(p_rlevelInfo,p_rsolverPrimary,rprojPrimary,&
                                    p_rsmoother, p_rsmoother, null())
      
      ! Set up a SSOR smoother for the secondary system.
      call linsol_initSSOR (p_rsmoother)
      call linsol_convertToSmoother(p_rsmoother, 4, 1.0_DP)
      call linsol_addMultiGridLevel(p_rlevelInfo,p_rsolverSecondary,rprojSecondary,&
                                    p_rsmoother, p_rsmoother, null())

    end do
    
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
    allocate(Rmatrices(NLMIN:NLMAX))
    ! Primary System
    do i = NLMIN, NLMAX
      call lsysbl_duplicateMatrix (Rlevels(i)%rmatrixPrimary,&
          Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    end do
    call linsol_setMatrices(p_RsolverPrimary,Rmatrices(NLMIN:NLMAX))
    do i=NLMIN,NLMAX
      call lsysbl_releaseMatrix (Rmatrices(i))
    end do

    ! Secondary System
    do i = NLMIN, NLMAX
      call lsysbl_duplicateMatrix (Rlevels(i)%rmatrixSecondary,&
          Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    end do
    call linsol_setMatrices(p_RsolverSecondary,Rmatrices(NLMIN:NLMAX))
    do i=NLMIN,NLMAX
      call lsysbl_releaseMatrix (Rmatrices(i))
    end do

    deallocate(Rmatrices)
    
    ! Initialise only the structure of the solvers - the data is initialised
    ! in the non-linear loop below.
    call linsol_initStructure (p_rsolverPrimary, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) stop
    call linsol_initStructure (p_rsolverSecondary, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) stop


    ! Now prepare the trilinear forms for the convection in the secondary system.
    ! u1 * d_x * T
    rtriform1%itermCount = 1
    rtriform1%ballCoeffConstant = .true.
    rtriform1%BconstantCoeff = .true.
    rtriform1%Dcoefficients = 1.0_DP
    rtriform1%Idescriptors(1,1) = DER_FUNC
    rtriform1%Idescriptors(2,1) = DER_DERIV_X
    rtriform1%Idescriptors(3,1) = DER_FUNC
    ! u2 * d_y * T
    rtriform2%itermCount = 1
    rtriform2%ballCoeffConstant = .true.
    rtriform2%BconstantCoeff = .true.
    rtriform2%Dcoefficients = 1.0_DP
    rtriform2%Idescriptors(1,1) = DER_FUNC
    rtriform2%Idescriptors(2,1) = DER_DERIV_Y
    rtriform2%Idescriptors(3,1) = DER_FUNC
    
    ! And prepare the streamline diffusion object for the convection in the
    ! primary system.
    rsd%dnu = dnu

    ! Start the timer
    call stat_clearTimer(rtimer)
    call stat_startTimer(rtimer)

    ! Now comes the fun part - the non-linear loop.
    do NLiter = 0, 50
    
      ! The very first thing is we need to update the non-linear parts of
      ! the system matrices, but to do this, we need the solution vectors
      ! on all levels - so let's restrict the solution vectors.
      do i = NLMAX, NLMIN+1, -1
        
        ! Restrict primary solution vector.
        call mlprj_performInterpolation(rprojPrimary, Rlevels(i-1)%rvecSolPrimary,&
                         Rlevels(i)%rvecSolPrimary, rvecTmpPrimary%RvectorBlock(1))

        ! Restrict secondary solution vector.
        call mlprj_performInterpolation(rprojSecondary, Rlevels(i-1)%rvecSolSecondary,&
                         Rlevels(i)%rvecSolSecondary, rvecTmpSecondary%RvectorBlock(1))
        
        ! And filter the restricted vectors.
        call vecfil_discreteBCsol (Rlevels(i-1)%rvecSolPrimary)
        call vecfil_discreteBCsol (Rlevels(i-1)%rvecSolSecondary)

      end do
      
      ! Now we can update the matrices - so go for it.
      do i = NLMIN, NLMAX
      
        ! Copy A-matrix into the X-velocity block of the primary system.
        call lsyssc_duplicateMatrix(Rlevels(i)%rmatrixA,&
            Rlevels(i)%rmatrixPrimary%RmatrixBlock(1,1),&
            LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPY)
        
        ! Assemble the convection of the primary system using streamline diffusion.
        call conv_streamlineDiffusion2d(Rlevels(i)%rvecSolPrimary, &
            Rlevels(i)%rvecSolPrimary, 1.0_DP, 0.0_DP, rsd, CONV_MODMATRIX, &
            Rlevels(i)%rmatrixPrimary%RmatrixBlock(1,1))
        
        ! And filter the primary system matrix.
        call matfil_discreteBC (Rlevels(i)%rmatrixPrimary)
        
        ! The Y-velocity block is automatically updated, as it is just a shared
        ! copy of the X-velocity block.
        
        ! Copy N-matrix into the secondary system block matrix.
        call lsyssc_duplicateMatrix(Rlevels(i)%rmatrixN,&
            Rlevels(i)%rmatrixSecondary%RmatrixBlock(1,1),&
            LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPY)

        ! Assemble the convection of the secondary system.
        ! Assemble u_1 * d_x T
        call trilf_buildMatrixScalar(rtriform1,.false.,&
             Rlevels(i)%rmatrixSecondary%RmatrixBlock(1,1), &
             Rlevels(i)%rvecSolPrimary%RvectorBlock(1))
        ! Assemble u_2 * d_y T
        call trilf_buildMatrixScalar(rtriform2,.false.,&
             Rlevels(i)%rmatrixSecondary%RmatrixBlock(1,1), &
             Rlevels(i)%rvecSolPrimary%RvectorBlock(2))
        
        ! And filter the secondary system matrix.
        call matfil_discreteBC (Rlevels(i)%rmatrixSecondary)

      end do
      
      ! Okay, the system matrices are up-to-date. Now we need to calculate
      ! the non-linear defect. But as we do not have a global system matrix,
      ! we need to calculate the defect of both the primary and secondary
      ! systems.
      ! So let's begin with the secondary system.
      call lsysbl_copyVector(rrhsSecondary, rvecDefSecondary)
      call lsysbl_blockMatVec(Rlevels(NLMAX)%rmatrixSecondary, &
          Rlevels(NLMAX)%rvecSolSecondary, rvecDefSecondary, -1.0_DP, 1.0_DP)
      
      ! Filter the defect vector.
      call vecfil_discreteBCdef(rvecDefSecondary)
      
      ! Calculate the secondary residual.
      dnlResSecondary = lsysbl_vectorNorm(rvecDefSecondary, LINALG_NORMEUCLID)
      
      ! Calculate the defect of the primary system.
      call lsysbl_copyVector(rrhsPrimary, rvecDefPrimary)
      call lsysbl_blockMatVec(Rlevels(NLMAX)%rmatrixPrimary, &
          Rlevels(NLMAX)%rvecSolPrimary, rvecDefPrimary, -1.0_DP, 1.0_DP)
      
      ! Don't forget the M-matrices which couple our primary system to
      ! the secondary system.
      call lsyssc_scalarMatVec(Rlevels(NLMAX)%rmatrixM,&
          Rlevels(NLMAX)%rvecSolSecondary%RvectorBlock(1),&
          rvecDefPrimary%RvectorBlock(1), -dgrav1, 1.0_DP)
      call lsyssc_scalarMatVec(Rlevels(NLMAX)%rmatrixM,&
          Rlevels(NLMAX)%rvecSolSecondary%RvectorBlock(1),&
          rvecDefPrimary%RvectorBlock(2), -dgrav2, 1.0_DP)
      
      ! Filter the defect vector
      call vecfil_discreteBCdef(rvecDefPrimary)
      
      ! Calculate the primary residual.
      dnlResPrimary = lsysbl_vectorNorm(rvecDefPrimary, LINALG_NORMEUCLID)
      
      ! And calculate the total residual.
      dnlResTotal = sqrt(dnlResPrimary**2 + dnlResSecondary**2)
      
      ! Print the residual
      !CALL output_separator(OU_SEP_MINUS)
      call output_line('NL-Iteration: ' // trim(sys_si(NLiter,2)) // &
                       '   |RES| = ' // trim(sys_sdEP(dnlResTotal,20,12)))
      !CALL output_separator(OU_SEP_MINUS)
      
      ! Is this precise enough?
      if(dnlResTotal .le. 1E-8_DP) exit
      
      ! Okay, we need to perform another non-linear iteration.
      
      !IF(dnlResSecondary .GT. 1E-10_DP) THEN
      if(.true.) then
      
        ! So let's first solve the secondary system to get a new temperature.

        ! Initialise the data of the secondary solver.
        call linsol_initData(p_rsolverSecondary, ierror)
        if(ierror .ne. LINSOL_ERR_NOERROR) stop
        
        ! And call the secondary solver
        !CALL output_line('Solving secondary system...')
        !CALL output_line('---------------------------')
        call linsol_precondDefect(p_rsolverSecondary, rvecDefSecondary)
        bError = (p_rsolverSecondary%iresult .ne. 0)
        
        ! Release the secondary solver data.
        call linsol_doneData(p_rsolverSecondary)
        
        ! Did the secondary solver break down?
        if(bError) then

          ! Print an error
          call output_separator(OU_SEP_STAR)
          call output_line('NL-Iteration: ERROR: secondary solver broke down')
          call output_separator(OU_SEP_STAR)
          
          ! Exit the non-linear loop
          exit
        
        end if
        
        ! Now update the secondary solution vector - which is the temperature.
        call lsysbl_vectorLinearComb(rvecDefSecondary, &
            Rlevels(NLMAX)%rvecSolSecondary, dnlDampSecondary, 1.0_DP)
        
        ! Now recalculate the defect of the primary system as we have updated
        ! the temperature.
        call lsysbl_copyVector(rrhsPrimary, rvecDefPrimary)
        call lsysbl_blockMatVec(Rlevels(NLMAX)%rmatrixPrimary, &
            Rlevels(NLMAX)%rvecSolPrimary, rvecDefPrimary, -1.0_DP, 1.0_DP)
        call lsyssc_scalarMatVec(Rlevels(NLMAX)%rmatrixM,&
            Rlevels(NLMAX)%rvecSolSecondary%RvectorBlock(1),&
            rvecDefPrimary%RvectorBlock(1), -dgrav1, 1.0_DP)
        call lsyssc_scalarMatVec(Rlevels(NLMAX)%rmatrixM,&
            Rlevels(NLMAX)%rvecSolSecondary%RvectorBlock(1),&
            rvecDefPrimary%RvectorBlock(2), -dgrav2, 1.0_DP)
        
        ! Filter the defect vector
        call vecfil_discreteBCdef(rvecDefPrimary)

        ! Calculate the primary residual.
        dnlResPrimary = lsysbl_vectorNorm(rvecDefPrimary, LINALG_NORMEUCLID)
      
      end if
      
      !IF (dnlResPrimary .GT. 1E-10_DP) THEN
      if(.true.) then

        ! Initialise the data of the primary solver.
        call linsol_initData(p_rsolverPrimary, ierror)
        if(ierror .ne. LINSOL_ERR_NOERROR) stop
        
        ! And call the primary solver
        !CALL output_line('Solving primary system...')
        !CALL output_line('---------------------------')
        call linsol_precondDefect(p_rsolverPrimary, rvecDefPrimary)
        bError = (p_rsolverPrimary%iresult .ne. 0)
        
        ! Release the primary solver data.
        call linsol_doneData(p_rsolverPrimary)
        
        ! Did the primary solver break down?
        if(bError) then

          ! Print an error
          call output_separator(OU_SEP_STAR)
          call output_line('NL-Iteration: ERROR: Primary solver broke down')
          call output_separator(OU_SEP_STAR)
          
          ! Exit the non-linear loop
          exit
        
        end if
        
        ! Now update the primary solution vector.
        call lsysbl_vectorLinearComb(rvecDefPrimary, &
            Rlevels(NLMAX)%rvecSolPrimary, dnlDampPrimary, 1.0_DP)
      
      end if
      
      ! That's it for this non-linear iteration...
          
    end do ! NLiter
    
    ! Stop the timer and print the elapsed time
    call stat_stopTimer(rtimer)
    call output_line('Total time for non-linear loop: ' // stat_sgetTime_byTimer(rtimer))

    ! We can now start the postprocessing.
    ! Start UCD export to GMV file:
    call ucd_startGMV (rexport,UCD_FLAG_STANDARD,&
        Rlevels(NLMAX)%rtriangulation,'gmv/u_mit_2.gmv')
    
    ! Project X- and Y-velocity to the vertices of the mesh.
    nullify(p_Ddata)
    nullify(p_Ddata2)
    call spdp_projectToVertices(Rlevels(NLMAX)%rvecSolPrimary%RvectorBlock(1), p_Ddata)
    call spdp_projectToVertices(Rlevels(NLMAX)%rvecSolPrimary%RvectorBlock(2), p_Ddata2)

    ! In case we use the VTK exporter, which supports vector output, we will
    ! pass the X- and Y-velocity at once to the ucd module.
    call ucd_addVarVertBasedVec(rexport,'velocity',p_Ddata,p_Ddata2)
    
    ! Release the temporary memory.
    deallocate(p_Ddata)
    deallocate(p_Ddata2)
        
    ! Write pressure
    call lsyssc_getbase_double (Rlevels(NLMAX)%rvecSolPrimary%RvectorBlock(3),p_Ddata)
    call ucd_addVariableElementBased (rexport,'pressure',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write temperature
    call lsyssc_getbase_double (Rlevels(NLMAX)%rvecSolSecondary%RvectorBlock(1),p_Ddata)
    call ucd_addVariableVertexBased (rexport,'temperature',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that's it.
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! We are finished - but not completely!
    ! Now, clean up so that all the memory is available again.
    !
    ! Release solver data and structure
   ! CALL linsol_doneData (p_rsolverNode)
    call linsol_doneStructure (p_rsolverSecondary)
    call linsol_doneStructure (p_rsolverPrimary)
    
    ! Release the interlevel projection structure
    call mlprj_doneProjection (rprojSecondary)
    call mlprj_doneProjection (rprojPrimary)

    ! Release the solver node and all subnodes attached to it (if at all):
    call linsol_releaseSolver (p_rsolverSecondary)
    call linsol_releaseSolver (p_rsolverPrimary)
    
    ! Release the block matrix/vectors
    call lsysbl_releaseVector (rvecDefPrimary)
    call lsysbl_releaseVector (rvecDefSecondary)
    call lsysbl_releaseVector (rvecTmpPrimary)
    call lsysbl_releaseVector (rvecTmpSecondary)
    call lsysbl_releaseVector (rrhsPrimary)
    call lsysbl_releaseVector (rrhsSecondary)
    do i = NLMAX, NLMIN, -1
      call lsysbl_releaseVector (Rlevels(i)%rvecSolSecondary)
      call lsysbl_releaseVector (Rlevels(i)%rvecSolPrimary)
      call lsysbl_releaseMatrix (Rlevels(i)%rmatrixSecondary)
      call lsysbl_releaseMatrix (Rlevels(i)%rmatrixPrimary)
    end do
    
    ! Release scalar matrices
    do i = NLMAX, NLMIN, -1
      call lsyssc_releaseMatrix (Rlevels(i)%rmatrixM)
      call lsyssc_releaseMatrix (Rlevels(i)%rmatrixN)
      call lsyssc_releaseMatrix (Rlevels(i)%rmatrixB2)
      call lsyssc_releaseMatrix (Rlevels(i)%rmatrixB1)
      call lsyssc_releaseMatrix (Rlevels(i)%rmatrixA)
    end do
    
    ! Release our discrete version of the boundary conditions
    do i = NLMAX, NLMIN, -1
      call bcasm_releaseDiscreteBC (Rlevels(i)%rdiscreteBCsecondary)
      call bcasm_releaseDiscreteBC (Rlevels(i)%rdiscreteBCprimary)
    end do

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    do i = NLMAX, NLMIN, -1
      call spdiscr_releaseBlockDiscr(Rlevels(i)%rdiscrSecondary)
      call spdiscr_releaseBlockDiscr(Rlevels(i)%rdiscrPrimary)
    end do
    
    ! Release the triangulation.
    do i = NLMAX, NLMIN, -1
      call tria_done (Rlevels(i)%rtriangulation)
    end do
    
    deallocate(Rlevels)
    
    ! Finally release the domain, that's it.
    call boundary_release (rboundary)

  end subroutine

end module
