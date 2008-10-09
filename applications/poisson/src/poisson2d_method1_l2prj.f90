!##############################################################################
!# ****************************************************************************
!# <name> poisson2d_method1_l2prj </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple Poisson
!# problem with constant coefficients on a simple domain.
!#
!# This module is based on poisson2d_method1_mg, but using the EB30 element
!# (Q1~ with bubble) and L2-projection for multi-level operations and
!# post-processing, i.e. prolongation, restiction and the spatial projection
!# for the UCD output.
!#
!# </purpose>
!##############################################################################

MODULE poisson2d_method1_l2prj

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
  USE stdoperators
  USE multileveloperators
    
  USE poisson2d_callback
  
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
    
    ! A mass matrix for that specific level. This matrix is needed by the
    ! L2-projection multi-level operators and is defined on each level except
    ! for the coarse-most one.
    TYPE(t_matrixScalar) :: rmatrixMass
    
    ! A 2-Level-Mass matrix for that specific level. This matrix is needed by
    ! L2-projection multi-level operators and is defined on each level except
    ! for the coarse-most one.
    TYPE(t_matrixScalar) :: rmatrix2Lvl

    ! An interlevel-projection structure for that specific level.
    TYPE(t_interlevelProjectionBlock) :: rprojection

    ! A variable describing the discrete boundary conditions.    
    TYPE(t_discreteBC) :: rdiscreteBC
  
  END TYPE
  
!</typeblock>

!</types>

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE poisson2d_1_l2prj
  
!<description>
  ! This is an all-in-one poisson solver for directly solving a Poisson
  ! problem without making use of special features like collections
  ! and so on. The routine performs the following tasks:
  !
  !  1.) Read in parametrisation
  !  2.) Read in triangulation
  !  3.) Set up RHS
  !  4.) Set up system matrix
  !  5.) Set up 2-Level-projection operators for multi-grid
  !  6.) Create solver structure
  !  7.) Solve the problem
  !  8.) Perform L2-Projection of the solution to Q1 space
  !  9.) Write solution to GMV file
  ! 10.) Release all variables, finish
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
    
    ! A bilinear and linear form describing the analytic problem to solve
    TYPE(t_bilinearForm) :: rform
    TYPE(t_linearForm) :: rlinform
    
    ! A couple of block vectors. These will be filled
    ! with data for the linear solver.
    TYPE(t_vectorBlock) :: rvectorBlock,rrhsBlock,rtempBlock

    ! A variable that is used to specify a region on the boundary.
    TYPE(t_boundaryRegion) :: rboundaryRegion

    ! A solver node that accepts parameters for the linear solver    
    TYPE(t_linsolNode), POINTER :: p_rsolverNode,p_rcoarseGridSolver,p_rsmoother

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    TYPE(t_matrixBlock), DIMENSION(:), ALLOCATABLE :: Rmatrices

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    TYPE(t_filterChain), DIMENSION(1), TARGET :: RfilterChain
    TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain

    ! One level of multigrid
    TYPE(t_linsolMGLevelInfo), POINTER :: p_rlevelInfo
    
    ! NLMIN receives the level of the coarse grid.
    INTEGER :: NLMIN

    ! NLMAX receives the level where we want to solve.
    INTEGER :: NLMAX
    
    ! Error indicator during initialisation of the solver
    INTEGER :: ierror
    
    ! One spatial discretisation for the L2-projection of the solution
    TYPE(t_spatialDiscretisation) :: rdiscrQ1
    
    ! Two scalar matrices for the L2-projection
    TYPE(t_matrixScalar) :: rmatrixMassPrj, rlumpedMassPrj
    
    ! Three scalar vectors for the L2-projection
    TYPE(t_vectorScalar) :: rvecSolQ1, rvecRhsQ1, rvecDefQ1
    
    ! A real for the defect-correction-loop
    REAL(DP) :: dres
    
    ! Error of FE function to reference function
    REAL(DP) :: derror
    
    ! Output block for UCD output to GMV file
    TYPE(t_ucdExport) :: rexport
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata

    ! A simple counter variable
    INTEGER :: i
    
    ! Ok, let's start. 
    !
    ! We want to solve our Poisson problem on level...
    NLMIN = 2
    NLMAX = 6
    
    ! Allocate memory for all levels
    ALLOCATE(Rlevels(NLMIN:NLMAX))
    
    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    CALL boundary_read_prm(rboundary, './pre/QUAD.prm')
        
    ! Now read in the basic triangulation into our coarse level.
    CALL tria_readTriFile2D (Rlevels(NLMIN)%rtriangulation, &
                             './pre/QUAD.tri', rboundary)
    
    ! Refine it.
    CALL tria_quickRefine2LevelOrdering (NLMIN-1,&
        Rlevels(NLMIN)%rtriangulation,rboundary)
    
    ! And create information about adjacencies and everything one needs
    ! from a triangulation.
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
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we only have one block.
    ! Do this for all levels
    DO i = NLMIN, NLMAX
      CALL spdiscr_initBlockDiscr (Rlevels(i)%rdiscretisation, 1, &
                                   Rlevels(i)%rtriangulation, rboundary)
    END DO
    
    ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
    ! structures for every component of the solution vector.
    ! Initialise the first element of the list to specify the element
    ! and cubature rule for this solution component:
    DO i = NLMIN, NLMAX
      CALL spdiscr_initDiscr_simple (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1), &
          EL_EB30,CUB_G3X3,Rlevels(i)%rtriangulation, rboundary)
    END DO
                 
    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! We create a scalar matrix, based on the discretisation structure
    ! for our one and only solution component.
    DO i = NLMIN, NLMAX

      ! Initialise the block matrix with default values based on
      ! the discretisation.
      CALL lsysbl_createMatBlockByDiscr (&
          Rlevels(i)%rdiscretisation,Rlevels(i)%rmatrix)    

      ! Now as the discretisation is set up, we can start to generate
      ! the structure of the system matrix which is to solve.
      ! We create that directly in the block (1,1) of the block matrix
      ! using the discretisation structure of the first block.
      CALL bilf_createMatrixStructure ( &
           Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
           LSYSSC_MATRIX9,Rlevels(i)%rmatrix%RmatrixBlock(1,1))
      
      ! Update the structural information of the block matrix, as we manually
      ! changed one of the submatrices:
      CALL lsysbl_updateMatStrucInfo (Rlevels(i)%rmatrix)

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
      CALL bilf_buildMatrixScalar (rform,.TRUE.,&
           Rlevels(i)%rmatrix%RmatrixBlock(1,1),coeff_Laplace_2D)
    
    END DO
      
    ! Although we could manually create the solution/RHS vector,
    ! the easiest way to set up the vector structure is
    ! to create it by using our matrix as template:
    CALL lsysbl_createVecBlockIndMat (Rlevels(NLMAX)%rmatrix,rrhsBlock, .FALSE.)

    ! The vector structure is ready but the entries are missing. 
    ! So the next thing is to calculate the content of that vector.
    !
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    
    ! ... and then discretise the RHS to get a discrete version of it.
    ! Again we simply create a scalar vector based on the one and only
    ! discretisation structure.
    ! This scalar vector will later be used as the one and only first
    ! component in a block vector.
    CALL linf_buildVectorScalar (&
        Rlevels(NLMAX)%rdiscretisation%RspatialDiscr(1),&
        rlinform,.TRUE.,rrhsBlock%RvectorBlock(1),coeff_RHS_2D)
    
    DO i = NLMIN, NLMAX
    
      ! Initialise the discrete BC structure
      CALL bcasm_initDiscreteBC(Rlevels(i)%rdiscreteBC)

      ! On edge 1 of boundary component 1 add Dirichlet boundary conditions.      
      CALL boundary_createRegion(rboundary,1,1,rboundaryRegion)
      CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                        rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                        getBoundaryValues_2D)
                               
      ! Now to the edge 2 of boundary component 1 the domain. 
      CALL boundary_createRegion(rboundary,1,2,rboundaryRegion)
      CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                        rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                        getBoundaryValues_2D)
                               
      ! Edge 3 of boundary component 1.
      CALL boundary_createRegion(rboundary,1,3,rboundaryRegion)
      CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                        rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                        getBoundaryValues_2D)
      
      ! Edge 4 of boundary component 1. That's it.
      CALL boundary_createRegion(rboundary,1,4,rboundaryRegion)
      CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                        rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                        getBoundaryValues_2D)
      
      ! Hang the pointer into the matrix. That way, these
      ! boundary conditions are always connected to that matrix.
      Rlevels(i)%rmatrix%p_rdiscreteBC => Rlevels(i)%rdiscreteBC
  
      ! Also implement the boundary conditions into the matrix.
      CALL matfil_discreteBC (Rlevels(i)%rmatrix)
      
    END DO

    ! Our right-hand-side also needs to know the boundary conditions.
    rrhsBlock%p_rdiscreteBC => Rlevels(NLMAX)%rdiscreteBC
    
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
    
    ! During the linear solver, the boundary conditions are also
    ! frequently imposed to the vectors. But as the linear solver
    ! does not work with the actual solution vectors but with
    ! defect vectors instead.
    ! So, set up a filter chain that filters the defect vector
    ! during the solution process to implement discrete boundary conditions.
    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up the L2-projection for Multigrid
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Now we need to set up the (2-Level) mass matrices for all levels except
    ! for the coarse-most one.
    DO i = NLMIN+1, NLMAX
    
      ! Since the structure of the mass matrix is equal to the one of the
      ! Laplace matrix on the current level, we will simply create a shared
      ! copy of the structure.
      CALL lsyssc_duplicateMatrix (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
          Rlevels(i)%rmatrixMass, LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      
      ! Assemble the mass matrix
      CALL stdop_assembleSimpleMatrix(Rlevels(i)%rmatrixMass,&
                                      DER_FUNC, DER_FUNC)

      ! Now create the matrix structure of the 2-Level mass matrix.
      CALL mlop_create2LvlMatrixStruct(&
          Rlevels(i-1)%rdiscretisation%RspatialDiscr(1),&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
          LSYSSC_MATRIX9, Rlevels(i)%rmatrix2Lvl)
      
      ! And assemble the entries of the 2-Level mass matrix:
      CALL mlop_build2LvlMassMatrix (&
          Rlevels(i-1)%rdiscretisation%RspatialDiscr(1),&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
          .TRUE., Rlevels(i)%rmatrix2Lvl)
      
      ! Now set up an interlevel projecton structure for this level
      ! based on the Laplace matrix on this level.
      CALL mlprj_initProjectionMat (Rlevels(i)%rprojection,&
                                    Rlevels(i)%rmatrix)
      
      ! And initialise the L2-projection
      CALL mlprj_initL2Projection (&
          Rlevels(i)%rprojection%RscalarProjection(1,1),&
          Rlevels(i)%rmatrix2Lvl, Rlevels(i)%rmatrixMass)
      
      Rlevels(i)%rprojection%RscalarProjection(1,1)%depsL2 = 1e-20_DP
      
    END DO

    ! And set up an interlevel projecton structure for the coarse-most level.
    CALL mlprj_initProjectionMat (Rlevels(NLMIN)%rprojection,&
                                  Rlevels(NLMIN)%rmatrix)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! L2-projection for Multigrid is set up now
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Now we have to build up the level information for multigrid.
    !
    ! Create a Multigrid-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    p_RfilterChain => RfilterChain
    CALL linsol_initMultigrid (p_rsolverNode,p_RfilterChain)
    
    ! Set up a coarse grid solver.
    CALL linsol_initUMFPACK4 (p_rcoarsegridSolver)

    ! Add the coarse grid level.
    CALL linsol_addMultiGridLevel(p_rlevelInfo, p_rsolverNode,&
        Rlevels(NLMIN)%rprojection, NULL(), NULL(), p_rcoarseGridSolver)

    ! Now set up the other levels...
    DO i = NLMIN+1, NLMAX
    
      ! Create a Jacobi smoother
      CALL linsol_initJacobi(p_rsmoother)

      ! We will use 4 smoothing steps with damping parameter 0.7
      CALL linsol_convertToSmoother(p_rsmoother, 4, 0.7_DP)

      ! And add this multi-grid level. We will use the same smoother
      ! for pre- and post-smoothing.
      CALL linsol_addMultiGridLevel(p_rlevelInfo,p_rsolverNode,&
          Rlevels(i)%rprojection, p_rsmoother, p_rsmoother, NULL())
      
    END DO
    
    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2
    
    ! Set the coarse grid correction type to calculate an optimal
    ! damping parameter in repsect to the energy norm, as we apply
    ! the multigrid solver on a non-conforming element - this improves
    ! the multigrid convergence rate a bit.
    p_rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection%ccorrectionType = &
      CGCOR_SCALARENERGYMIN
    
    ! Attach the system matrices to the solver.
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
    
    CALL linsol_setMatrices(p_RsolverNode,Rmatrices(NLMIN:NLMAX))

    ! We can release Rmatrices immediately -- as long as we don't
    ! release Rlevels(i)%rmatrix!
    DO i=NLMIN,NLMAX
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
      
    ! That's it, rvectorBlock now contains our solution.

    ! Calculate the error to the reference function.
    CALL pperr_scalar (rvectorBlock%RvectorBlock(1),PPERR_L2ERROR,derror,&
                       getReferenceFunction_2D)
    CALL output_line ('L2-error: ' // sys_sdEL(derror,10) )

    CALL pperr_scalar (rvectorBlock%RvectorBlock(1),PPERR_H1ERROR,derror,&
                       getReferenceFunction_2D)
    CALL output_line ('H1-error: ' // sys_sdEL(derror,10) )
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Perform L2-projection of solution into Q1 space
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    CALL output_lbrk()
    CALL output_line('Performing L2-projection of solution to Q1 space')
    CALL output_lbrk()

    ! We now have the solution vector, but unfortunately, it's a Q1~ solution
    ! vector and what we need are the function values in the vertices of the
    ! mesh. Instead of calling the spatial-projection module to interpolate
    ! the solution we will perform a "true" L2-projection of the Q1~ solution
    ! into the Q1 space, i.e. we want to "solve"
    !
    !                   v_h = u_h                                    (1)
    !
    ! with v_h in Q1 and u_h being our actual discrete Q1~ solution of the
    ! Poisson equation. So as a weak formulation of (1) we get
    !
    !          (v_h, phi_i) = (u_h, phi_i)  for all 1 <= i <= m      (2)
    !
    ! with (phi_1,...,phi_m) being the basis functions of Q1.
    ! Furthermore, we know that
    !
    !                          n
    !                   u_h = sum (x_k * psi_k)                      (3)
    !                         k=1
    !
    ! with (psi_1,...,psi_n) being the basis functions of Q1~ and
    ! x = (x_1,...,x_n) being our actual solution vector rvectorBlock.
    ! Now the discrete Q1 function v_h we search for also has a coefficient
    ! vector y = (y_1,...,y_m) and can be written as
    !
    !                          m
    !                   v_h = sum (y_j * phi_j)                      (4)
    !                         j=1
    !
    ! Inserting (3) and (4) into (2) and rewriting the result as a
    ! linear system we get
    !
    !                 M * y = N * x                                  (5)
    !
    !     for all 1 <= i,j <= m ; 1 <= k <= n :
    !
    !                 M_ij := (phi_j, phi_i)
    !                 N_ik := (psi_k, phi_i)
    !
    ! So there are 2 matrices we need to build: the Q1 mass matrix (M) and a
    ! mass matrix with Q1~ being its trial and Q1 being its test space (N).
    ! Afterwards we solve (5) and get the coefficient vector y of our
    ! Q1 solution v_h...
    !
    ! The first thing that we need is a Q1 discretisation on the fine mesh.
    CALL spdiscr_initDiscr_simple(rdiscrQ1, EL_Q1, CUB_G3X3, &
                                  Rlevels(NLMAX)%rtriangulation, rboundary)

    ! Now create the the matrix structure of N.
    ! The trial space is EB30 and the test space is Q1:
    CALL bilf_createMatrixStructure (&
        Rlevels(NLMAX)%rdiscretisation%RspatialDiscr(1),&
        LSYSSC_MATRIX9, rmatrixMassPrj, rdiscrQ1)

    ! And assemble the mass matrix entries of N:
    CALL stdop_assembleSimpleMatrix(rmatrixMassPrj, DER_FUNC, DER_FUNC)

    ! We need to create the Q1 coefficient vector y:
    CALL lsyssc_createVecByDiscr (rdiscrQ1, rvecSolQ1, .TRUE.)
    
    ! And we also need a Q1 rhs vector r which recieves N*x
    CALL lsyssc_createVecByDiscr (rdiscrQ1, rvecRhsQ1, .FALSE.)
    
    ! Calculate r := N*x
    CALL lsyssc_scalarMatVec(rmatrixMassPrj, rvectorBlock%rvectorBlock(1),&
                             rvecRhsQ1, 1.0_DP, 0.0_DP)
    
    ! At this point we won't need the matrix N anymore, as we just needed
    ! it to get a rhs vector for our Q1 mass system - so we'll release it now.
    CALL lsyssc_releaseMatrix(rmatrixMassPrj)

    ! As rmatrixMassPrj is free now, we will use it to store M
    CALL bilf_createMatrixStructure (rdiscrQ1,LSYSSC_MATRIX9,rmatrixMassPrj)
    CALL stdop_assembleSimpleMatrix(rmatrixMassPrj, DER_FUNC, DER_FUNC)
    
    ! We now have assembled the Q1 mass matrix M, created a Q1 solution
    ! vector y and calculated the right-hand-side vector r := N*x.
    ! Basically, we could now create a linear solver as CG or whatever, but
    ! this is not necessary as a system with a mass matrix is quite easy
    ! to solve (in contrast to a system with a Laplace matrix).
    ! We will solve the linear system by a Defect-Correction-Approach
    ! using the inverse "lumped" mass matrix as a preconditioner:
    !
    !                          -1
    !          y_k+1 := y_k + L   * (r - M * y_k)
    !
    !                    m
    !           L_ii := sum M_ij         L_ij := 0 for i /= j
    !                   j=1
    !
    ! This simple Defect-Correction-Loop converges quite fast (at least as
    ! long as L is regular, that is ^_^) so there's no need to seek for
    ! solvers as CG or even multigrid here...
    !
    ! Remark:
    ! Basically, we could also discretise the boundary conditions of our
    ! Poisson equation in Q1 and implement them into the mass matrix M.
    ! However, this is not crucial as the mass matrix M is regular even
    ! without any boundary conditions and the violation of the boundary
    ! conditions in our resulting Q1 solution will be tolerable...
    !
    ! So let's create a "lumped" mass matrix first...
    
    ! Create a copy of M
    CALL lsyssc_duplicateMatrix(rmatrixMassPrj, rlumpedMassPrj,&
                                LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)
    
    ! And lump the copy to get L
    CALL lsyssc_lumpMatrixScalar(rlumpedMassPrj, LSYSSC_LUMP_DIAG)
    
    ! Furthermore we need a defect vector
    CALL lsyssc_createVecIndMat(rmatrixMassPrj, rvecDefQ1)
    
    ! And let's start the Defect-Correction-Loop
    DO i = 1, 50
    
      ! Calculate current defect d_i := r - M * y_i
      CALL lsyssc_copyVector(rvecRhsQ1, rvecDefQ1)
      CALL lsyssc_scalarMatVec(rmatrixMassPrj, rvecSolQ1, rvecDefQ1,&
                               -1.0_DP, 1.0_DP)
      
      ! Calculate the current residual = || d_i ||_2
      dres = lsyssc_vectorNorm(rvecDefQ1, LINALG_NORML2)
      
      ! Print the residual to screen
      CALL output_line('L2prj: Iteration ' // TRIM(sys_siL(i,10)) // &
                       ',  !!RES!! = ' // TRIM(sys_sdEL(dres,15)))
      
      ! Is our L2-projection good enough?
      IF (dres .LE. 1e-7_DP) EXIT
      
      ! Otherwise multiply the defect by the inverse of L
      CALL lsyssc_invertedDiagMatVec(rlumpedMassPrj, rvecDefQ1, 1.0_DP,&
                                     rvecDefQ1)
      
      ! And add it onto y_i
      CALL lsyssc_vectorLinearComb(rvecDefQ1, rvecSolQ1, 1.0_DP, 1.0_DP)
    
    END DO

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! L2-projection of solution into Q1 space done
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Start UCD export to GMV file:
    CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,&
        Rlevels(NLMAX)%rtriangulation,'gmv/u2d_1_l2prj.gmv')
    
    ! Add our Q1-solution to the UCD exporter:
    CALL lsyssc_getbase_double (rvecSolQ1,p_Ddata)
    CALL ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that's it.
    CALL ucd_write (rexport)
    CALL ucd_release (rexport)
    
    ! We are finished - but not completely!
    ! Now, clean up so that all the memory is available again.
    !
    ! Release the three Q1 vectors
    CALL lsyssc_releaseVector(rvecDefQ1)
    CALL lsyssc_releaseVector(rvecRhsQ1)
    CALL lsyssc_releaseVector(rvecSolQ1)
    
    ! Release the Q1 mass matrix and its lumped version
    CALL lsyssc_releaseMatrix(rlumpedMassPrj)
    CALL lsyssc_releaseMatrix(rmatrixMassPrj)
    
    ! And release the Q1 discretisation
    CALL spdiscr_releaseDiscr(rdiscrQ1)
    
    ! That was all we used in the L2-projection of the solution and
    ! that has not been already released before.
    
    ! Release solver data and structure
    CALL linsol_doneData (p_rsolverNode)
    CALL linsol_doneStructure (p_rsolverNode)

    ! Release the solver node and all subnodes attached to it (if at all):
    CALL linsol_releaseSolver (p_rsolverNode)
    
    ! Release the L2-projection
    DO i = NLMAX, NLMIN+1, -1
    
      ! Release the projection structure itself
      CALL mlprj_doneProjection(Rlevels(i)%rprojection)
      
      ! Release the 2-Level-Mass matrix
      CALL lsyssc_releaseMatrix (Rlevels(i)%rmatrix2Lvl)
      
      ! Release the mass matrix
      CALL lsyssc_releaseMatrix (Rlevels(i)%rmatrixMass)

    END DO
    
    ! Release the projection structure on the coarse mesh
    CALL mlprj_doneProjection(Rlevels(NLMIN)%rprojection)
        
    ! Release the block matrix/vectors
    CALL lsysbl_releaseVector (rtempBlock)
    CALL lsysbl_releaseVector (rvectorBlock)
    CALL lsysbl_releaseVector (rrhsBlock)
    DO i = NLMAX, NLMIN, -1
      CALL lsysbl_releaseMatrix (Rlevels(i)%rmatrix)
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
