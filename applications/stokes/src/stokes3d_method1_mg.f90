!##############################################################################
!# ****************************************************************************
!# <name> stokes3d_method1_mg </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstation program how to solve a Stokes
!# problem on a simple domain.
!#
!# The routine uses the simple-VANKA smoother for 3D saddle point problems,
!# Jacobi-Type, for a multigrid solver.
!# </purpose>
!##############################################################################

MODULE stokes3d_method1_mg

  USE fsystem
  USE storage
  USE linearsolver
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
  USE meshregion
  
  USE stokes3d_callback
  USE dom3d_cube
  
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

    ! B1-matrix for that specific level. 
    TYPE(t_matrixScalar) :: rmatrixB1

    ! B2-matrix for that specific level. 
    TYPE(t_matrixScalar) :: rmatrixB2
    
    ! B3-matrix for that specific level.
    TYPE(t_matrixScalar) :: rmatrixB3

    ! A variable describing the discrete boundary conditions.    
    TYPE(t_discreteBC) :: rdiscreteBC
  
  END TYPE
  
!</typeblock>

!</types>

CONTAINS
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE stokes3d_1_mg
  
 
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

    ! An object for saving the boundary mesh region
    TYPE(t_meshregion) :: rmeshRegion

    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    TYPE(t_blockDiscretisation) :: rprjDiscretisation
    
    ! A bilinear and linear form describing the analytic problem to solve
    TYPE(t_bilinearForm) :: rform
    TYPE(t_linearForm) :: rlinform

    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    TYPE(t_vectorBlock) :: rvector,rrhs,rtempBlock,rprjVector
    
    ! A set of variables describing the analytic and discrete boundary
    ! conditions.    
    TYPE(t_discreteBC), TARGET :: rprjDiscreteBC

    ! A solver node that accepts parameters for the linear solver    
    TYPE(t_linsolNode), POINTER :: p_rsolverNode,p_rpreconditioner,&
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
    
    ! Error indicator during initialisation of the solver
    INTEGER :: ierror    
    
    ! Output block for UCD output to GMV file
    TYPE(t_ucdExport) :: rexport
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_Ddata2,p_Ddata3

    ! A counter variable
    INTEGER :: i

    ! Ok, let's start. 
    !
    ! We want to solve our Stokes problem on level...
    NLMIN = 2
    NLMAX = 4
    
    ! Viscosity parameter:
    dnu = 1.0_DP

    ! Allocate memory for all levels
    ALLOCATE(Rlevels(NLMIN:NLMAX))

    ! At first read in the basic triangulation.
    CALL tria_readTriFile3D (Rlevels(NLMIN)%rtriangulation, &
                             './pre/CUBE.tri')
    
    ! Refine the mesh up to the minimum level
    CALL tria_quickRefine2LevelOrdering (NLMIN-1,&
        Rlevels(NLMIN)%rtriangulation)
    
    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    CALL tria_initStandardMeshFromRaw (Rlevels(NLMIN)%rtriangulation)
    
    ! Now refine the grid for the fine levels.
    DO i = NLMIN+1, NLMAX

      ! Refine the grid using the 2-Level-Ordering algorithm
      CALL tria_refine2LevelOrdering(Rlevels(i-1)%rtriangulation,&
          Rlevels(i)%rtriangulation)
      
      ! Create a standard mesh
      CALL tria_initStandardMeshFromRaw(Rlevels(i)%rtriangulation)
    
    END DO

    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies 4 blocks in the
    ! solution vector.
    DO i = NLMIN, NLMAX
      CALL spdiscr_initBlockDiscr (Rlevels(i)%rdiscretisation, 4, &
                                   Rlevels(i)%rtriangulation)
    END DO

    ! rdiscretisation%RspatialDiscr is a list of scalar 
    ! discretisation structures for every component of the solution vector.
    ! We have a solution vector with three components:
    !  Component 1 = X-velocity
    !  Component 2 = Y-velocity
    !  Component 3 = Z-velocity
    !  Component 4 = Pressure
    DO i = NLMIN, NLMAX
      ! For simplicity, we set up one discretisation structure for the 
      ! velocity...
      CALL spdiscr_initDiscr_simple (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
          EL_EM30_3D, CUB_G2_3D, Rlevels(i)%rtriangulation)
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd and 3rd component (Y-/Z-velocity). This needs no
      ! additional memory, as all three structures will share the same dynamic
      ! information afterwards.
      CALL spdiscr_duplicateDiscrSc (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
          Rlevels(i)%rdiscretisation%RspatialDiscr(2))
      CALL spdiscr_duplicateDiscrSc (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
          Rlevels(i)%rdiscretisation%RspatialDiscr(3))

      ! For the pressure (4th component), we set up a separate discretisation 
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      CALL spdiscr_deriveSimpleDiscrSc (Rlevels(i)%rdiscretisation%RspatialDiscr(1), &
          EL_Q0_3D, CUB_G3_3D, Rlevels(i)%rdiscretisation%RspatialDiscr(4))
    
    END DO

    DO i = NLMIN, NLMAX
    
      ! Initialise the block matrix with default values based on
      ! the discretisation.
      CALL lsysbl_createMatBlockByDiscr (Rlevels(i)%rdiscretisation,&
                                         Rlevels(i)%rmatrix)    
      
      ! Inform the matrix that we build a saddle-point problem.
      ! Normally, imatrixSpec has the value LSYSBS_MSPEC_GENERAL,
      ! but probably some solvers can use the special structure later.
      Rlevels(i)%rmatrix%imatrixSpec = LSYSBS_MSPEC_SADDLEPOINT
      
      ! Now as the discretisation is set up, we can start to generate
      ! the structure of the system matrix which is to solve.
      ! We create that directly in the block (1,1) of the block matrix
      ! using the discretisation structure of the first block.
      !
      ! Create the matrix structure of the X-velocity.
      CALL bilf_createMatrixStructure (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
          LSYSSC_MATRIX9, Rlevels(i)%rmatrix%RmatrixBlock(1,1))

      ! In the Stokes problem, the matrices for the Y-velocity and Z-velocity
      ! are identical to the matrix for the X-velocity; all three are
      ! Laplace-matrices!
      ! Therefore, we can simply make a copy of the matrix for the X-velocity.
      ! This we do later after the entries are created.
      !
      ! In the global system, there are three coupling matrices B1, B2 and B3.
      ! All three have the same structure.
      !
      !    / A              B1 \
      !    |      A         B2 |
      !    |           A    B3 |
      !    \ B1^T B2^T B3^T    /
      !
      ! Create the matrices structure of the pressure using the 4th
      ! spatial discretisation structure in p_rdiscretisation%RspatialDiscr.
      CALL bilf_createMatrixStructure (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(4),&
          LSYSSC_MATRIX9, Rlevels(i)%rmatrixB1,&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1))
                
      ! Duplicate the B1 matrix structure to the B2/B3 matrix, so use
      ! lsyssc_duplicateMatrix to create B2/B3. Share the matrix 
      ! structure between B1 and B2/B3 (B1 is the parent and B2/B3 the children). 
      ! Don't create a content array yet, it will be created by 
      ! the assembly routines later.
      CALL lsyssc_duplicateMatrix (Rlevels(i)%rmatrixB1, Rlevels(i)%rmatrixB2,&
                                   LSYSSC_DUP_COPY, LSYSSC_DUP_REMOVE)
      CALL lsyssc_duplicateMatrix (Rlevels(i)%rmatrixB1, Rlevels(i)%rmatrixB3,&
                                   LSYSSC_DUP_COPY, LSYSSC_DUP_REMOVE)
                                       
      ! And now to the entries of the matrix. For assembling of the entries,
      ! we need a bilinear form, which first has to be set up manually.
      ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
      ! scalar system matrix in 3D.
      rform%itermCount = 3
      rform%Idescriptors(1,1) = DER_DERIV3D_X
      rform%Idescriptors(2,1) = DER_DERIV3D_X
      rform%Idescriptors(1,2) = DER_DERIV3D_Y
      rform%Idescriptors(2,2) = DER_DERIV3D_Y
      rform%Idescriptors(1,3) = DER_DERIV3D_Z
      rform%Idescriptors(2,3) = DER_DERIV3D_Z

      ! In the standard case, we have constant coefficients:
      rform%ballCoeffConstant = .TRUE.
      rform%BconstantCoeff = .TRUE.
      rform%Dcoefficients(1)  = dnu
      rform%Dcoefficients(2)  = dnu
      rform%Dcoefficients(3)  = dnu

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
      CALL bilf_buildMatrixScalar (rform,.TRUE.,&
          Rlevels(i)%rmatrix%RmatrixBlock(1,1), coeff_Stokes_3D)
      
      ! Duplicate the matrix to the Y-/Z-velocity matrix, share structure and
      ! content between them (as the matrices are the same).
      CALL lsyssc_duplicateMatrix (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
          Rlevels(i)%rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      CALL lsyssc_duplicateMatrix (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
          Rlevels(i)%rmatrix%RmatrixBlock(3,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      
      ! Manually change the discretisation structure of the Y-/Z-velocity 
      ! matrix to the Y-/Z-discretisation structure.
      ! Ok, we use the same discretisation structure for both, X- and Y-velocity,
      ! so this is not really necessary - we do this for sure...
      CALL lsyssc_assignDiscretDirectMat (Rlevels(i)%rmatrix%RmatrixBlock(2,2),&
          Rlevels(i)%rdiscretisation%RspatialDiscr(2))
      CALL lsyssc_assignDiscretDirectMat (Rlevels(i)%rmatrix%RmatrixBlock(3,3),&
          Rlevels(i)%rdiscretisation%RspatialDiscr(3))
                                  
      ! Build the first pressure matrix B1.
      ! Again first set up the bilinear form, then call the matrix assembly.
      rform%itermCount = 1
      rform%Idescriptors(1,1) = DER_FUNC3D
      rform%Idescriptors(2,1) = DER_DERIV3D_X

      ! In the standard case, we have constant coefficients:
      rform%ballCoeffConstant = .TRUE.
      rform%BconstantCoeff = .TRUE.
      rform%Dcoefficients(1)  = -1.0_DP
      
      CALL bilf_buildMatrixScalar (rform,.TRUE.,Rlevels(i)%rmatrixB1,&
                                   coeff_Pressure_3D)

      ! Build the second pressure matrix B2.
      ! Again first set up the bilinear form, then call the matrix assembly.
      rform%itermCount = 1
      rform%Idescriptors(1,1) = DER_FUNC3D
      rform%Idescriptors(2,1) = DER_DERIV3D_Y

      ! In the standard case, we have constant coefficients:
      rform%ballCoeffConstant = .TRUE.
      rform%BconstantCoeff = .TRUE.
      rform%Dcoefficients(1)  = -1.0_DP
      
      CALL bilf_buildMatrixScalar (rform,.TRUE.,Rlevels(i)%rmatrixB2,&
                                   coeff_Pressure_3D)
                                  
      ! Build the third pressure matrix B3.
      ! Again first set up the bilinear form, then call the matrix assembly.
      rform%itermCount = 1
      rform%Idescriptors(1,1) = DER_FUNC3D
      rform%Idescriptors(2,1) = DER_DERIV3D_Z

      ! In the standard case, we have constant coefficients:
      rform%ballCoeffConstant = .TRUE.
      rform%BconstantCoeff = .TRUE.
      rform%Dcoefficients(1)  = -1.0_DP
      
      CALL bilf_buildMatrixScalar (rform,.TRUE.,Rlevels(i)%rmatrixB3,&
                                   coeff_Pressure_3D)

      ! The B1/B2/B3 matrices exist up to now only in our local problem structure.
      ! Put a copy of them into the block matrix.
      !
      ! Note that we share the structure of B1/B2/B3 with those B1/B2/B3 of the
      ! block matrix, while we create copies of the entries. The reason is
      ! that these matrices are modified for boundary conditions later.
      CALL lsyssc_duplicateMatrix (Rlevels(i)%rmatrixB1, &
          Rlevels(i)%rmatrix%RmatrixBlock(1,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

      CALL lsyssc_duplicateMatrix (Rlevels(i)%rmatrixB2, &
          Rlevels(i)%rmatrix%RmatrixBlock(2,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      
      CALL lsyssc_duplicateMatrix (Rlevels(i)%rmatrixB3, &
          Rlevels(i)%rmatrix%RmatrixBlock(3,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

      ! Furthermore, put B1^T, B2^T and B3^T to the block matrix.
      CALL lsyssc_transposeMatrix (Rlevels(i)%rmatrixB1, &
          Rlevels(i)%rmatrix%RmatrixBlock(4,1),LSYSSC_TR_VIRTUAL)

      CALL lsyssc_transposeMatrix (Rlevels(i)%rmatrixB2, &
          Rlevels(i)%rmatrix%RmatrixBlock(4,2),LSYSSC_TR_VIRTUAL)

      CALL lsyssc_transposeMatrix (Rlevels(i)%rmatrixB3, &
          Rlevels(i)%rmatrix%RmatrixBlock(4,3),LSYSSC_TR_VIRTUAL)

      ! Update the structural information of the block matrix, as we manually
      ! changed the submatrices:
      CALL lsysbl_updateMatStrucInfo (Rlevels(i)%rmatrix)
    
    END DO

    ! Although we could manually create the solution/RHS vector,
    ! the easiest way to set up the vector structure is
    ! to create it by using our matrix as template:
    CALL lsysbl_createVecBlockIndMat (Rlevels(NLMAX)%rmatrix,rrhs,.FALSE.)
    CALL lsysbl_createVecBlockIndMat (Rlevels(NLMAX)%rmatrix,rvector,.FALSE.)

    ! The vector structure is ready but the entries are missing. 
    ! So the next thing is to calculate the content of that vector.
    !
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC3D
    
    ! ... and then discretise the RHS to the first three subvectors of
    ! the block vector using the discretisation structure of the 
    ! corresponding block.
    !
    ! Note that the vector is unsorted after calling this routine!
    CALL linf_buildVectorScalar (&
        Rlevels(NLMAX)%rdiscretisation%RspatialDiscr(1),&
        rlinform,.TRUE.,rrhs%RvectorBlock(1),coeff_RHS_X_3D)

    CALL linf_buildVectorScalar (&
        Rlevels(NLMAX)%rdiscretisation%RspatialDiscr(2),&
        rlinform,.TRUE.,rrhs%RvectorBlock(2),coeff_RHS_X_3D)

    CALL linf_buildVectorScalar (&
        Rlevels(NLMAX)%rdiscretisation%RspatialDiscr(3),&
        rlinform,.TRUE.,rrhs%RvectorBlock(3),coeff_RHS_X_3D)
                                
    ! The fourth subvector must be zero - as it represents the RHS of
    ! the equation "div(u) = 0".
    CALL lsyssc_clearVector(rrhs%RvectorBlock(4))
                                
    ! Clear the solution vector on the finest level.
    CALL lsysbl_clearVector(rvector)

    ! Now set up the boundary conditions
    DO i = NLMIN, NLMAX

      ! Now we need to implement the boundary conditions. To do this, we
      ! first need to create a mesh region describing the mesh's boundary.
      ! We want to prescribe Dirichlet on the cube's boundary, except for
      ! the face where the X-coordinate is 1. 
      ! We could now manually create a mesh region based on the triangulation's
      ! nodal-property array and then kick out everything that belongs to the
      ! right face. But we will use the dom3d_cube module, which performs
      ! this task for us.
      CALL dom3d_cube_calcMeshRegion(rmeshRegion, Rlevels(i)%rtriangulation, &
                                     DOM3D_CUBE_REG_STOKES)

      ! Initialise the structure for discrete boundary conditions.
      CALL bcasm_initDiscreteBC (Rlevels(i)%rdiscreteBC)

      ! Set Dirichlet BCs for X-velocity:
      CALL bcasm_newDirichletBConMR(Rlevels(i)%rdiscretisation, 1, &
          Rlevels(i)%rdiscreteBC, rmeshRegion, &
          getBoundaryValuesMR_3D)

      ! Set Dirichlet BCs for Y-velocity:
      CALL bcasm_newDirichletBConMR(Rlevels(i)%rdiscretisation, 2, &
          Rlevels(i)%rdiscreteBC, rmeshRegion, &
          getBoundaryValuesMR_3D)

      ! Set Dirichlet BCs for Z-velocity:
      CALL bcasm_newDirichletBConMR(Rlevels(i)%rdiscretisation, 3, &
          Rlevels(i)%rdiscreteBC, rmeshRegion, &
          getBoundaryValuesMR_3D)

      ! Hang the pointer into the vector and matrix. That way, these
      ! boundary conditions are always connected to that matrix and that
      ! vector.
      Rlevels(i)%rmatrix%p_rdiscreteBC => Rlevels(i)%rdiscreteBC
      
      ! Next step is to implement boundary conditions into the matrix.
      ! This is done using a matrix filter for discrete boundary conditions.
      ! The discrete boundary conditions are already attached to the
      ! matrix. Call the appropriate matrix filter that modifies the matrix
      ! according to the boundary conditions.
      CALL matfil_discreteBC (Rlevels(i)%rmatrix)
      
      ! Don't forget to release the mesh region
      CALL mshreg_done(rmeshRegion)
    
    END DO

    ! Also implement the discrete boundary conditions on the finest level
    ! onto our right-hand-side and solution vectors.
    rrhs%p_rdiscreteBC => Rlevels(NLMAX)%rdiscreteBC
    rvector%p_rdiscreteBC => Rlevels(NLMAX)%rdiscreteBC
    CALL vecfil_discreteBCrhs (rrhs)
    CALL vecfil_discreteBCsol (rvector)

    ! Create a temporary vector we need that for some preparation.
    CALL lsysbl_createVecBlockIndirect (rrhs, rtempBlock, .FALSE.)

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
    CALL linsol_initMultigrid (p_rsolverNode,p_RfilterChain)
    
    ! Set up a BiCGStab solver with VANKA preconditioning as coarse grid solver:
    CALL linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_3DNAVST)
    CALL linsol_initBiCGStab (p_rcoarseGridSolver,p_rpreconditioner,p_RfilterChain)
    
    ! Add the coarse grid level.
    CALL linsol_addMultiGridLevel(p_rlevelInfo,p_rsolverNode,rprojection,&
                                  NULL(), NULL(), p_rcoarseGridSolver)

    ! Now set up the other levels...
    DO i = NLMIN+1, NLMAX
    
      ! Set up the general VANKA smoother.
      CALL linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_3DNAVST)
      
      ! We will use 4 smoothing steps with damping parameter 0.7
      CALL linsol_convertToSmoother(p_rsmoother, 4, 0.7_DP)
      
      ! And add this multi-grid level. We will use the same smoother
      ! for pre- and post-smoothing.
      CALL linsol_addMultiGridLevel(p_rlevelInfo,p_rsolverNode,rprojection,&
                                    p_rsmoother, p_rsmoother, NULL())
      
    END DO
    
    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2

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
    
    CALL spdiscr_duplicateBlockDiscr (Rlevels(NLMAX)%rdiscretisation,&
                                      rprjDiscretisation)
    
    CALL spdiscr_deriveSimpleDiscrSc (&
        Rlevels(NLMAX)%rdiscretisation%RspatialDiscr(1), &
        EL_Q1_3D, CUB_G2_3D, rprjDiscretisation%RspatialDiscr(1))
    CALL spdiscr_deriveSimpleDiscrSc (&
        Rlevels(NLMAX)%rdiscretisation%RspatialDiscr(2), &
        EL_Q1_3D, CUB_G2_3D, rprjDiscretisation%RspatialDiscr(2))
    CALL spdiscr_deriveSimpleDiscrSc (&
        Rlevels(NLMAX)%rdiscretisation%RspatialDiscr(3), &
        EL_Q1_3D, CUB_G2_3D, rprjDiscretisation%RspatialDiscr(3))
                 
    ! The pressure discretisation substructure stays the old.
    !
    ! Now set up a new solution vector based on this discretisation,
    ! allocate memory.
    CALL lsysbl_createVecBlockByDiscr (rprjDiscretisation,rprjVector,.FALSE.)
    
    ! Then take our original solution vector and convert it according to the
    ! new discretisation:
    CALL spdp_projectSolution (rvector,rprjVector)
    
    ! Discretise the boundary conditions according to the Q1/Q1/Q1/Q0 
    ! discretisation:
    
    ! Get the boundary mesh region on the finest level
    CALL dom3d_cube_calcMeshRegion(rmeshRegion, Rlevels(NLMAX)%rtriangulation, &
                                   DOM3D_CUBE_REG_STOKES)

    ! And prescribe Dirichlet boundary conditions
    CALL bcasm_initDiscreteBC (rprjDiscreteBC)
    CALL bcasm_newDirichletBConMR(rprjDiscretisation, 1, rprjDiscreteBC, &
                       rmeshRegion, getBoundaryValuesMR_3D)
    CALL bcasm_newDirichletBConMR(rprjDiscretisation, 2, rprjDiscreteBC, &
                       rmeshRegion, getBoundaryValuesMR_3D)
    CALL bcasm_newDirichletBConMR(rprjDiscretisation, 3, rprjDiscreteBC, &
                       rmeshRegion, getBoundaryValuesMR_3D)

    ! Now we don't need the mesh region anymore, so release it
    CALL mshreg_done(rmeshRegion)

    ! Connect the vector to the BC's
    rprjVector%p_rdiscreteBC => rprjDiscreteBC
    
    ! Send the vector to the boundary-condition implementation filter.
    ! This modifies the vector according to the attached discrete boundary
    ! conditions.
    CALL vecfil_discreteBCsol (rprjVector)
    
    ! Now we have a Q1/Q1/Q0 solution in rprjVector.
    ! We can now start the postprocessing. 
    ! Start UCD export to GMV file:
    CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,&
        Rlevels(NLMAX)%rtriangulation,'gmv/u3d_1_mg.gmv')

    ! Write velocity field
    CALL lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
    CALL lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
    CALL lsyssc_getbase_double (rprjVector%RvectorBlock(3),p_Ddata3)
    
    ! In case we use the VTK exporter, which supports vector output, we will
    ! pass the X-,Y- and Z-velocity at once to the ucd module.
    CALL ucd_addVarVertBasedVec(rexport,'velocity',p_Ddata,p_Ddata2,p_Ddata3)

    ! If we use the GMV exporter, we might replace the line above by the
    ! following two lines:
    !CALL ucd_addVariableVertexBased (rexport,'X-vel',UCD_VAR_XVELOCITY, p_Ddata)
    !CALL ucd_addVariableVertexBased (rexport,'Y-vel',UCD_VAR_YVELOCITY, p_Ddata2)
    !CALL ucd_addVariableVertexBased (rexport,'Z-vel',UCD_VAR_ZVELOCITY, p_Ddata3)
        
    ! Write pressure
    CALL lsyssc_getbase_double (rprjVector%RvectorBlock(4),p_Ddata)
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
    
    ! Release the interlevel projection structure
    CALL mlprj_doneProjection (rprojection)

    ! Release the solver node and all subnodes attached to it (if at all):
    CALL linsol_releaseSolver (p_rsolverNode)
    
    ! Release the block matrix/vectors
    CALL lsysbl_releaseVector (rprjVector)
    CALL lsysbl_releaseVector (rtempBlock)
    CALL lsysbl_releaseVector (rvector)
    CALL lsysbl_releaseVector (rrhs)
    DO i = NLMAX, NLMIN, -1
      CALL lsysbl_releaseMatrix (Rlevels(i)%rmatrix)
    END DO
    
    ! Release B1 and B2 matrix
    DO i = NLMAX, NLMIN, -1
      CALL lsyssc_releaseMatrix (Rlevels(i)%rmatrixB3)
      CALL lsyssc_releaseMatrix (Rlevels(i)%rmatrixB2)
      CALL lsyssc_releaseMatrix (Rlevels(i)%rmatrixB1)
    END DO
    
    ! Release our discrete version of the boundary conditions
    CALL bcasm_releaseDiscreteBC (rprjDiscreteBC)
    DO i = NLMAX, NLMIN, -1
      CALL bcasm_releaseDiscreteBC (Rlevels(i)%rdiscreteBC)
    END DO

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    CALL spdiscr_releaseBlockDiscr(rprjDiscretisation)
    DO i = NLMAX, NLMIN, -1
      CALL spdiscr_releaseBlockDiscr(Rlevels(i)%rdiscretisation)
    END DO
    
    ! Release the triangulation. 
    DO i = NLMAX, NLMIN, -1
      CALL tria_done (Rlevels(i)%rtriangulation)
    END DO
    
    DEALLOCATE(Rlevels)
    
  END SUBROUTINE

END MODULE