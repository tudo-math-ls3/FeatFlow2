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

module stokes3d_method1_mg

  use fsystem
  use genoutput
  use storage
  use cubature
  use derivatives
  use matrixfilters
  use vectorfilters
  use discretebc
  use bcassembly
  use triangulation
  use element
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use scalarpde
  use linearsystemscalar
  use linearsystemblock
  use multilevelprojection
  use bilinearformevaluation
  use linearformevaluation
  use filtersupport
  use linearsolver
  use ucd
  use meshregion
  
  use stokes3d_callback
  use dom3d_cube
  
  implicit none

!<types>

!<typeblock description="Type block defining all information about one level">

  type t_level
  
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! Cubature info structure which encapsules the cubature formula
    type(t_scalarCubatureInfo) :: rcubatureInfo

    ! An object specifying the discretisation (structure of the
    ! solution, trial/test functions,...)
    type(t_blockDiscretisation) :: rdiscretisation
    
    ! A system matrix for that specific level. The matrix will receive the
    ! discrete Laplace operator.
    type(t_matrixBlock) :: rmatrix

    ! B1-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixB1

    ! B2-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixB2
    
    ! B3-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixB3

    ! A variable describing the discrete boundary conditions.
    type(t_discreteBC) :: rdiscreteBC
  
    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    type(t_filterChain), dimension(1) :: RfilterChain
    
    ! Number of filters in the filter chain
    integer :: nfilters

  end type
  
!</typeblock>

!</types>

contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine stokes3d_1_mg
  
 
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
    ! We need a couple of variables for this problem. Let us see...
    !
    ! An array of problem levels for the multigrid solver
    type(t_level), dimension(:), pointer :: Rlevels

    ! An object for saving the boundary mesh region
    type(t_meshregion) :: rmeshRegion

    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    type(t_blockDiscretisation) :: rprjDiscretisation
    
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform

    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_vectorBlock) :: rvector,rrhs,rtempBlock,rprjVector
    
    ! A set of variables describing the analytic and discrete boundary
    ! conditions.
    type(t_discreteBC), target :: rprjDiscreteBC

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode,p_rpreconditioner,&
                                   p_rcoarseGridSolver,p_rsmoother

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(:), pointer :: Rmatrices

    ! One level of multigrid
    type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo
    
    ! NLMIN receives the level of the coarse grid.
    integer :: NLMIN

    ! NLMAX receives the level where we want to solve.
    integer :: NLMAX
    
    ! Viscosity parameter nu = 1/Re
    real(DP) :: dnu
    
    ! Error indicator during initialisation of the solver
    integer :: ierror
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir
    real(DP), dimension(:), pointer :: p_Ddata,p_Ddata2,p_Ddata3

    ! A counter variable
    integer :: i

    ! Path to the mesh
    character(len=SYS_STRLEN) :: spredir
    
    ! Ok, let us start.
    !
    ! We want to solve our Stokes problem on level...
    NLMIN = 2
    NLMAX = 4
    
    ! Viscosity parameter:
    dnu = 1.0_DP

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Read the mesh, refine
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Allocate memory for all levels
    allocate(Rlevels(NLMIN:NLMAX))

    ! Get the path $PREDIR from the environment, where to read .prm/.tri files
    ! from. If that does not exist, write to the directory "./pre".
    if (.not. sys_getenv_string("PREDIR", spredir)) spredir = "./pre"

    ! At first read in the basic triangulation.
    call tria_readTriFile3D (Rlevels(NLMIN)%rtriangulation, &
                             trim(spredir)//"/CUBE.tri")
    
    ! Refine the mesh up to the minimum level
    call tria_quickRefine2LevelOrdering (NLMIN-1,&
        Rlevels(NLMIN)%rtriangulation)
    
    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    call tria_initStandardMeshFromRaw (Rlevels(NLMIN)%rtriangulation)
    
    ! Now refine the grid for the fine levels.
    do i = NLMIN+1, NLMAX

      ! Refine the grid using the 2-Level-Ordering algorithm
      call tria_refine2LevelOrdering(Rlevels(i-1)%rtriangulation,&
          Rlevels(i)%rtriangulation)
      
      ! Create a standard mesh
      call tria_initStandardMeshFromRaw(Rlevels(i)%rtriangulation)
    
    end do

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up a discretisation structure which tells the code which
    ! finite element to use
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies 4 blocks in the
    ! solution vector.
    do i = NLMIN, NLMAX
      call spdiscr_initBlockDiscr (Rlevels(i)%rdiscretisation, 4, &
                                   Rlevels(i)%rtriangulation)
    end do

    ! rdiscretisation%RspatialDiscr is a list of scalar
    ! discretisation structures for every component of the solution vector.
    ! We have a solution vector with three components:
    !  Component 1 = X-velocity
    !  Component 2 = Y-velocity
    !  Component 3 = Z-velocity
    !  Component 4 = Pressure
    do i = NLMIN, NLMAX
      ! For simplicity, we set up one discretisation structure for the
      ! velocity...
      call spdiscr_initDiscr_simple (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
          EL_EM30_3D, Rlevels(i)%rtriangulation)
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd and 3rd component (Y-/Z-velocity). This needs no
      ! additional memory, as all three structures will share the same dynamic
      ! information afterwards.
      call spdiscr_duplicateDiscrSc (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
          Rlevels(i)%rdiscretisation%RspatialDiscr(2))
      call spdiscr_duplicateDiscrSc (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
          Rlevels(i)%rdiscretisation%RspatialDiscr(3))

      ! For the pressure (4th component), we set up a separate discretisation
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      call spdiscr_deriveSimpleDiscrSc (Rlevels(i)%rdiscretisation%RspatialDiscr(1), &
          EL_Q0_3D, Rlevels(i)%rdiscretisation%RspatialDiscr(4))
    
    end do

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up an cubature info structure to tell the code which cubature
    ! formula to use
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
                 
    do i = NLMIN, NLMAX
      ! Create an cubature information structure which tells the code
      ! the cubature formula to use. Standard: Gauss 3x3.
      call spdiscr_createDefCubStructure(&  
          Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
          Rlevels(i)%rcubatureInfo,CUB_GEN_AUTO_G3)
    end do

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Create a 3x3 block matrix with the operator
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    do i = NLMIN, NLMAX
    
      ! Initialise the block matrix with default values based on
      ! the discretisation.
      call lsysbl_createMatBlockByDiscr (Rlevels(i)%rdiscretisation,&
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
      call bilf_createMatrixStructure (&
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
      call bilf_createMatrixStructure (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(4),&
          LSYSSC_MATRIX9, Rlevels(i)%rmatrixB1,&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1))
                
      ! Duplicate the B1 matrix structure to the B2/B3 matrix, so use
      ! lsyssc_duplicateMatrix to create B2/B3. Share the matrix
      ! structure between B1 and B2/B3 (B1 is the parent and B2/B3 the children).
      ! Do not create a content array yet, it will be created by
      ! the assembly routines later.
      call lsyssc_duplicateMatrix (Rlevels(i)%rmatrixB1, Rlevels(i)%rmatrixB2,&
          LSYSSC_DUP_COPY, LSYSSC_DUP_REMOVE)
      call lsyssc_duplicateMatrix (Rlevels(i)%rmatrixB1, Rlevels(i)%rmatrixB3,&
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
      rform%ballCoeffConstant = .true.
      rform%Dcoefficients(1)  = dnu
      rform%Dcoefficients(2)  = dnu
      rform%Dcoefficients(3)  = dnu

      ! Now we can build the matrix entries.
      ! We specify the callback function coeff_Laplace for the coefficients.
      ! As long as we use constant coefficients, this routine is not used.
      ! By specifying ballCoeffConstant = .FALSE. above,
      ! the framework will call the callback routine to get analytical data.
      !
      ! We pass our collection structure as well to this routine,
      ! so the callback routine has access to everything what is
      ! in the collection.
      !
      ! Build the X-velocity matrix:
      call bilf_buildMatrixScalar (rform,.true.,&
          Rlevels(i)%rmatrix%RmatrixBlock(1,1), Rlevels(i)%rcubatureInfo,coeff_Stokes_3D)
      
      ! Duplicate the matrix to the Y-/Z-velocity matrix, share structure and
      ! content between them (as the matrices are the same).
      call lsyssc_duplicateMatrix (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
          Rlevels(i)%rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      call lsyssc_duplicateMatrix (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
          Rlevels(i)%rmatrix%RmatrixBlock(3,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      
      ! Manually change the discretisation structure of the Y-/Z-velocity
      ! matrix to the Y-/Z-discretisation structure.
      ! Ok, we use the same discretisation structure for both, X- and Y-velocity,
      ! so this is not really necessary - we do this for sure...
      call lsyssc_assignDiscretisation (Rlevels(i)%rmatrix%RmatrixBlock(2,2),&
          Rlevels(i)%rdiscretisation%RspatialDiscr(2))
      call lsyssc_assignDiscretisation (Rlevels(i)%rmatrix%RmatrixBlock(3,3),&
          Rlevels(i)%rdiscretisation%RspatialDiscr(3))
                                  
      ! Build the first pressure matrix B1.
      ! Again first set up the bilinear form, then call the matrix assembly.
      rform%itermCount = 1
      rform%Idescriptors(1,1) = DER_FUNC3D
      rform%Idescriptors(2,1) = DER_DERIV3D_X

      ! In the standard case, we have constant coefficients:
      rform%ballCoeffConstant = .true.
      rform%Dcoefficients(1)  = -1.0_DP
      
      call bilf_buildMatrixScalar (rform,.true.,Rlevels(i)%rmatrixB1,&
          Rlevels(i)%rcubatureInfo,coeff_Pressure_3D)

      ! Build the second pressure matrix B2.
      ! Again first set up the bilinear form, then call the matrix assembly.
      rform%itermCount = 1
      rform%Idescriptors(1,1) = DER_FUNC3D
      rform%Idescriptors(2,1) = DER_DERIV3D_Y

      ! In the standard case, we have constant coefficients:
      rform%ballCoeffConstant = .true.
      rform%Dcoefficients(1)  = -1.0_DP
      
      call bilf_buildMatrixScalar (rform,.true.,Rlevels(i)%rmatrixB2,&
          Rlevels(i)%rcubatureInfo,coeff_Pressure_3D)
                                  
      ! Build the third pressure matrix B3.
      ! Again first set up the bilinear form, then call the matrix assembly.
      rform%itermCount = 1
      rform%Idescriptors(1,1) = DER_FUNC3D
      rform%Idescriptors(2,1) = DER_DERIV3D_Z

      ! In the standard case, we have constant coefficients:
      rform%ballCoeffConstant = .true.
      rform%Dcoefficients(1)  = -1.0_DP
      
      call bilf_buildMatrixScalar (rform,.true.,Rlevels(i)%rmatrixB3,&
          Rlevels(i)%rcubatureInfo,coeff_Pressure_3D)

      ! The B1/B2/B3 matrices exist up to now only in our local problem structure.
      ! Put a copy of them into the block matrix.
      !
      ! Note that we share the structure of B1/B2/B3 with those B1/B2/B3 of the
      ! block matrix, while we create copies of the entries. The reason is
      ! that these matrices are modified for boundary conditions later.
      call lsyssc_duplicateMatrix (Rlevels(i)%rmatrixB1, &
          Rlevels(i)%rmatrix%RmatrixBlock(1,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

      call lsyssc_duplicateMatrix (Rlevels(i)%rmatrixB2, &
          Rlevels(i)%rmatrix%RmatrixBlock(2,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      
      call lsyssc_duplicateMatrix (Rlevels(i)%rmatrixB3, &
          Rlevels(i)%rmatrix%RmatrixBlock(3,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

      ! Furthermore, put B1^T, B2^T and B3^T to the block matrix.
      call lsyssc_transposeMatrix (Rlevels(i)%rmatrixB1, &
          Rlevels(i)%rmatrix%RmatrixBlock(4,1),LSYSSC_TR_VIRTUAL)

      call lsyssc_transposeMatrix (Rlevels(i)%rmatrixB2, &
          Rlevels(i)%rmatrix%RmatrixBlock(4,2),LSYSSC_TR_VIRTUAL)

      call lsyssc_transposeMatrix (Rlevels(i)%rmatrixB3, &
          Rlevels(i)%rmatrix%RmatrixBlock(4,3),LSYSSC_TR_VIRTUAL)

    end do

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Create RHS and solution vectors
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Next step: Create a RHS vector, a solution vector and a temporary
    ! vector. All are filled with zero.
    call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,rrhs,.true.)
    call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,rvector,.true.)
    call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,rtempBlock,.true.)

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
    call linf_buildVectorScalar (rlinform,.true.,rrhs%RvectorBlock(1),&
        Rlevels(NLMAX)%rcubatureInfo,coeff_RHS_X_3D)

    call linf_buildVectorScalar (rlinform,.true.,rrhs%RvectorBlock(2),&
        Rlevels(NLMAX)%rcubatureInfo,coeff_RHS_Y_3D)
                                
    call linf_buildVectorScalar (rlinform,.true.,rrhs%RvectorBlock(3),&
        Rlevels(NLMAX)%rcubatureInfo,coeff_RHS_Z_3D)

    ! The fourth subvector must be zero - as it represents the RHS of
    ! the equation "div(u) = 0".
    call lsyssc_clearVector(rrhs%RvectorBlock(4))
                                
    ! Clear the solution vector on the finest level.
    call lsysbl_clearVector(rvector)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Assembly of matrices/vectors finished
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Discretise the boundary conditions and apply them to the matrix/RHS/sol.
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Now set up the boundary conditions
    do i = NLMIN, NLMAX

      ! Now we need to implement the boundary conditions. To do this, we
      ! first need to create a mesh region describing the mesh`s boundary.
      ! We want to prescribe Dirichlet on the cube`s boundary, except for
      ! the face where the X-coordinate is 1.
      ! We could now manually create a mesh region based on the triangulation`s
      ! nodal-property array and then kick out everything that belongs to the
      ! right face. But we will use the dom3d_cube module, which performs
      ! this task for us.
      call dom3d_cube_calcMeshRegion(rmeshRegion, Rlevels(i)%rtriangulation, &
                                     DOM3D_CUBE_REG_STOKES)

      ! Initialise the structure for discrete boundary conditions.
      call bcasm_initDiscreteBC (Rlevels(i)%rdiscreteBC)

      ! Set Dirichlet BCs for X-velocity:
      call bcasm_newDirichletBConMR(Rlevels(i)%rdiscretisation, 1, &
          Rlevels(i)%rdiscreteBC, rmeshRegion, &
          getBoundaryValuesMR_3D)

      ! Set Dirichlet BCs for Y-velocity:
      call bcasm_newDirichletBConMR(Rlevels(i)%rdiscretisation, 2, &
          Rlevels(i)%rdiscreteBC, rmeshRegion, &
          getBoundaryValuesMR_3D)

      ! Set Dirichlet BCs for Z-velocity:
      call bcasm_newDirichletBConMR(Rlevels(i)%rdiscretisation, 3, &
          Rlevels(i)%rdiscreteBC, rmeshRegion, &
          getBoundaryValuesMR_3D)

      ! The pressure does not need boundary conditions.
            
      ! Do not forget to release the mesh region
      call mshreg_done(rmeshRegion)

      ! Next step is to implement boundary conditions into the matrix.
      ! This is done using a matrix filter for discrete boundary conditions.
      ! The discrete boundary conditions are already attached to the
      ! matrix. Call the appropriate matrix filter that modifies the matrix
      ! according to the boundary conditions.
      call matfil_discreteBC (Rlevels(i)%rmatrix,Rlevels(i)%rdiscreteBC)
    
      ! Create a filter chain for the solver that implements boundary conditions
      ! during the solution process.
      call filter_initFilterChain (Rlevels(i)%RfilterChain,Rlevels(i)%nfilters)
      call filter_newFilterDiscBCDef (&
          Rlevels(i)%RfilterChain,Rlevels(i)%nfilters,Rlevels(i)%rdiscreteBC)

    end do

    ! Also implement the discrete boundary conditions on the finest level
    ! onto our right-hand-side and solution vectors.

    call vecfil_discreteBCrhs (rrhs,Rlevels(NLMAX)%rdiscreteBC)
    call vecfil_discreteBCsol (rvector,Rlevels(NLMAX)%rdiscreteBC)
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up a linear solver
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Now we have to build up the level information for multigrid.
    !
    ! Create a Multigrid-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    call linsol_initMultigrid2 (p_rsolverNode,NLMAX-NLMIN+1)
    
    ! Set up a BiCGStab solver with VANKA preconditioning as coarse grid solver:
    call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_3DNAVST)
    call linsol_initBiCGStab (p_rcoarseGridSolver,p_rpreconditioner,&
        Rlevels(NLMIN)%RfilterChain)

    ! The coarse grid in multigrid is always grid 1!
    call linsol_getMultigrid2Level (p_rsolverNode,1,p_rlevelInfo)
    p_rlevelInfo%p_rcoarseGridSolver => p_rcoarseGridSolver
    p_rlevelInfo%p_rfilterChain => Rlevels(NLMIN)%RfilterChain

    ! Now set up the other levels...
    do i = NLMIN+1, NLMAX
    
      ! Set up the general VANKA smoother.
      call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_3DNAVST)
      
      ! We will use 4 smoothing steps with damping parameter 0.7
      call linsol_convertToSmoother(p_rsmoother, 4, 0.7_DP)
      
      ! And add this multi-grid level. We will use the same smoother
      ! for pre- and post-smoothing.
      call linsol_getMultigrid2Level (p_rsolverNode,i-NLMIN+1,p_rlevelInfo)
      p_rlevelInfo%p_rpresmoother => p_rsmoother
      p_rlevelInfo%p_rpostsmoother => p_rsmoother
      p_rlevelInfo%p_rfilterChain => Rlevels(i)%RfilterChain
      
    end do
    
    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2
    
    ! Attach the system matrices to the solver.
    !
    ! We copy our matrices to a big matrix array and transfer that
    ! to the setMatrices routines. This intitialises then the matrices
    ! on all levels according to that array. Note that this does not
    ! allocate new memory, we create only "links" to existing matrices
    ! into Rmatrices(:)!
    allocate(Rmatrices(NLMIN:NLMAX))
    do i = NLMIN, NLMAX
      call lsysbl_duplicateMatrix (Rlevels(i)%rmatrix,&
          Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    end do
    
    call linsol_setMatrices(p_RsolverNode,Rmatrices(NLMIN:NLMAX))

    ! We can release Rmatrices immediately -- as long as we do not
    ! release Rlevels(i)%rmatrix!
    do i=NLMIN,NLMAX
      call lsysbl_releaseMatrix (Rmatrices(i))
    end do
    deallocate(Rmatrices)

    ! Initialise structure/data of the solver. This allows the
    ! solver to allocate memory / perform some precalculation
    ! to the problem.
    call linsol_initStructure (p_rsolverNode, ierror)
    
    if (ierror .ne. LINSOL_ERR_NOERROR) then
      call output_line("Matrix structure invalid!",OU_CLASS_ERROR)
      call sys_halt()
    end if

    call linsol_initData (p_rsolverNode, ierror)
    
    if (ierror .ne. LINSOL_ERR_NOERROR) then
      call output_line("Matrix singular!",OU_CLASS_ERROR)
      call sys_halt()
    end if
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Solve the system
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Finally solve the system. As we want to solve Ax=b with
    ! b being the real RHS and x being the real solution vector,
    ! we use linsol_solveAdaptively. If b is a defect
    ! RHS and x a defect update to be added to a solution vector,
    ! we would have to use linsol_precondDefect instead.
    call linsol_solveAdaptively (p_rsolverNode,rvector,rrhs,rtempBlock)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Postprocessing of the solution
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! The solution vector is probably not in the way Paraview likes it!
    ! Paraview for example does not understand Q1~ vectors!
    ! Therefore, we first have to convert the vector to a form that
    ! Paraview understands.
    ! Paraview understands only Q1 solutions! So the task is now to create
    ! a Q1 solution from p_rvector and write that out.
    !
    ! For this purpose, first create a "derived" simple discretisation
    ! structure based on Q1 by copying the main guiding block discretisation
    ! structure and modifying the discretisation structures of the
    ! two velocity subvectors:
    
    call spdiscr_duplicateBlockDiscr (Rlevels(NLMAX)%rdiscretisation,&
                                      rprjDiscretisation)
    
    call spdiscr_deriveSimpleDiscrSc (&
        Rlevels(NLMAX)%rdiscretisation%RspatialDiscr(1), &
        EL_Q1_3D, rprjDiscretisation%RspatialDiscr(1))
    call spdiscr_deriveSimpleDiscrSc (&
        Rlevels(NLMAX)%rdiscretisation%RspatialDiscr(2), &
        EL_Q1_3D, rprjDiscretisation%RspatialDiscr(2))
    call spdiscr_deriveSimpleDiscrSc (&
        Rlevels(NLMAX)%rdiscretisation%RspatialDiscr(3), &
        EL_Q1_3D, rprjDiscretisation%RspatialDiscr(3))
                 
    ! The pressure discretisation substructure stays the old.
    !
    ! Now set up a new solution vector based on this discretisation,
    ! allocate memory.
    call lsysbl_createVecBlockByDiscr (rprjDiscretisation,rprjVector,.false.)
    
    ! Then take our original solution vector and convert it according to the
    ! new discretisation:
    call spdp_projectSolution (rvector,rprjVector)
    
    ! Discretise the boundary conditions according to the Q1/Q1/Q1/Q0
    ! discretisation:
    
    ! Get the boundary mesh region on the finest level
    call dom3d_cube_calcMeshRegion(rmeshRegion, Rlevels(NLMAX)%rtriangulation, &
                                   DOM3D_CUBE_REG_STOKES)

    ! And prescribe Dirichlet boundary conditions
    call bcasm_initDiscreteBC (rprjDiscreteBC)
    call bcasm_newDirichletBConMR(rprjDiscretisation, 1, rprjDiscreteBC, &
                       rmeshRegion, getBoundaryValuesMR_3D)
    call bcasm_newDirichletBConMR(rprjDiscretisation, 2, rprjDiscreteBC, &
                       rmeshRegion, getBoundaryValuesMR_3D)
    call bcasm_newDirichletBConMR(rprjDiscretisation, 3, rprjDiscreteBC, &
                       rmeshRegion, getBoundaryValuesMR_3D)

    ! Now we do not need the mesh region anymore, so release it
    call mshreg_done(rmeshRegion)

    
    ! Send the vector to the boundary-condition implementation filter.
    ! This modifies the vector according to the discrete boundary
    ! conditions.
    call vecfil_discreteBCsol (rprjVector,rprjDiscreteBC)
    
    ! Get the path for writing postprocessing files from the environment variable
    ! $UCDDIR. If that does not exist, write to the directory "./gmv".
    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = "./gmv"

    ! Now we have a Q1/Q1/Q0 solution in rprjVector.
    ! We can now start the postprocessing.
    ! Start UCD export to VTK file:
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,&
        Rlevels(NLMAX)%rtriangulation,trim(sucddir)//"/u3d_1_mg.vtk")

    ! Write velocity field
    call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
    call lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
    call lsyssc_getbase_double (rprjVector%RvectorBlock(3),p_Ddata3)
    
    ! In case we use the VTK exporter, which supports vector output, we will
    ! pass the X-,Y- and Z-velocity at once to the ucd module.
    call ucd_addVarVertBasedVec(rexport,"velocity",p_Ddata,p_Ddata2,p_Ddata3)

    ! If we use the GMV exporter, we might replace the line above by the
    ! following two lines:
    !CALL ucd_addVariableVertexBased (rexport,"X-vel",UCD_VAR_XVELOCITY, p_Ddata)
    !CALL ucd_addVariableVertexBased (rexport,"Y-vel",UCD_VAR_YVELOCITY, p_Ddata2)
    !CALL ucd_addVariableVertexBased (rexport,"Z-vel",UCD_VAR_ZVELOCITY, p_Ddata3)
        
    ! Write pressure
    call lsyssc_getbase_double (rprjVector%RvectorBlock(4),p_Ddata)
    call ucd_addVariableElementBased (rexport,"pressure",UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that is it.
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! We are finished - but not completely!
    ! Now, clean up so that all the memory is available again.
    !
    ! Release solver data and structure
    call linsol_doneData (p_rsolverNode)
    call linsol_doneStructure (p_rsolverNode)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    call linsol_releaseSolver (p_rsolverNode)
    
    ! Release the filter chain
    do i=NLMIN,NLMAX
      call filter_doneFilterChain (Rlevels(i)%RfilterChain,Rlevels(i)%nfilters)
    end do

    ! Release the block matrix/vectors
    call lsysbl_releaseVector (rprjVector)
    call lsysbl_releaseVector (rtempBlock)
    call lsysbl_releaseVector (rvector)
    call lsysbl_releaseVector (rrhs)
    do i = NLMAX, NLMIN, -1
      call lsysbl_releaseMatrix (Rlevels(i)%rmatrix)
    end do
    
    ! Release B1 and B2 matrix
    do i = NLMAX, NLMIN, -1
      call lsyssc_releaseMatrix (Rlevels(i)%rmatrixB3)
      call lsyssc_releaseMatrix (Rlevels(i)%rmatrixB2)
      call lsyssc_releaseMatrix (Rlevels(i)%rmatrixB1)
    end do
    
    ! Release the cubature info structures.
    do i = NLMAX, NLMIN, -1
      call spdiscr_releaseCubStructure(Rlevels(i)%rcubatureInfo)
    end do

    ! Release our discrete version of the boundary conditions
    call bcasm_releaseDiscreteBC (rprjDiscreteBC)
    do i = NLMAX, NLMIN, -1
      call bcasm_releaseDiscreteBC (Rlevels(i)%rdiscreteBC)
    end do

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    call spdiscr_releaseBlockDiscr(rprjDiscretisation)
    do i = NLMAX, NLMIN, -1
      call spdiscr_releaseBlockDiscr(Rlevels(i)%rdiscretisation)
    end do
    
    ! Release the triangulation.
    do i = NLMAX, NLMIN, -1
      call tria_done (Rlevels(i)%rtriangulation)
    end do
    
    deallocate(Rlevels)
    
  end subroutine

end module
