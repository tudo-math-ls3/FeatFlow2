!##############################################################################
!# ****************************************************************************
!# <name> navst2d_method1_mg </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstation program how to solve a (Navier-)Stokes
!# problem on a simple domain.
!#
!# The routine uses a SP-SOR smoother for a multigrid solver.
!#
!# This module is based on the stokes2d_method1_mg example, but in contrast
!# this module uses the bench1 domain and, optionally, solves a Navier-Stokes
!# system instead of a Stokes system.
!#
!# This example can be seen as the most simple stationary 2D Navier-Stokes
!# solver without any fancy features as optimal damping parameter control or
!# Newton-iteration.
!# </purpose>
!##############################################################################

module navst2d_method1_mg

  use fsystem
  use storage
  use genoutput
  use boundary
  use cubature
  use derivatives
  use matrixfilters
  use vectorfilters
  use linearalgebra
  use discretebc
  use bcassembly
  use triangulation
  use element
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use spdiscprojection
  use scalarpde
  use bilinearformevaluation
  use trilinearformevaluation
  use linearformevaluation
  use discretebc
  use filtersupport
  use coarsegridcorrection
  use multilevelprojection
  use linearsolver
  use ucd
  use convection
  
  use stokes2d_callback
  
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
    
    ! A system matrix for that specific level.
    type(t_matrixBlock) :: rmatrix

    ! The solution vector for that specific level.
    type(t_vectorBlock) :: rvecSol

    ! The Laplace matrix for that level.
    type(t_matrixScalar) :: rmatrixStokes

    ! B1-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixB1

    ! B2-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixB2

    ! A variable describing the discrete boundary conditions.
    type(t_discreteBC) :: rdiscreteBC
  
  end type
  
!</typeblock>

!</types>

contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine navst2d_1_mg
  
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

    ! An object for saving the domain:
    type(t_boundary) :: rboundary
    
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform
    type(t_trilinearForm) :: rtriform1, rtriform2
    type(t_linearForm) :: rlinform
    
    ! A streamline-diffusion structure
    type(t_convStreamlineDiffusion) :: rsd
    
    ! An upwind structure
    type(t_convUpwind) :: rupwind

    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_vectorBlock) :: rvecDef,rrhs,rtempBlock
    
    ! A set of variables describing the analytic and discrete boundary
    ! conditions.
    type(t_boundaryRegion) :: rboundaryRegion

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode,p_rcoarseGridSolver,p_rsmoother

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(:), pointer :: Rmatrices

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    type(t_filterChain), dimension(1), target :: RfilterChain
    
    ! An interlevel projection structure for changing levels
    type(t_interlevelProjectionBlock) :: rprojection

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
    real(DP), dimension(:), pointer :: p_Du1,p_Du2,p_Dp

    ! A counter variable
    integer :: i,nl,niterMaxNL
    
    ! Stokes or Navier-Stokes?
    logical :: bNavier
    
    ! What kind of stabilisation?
    integer :: iConvAsm
    
    ! Streamline-diffusion parameter
    real(DP) :: dupsam
    
    ! Residual for the non-linear iteration
    real(DP) :: dnlres, dnlresInit
    
    ! Damping parameter
    real(DP) :: dnlDamping

    ! Path to the mesh
    character(len=SYS_STRLEN) :: spredir
    
    ! We want to solve our (Navier-)Stokes problem on level...
    NLMIN = 1
    NLMAX = 4
    
    ! Viscosity parameter: nu = 1 / RE
    dnu = 1.0_DP / 1000.0_DP
    
    ! Do we want to assemble a Navier-Stokes system?
    bNavier = .true.
    
    ! If we want to solve Navier-Stokes, how do we assemble the
    ! convective term?
    ! 0 => trilinearform (no stabilisation)
    ! 1 => streamline-diffusion
    ! 2 => upwind
    iConvAsm = 1
    
    ! In the case of streamline-diffusion or upwind, what is the parameter?
    dupsam = 1.0_DP
    
    ! Maximum number of non-linear loop iterations
    niterMaxNL = 20
    
    ! Fixed damping parameter for non-linear loop
    dnlDamping = 1.0_DP
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Read the domain, read the mesh, refine
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Allocate memory for all levels
    allocate(Rlevels(NLMIN:NLMAX))

    ! Get the path $PREDIR from the environment, where to read .prm/.tri files
    ! from. If that does not exist, write to the directory "./pre".
    if (.not. sys_getenv_string("PREDIR", spredir)) spredir = './pre'

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call boundary_read_prm(rboundary, trim(spredir)//'/bench1.prm')
        
    ! Now read in the basic triangulation.
    call tria_readTriFile2D (Rlevels(NLMIN)%rtriangulation, &
                             trim(spredir)//'/bench1.tri', rboundary)
    
    ! Refine the mesh up to the minimum level
    call tria_quickRefine2LevelOrdering (NLMIN-1,&
        Rlevels(NLMIN)%rtriangulation,rboundary)
    
    ! And create information about adjacencies and everything one needs
    ! from a triangulation.
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

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up a discretisation structure which tells the code which
    ! finite element to use
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies 3 blocks in the
    ! solution vector.
    do i = NLMIN, NLMAX
      call spdiscr_initBlockDiscr (Rlevels(i)%rdiscretisation, 3, &
                                   Rlevels(i)%rtriangulation, rboundary)
    end do

    ! rdiscretisation%RspatialDiscr is a list of scalar
    ! discretisation structures for every component of the solution vector.
    ! We have a solution vector with three components:
    !  Component 1 = X-velocity
    !  Component 2 = Y-velocity
    !  Component 3 = Pressure
    do i = NLMIN, NLMAX
      ! For simplicity, we set up one discretisation structure for the
      ! velocity...
      call spdiscr_initDiscr_simple (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
          EL_EM30, Rlevels(i)%rtriangulation, rboundary)
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory,
      ! as both structures will share the same dynamic information afterwards.
      call spdiscr_duplicateDiscrSc (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
          Rlevels(i)%rdiscretisation%RspatialDiscr(2))

      ! For the pressure (3rd component), we set up a separate discretisation
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      call spdiscr_deriveSimpleDiscrSc (Rlevels(i)%rdiscretisation%RspatialDiscr(1), &
          EL_Q0, Rlevels(i)%rdiscretisation%RspatialDiscr(3))
    
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
          LSYSSC_MATRIX9, Rlevels(i)%rmatrixStokes)

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
      ! spatial discretisation structure in p_rdiscretisation%RspatialDiscr.
      call bilf_createMatrixStructure (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(3),&
          LSYSSC_MATRIX9, Rlevels(i)%rmatrixB1,&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1))
                
      ! Duplicate the B1 matrix structure to the B2 matrix, so use
      ! lsyssc_duplicateMatrix to create B2. Share the matrix
      ! structure between B1 and B2 (B1 is the parent and B2 the child).
      ! Do not create a content array yet, it will be created by
      ! the assembly routines later.
      call lsyssc_duplicateMatrix (Rlevels(i)%rmatrixB1, Rlevels(i)%rmatrixB2,&
                                   LSYSSC_DUP_COPY, LSYSSC_DUP_REMOVE)
                                       
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
      rform%ballCoeffConstant = .true.
      rform%BconstantCoeff = .true.
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
      call bilf_buildMatrixScalar (rform,.true.,&
          Rlevels(i)%rmatrixStokes, Rlevels(i)%rcubatureInfo, coeff_Stokes_2D)

      ! Duplicate the Laplace matrix to the X-velocity matrix and share
      ! the structure. If we want to assemble a Stokes system, we can also
      ! share the content with the Laplace matrix, in the case of a
      ! Navier-Stokes system the velocity blocks will also recieve the
      ! convective term (u * \nabla)u, so the content needs to be copied.
      if (bNavier) then
        call lsyssc_duplicateMatrix (Rlevels(i)%rmatrixStokes,&
            Rlevels(i)%rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      else
        call lsyssc_duplicateMatrix (Rlevels(i)%rmatrixStokes,&
            Rlevels(i)%rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      end if
      
      ! Duplicate the matrix to the Y-velocity matrix, share structure and
      ! content between them (as the matrices are the same).
      call lsyssc_duplicateMatrix (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
          Rlevels(i)%rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      
      ! Manually change the discretisation structure of the Y-velocity
      ! matrix to the Y-discretisation structure.
      ! Ok, we use the same discretisation structure for both, X- and Y-velocity,
      ! so this is not really necessary - we do this for sure...
      call lsyssc_assignDiscrDirectMat (Rlevels(i)%rmatrix%RmatrixBlock(2,2),&
          Rlevels(i)%rdiscretisation%RspatialDiscr(2))
                                  
      ! Build the first pressure matrix B1.
      ! Again first set up the bilinear form, then call the matrix assembly.
      rform%itermCount = 1
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_DERIV_X

      ! In the standard case, we have constant coefficients:
      rform%ballCoeffConstant = .true.
      rform%BconstantCoeff = .true.
      rform%Dcoefficients(1)  = -1.0_DP
      
      call bilf_buildMatrixScalar (rform,.true.,Rlevels(i)%rmatrixB1,&
          Rlevels(i)%rcubatureInfo,coeff_Pressure_2D)

      ! Build the second pressure matrix B2.
      ! Again first set up the bilinear form, then call the matrix assembly.
      rform%itermCount = 1
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_DERIV_Y

      ! In the standard case, we have constant coefficients:
      rform%ballCoeffConstant = .true.
      rform%BconstantCoeff = .true.
      rform%Dcoefficients(1)  = -1.0_DP
      
      call bilf_buildMatrixScalar (rform,.true.,Rlevels(i)%rmatrixB2,&
          Rlevels(i)%rcubatureInfo,coeff_Pressure_2D)
                                  
      ! The B1/B2 matrices exist up to now only in our local problem structure.
      ! Put a copy of them into the block matrix.
      !
      ! Note that we share the structure of B1/B2 with those B1/B2 of the
      ! block matrix, while we create copies of the entries. The reason is
      ! that these matrices are modified for boundary conditions later.
      call lsyssc_duplicateMatrix (Rlevels(i)%rmatrixB1, &
          Rlevels(i)%rmatrix%RmatrixBlock(1,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

      call lsyssc_duplicateMatrix (Rlevels(i)%rmatrixB2, &
          Rlevels(i)%rmatrix%RmatrixBlock(2,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      
      ! Furthermore, put B1^T and B2^T to the block matrix.
      call lsyssc_transposeMatrix (Rlevels(i)%rmatrixB1, &
          Rlevels(i)%rmatrix%RmatrixBlock(3,1),LSYSSC_TR_ALL)

      call lsyssc_transposeMatrix (Rlevels(i)%rmatrixB2, &
          Rlevels(i)%rmatrix%RmatrixBlock(3,2),LSYSSC_TR_ALL)

      ! Create the solution vector for this level
      call lsysbl_createVectorBlock (Rlevels(i)%rdiscretisation,Rlevels(i)%rvecSol,.true.)
    
    end do

    ! Prepare the trilinearform for the convection operator (u * \nabla) u
    ! The trilinearform is evaluated in the nonlinear defect correction loop
    ! below.
    ! First convective term: u1 * d_x * u
    rtriform1%itermCount = 1
    rtriform1%ballCoeffConstant = .true.
    rtriform1%BconstantCoeff = .true.
    rtriform1%Dcoefficients = 1.0_DP
    rtriform1%Idescriptors(1,1) = DER_FUNC
    rtriform1%Idescriptors(2,1) = DER_DERIV_X
    rtriform1%Idescriptors(3,1) = DER_FUNC
    ! Second convective term: u2 * d_y * u
    rtriform2%itermCount = 1
    rtriform2%ballCoeffConstant = .true.
    rtriform2%BconstantCoeff = .true.
    rtriform2%Dcoefficients = 1.0_DP
    rtriform2%Idescriptors(1,1) = DER_FUNC
    rtriform2%Idescriptors(2,1) = DER_DERIV_Y
    rtriform2%Idescriptors(3,1) = DER_FUNC

    ! Prepare the streamline-diffusion structure
    rsd%dnu = dnu
    rsd%dupsam = dupsam
    
    ! Prepare the upwind structure
    rupwind%dnu = dnu
    rupwind%dupsam = dupsam

    ! Next step: Create a RHS vector, a solution vector and a temporary
    ! vector. All are filled with zero.
    call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,rrhs,.true.)
    call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,rvecDef,.true.)
    call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,rtempBlock,.true.)

    ! The vector structure is ready but the entries are missing.
    ! So the next thing is to calculate the content of that vector.
    !
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    
    ! ... and then discretise the RHS to the first two subvectors of
    ! the block vector using the discretisation structure of the
    ! corresponding blocks.
    !
    ! Note that the vector is unsorted after calling this routine!
    call linf_buildVectorScalar (&
        rlinform,.true.,rrhs%RvectorBlock(1),Rlevels(NLMAX)%rcubatureInfo,coeff_RHS_X_2D)

    call linf_buildVectorScalar (&
        rlinform,.true.,rrhs%RvectorBlock(2),Rlevels(NLMAX)%rcubatureInfo,coeff_RHS_Y_2D)
                                
    ! The third subvector must be zero - as it represents the RHS of
    ! the equation "div(u) = 0".
    call lsyssc_clearVector(rrhs%RvectorBlock(3))
                                
    ! Clear the solution vector on the finest level.
    call lsysbl_clearVector(Rlevels(NLMAX)%rvecSol)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Assembly of matrices/vectors finished
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Discretise the boundary conditions and apply them to the matrix/RHS/sol.
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    do i = NLMIN, NLMAX

      ! For implementing boundary conditions, we use a `filter technique with
      ! discretised boundary conditions`. This means, we first have to calculate
      ! a discrete version of the analytic BC, which we can implement into the
      ! solution/RHS vectors using the corresponding filter.
      !
      ! Create a t_discreteBC structure where we store all discretised boundary
      ! conditions.
      call bcasm_initDiscreteBC(Rlevels(i)%rdiscreteBC)
      
      ! We first set up the boundary conditions for the X-velocity, then those
      ! of the Y-velocity.
      !
      ! We 'know' already (from the problem definition) that we have four boundary
      ! segments in the domain. Each of these, we want to use for enforcing
      ! some kind of boundary condition.
      !
      ! We ask the bondary routines to create a 'boundary region' - which is
      ! simply a part of the boundary corresponding to a boundary segment.
      ! A boundary region roughly contains the type, the min/max parameter value
      ! and whether the endpoints are inside the region or not.
      call boundary_createRegion(rboundary,1,1,rboundaryRegion)
      
      ! The endpoint of this segment should also be Dirichlet. We set this by
      ! changing the region properties in rboundaryRegion.
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
      
      ! We use this boundary region and specify that we want to have Dirichlet
      ! boundary there. The following call does the following:
      ! - Create Dirichlet boundary conditions on the region rboundaryRegion.
      !   We specify icomponent='1' to indicate that we set up the
      !   Dirichlet BC`s for the first (here: one and only) component in the
      !   solution vector.
      ! - Discretise the boundary condition so that the BC`s can be applied
      !   to matrices and vectors
      ! - Add the calculated discrete BC`s to Rlevels(i)%rdiscreteBC for later use.
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                        rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                        getBoundaryValues_2D)
                               
      ! Edge 2 is Neumann boundary, so it is commented out.
      ! CALL boundary_createRegion(rboundary,1,2,rboundaryRegion)
      ! CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
      !                                    rboundaryRegion,Rlevels(i)%rdiscreteBC,&
      !                                    getBoundaryValues_2D)
                               
      ! Edge 3 of boundary component 1.
      call boundary_createRegion(rboundary,1,3,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                        rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                        getBoundaryValues_2D)
      
      ! Edge 4 of boundary component 1.
      call boundary_createRegion(rboundary,1,4,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                        rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                        getBoundaryValues_2D)

      ! Circle boundary component. That is it.
      call boundary_createRegion(rboundary,2,1,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                        rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                        getBoundaryValues_2D)

      ! Now continue with defining the boundary conditions of the Y-velocity:
      !
      ! Define edge 1.
      call boundary_createRegion(rboundary,1,1,rboundaryRegion)
      
      ! Edge with start- and endpoint.
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
      
      ! As we define the Y-velocity, we now set icomponent=2 in the following call.
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,2,&
                                        rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                        getBoundaryValues_2D)
                               
      ! Edge 2 is Neumann boundary, so it is commented out.
      ! CALL boundary_createRegion(rboundary,1,2,rboundaryRegion)
      ! CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,2,&
      !                                    rboundaryRegion,Rlevels(i)%rdiscreteBC,&
      !                                    getBoundaryValues_2D)
                               
      ! Edge 3 of boundary component 1.
      call boundary_createRegion(rboundary,1,3,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,2,&
                                        rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                        getBoundaryValues_2D)
      
      ! Edge 4 of boundary component 1.
      call boundary_createRegion(rboundary,1,4,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,2,&
                                        rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                        getBoundaryValues_2D)

      ! Circle boundary component. That is it.
      call boundary_createRegion(rboundary,2,1,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,2,&
                                        rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                        getBoundaryValues_2D)

      ! The pressure does not need boundary conditions.
      ! Hang the pointer into the vector and matrix. That way, these
      ! boundary conditions are always connected to that matrix and that
      ! vector.
      call lsysbl_assignDiscreteBC(Rlevels(i)%rmatrix,Rlevels(i)%rdiscreteBC)
      call lsysbl_assignDiscreteBC(Rlevels(i)%rvecSol,Rlevels(i)%rdiscreteBC)

      ! Next step is to implement boundary conditions into the matrix.
      ! This is done using a matrix filter for discrete boundary conditions.
      ! The discrete boundary conditions are already attached to the
      ! matrix. Call the appropriate matrix filter that modifies the matrix
      ! according to the boundary conditions.
      call matfil_discreteBC (Rlevels(i)%rmatrix)
    
    end do

    ! Also implement the discrete boundary conditions on the finest level
    ! onto our right-hand-side and solution vector.
    call lsysbl_assignDiscreteBC(rrhs,Rlevels(NLMAX)%rdiscreteBC)
    call lsysbl_assignDiscreteBC(rvecDef,Rlevels(NLMAX)%rdiscreteBC)
    call vecfil_discreteBCrhs (rrhs)
    call vecfil_discreteBCsol (Rlevels(NLMAX)%rvecSol)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up a linear solver
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

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
    call mlprj_initProjectionMat (rprojection,Rlevels(NLMAX)%rmatrix)

    ! Create a Multigrid-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    call linsol_initMultigrid2 (p_rsolverNode,NLMAX-NLMIN+1,RfilterChain)

    ! As we will use multigrid as a preconditioner for the non-linear loop,
    ! we set the maximum allowed iterations to 10 and the relative convergence
    ! tolerance to 0.01.
    p_rsolverNode%nmaxIterations = 10
    p_rsolverNode%depsRel=1E-2

    ! Set the output level of multigrid to 2 for some output
    p_rsolverNode%ioutputLevel = 2
    
    ! Set up a UMFPACK as coarse grid solver:
    call linsol_initUMFPACK4 (p_rcoarseGridSolver)
    
    ! Set the output level of the coarse grid solver to -1, so that it
    ! does not print any residuals or warning messages...
    p_rcoarseGridSolver%ioutputLevel = -1
    
    ! The coarse grid in multigrid is always grid 1!
    call linsol_getMultigrid2Level (p_rsolverNode,1,p_rlevelInfo)
    p_rlevelInfo%p_rcoarseGridSolver => p_rcoarseGridSolver

    ! Now set up the other levels...
    do i = NLMIN+1, NLMAX
    
      ! Set up the SP-SOR smoother.
      call linsol_initSPSOR (p_rsmoother,LINSOL_SPSOR_NAVST2D)
      
      ! We will use 4 smoothing steps with damping parameter 1.0
      call linsol_convertToSmoother(p_rsmoother, 8, 1.0_DP)
      
      ! And add this multi-grid level. We will use the same smoother
      ! for pre- and post-smoothing.
      call linsol_getMultigrid2Level (p_rsolverNode,i-NLMIN+1,p_rlevelInfo)
      p_rlevelInfo%p_rpresmoother => p_rsmoother
      p_rlevelInfo%p_rpostsmoother => p_rsmoother
      
    end do
    
    ! Attach the system matrix to the solver.
    !
    ! We copy our matrices to a big matrix array and transfer that
    ! to the setMatrices routines. This intitialises then the matrices
    ! on all levels according to that array. Note that this does not
    ! allocate new memory, we create only 'links' to existing matrices
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
    
    if(.not. bNavier) then
      ! Initialise solver data
      call linsol_initData (p_rsolverNode, ierror)
      if (ierror .ne. LINSOL_ERR_NOERROR) then
        call output_line("Matrix singular!",OU_CLASS_ERROR)
        call sys_halt()
      end if
    end if

    ! Okay, everything is set up - so we can start our nonlinear iteration
      
    ! First, calculate the initial non-linear defect
    call lsysbl_copyVector(rrhs, rvecDef)
    call lsysbl_blockMatVec(Rlevels(NLMAX)%rmatrix, &
             Rlevels(NLMAX)%rvecSol, rvecDef, -1.0_DP, 1.0_DP)
    call vecfil_discreteBCdef (rvecDef)
    dnlresInit = lsysbl_vectorNorm(rvecDef, LINALG_NORML2)
    
    ! Print the defect
    call output_separator(OU_SEP_MINUS)
    call output_line('NL-Iteration:    0 |RES| = ' // &
                     trim(sys_sdEP(dnlresInit,20,12)))
    call output_separator(OU_SEP_MINUS)

    ! Make sure the inital defect is not zero, as we need to divide
    ! by it later...
    if (dnlresInit .le. SYS_EPSREAL_DP) dnlresInit = 1.0_DP

    ! Start the non-linear defect-correction loop
    do nl = 1, niterMaxNL

      ! If we want to solve a Navier-Stokes system, we need to
      ! initialise the solver data here:
      if (bNavier) then
      
        ! Initialize solver data
        call linsol_initData (p_rsolverNode, ierror)
        if (ierror .ne. LINSOL_ERR_NOERROR) stop
        
      end if

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Solve the linear sussystem
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      
      ! Call the preconditioner - which is our linear solver
      call linsol_precondDefect(p_rsolverNode, rvecDef)
      
      ! Check the solver result
      if (p_rsolverNode%iresult .ne. 0) then
      
        ! Release the solver data in case of a Navier-Stokes system.
        if (bNavier) then
          call linsol_doneData(p_rsolverNode)
        end if

        ! Print an error
        call output_separator(OU_SEP_MINUS)
        call output_line('NL-Iteration: ERROR: linear solver broke down')
        call output_separator(OU_SEP_MINUS)
        
        ! Exit the non-linear loop
        exit
        
      end if
      
      ! Add the preconditioned defect onto the solution
      call lsysbl_vectorLinearComb(rvecDef, Rlevels(NLMAX)%rvecSol, &
                                   dnlDamping, 1.0_DP)

      ! In the case of a Navier-Stokes system, we need to reassemble
      ! the system matrices on all levels.
      if (bNavier) then
      
        ! Release solver data
        call linsol_doneData(p_rsolverNode)
        
        ! Restrict the new solution vector to all levels
        do i = NLMAX, NLMIN+1, -1
        
          call mlprj_performInterpolation(rprojection, Rlevels(i-1)%rvecSol,&
                             Rlevels(i)%rvecSol, rtempBlock%RvectorBlock(1))
          
          ! And filter the restricted vector.
          ! Note: We do not need to filter the solution on the finest level
          call vecfil_discreteBCsol (Rlevels(i-1)%rvecSol)
      
        end do
        
        ! Update the matrices on all levels
        do i = NLMIN, NLMAX
        
          ! Copy Stokes matrix to X-velocity block
          call lsyssc_duplicateMatrix(Rlevels(i)%rmatrixStokes,&
                            Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
                            LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPY)
          
          ! Now how do we assemble the convective term?
          select case (iConvAsm)
          case (1)
            ! Use streamline-diffusion
            call conv_streamlineDiffusion2d(Rlevels(i)%rvecSol, &
                 Rlevels(i)%rvecSol, 1.0_DP, 0.0_DP, rsd, CONV_MODMATRIX, &
                 Rlevels(i)%rmatrix%RmatrixBlock(1,1), &
                 rcubatureInfo=Rlevels(i)%rcubatureInfo)
          
          case (2)
            ! Use upwind
            call conv_upwind2d (Rlevels(i)%rvecSol, Rlevels(i)%rvecSol, &
                 1.0_DP, 0.0_DP, rupwind, CONV_MODMATRIX, &
                 Rlevels(i)%rmatrix%RmatrixBlock(1,1))
          
          case default
            ! Use trilinearform
            ! Assemble u_1 * d_x u
            call trilf_buildMatrixScalar(rtriform1,.false.,&
                 Rlevels(i)%rmatrix%RmatrixBlock(1,1), &
                 Rlevels(i)%rvecSol%RvectorBlock(1))
            ! Assemble u_2 * d_y u
            call trilf_buildMatrixScalar(rtriform2,.false.,&
                 Rlevels(i)%rmatrix%RmatrixBlock(1,1), &
                 Rlevels(i)%rvecSol%RvectorBlock(2))
          
          end select
          
          ! And filter the matrix
          call matfil_discreteBC (Rlevels(i)%rmatrix)
          
          ! The other velocity blocks are automatically updated, since they
          ! are just a shared copy of the X-velocity block
          
        end do
      
      end if ! bNavier

      ! Calculate non-linear defect:
      call lsysbl_copyVector(rrhs, rvecDef)
      call lsysbl_blockMatVec(Rlevels(NLMAX)%rmatrix, &
               Rlevels(NLMAX)%rvecSol, rvecDef, -1.0_DP, 1.0_DP)
      
      ! Filter the defect vector
      call vecfil_discreteBCdef (rvecDef)
      
      ! Calculate residual
      dnlres = lsysbl_vectorNorm(rvecDef, LINALG_NORML2)
      call output_separator(OU_SEP_MINUS)
      call output_line('NL-Iteration: ' // trim(sys_si(nl,4)) // &
                       ' |RES| = ' // trim(sys_sdEP(dnlres,20,12)))
      call output_separator(OU_SEP_MINUS)
      
      ! Are we finished?
      if ((dnlres .le. 1E-8_DP) .and. ((dnlres / dnlresInit) .le. 1E-5_DP)) exit
      
      ! Proceed with next non-linear iteration
      
    end do ! nl
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Postprocessing of the solution
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Project the velocity onto the mesh`s vertices and the pressure onto the
    ! mesh`s cells.
    nullify(p_Du1)
    nullify(p_Du2)
    nullify(p_Dp)
    call spdp_projectToVertices(Rlevels(NLMAX)%rvecSol%RvectorBlock(1),p_Du1)
    call spdp_projectToVertices(Rlevels(NLMAX)%rvecSol%RvectorBlock(2),p_Du2)
    call spdp_projectToCells(Rlevels(NLMAX)%rvecSol%RvectorBlock(3),p_Dp)
    
    ! Get the path for writing postprocessing files from the environment variable
    ! $UCDDIR. If that does not exist, write to the directory "./gmv".
    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = './gmv'

    ! Start UCD export to VTK file:
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,&
        Rlevels(NLMAX)%rtriangulation,trim(sucddir)//'/u2d_navst_mg.vtk')

    ! Write velocity field
    call ucd_addVarVertBasedVec(rexport,'velocity',p_Du1,p_Du2)
        
    ! Write pressure
    call ucd_addVariableElementBased (rexport,'pressure',UCD_VAR_STANDARD, p_Dp)
    
    ! Write the file to disc, that is it.
    call ucd_write (rexport)
    call ucd_release (rexport)
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Clean up
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! And release the memory
    deallocate(p_Dp)
    deallocate(p_Du2)
    deallocate(p_Du1)

    ! We are finished - but not completely!
    ! Now, clean up so that all the memory is available again.
    !
    ! If we have solved a Stokes system, then we need to release the
    ! solver data - in case of a Naver-Stokes system this was already
    ! done inside the non-linear loop
    if (.not. bNavier) then
      ! Release solver data
      call linsol_doneData(p_rsolverNode)
    end if

    ! Release solver data and structure
    call linsol_doneStructure (p_rsolverNode)
    
    ! Release the interlevel projection structure
    call mlprj_doneProjection (rprojection)

    ! Release the solver node and all subnodes attached to it (if at all):
    call linsol_releaseSolver (p_rsolverNode)
    
    ! Release the block matrix/vectors
    call lsysbl_releaseVector (rtempBlock)
    call lsysbl_releaseVector (rvecDef)
    call lsysbl_releaseVector (rrhs)
    do i = NLMAX, NLMIN, -1
      call lsysbl_releaseVector (Rlevels(i)%rvecSol)
      call lsysbl_releaseMatrix (Rlevels(i)%rmatrix)
    end do
    
    ! Release B1 and B2 matrix
    do i = NLMAX, NLMIN, -1
      call lsyssc_releaseMatrix (Rlevels(i)%rmatrixB2)
      call lsyssc_releaseMatrix (Rlevels(i)%rmatrixB1)
      call lsyssc_releaseMatrix (Rlevels(i)%rmatrixStokes)
    end do
    
    ! Release the cubature info structures.
    do i = NLMAX, NLMIN, -1
      call spdiscr_releaseCubStructure(Rlevels(i)%rcubatureInfo)
    end do

    ! Release our discrete version of the boundary conditions
    do i = NLMAX, NLMIN, -1
      call bcasm_releaseDiscreteBC (Rlevels(i)%rdiscreteBC)
    end do

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    do i = NLMAX, NLMIN, -1
      call spdiscr_releaseBlockDiscr(Rlevels(i)%rdiscretisation)
    end do
    
    ! Release the triangulation.
    do i = NLMAX, NLMIN, -1
      call tria_done (Rlevels(i)%rtriangulation)
    end do
    
    deallocate(Rlevels)
    
    ! Finally release the domain, that is it.
    call boundary_release (rboundary)

  end subroutine

end module
