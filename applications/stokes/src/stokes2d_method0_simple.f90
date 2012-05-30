!##############################################################################
!# ****************************************************************************
!# <name> stokes2d_method0_simple </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstation program how to solve a Stokes
!# problem on a simple domain.
!#
!# The routine uses the BiCGStab solver with a simple-VANKA preconditioner
!# for 2D saddle point problems, Jacobi-Type.
!# </purpose>
!##############################################################################

module stokes2d_method0_simple

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
  use linearformevaluation
  use discretebc
  use filtersupport
  use coarsegridcorrection
  use linearsolver
  use ucd
  use collection, only: t_collection
  use pprocerror
  use genoutput
  
  use stokes2d_callback
  
  implicit none

contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine stokes2d_0_simple
  
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
    ! An object for saving the domain:
    type(t_boundary) :: rboundary

    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    type(t_blockDiscretisation) :: rdiscretisation,rprjDiscretisation
    
    ! Cubature info structure which encapsules the cubature formula
    type(t_scalarCubatureInfo) :: rcubatureInfo

    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform

    ! A scalar matrix and vector. The vector accepts the RHS of the problem
    ! in scalar form.
    type(t_matrixScalar) :: rmatrixB1, rmatrixB2

    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_matrixBlock) :: rmatrix
    type(t_vectorBlock) :: rvector,rrhs,rtempBlock,rprjVector
    
    ! A set of variables describing the analytic and discrete boundary
    ! conditions.
    type(t_boundaryRegion) :: rboundaryRegion
    type(t_discreteBC), target :: rdiscreteBC, rprjDiscreteBC

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode,p_rpreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(1) :: Rmatrices

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    type(t_filterChain), dimension(1), target :: RfilterChain
    
    ! NLMAX receives the level where we want to solve.
    integer :: NLMAX
    
    ! Viscosity parameter nu = 1/Re
    real(DP) :: dnu
    
    ! Error indicator during initialisation of the solver
    integer :: ierror
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir
    real(DP), dimension(:), pointer :: p_Ddata,p_Ddata2

    ! Path to the mesh
    character(len=SYS_STRLEN) :: spredir

    ! A collection structure for post-processing
    type(t_collection) :: rcollection

    ! Error data structures for post-processing
    type(t_errorScVec) :: rerrorU, rerrorP

    ! Error arrays for post-processing
    real(DP), dimension(2), target :: DerrorUL2, DerrorUH1
    real(DP), dimension(1), target :: DerrorPL2

    ! Ok, let us start.
    !
    ! We want to solve our Poisson problem on level...
    ! As we do not use a multigrid solver here, we will set the level
    ! to 5 instead of 7.
    NLMAX = 5
    
    ! Viscosity parameter:
    dnu = 1.0_DP

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Read the the mesh, refine
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Get the path $PREDIR from the environment, where to read .prm/.tri files
    ! from. If that does not exist, write to the directory "./pre".
    if (.not. sys_getenv_string("PREDIR", spredir)) spredir = './pre'

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call boundary_read_prm(rboundary, trim(spredir)//'/QUAD.prm')
        
    ! Now read in the basic triangulation.
    call tria_readTriFile2D (rtriangulation, trim(spredir)//'/QUAD.tri', rboundary)
    
    ! Refine the mesh up to the minimum level
    call tria_quickRefine2LevelOrdering(NLMAX-1,rtriangulation,rboundary)
    
    ! Create information about adjacencies and everything one needs from
    ! a triangulation. Afterwards, we have the coarse mesh.
    call tria_initStandardMeshFromRaw (rtriangulation,rboundary)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up a discretisation structure which tells the code which
    ! finite element to use
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies 3 blocks in the
    ! solution vector.
    call spdiscr_initBlockDiscr (rdiscretisation,3,rtriangulation, rboundary)

    ! rdiscretisation%RspatialDiscr is a list of scalar
    ! discretisation structures for every component of the solution vector.
    ! We have a solution vector with three components:
    !  Component 1 = X-velocity
    !  Component 2 = Y-velocity
    !  Component 3 = Pressure
    ! For simplicity, we set up one discretisation structure for the
    ! velocity...
    call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1),&
        EL_EM30, rtriangulation, rboundary)
                
    ! ...and copy this structure also to the discretisation structure
    ! of the 2nd component (Y-velocity). This needs no additional memory,
    ! as both structures will share the same dynamic information afterwards.
    call spdiscr_duplicateDiscrSc(rdiscretisation%RspatialDiscr(1),&
        rdiscretisation%RspatialDiscr(2))

    ! For the pressure (3rd component), we set up a separate discretisation
    ! structure, as this uses different finite elements for trial and test
    ! functions.
    call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(1), &
        EL_Q0, rdiscretisation%RspatialDiscr(3))

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up an cubature info structure to tell the code which cubature
    ! formula to use
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
                 
    ! Create an assembly information structure which tells the code
    ! the cubature formula to use. Standard: Gauss 3x3.
    call spdiscr_createDefCubStructure(&  
        rdiscretisation%RspatialDiscr(1),rcubatureInfo,CUB_GEN_AUTO_G3)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Create a 3x3 block matrix with the operator
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Initialise the block matrix with default values based on
    ! the discretisation.
    call lsysbl_createMatBlockByDiscr (rdiscretisation,rmatrix)
    
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
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
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
    ! spatial discretisation structure in p_rdiscretisation%RspatialDiscr.
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(3),&
        LSYSSC_MATRIX9, rmatrixB1,rdiscretisation%RspatialDiscr(1))
              
    ! Duplicate the B1 matrix structure to the B2 matrix, so use
    ! lsyssc_duplicateMatrix to create B2. Share the matrix
    ! structure between B1 and B2 (B1 is the parent and B2 the child).
    ! Do not create a content array yet, it will be created by
    ! the assembly routines later.
    call lsyssc_duplicateMatrix (rmatrixB1, rmatrixB2, LSYSSC_DUP_COPY,&
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
    call bilf_buildMatrixScalar (rform,.true.,rmatrix%RmatrixBlock(1,1),&
        rcubatureInfo,coeff_Stokes_2D)
    
    ! Duplicate the matrix to the Y-velocity matrix, share structure and
    ! content between them (as the matrices are the same).
    call lsyssc_duplicateMatrix (rmatrix%RmatrixBlock(1,1),&
        rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    
    ! Manually change the discretisation structure of the Y-velocity
    ! matrix to the Y-discretisation structure.
    ! Ok, we use the same discretisation structure for both, X- and Y-velocity,
    ! so this is not really necessary - we do this for sure...
    call lsyssc_assignDiscrDirectMat (rmatrix%RmatrixBlock(2,2),&
        rdiscretisation%RspatialDiscr(2))

    ! Build the first pressure matrix B1.
    ! Again first set up the bilinear form, then call the matrix assembly.
    rform%itermCount = 1
    rform%Idescriptors(1,1) = DER_FUNC
    rform%Idescriptors(2,1) = DER_DERIV_X

    ! In the standard case, we have constant coefficients:
    rform%ballCoeffConstant = .true.
    rform%BconstantCoeff = .true.
    rform%Dcoefficients(1)  = -1.0_DP
    
    call bilf_buildMatrixScalar (rform,.true.,rmatrixB1,&
        rcubatureInfo,coeff_Pressure_2D)

    ! Build the second pressure matrix B2.
    ! Again first set up the bilinear form, then call the matrix assembly.
    rform%itermCount = 1
    rform%Idescriptors(1,1) = DER_FUNC
    rform%Idescriptors(2,1) = DER_DERIV_Y

    ! In the standard case, we have constant coefficients:
    rform%ballCoeffConstant = .true.
    rform%BconstantCoeff = .true.
    rform%Dcoefficients(1)  = -1.0_DP
    
    call bilf_buildMatrixScalar (rform,.true.,rmatrixB2,&
        rcubatureInfo,coeff_Pressure_2D)
                                
    ! The B1/B2 matrices exist up to now only in our local problem structure.
    ! Put a copy of them into the block matrix.
    !
    ! Note that we share the structure of B1/B2 with those B1/B2 of the
    ! block matrix, while we create copies of the entries. The reason is
    ! that these matrices are modified for boundary conditions later.
    call lsyssc_duplicateMatrix (rmatrixB1, rmatrix%RmatrixBlock(1,3),&
                                 LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

    call lsyssc_duplicateMatrix (rmatrixB2, rmatrix%RmatrixBlock(2,3),&
                                 LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
    
    ! Furthermore, put B1^T and B2^T to the block matrix.
    call lsyssc_transposeMatrix (rmatrixB1, rmatrix%RmatrixBlock(3,1),&
                                 LSYSSC_TR_VIRTUAL)

    call lsyssc_transposeMatrix (rmatrixB2, rmatrix%RmatrixBlock(3,2),&
                                 LSYSSC_TR_VIRTUAL)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Create RHS and solution vectors
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! Next step: Create a RHS vector, a solution vector and a temporary
    ! vector. All are filled with zero.
    call lsysbl_createVectorBlock (rdiscretisation,rrhs,.true.)
    call lsysbl_createVectorBlock (rdiscretisation,rvector,.true.)
    call lsysbl_createVectorBlock (rdiscretisation,rtempBlock,.true.)

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
        rlinform,.true.,rrhs%RvectorBlock(1),rcubatureInfo,coeff_RHS_X_2D)

    call linf_buildVectorScalar (&
        rlinform,.true.,rrhs%RvectorBlock(2),rcubatureInfo,coeff_RHS_Y_2D)
                                
    ! The third subvector must be zero - as it represents the RHS of
    ! the equation "div(u) = 0".
    call lsyssc_clearVector(rrhs%RvectorBlock(3))
                                
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Assembly of matrices/vectors finished
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Discretise the boundary conditions and apply them to the matrix/RHS/sol.
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! For implementing boundary conditions, we use a `filter technique with
    ! discretised boundary conditions`. This means, we first have to calculate
    ! a discrete version of the analytic BC, which we can implement into the
    ! solution/RHS vectors using the corresponding filter.
    !
    ! Create a t_discreteBC structure where we store all discretised boundary
    ! conditions.
    call bcasm_initDiscreteBC(rdiscreteBC)
    
    ! We first set up the boundary conditions for the X-velocity, then those
    ! of the Y-velocity.
    !
    ! We 'know' already (from the problem definition) that we have four boundary
    ! segments in the domain. Each of these, we want to use for enforcing
    ! some kind of boundary condition.
    !
    ! We ask the boundary routines to create a 'boundary region' - which is
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
    ! - Add the calculated discrete BC`s to rdiscreteBC for later use.
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
                             
    ! Edge 2 is Neumann boundary, so it is commented out.
    ! CALL boundary_createRegion(rboundary,1,2,rboundaryRegion)
    ! CALL bcasm_newDirichletBConRealBD (rdiscretisation,1,&
    !                                    rboundaryRegion,rdiscreteBC,&
    !                                    getBoundaryValues_2D)
                             
    ! Edge 3 of boundary component 1.
    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
    
    ! Edge 4 of boundary component 1. That is it.
    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)

    ! Now continue with defining the boundary conditions of the Y-velocity:
    !
    ! Define edge 1.
    call boundary_createRegion(rboundary,1,1,rboundaryRegion)
    
    ! Edge with start- and endpoint.
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    
    ! As we define the Y-velocity, we now set icomponent=2 in the following call.
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
                             
    ! Edge 2 is Neumann boundary, so it is commented out.
    ! CALL boundary_createRegion(rboundary,1,2,rboundaryRegion)
    ! CALL bcasm_newDirichletBConRealBD (rdiscretisation,2,&
    !                                    rboundaryRegion,rdiscreteBC,&
    !                                    getBoundaryValues_2D)
                             
    ! Edge 3 of boundary component 1.
    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
    
    ! Edge 4 of boundary component 1. That is it.
    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)

    ! The pressure does not need boundary conditions.

    ! Assign the BC`s to the vectors and the matrix. That way, these
    ! boundary conditions are always connected to that matrix and that
    ! vector.
    call lsysbl_assignDiscreteBC(rmatrix,rdiscreteBC)
    call lsysbl_assignDiscreteBC(rrhs,rdiscreteBC)
    call lsysbl_assignDiscreteBC(rvector,rdiscreteBC)
    call lsysbl_assignDiscreteBC(rtempBlock,rdiscreteBC)
    
    ! Next step is to implement boundary conditions into the RHS,
    ! solution and matrix. This is done using a vector/matrix filter
    ! for discrete boundary conditions.
    ! The discrete boundary conditions are already attached to the
    ! vectors/matrix. Call the appropriate vector/matrix filter that
    ! modifies the vectors/matrix according to the boundary conditions.
    call vecfil_discreteBCrhs (rrhs)
    call vecfil_discreteBCsol (rvector)
    call matfil_discreteBC (rmatrix)

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

    ! Create a BiCGStab-solver with VANKA preconditioner.
    ! Attach the above filter chain to the solver, so that the solver
    ! automatically filters the vector during the solution process.
    nullify(p_rpreconditioner)
    call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DNAVST)
    call linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,RfilterChain)

    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2

    ! We will allow the solver to perform 200 iterations
    p_rsolverNode%nmaxIterations = 200

    ! Attach the system matrix to the solver.
    ! First create an array with the matrix data (on all levels, but we
    ! only have one level here), then call the initialisation
    ! routine to attach all these matrices.
    ! Remark: Do not make a call like
    !    CALL linsol_setMatrices(p_RsolverNode,(/p_rmatrix/))
    ! This does not work on all compilers, since the compiler would have
    ! to create a temp array on the stack - which does not always work!
    Rmatrices = (/rmatrix/)
    call linsol_setMatrices(p_rsolverNode,Rmatrices)
    
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
    ! For this purpose, first create a 'derived' simple discretisation
    ! structure based on Q1 by copying the main guiding block discretisation
    ! structure and modifying the discretisation structures of the
    ! two velocity subvectors:
    
    call spdiscr_duplicateBlockDiscr (rdiscretisation,rprjDiscretisation)
    
    call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(1), &
                 EL_Q1, rprjDiscretisation%RspatialDiscr(1))

    call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(2), &
                 EL_Q1, rprjDiscretisation%RspatialDiscr(2))
                 
    ! The pressure discretisation substructure stays the old.
    !
    ! Now set up a new solution vector based on this discretisation,
    ! allocate memory.
    call lsysbl_createVecBlockByDiscr (rprjDiscretisation,rprjVector,.false.)
    
    ! Then take our original solution vector and convert it according to the
    ! new discretisation:
    call spdp_projectSolution (rvector,rprjVector)
    
    ! Discretise the boundary conditions according to the Q1/Q1/Q0
    ! discretisation.
    !
    ! Create a t_discreteBC structure where we store all discretised boundary
    ! conditions.
    call bcasm_initDiscreteBC(rprjDiscreteBC)
    !
    ! Edge 1 of boundary component 1, X-velocity.
    call boundary_createRegion(rboundary,1,1,rboundaryRegion)

    ! Edge with start- and endpoint.
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    
    call bcasm_newDirichletBConRealBD (rprjDiscretisation,1,&
                                       rboundaryRegion,rprjDiscreteBC,&
                                       getBoundaryValues_2D)
                             
    ! Edge 2 is Neumann boundary, so it is commented out.
    ! CALL boundary_createRegion(rboundary,1,2,rboundaryRegion)
    ! CALL bcasm_newDirichletBConRealBD (rprjDiscretisation,1,&
    !                                    rboundaryRegion,rprjDiscreteBC,&
    !                                    getBoundaryValues_2D)
                             
    ! Edge 3 of boundary component 1.
    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rprjDiscretisation,1,&
                                       rboundaryRegion,rprjDiscreteBC,&
                                       getBoundaryValues_2D)
    
    ! Edge 4 of boundary component 1. That is it.
    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rprjDiscretisation,1,&
                                       rboundaryRegion,rprjDiscreteBC,&
                                       getBoundaryValues_2D)

    ! Edge 1 of boundary component 1, Y-velocity.
    call boundary_createRegion(rboundary,1,1,rboundaryRegion)
  
    ! Edge with start- and endpoint.
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    
    ! As we define the Y-velocity, we now set icomponent=2 in the following call.
    call bcasm_newDirichletBConRealBD (rprjDiscretisation,2,&
                                       rboundaryRegion,rprjDiscreteBC,&
                                       getBoundaryValues_2D)
                             
    ! Edge 2 is Neumann boundary, so it is commented out.
    ! CALL boundary_createRegion(rboundary,1,2,rboundaryRegion)
    ! CALL bcasm_newDirichletBConRealBD (rprjDiscretisation,2,&
    !                                    rboundaryRegion,rprjDiscreteBC,&
    !                                    getBoundaryValues_2D)
                             
    ! Edge 3 of boundary component 1.
    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rprjDiscretisation,2,&
                                       rboundaryRegion,rprjDiscreteBC,&
                                       getBoundaryValues_2D)
    
    ! Edge 4 of boundary component 1. That is it.
    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rprjDiscretisation,2,&
                                       rboundaryRegion,rprjDiscreteBC,&
                                       getBoundaryValues_2D)

    ! Hang the pointer into the vector.
    rprjVector%p_rdiscreteBC => rprjDiscreteBC

    ! Send the vector to the boundary-condition implementation filter.
    ! This modifies the vector according to the discrete boundary
    ! conditions.
    call vecfil_discreteBCsol (rprjVector)
    
    ! Get the path for writing postprocessing files from the environment variable
    ! $UCDDIR. If that does not exist, write to the directory "./gmv".
    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = './gmv'

    ! Now we have a Q1/Q1/Q0 solution in rprjVector.
    ! We can now start the postprocessing.
    ! Start UCD export to VTK file:
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
        trim(sucddir)//'/u2d_0_simple.vtk')

    ! Write velocity field
    call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
    call lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
    
    ! In case we use the VTK exporter, which supports vector output, we will
    ! pass the X- and Y-velocity at once to the ucd module.
    call ucd_addVarVertBasedVec(rexport,'velocity',p_Ddata,p_Ddata2)

    ! If we use the GMV exporter, we might replace the line above by the
    ! following two lines:
    !CALL ucd_addVariableVertexBased (rexport,'X-vel',UCD_VAR_XVELOCITY, p_Ddata)
    !CALL ucd_addVariableVertexBased (rexport,'Y-vel',UCD_VAR_YVELOCITY, p_Ddata2)
        
    ! Write pressure
    call lsyssc_getbase_double (rprjVector%RvectorBlock(3),p_Ddata)
    call ucd_addVariableElementBased (rexport,'pressure',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that is it.
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! Store the viscosity parameter nu in the collection's quick access array
    rcollection%DquickAccess(1) = dnu

    ! Set up the error structure for velocity
    rerrorU%p_RvecCoeff => rvector%RvectorBlock(1:2)
    rerrorU%p_DerrorL2 => DerrorUL2
    rerrorU%p_DerrorH1 => DerrorUH1

    ! Set up the error structure for pressure
    rerrorP%p_RvecCoeff => rvector%RvectorBlock(3:3)
    rerrorP%p_DerrorL2 => DerrorPL2

    ! Calculate errors of velocity and pressure against analytic solutions.
    call pperr_scalarVec(rerrorU, funcVelocity2D, rcollection, rcubatureInfo)
    call pperr_scalarVec(rerrorP, funcPressure2D, rcollection, rcubatureInfo)

    ! Print the errors.
    call output_lbrk()
    call output_line('|u - u_h|_L2 = ' // trim(sys_sdEL(DerrorUL2(1), 10)) &
                                // ' ' // trim(sys_sdEL(DerrorUL2(2), 10)))
    call output_line('|u - u_h|_H1 = ' // trim(sys_sdEL(DerrorUH1(1), 10)) &
                                // ' ' // trim(sys_sdEL(DerrorUH1(2), 10)))
    call output_line('|p - p_h|_L2 = ' // trim(sys_sdEL(derrorPL2(1), 10)))

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Clean up
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! We are finished - but not completely!
    ! Now, clean up so that all the memory is available again.
    !
    ! Release solver data and structure
    call linsol_doneData (p_rsolverNode)
    call linsol_doneStructure (p_rsolverNode)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    call linsol_releaseSolver (p_rsolverNode)
    
    ! Release the block matrix/vectors
    call lsysbl_releaseVector (rprjVector)
    call lsysbl_releaseVector (rtempBlock)
    call lsysbl_releaseVector (rvector)
    call lsysbl_releaseVector (rrhs)
    call lsysbl_releaseMatrix (rmatrix)
    
    ! Release B1 and B2 matrix
    call lsyssc_releaseMatrix (rmatrixB2)
    call lsyssc_releaseMatrix (rmatrixB1)
    
    ! Release the cubature info structure.
    call spdiscr_releaseCubStructure(rcubatureInfo)

    ! Release our discrete version of the boundary conditions
    call bcasm_releaseDiscreteBC (rprjDiscreteBC)
    call bcasm_releaseDiscreteBC (rdiscreteBC)

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    call spdiscr_releaseBlockDiscr(rprjDiscretisation)
    call spdiscr_releaseBlockDiscr(rdiscretisation)
    
    ! Release the triangulation.
    call tria_done (rtriangulation)
    
    ! Finally release the domain, that is it.
    call boundary_release (rboundary)

  end subroutine

end module
