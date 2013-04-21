!##############################################################################
!# ****************************************************************************
!# <name> stokes3d_method0_simple </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstation program how to solve a Stokes
!# problem on a simple domain.
!#
!# The routine uses the BiCGStab solver with a simple-VANKA preconditioner
!# for 3D saddle point problems, Jacobi-Type on a single grid.
!# </purpose>
!##############################################################################

module stokes3d_method0_simple

  use fsystem
  use storage
  use genoutput
  use boundary
  use cubature
  use derivatives
  use matrixfilters
  use vectorfilters
  use linearalgebra
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
  use meshregion
  
  use stokes3d_callback
  use dom3d_cube
  
  implicit none

contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine stokes3d_0_simple
  
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
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! An object for saving the boundary mesh region
    type(t_meshregion) :: rmeshRegion

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
    type(t_matrixScalar) :: rmatrixB1, rmatrixB2, rmatrixB3

    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_matrixBlock) :: rmatrix
    type(t_vectorBlock) :: rvector,rrhs,rtempBlock,rprjVector
    
    ! A set of variables describing the analytic and discrete boundary
    ! conditions.
    type(t_discreteBC), target :: rdiscreteBC, rprjDiscreteBC

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode,p_rpreconditioner

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    type(t_filterChain), dimension(1), target :: RfilterChain
    
    ! Number of filters in the filter chain
    integer :: nfilters
    
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
    
    ! Path to the mesh
    character(len=SYS_STRLEN) :: spredir
    
    ! Ok, let us start.
    !
    ! We want to solve our Poisson problem on level...
    ! As we do not use a multigrid solver here, we will set the level
    ! to 4 instead of 7.
    NLMAX = 4
    
    ! Viscosity parameter:
    dnu = 1.0_DP

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Read the domain, read the mesh, refine
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Get the path $PREDIR from the environment, where to read .prm/.tri files
    ! from. If that does not exist, write to the directory "./pre".
    if (.not. sys_getenv_string("PREDIR", spredir)) spredir = "./pre"

    ! At first, read in the basic triangulation.
    call tria_readTriFile3D (rtriangulation, trim(spredir)//"/CUBE.tri")
    
    ! Refine the mesh up to the minimum level
    call tria_quickRefine2LevelOrdering(NLMAX-1,rtriangulation)
    
    ! Create information about adjacencies and everything one needs from
    ! a triangulation. Afterwards, we have the coarse mesh.
    call tria_initStandardMeshFromRaw (rtriangulation)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up a discretisation structure which tells the code which
    ! finite element to use
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies 4 blocks in the
    ! solution vector.
    call spdiscr_initBlockDiscr (rdiscretisation,4,rtriangulation)

    ! rdiscretisation%RspatialDiscr is a list of scalar
    ! discretisation structures for every component of the solution vector.
    ! We have a solution vector with three components:
    !  Component 1 = X-velocity
    !  Component 2 = Y-velocity
    !  Component 3 = Z-velocity
    !  Component 4 = Pressure
    ! For simplicity, we set up one discretisation structure for the
    ! velocity...
    call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1),&
                                   EL_EM30_3D, rtriangulation)
                
    ! ...and copy this structure also to the discretisation structure
    ! of the 2nd and 3rd component (Y-/Z-velocity). This needs no additional
    ! memory, as all three structures will share the same dynamic
    ! information afterwards.
    call spdiscr_duplicateDiscrSc(rdiscretisation%RspatialDiscr(1),&
        rdiscretisation%RspatialDiscr(2))
    call spdiscr_duplicateDiscrSc(rdiscretisation%RspatialDiscr(1),&
        rdiscretisation%RspatialDiscr(3))

    ! For the pressure (4th component), we set up a separate discretisation
    ! structure, as this uses different finite elements for trial and test
    ! functions.
    call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(1), &
        EL_Q0_3D, rdiscretisation%RspatialDiscr(4))

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

    ! In the Stokes problem, the matrices for the Y-/Z-velocity are identical
    ! to the matrix for the X-velocity; all are Laplace-matrices!
    ! Therefore, we can simply make a copy of the matrix for the X-velocity.
    ! This we do later after the entries are created.
    !
    ! In the global system, there are three coupling matrices B1, B2 and B§.
    ! Both have the same structure.
    !
    !    / A              B1 \
    !    |      A         B2 |
    !    |           A    B3 |
    !    \ B1^T B2^T B3^T    /
    !
    ! Create the matrices structure of the pressure using the 4th
    ! spatial discretisation structure in p_rdiscretisation%RspatialDiscr.
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(4),&
        LSYSSC_MATRIX9, rmatrixB1,rdiscretisation%RspatialDiscr(1))
              
    ! Duplicate the B1 matrix structure to the B2/B3 matrix, so use
    ! lsyssc_duplicateMatrix to create B2/B3. Share the matrix
    ! structure between B1, B2 and B3 (B1 is the parent and B2/B3 the child).
    ! Do not create a content array yet, it will be created by
    ! the assembly routines later.
    call lsyssc_duplicateMatrix (rmatrixB1, rmatrixB2, &
        LSYSSC_DUP_COPY,LSYSSC_DUP_REMOVE)
    call lsyssc_duplicateMatrix (rmatrixB1, rmatrixB3, &
        LSYSSC_DUP_COPY,LSYSSC_DUP_REMOVE)
                                     
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
    call bilf_buildMatrixScalar (rform,.true.,rmatrix%RmatrixBlock(1,1),&
        rcubatureInfo,coeff_Stokes_3D)
    
    ! Duplicate the matrix to the Y-/Z-velocity matrix, share structure and
    ! content between them (as the matrices are the same).
    call lsyssc_duplicateMatrix (rmatrix%RmatrixBlock(1,1),&
        rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

    call lsyssc_duplicateMatrix (rmatrix%RmatrixBlock(1,1),&
        rmatrix%RmatrixBlock(3,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    
    ! Manually change the discretisation structure of the Y-/Z-velocity
    ! matrix to the Y-/Z-discretisation structure.
    ! Ok, we use the same discretisation structure for both, X-,Y- and
    ! Z-velocity, so this is not really necessary - we do this for sure...
    call lsyssc_assignDiscretisation (rmatrix%RmatrixBlock(2,2),&
        rdiscretisation%RspatialDiscr(2))
    call lsyssc_assignDiscretisation (rmatrix%RmatrixBlock(3,3),&
        rdiscretisation%RspatialDiscr(3))
                                
    ! Build the first pressure matrix B1.
    ! Again first set up the bilinear form, then call the matrix assembly.
    rform%itermCount = 1
    rform%Idescriptors(1,1) = DER_FUNC3D
    rform%Idescriptors(2,1) = DER_DERIV3D_X

    ! In the standard case, we have constant coefficients:
    rform%ballCoeffConstant = .true.
    rform%Dcoefficients(1)  = -1.0_DP
    
    call bilf_buildMatrixScalar (rform,.true.,rmatrixB1,&
        rcubatureInfo,coeff_Pressure_3D)

    ! Build the second pressure matrix B2.
    ! Again first set up the bilinear form, then call the matrix assembly.
    rform%itermCount = 1
    rform%Idescriptors(1,1) = DER_FUNC3D
    rform%Idescriptors(2,1) = DER_DERIV3D_Y

    ! In the standard case, we have constant coefficients:
    rform%ballCoeffConstant = .true.
    rform%Dcoefficients(1)  = -1.0_DP
    
    call bilf_buildMatrixScalar (rform,.true.,rmatrixB2,&
        rcubatureInfo,coeff_Pressure_3D)
                                
    ! Build the second pressure matrix B3.
    ! Again first set up the bilinear form, then call the matrix assembly.
    rform%itermCount = 1
    rform%Idescriptors(1,1) = DER_FUNC3D
    rform%Idescriptors(2,1) = DER_DERIV3D_Z

    ! In the standard case, we have constant coefficients:
    rform%ballCoeffConstant = .true.
    rform%Dcoefficients(1)  = -1.0_DP
    
    call bilf_buildMatrixScalar (rform,.true.,rmatrixB3,&
        rcubatureInfo,coeff_Pressure_3D)

    ! The B1/B2 matrices exist up to now only in our local problem structure.
    ! Put a copy of them into the block matrix.
    !
    ! Note that we share the structure of B1/B2/B3 with those B1/B2/B3 of the
    ! block matrix, while we create copies of the entries. The reason is
    ! that these matrices are modified for boundary conditions later.
    call lsyssc_duplicateMatrix (rmatrixB1, rmatrix%RmatrixBlock(1,4),&
                                 LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

    call lsyssc_duplicateMatrix (rmatrixB2, rmatrix%RmatrixBlock(2,4),&
                                 LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

    call lsyssc_duplicateMatrix (rmatrixB3, rmatrix%RmatrixBlock(3,4),&
                                 LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
    
    ! Furthermore, put B1^T and B2^T to the block matrix.
    call lsyssc_transposeMatrix (rmatrixB1, rmatrix%RmatrixBlock(4,1),&
                                 LSYSSC_TR_VIRTUAL)

    call lsyssc_transposeMatrix (rmatrixB2, rmatrix%RmatrixBlock(4,2),&
                                 LSYSSC_TR_VIRTUAL)

    call lsyssc_transposeMatrix (rmatrixB3, rmatrix%RmatrixBlock(4,3),&
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
    rlinform%Idescriptors(1) = DER_FUNC3D
    
    ! ... and then discretise the RHS to the first three subvectors of
    ! the block vector using the discretisation structure of the
    ! corresponding block.
    call linf_buildVectorScalar (&
        rlinform,.true.,rrhs%RvectorBlock(1),rcubatureInfo,coeff_RHS_X_3D)
                  
    call linf_buildVectorScalar (&
        rlinform,.true.,rrhs%RvectorBlock(2),rcubatureInfo,coeff_RHS_Y_3D)
                  
    call linf_buildVectorScalar (&
        rlinform,.true.,rrhs%RvectorBlock(3),rcubatureInfo,coeff_RHS_Z_3D)
                                
    ! The fourth subvector must be zero - as it represents the RHS of
    ! the equation "div(u) = 0".
    call lsyssc_clearVector(rrhs%RvectorBlock(4))
                                
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Assembly of matrices/vectors finished
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Discretise the boundary conditions and apply them to the matrix/RHS/sol.
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Clear the solution vector on the finest level.
    call lsysbl_clearVector(rvector)

    ! Now we need to implement the boundary conditions. To do this, we
    ! first need to create a mesh region describing the mesh`s boundary.
    ! We want to prescribe Dirichlet on the cube`s boundary, except for
    ! the face where the X-coordinate is 1.
    ! We could now manually create a mesh region based on the triangulation`s
    ! nodal-property array and then kick out everything that belongs to the
    ! right face. But we will use the dom3d_cube module, which performs
    ! this task for us.
    call dom3d_cube_calcMeshRegion(rmeshRegion, rtriangulation, &
                                   DOM3D_CUBE_REG_STOKES)
    
    ! Initialise the structure for discrete boundary conditions.
    call bcasm_initDiscreteBC (rdiscreteBC)
    
    ! Assemble Dirichlet BCs for X-velocity:
    call bcasm_newDirichletBConMR(rdiscretisation, 1, rdiscreteBC, &
                       rmeshRegion, getBoundaryValuesMR_3D)

    ! Assemble Dirichlet BCs for Y-velocity:
    call bcasm_newDirichletBConMR(rdiscretisation, 2, rdiscreteBC, &
                       rmeshRegion, getBoundaryValuesMR_3D)

    ! Assemble Dirichlet BCs for Z-velocity:
    call bcasm_newDirichletBConMR(rdiscretisation, 3, rdiscreteBC, &
                       rmeshRegion, getBoundaryValuesMR_3D)
    
    ! We do not need any boundary conditions for the pressure
    
    ! We will not release rmeshRegion yet, as we will need it for the
    ! postprocessing!
                             
    
    ! Next step is to implement boundary conditions into the RHS,
    ! solution and matrix. This is done using a vector/matrix filter
    ! for discrete boundary conditions.
    ! The discrete boundary conditions are already attached to the
    ! vectors/matrix. Call the appropriate vector/matrix filter that
    ! modifies the vectors/matrix according to the boundary conditions.
    call vecfil_discreteBCrhs (rrhs,rdiscreteBC)
    call vecfil_discreteBCsol (rvector,rdiscreteBC)
    call matfil_discreteBC (rmatrix,rdiscreteBC)

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
    call filter_clearFilterChain (RfilterChain,nfilters)
    call filter_newFilterDiscBCDef (RfilterChain,nfilters,rdiscreteBC)

    ! Create a BiCGStab-solver with VANKA preconditioner.
    ! Attach the above filter chain to the solver, so that the solver
    ! automatically filters the vector during the solution process.
    nullify(p_rpreconditioner)
    call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_3DNAVST)
    call linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,RfilterChain)

    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2

    ! We will allow the solver to perform 200 iterations
    p_rsolverNode%nmaxIterations = 200

    ! Attach the system matrix to the solver.
    call linsol_setMatrix(p_rsolverNode,rmatrix)
    
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
    call spdiscr_duplicateBlockDiscr (rdiscretisation,rprjDiscretisation)
    
    call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(1), &
        EL_Q1_3D, rprjDiscretisation%RspatialDiscr(1))
    call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(2), &
        EL_Q1_3D, rprjDiscretisation%RspatialDiscr(2))
    call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(3), &
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
    !
    ! Initialise the structure for discrete boundary conditions.
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
    ! This modifies the vector according to the attached discrete boundary
    ! conditions.
    call vecfil_discreteBCsol (rprjVector,rprjDiscreteBC)
    
    ! Get the path for writing postprocessing files from the environment variable
    ! $UCDDIR. If that does not exist, write to the directory "./gmv".
    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = "./gmv"

    ! Now we have a Q1/Q1/Q1/Q0 solution in rprjVector.
    ! We can now start the postprocessing.
    ! Start UCD export to VTK file:
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
        trim(sucddir)//"/u3d_0_simple.vtk")

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
    
    ! Release B1, B2 and B3 matrices
    call lsyssc_releaseMatrix (rmatrixB3)
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
    
  end subroutine

end module
