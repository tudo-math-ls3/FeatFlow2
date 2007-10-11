!##############################################################################
!# ****************************************************************************
!# <name> anisotropicdiffusion_method1 </name>
!# ****************************************************************************
!#
!# <purpose>
!# This program is a simple test program for discretising the equation
!#
!#              - div (A grad(u) = f
!#
!# on a 2D domain for a scalar function u and a matrix
!#
!#   A =  ( cos(t) -sin(t) )  ( a11 a12 )  (  cos(t) sin(t) ) 
!#        ( sin(t)  cos(t) )  ( a21 a22 )  ( -sin(t) cos(t) ) 
!#
!# All parameters (e.g. a11,...,a22) are read from a .DAT file.
!#
!# The example  discretises and solves this equation in a direct way, 
!# just listing all commands necessary  for initialisation, discretisation, 
!# solving (with Gauss elimination = UMFPACK) and cleanup.
!# </purpose>
!##############################################################################

MODULE anisotropicdiffusion_method1

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
  USE genoutput
  USE statistics
  USE paramlist
  USE collection
  USE meshmodification
  USE spdiscprojection
  
  USE jumpstabilisation
    
  USE anisotropicdiffusion_callback
  
  IMPLICIT NONE

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE anisotropicdiffusion1
  
!<description>
  ! This is an all-in-one poisson solver for directly solving a Poisson
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
    TYPE(t_blockDiscretisation) :: rdiscretisation
    
    ! A bilinear and linear form describing the analytic problem to solve
    TYPE(t_bilinearForm) :: rform
    TYPE(t_linearForm) :: rlinform
    
    ! A scalar matrix and vector. The vector accepts the RHS of the problem
    ! in scalar form.
    TYPE(t_matrixScalar) :: rmatrix
    TYPE(t_vectorScalar) :: rrhs

    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    TYPE(t_matrixBlock) :: rmatrixBlock
    TYPE(t_vectorBlock) :: rvectorBlock,rrhsBlock,rtempBlock

    ! A set of variables describing the analytic and discrete boundary
    ! conditions.    
    TYPE(t_boundaryConditions), POINTER :: p_rboundaryConditions
    TYPE(t_boundaryRegion) :: rboundaryRegion
    TYPE(t_bcRegion), POINTER :: p_rbcRegion
    TYPE(t_discreteBC), POINTER :: p_rdiscreteBC

    ! A solver node that accepts parameters for the linear solver    
    TYPE(t_linsolNode), POINTER :: p_rsolverNode

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    TYPE(t_matrixBlock), DIMENSION(1) :: Rmatrices

    ! NLMAX receives the level where we want to solve.
    INTEGER :: NLMAX
    
    ! Rotation angle, X-/Y-diffusion
    REAL(DP) :: dtheta,ddiffusionA11,ddiffusionA12
    REAL(DP) :: ddiffusionA21,ddiffusionA22
    
    ! Type of solution
    INTEGER :: isolution
    
    ! Stabilisation
    REAL(DP) :: dgamma
    INTEGER :: istabilisation
    
    ! Mesh distortion
    REAL(DP) :: dmeshDistortion
    
    ! Element type of the discretisation
    INTEGER :: ieltype
    
    ! Whether to convert to triangle mesh
    INTEGER :: iconvertToTriangleMesh
    
    ! Final diffusion matrix after rotation
    REAL(DP), DIMENSION(2,2) :: DdiffusionMatrix 
    
    ! PRM/TRI file
    CHARACTER(LEN=SYS_STRLEN) :: sstring,sfilePRM,sfileTRI
    
    ! Error indicator during initialisation of the solver
    INTEGER :: ierror   
    
    ! Error of FE function to reference function
    REAL(DP) :: derror
    
    ! Output block for UCD output to GMV file
    TYPE(t_ucdExport) :: rexport
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_DdataQ1
    TYPE(t_discreteBC), POINTER :: p_rdiscreteBCPostProc
    TYPE(t_vectorScalar) :: rvectorPostProc
    TYPE(t_vectorBlock) :: rvectorPostProcBlock
    TYPE(t_blockDiscretisation) :: rdiscretisationPostProc
    
    ! Structure for saving parameters from the DAT file
    TYPE (t_parlist) :: rparams
    
    ! Collection structure for setting up the RHS
    TYPE(t_collection) :: rcollection

    ! Ok, let's start. 
    !
    ! +------------------------------------------------------------------------
    ! | READ DAT FILE PARAMETERS
    ! +------------------------------------------------------------------------
    !
    ! Initialise the parameter structure and read the DAT file.
    CALL parlst_init(rparams)
    CALL parlst_readfromfile (rparams, './dat/anisotropicdiffusion.dat')
    
    ! Get the parameters...
    !
    ! PRM file
    CALL parlst_getvalue_string (rparams, '', &
                                 'sfilePRM', sstring)
    READ(sstring,*) sfilePRM
                                 
    ! TRI file
    CALL parlst_getvalue_string (rparams, '', &
                                 'sfileTRI', sstring)
    READ(sstring,*) sfileTRI

    ! Whether to convert to a triangle mesh
    CALL parlst_getvalue_int (rparams, '', &
      'iconvertToTriangleMesh', iconvertToTriangleMesh, 0)

    ! Mesh distortion
    CALL parlst_getvalue_double (rparams, '', &
      'dmeshDistortion', dmeshDistortion, 0.0_DP)

    ! Element type
    CALL parlst_getvalue_int (rparams, '', &
                              'ielementType', ieltype, 7)

    ! Type of stabilisation
    CALL parlst_getvalue_int (rparams, '', &
      'istabilisation', istabilisation, 0)
    CALL parlst_getvalue_double (rparams, '', &
      'dgamma', dgamma, 0.01_DP)

    ! Maximum level
    CALL parlst_getvalue_int (rparams, '', &
                                     'NLMAX', NLMAX, 7)

    ! Rotation angle
    CALL parlst_getvalue_double (rparams, '', 'dtheta', dtheta, 0.0_DP)
    
    ! in rad...
    dtheta = 2.0_DP*SYS_PI*dtheta/360.0_DP
    
    ! diffusion matrix
    CALL parlst_getvalue_double (rparams, '', &
        'ddiffusionA11', ddiffusionA11, 1.0_DP)
    CALL parlst_getvalue_double (rparams, '', &
        'ddiffusionA12', ddiffusionA12, 0.0_DP)
    CALL parlst_getvalue_double (rparams, '', &
        'ddiffusionA21', ddiffusionA21, 0.0_DP)
    CALL parlst_getvalue_double (rparams, '', &
        'ddiffusionA22', ddiffusionA22, 1.0_DP)
        
    ! solution
    CALL parlst_getvalue_int (rparams, '', &
                              'isolution', isolution, 0)
    
    ! Create the actual diffusion matrix:
    !
    ! A = (  cos t  -sin t )  ( d11 d12 ) (  cos t   sin t )
    !     (  sin t   cos t )  ( d21 d22 ) ( -sin t   cos t )
    !
    DdiffusionMatrix(1,1) = &
      (cos(dtheta)*ddiffusionA11-sin(dtheta)*ddiffusionA21)*cos(dtheta)-&
      (cos(dtheta)*ddiffusionA12-sin(dtheta)*ddiffusionA22)*sin(dtheta)

    DdiffusionMatrix(1,2) = &
        (cos(dtheta)*ddiffusionA11-sin(dtheta)*ddiffusionA21)*sin(dtheta)+&
        (cos(dtheta)*ddiffusionA12-sin(dtheta)*ddiffusionA22)*cos(dtheta)

    DdiffusionMatrix(2,1) = &
       (sin(dtheta)*ddiffusionA11+cos(dtheta)*ddiffusionA21)*cos(dtheta)-&
       (sin(dtheta)*ddiffusionA12+cos(dtheta)*ddiffusionA22)*sin(dtheta)

    DdiffusionMatrix(2,2) = &
        (sin(dtheta)*ddiffusionA11+cos(dtheta)*ddiffusionA21)*sin(dtheta)+&
        (sin(dtheta)*ddiffusionA12+cos(dtheta)*ddiffusionA22)*cos(dtheta)
    
    ! +------------------------------------------------------------------------
    ! | BOUNDARY AND TRIANGULATION
    ! +------------------------------------------------------------------------
    !
    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    ! Set p_rboundary to NULL to create a new structure on the heap.
    NULLIFY(p_rboundary)
    CALL boundary_read_prm(p_rboundary, sfilePRM)
        
    ! Now read in the basic triangulation.
    CALL tria_readTriFile2D (rtriangulation, sfileTRI, p_rboundary)
    
    ! Convert to triangle mesh?
    IF (iconvertToTriangleMesh .EQ. 1) THEN
      CALL tria_rawGridToTri(rtriangulation)
    END IF
    
    ! Refine it.
    CALL tria_quickRefine2LevelOrdering (NLMAX-1,rtriangulation,p_rboundary)

    ! Mesh distortion?
    IF (dmeshDistortion .NE. 0.0_DP) THEN
      CALL meshmod_disturbMesh (rtriangulation,dmeshDistortion)
    END IF
    
    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    CALL tria_initStandardMeshFromRaw (rtriangulation,p_rboundary)
    
    ! +------------------------------------------------------------------------
    ! | SETUP OF THE DISCRETISATION
    ! +------------------------------------------------------------------------
    !
    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we only have one block.
    CALL spdiscr_initBlockDiscr2D (rdiscretisation,1,&
                                   rtriangulation, p_rboundary)
    
    ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
    ! structures for every component of the solution vector.
    ! Initialise the first element of the list to specify the element
    ! and cubature rule for this solution component:
    SELECT CASE (ieltype)
    CASE (1)
      CALL spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscretisation(1), &
        EL_E001,SPDISC_CUB_AUTOMATIC,rtriangulation, p_rboundary)
    CASE (11)
      CALL spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscretisation(1), &
        EL_E011,SPDISC_CUB_AUTOMATIC,rtriangulation, p_rboundary)
    CASE (13)
      CALL spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscretisation(1), &
        EL_E013,SPDISC_CUB_AUTOMATIC,rtriangulation, p_rboundary)
    CASE (-30)
      CALL spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscretisation(1), &
        EL_EM30,SPDISC_CUB_AUTOMATIC,rtriangulation, p_rboundary)
    END SELECT
                 
    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! We create a scalar matrix, based on the discretisation structure
    ! for our one and only solution component.
    SELECT CASE (istabilisation)
    CASE (1)
      ! Jump stabilisation. Needs an extended matrix stencil.
      CALL bilf_createMatrixStructure (&
          rdiscretisation%RspatialDiscretisation(1),&
          LSYSSC_MATRIX9,rmatrix,BILF_MATC_EDGEBASED)
    CASE DEFAULT
      ! No stabilisation
      CALL bilf_createMatrixStructure (&
          rdiscretisation%RspatialDiscretisation(1),&
          LSYSSC_MATRIX9,rmatrix)
    END SELECT
    
    ! +------------------------------------------------------------------------
    ! | DISCRETISATION OF MATRICES/VECTORS
    ! +------------------------------------------------------------------------
    !
    ! And now to the entries of the matrix. For assembling of the entries,
    ! we need a bilinear form, which first has to be set up manually.
    ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
    ! scalar system matrix in 2D.
    
    rform%itermCount = 4
    rform%Idescriptors(1,1) = DER_DERIV_X
    rform%Idescriptors(2,1) = DER_DERIV_X
    
    rform%Idescriptors(1,2) = DER_DERIV_X
    rform%Idescriptors(2,2) = DER_DERIV_Y
    
    rform%Idescriptors(1,3) = DER_DERIV_Y
    rform%Idescriptors(2,3) = DER_DERIV_X
    
    rform%Idescriptors(1,4) = DER_DERIV_Y
    rform%Idescriptors(2,4) = DER_DERIV_Y

    ! In the standard case, we have constant coefficients:
    rform%ballCoeffConstant = .TRUE.
    rform%BconstantCoeff = .TRUE.
    rform%Dcoefficients(1)  = DdiffusionMatrix(1,1)
    rform%Dcoefficients(2)  = DdiffusionMatrix(1,2)
    rform%Dcoefficients(3)  = DdiffusionMatrix(2,1)
    rform%Dcoefficients(4)  = DdiffusionMatrix(2,2)

    ! Now we can build the matrix entries.
    ! We specify the callback function coeff_Laplace for the coefficients.
    ! As long as we use constant coefficients, this routine is not used.
    ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
    ! the framework will call the callback routine to get analytical
    ! data.
    SELECT CASE (istabilisation)
    CASE (1)
      ! Jump stabilisation. CReate the Laplace matrix and add the
      ! stabilisation
      CALL bilf_buildMatrixScalar (rform,.TRUE.,rmatrix,coeff_Laplace)
      CALL jstab_calcUEOJumpStabilisation (&
          rmatrix,dgamma,dgamma,1.0_DP,CUB_G2_1D,1.0_DP)
    CASE DEFAULT    
      ! No stabilisation. Create the Laplace matrix directly.
      CALL bilf_buildMatrixScalar (rform,.TRUE.,rmatrix,coeff_Laplace)
    END SELECT
    
    ! Create a collection; used for passing parameters to the RHS.
    CALL collct_init (rcollection)
    
    ! Put the parameters to the quick-access array
    rcollection%IquickAccess(1) = isolution
    rcollection%DquickAccess(1) = DdiffusionMatrix(1,1)
    rcollection%DquickAccess(2) = DdiffusionMatrix(1,2)
    rcollection%DquickAccess(3) = DdiffusionMatrix(2,1)
    rcollection%DquickAccess(4) = DdiffusionMatrix(2,2)
    
    ! The same has to be done for the right hand side of the problem.
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    
    ! ... and then discretise the RHS to get a discrete version of it.
    ! Again we simply create a scalar vector based on the one and only
    ! discretisation structure.
    ! This scalar vector will later be used as the one and only first
    ! component in a block vector.
    CALL linf_buildVectorScalar (rdiscretisation%RspatialDiscretisation(1),&
                                 rlinform,.TRUE.,rrhs,coeff_RHS,rcollection)
    
    ! The linear solver only works for block matrices/vectors - but above,
    ! we created scalar ones. So the next step is to make a 1x1 block
    ! system from the matrices/vectors above which the linear solver
    ! understands.
    CALL lsysbl_createMatFromScalar (rmatrix,rmatrixBlock,rdiscretisation)
    CALL lsysbl_createVecFromScalar (rrhs,rrhsBlock,rdiscretisation)
    
    ! +------------------------------------------------------------------------
    ! | DISCRETISATION AND IMPLEMENTATION OF BOUNDARY CONDITIONS
    ! +------------------------------------------------------------------------
    !
    ! Now we have the raw problem. What is missing is the definition of the boudary
    ! conditions.
    ! For implementing boundary conditions, we use a 'filter technique with
    ! discretised boundary conditions'. This means, we first have to calculate
    ! a discrete version of the analytic BC, which we can implement into the
    ! solution/RHS vectors using the corresponding filter.
    !
    ! At first, we need the analytic description of the boundary conditions.
    ! Initialise a structure for boundary conditions which accepts this:
    !
    ! Set p_rboundaryConditions to create a new structure on the heap.
    NULLIFY (p_rboundaryConditions)
    CALL bcond_initBC (p_rboundaryConditions,p_rboundary)
    
    SELECT CASE (isolution)
    CASE (0,1,3)
      ! We 'know' already (from the problem definition) that we have four boundary
      ! segments in the domain. Each of these, we want to use for inforcing
      ! some kind of boundary condition.
      !
      ! We ask the bondary routines to create a 'boundary region' - which is
      ! simply a part of the boundary corresponding to a boundary segment.
      ! A boundary region roughly contains the type, the min/max parameter value
      ! and whether the endpoints are inside the region or not.
      CALL boundary_createRegion(p_rboundary,1,1,rboundaryRegion)
      
      ! We use this boundary region and specify that we want to have Dirichlet
      ! boundary there. The following routine adds a new 'boundary condition region'
      ! for the first segment to the boundary condition structure.
      ! The region will be set up as 'Dirichlet boundary'.
      ! We specify icomponent='1' to indicate that we set up the
      ! Dirichlet BC's for the first (here: one and only) component in the solution
      ! vector.
      ! The routine also returns the created object in p_rbcRegion so that we can
      ! modify it - but accept it as it is, so we can ignore that.
      CALL bcond_newDirichletBConRealBD (p_rboundaryConditions,1,&
                                        rboundaryRegion,p_rbcRegion)
                               
      ! Now to the edge 2 of boundary component 1 the domain. We use the
      ! same two routines to add the boundary condition to p_rboundaryConditions.
      CALL boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
      CALL bcond_newDirichletBConRealBD (p_rboundaryConditions,1,&
                                        rboundaryRegion,p_rbcRegion)
                               
      ! Edge 3 of boundary component 1.
      CALL boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
      CALL bcond_newDirichletBConRealBD (p_rboundaryConditions,1,&
                                        rboundaryRegion,p_rbcRegion)
      
      ! Edge 4 of boundary component 1. That's it.
      CALL boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
      CALL bcond_newDirichletBConRealBD (p_rboundaryConditions,1,&
                                        rboundaryRegion,p_rbcRegion)
    CASE (2,4)
      ! Edge 1 of boundary component 1.
      CALL boundary_createRegion(p_rboundary,1,1,rboundaryRegion)
      CALL bcond_newDirichletBConRealBD (p_rboundaryConditions,1,&
                                        rboundaryRegion,p_rbcRegion)
                               
      ! Now to the edge 2 of boundary component 1 the domain. We use the
      ! same two routines to add the boundary condition to p_rboundaryConditions.
      CALL boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
      CALL bcond_newDirichletBConRealBD (p_rboundaryConditions,1,&
                                        rboundaryRegion,p_rbcRegion)
                               
      ! Edge 3 of boundary component 1.
      CALL boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
      CALL bcond_newDirichletBConRealBD (p_rboundaryConditions,1,&
                                        rboundaryRegion,p_rbcRegion)
      
      ! Edge 4 of boundary component 1. That's it.
      CALL boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
      CALL bcond_newDirichletBConRealBD (p_rboundaryConditions,1,&
                                        rboundaryRegion,p_rbcRegion)

      ! Edge 1 of boundary component 2.
      CALL boundary_createRegion(p_rboundary,2,1,rboundaryRegion)
      CALL bcond_newDirichletBConRealBD (p_rboundaryConditions,1,&
                                        rboundaryRegion,p_rbcRegion)
                               
      ! Edge 2 of boundary component 2.
      CALL boundary_createRegion(p_rboundary,2,2,rboundaryRegion)
      CALL bcond_newDirichletBConRealBD (p_rboundaryConditions,1,&
                                        rboundaryRegion,p_rbcRegion)
                               
      ! Edge 3 of boundary component 3.
      CALL boundary_createRegion(p_rboundary,2,3,rboundaryRegion)
      CALL bcond_newDirichletBConRealBD (p_rboundaryConditions,1,&
                                        rboundaryRegion,p_rbcRegion)
      
      ! Edge 4 of boundary component 4. That's it.
      CALL boundary_createRegion(p_rboundary,2,4,rboundaryRegion)
      CALL bcond_newDirichletBConRealBD (p_rboundaryConditions,1,&
                                        rboundaryRegion,p_rbcRegion)

    END SELECT
                             
    ! The boundary conditions are set up, but still the block discretisation
    ! does not know about it. So inform the discretisation which
    ! analytic boundary conditions to use:
    rdiscretisation%p_rboundaryConditions => p_rboundaryConditions

    ! For the discrete problem, we need a discrete version of the above
    ! boundary conditions. So we have to discretise them.
    ! The following routine gives back p_rdiscreteBC, a pointer to a
    ! discrete version of the boundary conditions. Remark that
    ! the pointer has to be nullified before calling the routine,
    ! otherwise, the routine tries to update the boundary conditions
    ! in p_rdiscreteBC!
    NULLIFY(p_rdiscreteBC)
    CALL bcasm_discretiseBC (rdiscretisation,p_rdiscreteBC,.FALSE., &
                             getBoundaryValues,rcollection)
                             
    ! Hang the pointer into the vector and matrix. That way, these
    ! boundary conditions are always connected to that matrix and that
    ! vector.
    rmatrixBlock%p_rdiscreteBC => p_rdiscreteBC
    rrhsBlock%p_rdiscreteBC => p_rdiscreteBC
                             
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
    CALL matfil_discreteBC (rmatrixBlock)
    
    ! +------------------------------------------------------------------------
    ! | INVOKE THE SOLVER
    ! +------------------------------------------------------------------------
    !
    ! Initialise the solver
    CALL linsol_initUMFPACK4 (p_rsolverNode)
    
    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2
    
    ! Attach the system matrix to the solver.
    ! First create an array with the matrix data (on all levels, but we
    ! only have one level here), then call the initialisation 
    ! routine to attach all these matrices.
    ! Remark: Don't make a call like
    !    CALL linsol_setMatrices(p_RsolverNode,(/p_rmatrix/))
    ! This doesn't work on all compilers, since the compiler would have
    ! to create a temp array on the stack - which does not always work!
    Rmatrices = (/rmatrixBlock/)
    CALL linsol_setMatrices(p_RsolverNode,Rmatrices)
    
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
    CALL linsol_solveAdaptively (p_rsolverNode,rvectorBlock,&
        rrhsBlock,rtempBlock)
    
    ! +------------------------------------------------------------------------
    ! | POSTPROCESSING
    ! +------------------------------------------------------------------------
    !
    ! That's it, rvectorBlock now contains our solution. We can now
    ! start the postprocessing. 

    CALL lsyssc_getbase_double (rvectorBlock%RvectorBlock(1),p_Ddata)

    ! Start UCD export to GMV file:
    SELECT CASE (ieltype)
    CASE (1,11)
      CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,'gmv/u1.gmv')
      CALL ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
    CASE (2)
      CALL ucd_startGMV (rexport,UCD_FLAG_ONCEREFINED,rtriangulation,'gmv/u1.gmv')
      CALL ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, &
          p_Ddata(1:rtriangulation%NVT),&
          p_Ddata(rtriangulation%NVT+1:rtriangulation%NVT+rtriangulation%NMT))
    CASE (13)
      CALL ucd_startGMV (rexport,UCD_FLAG_ONCEREFINED,rtriangulation,'gmv/u1.gmv')
      CALL ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, &
          p_Ddata(1:rtriangulation%NVT),&
          p_Ddata(rtriangulation%NVT+1:rtriangulation%NVT+rtriangulation%NMT),&
          p_Ddata(rtriangulation%NVT+rtriangulation%NMT+1:&
            rtriangulation%NVT+rtriangulation%NMT+rtriangulation%NEL))
    CASE (-30)
      !CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,'gmv/u1.gmv')
      CALL ucd_startGMV (rexport,UCD_FLAG_BULBQUADRATIC,rtriangulation,'gmv/u1.gmv')
      
      ! Project the solution to P1/Q1 -> rvectorPostProc
      CALL spdiscr_initBlockDiscr2D (rdiscretisationPostProc,1,&
          rtriangulation, p_rboundary, p_rboundaryConditions)
      CALL spdp_stdProjectionToP1Q1Scalar (rvectorBlock%RvectorBlock(1),&
          rvectorPostProc,rdiscretisationPostProc%RspatialDiscretisation(1))
          
      ! Discretise the boundary conditions and include them into the Q1
      ! solution vector rvectorPostProc using the corresponding vector filter.
      NULLIFY(p_rdiscreteBCPostProc)
      CALL bcasm_discretiseBC (&
          rdiscretisationPostProc,p_rdiscreteBCPostProc,.FALSE., &
          getBoundaryValues,rcollection)
      CALL lsysbl_createVecFromScalar(rvectorPostProc,rvectorPostProcBlock)
      CALL vecfil_discreteBCsol (rvectorPostProcBlock,1.0_DP,p_rdiscreteBCPostProc)
      CALL bcasm_releaseDiscreteBC (p_rdiscreteBCPostProc)
          
      ! Put the vector to the postprocessing file
      CALL lsyssc_getbase_double (rvectorPostProc,p_DdataQ1)
      CALL ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, &
         p_DdataQ1,p_Ddata)
      
      ! Release the allocated information
      CALL lsysbl_releaseVector (rvectorPostProcBlock)
      CALL lsyssc_releaseVector (rvectorPostProc)
      CALL spdiscr_releaseBlockDiscr (rdiscretisationPostProc)
    END SELECT
    
    ! Write the file to disc, that's it.
    CALL ucd_write (rexport)
    CALL ucd_release (rexport)
    
    ! +------------------------------------------------------------------------
    ! | CLEANUP
    ! +------------------------------------------------------------------------
    !
    ! We are finished - but not completely!
    ! Now, clean up so that all the memory is available again.
    !
    ! Release solver data and structure
    CALL linsol_doneData (p_rsolverNode)
    CALL linsol_doneStructure (p_rsolverNode)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    CALL linsol_releaseSolver (p_rsolverNode)
    
    ! Release the block matrix/vectors
    CALL lsysbl_releaseVector (rtempBlock)
    CALL lsysbl_releaseVector (rvectorBlock)
    CALL lsysbl_releaseVector (rrhsBlock)
    CALL lsysbl_releaseMatrix (rmatrixBlock)

    ! Release the scalar matrix/rhs vector which were used to create
    ! the block matrices/vectors. These must exist as long as the
    ! block matrices/vectors exist, as the block matrices/vectors are
    ! only 'copies' of the scalar ones, sharing the same handles!
    CALL lsyssc_releaseVector (rrhs)
    CALL lsyssc_releaseMatrix (rmatrix)
    
    ! Release our discrete version of the boundary conditions
    CALL bcasm_releaseDiscreteBC (p_rdiscreteBC)

    ! ...and also the corresponding analytic description.
    CALL bcond_doneBC (p_rboundaryConditions)
    
    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    CALL spdiscr_releaseBlockDiscr(rdiscretisation)
    
    ! Release the triangulation. 
    CALL tria_done (rtriangulation)
    
    ! Finally release the domain, that's it.
    CALL boundary_release (p_rboundary)

    ! Release the collection
    CALL collct_done (rcollection)

    ! Release the parameters from the DAT file.
    CALL parlst_done(rparams)
    
  END SUBROUTINE

END MODULE
