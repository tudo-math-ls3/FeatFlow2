!##############################################################################
!# ****************************************************************************
!# <name> poisson2d_method1_hadapt </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple Poisson
!# problem with constant coefficients on a simple domain.
!#
!# The triangulation is adaptively refined by the h-adaptivity module,
!# mixing triangle and quadrilateral elements. As solver, the Gauss elimination
!# UMFPACK solver is used.
!# </purpose>
!##############################################################################

module poisson2d_method1_hadapt

  use fsystem
  use genoutput
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use ucd
  use pprocerror
  use hadaptaux
  use hadaptivity
    
  use poisson2d_callback
  
  implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine poisson2d_1_hadapt
  
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
    type(t_boundary) :: rboundary
    
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    type(t_blockDiscretisation) :: rdiscretisation
    
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform
    
    ! A scalar matrix and vector. The vector accepts the RHS of the problem
    ! in scalar form.
    type(t_matrixScalar) :: rmatrix
    type(t_vectorScalar) :: rrhs

    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_matrixBlock) :: rmatrixBlock
    type(t_vectorBlock) :: rvectorBlock,rrhsBlock,rtempBlock

    ! A set of variables describing the analytic and discrete boundary
    ! conditions.    
    type(t_boundaryRegion) :: rboundaryRegion
    type(t_discreteBC), target :: rdiscreteBC

    ! A solver node that accepts parameters for the linear solver    
    type(t_linsolNode), pointer :: p_rsolverNode,p_rpreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(1) :: Rmatrices

    ! Adaptivity structure
    type(t_hadapt) :: rhadapt
    type(t_vectorScalar) :: rindicator

    ! NLMIN defines the pre-refinement level of the coarse mesh.
    integer :: NLMIN
    
    ! NLMAXhRefinement defines how much refinement levels are used by
    ! the adaptive refinement
    integer :: NLMAXhRefinement
    
    ! Maximum number of h-adaptivity cycles
    integer :: MAXhRefinementSteps

    ! Error indicator during initialisation of the solver
    integer :: ierror    
    
    ! Error of FE function to reference function
    real(DP) :: derror
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir
    real(DP), dimension(:), pointer :: p_Ddata

    ! Ok, let's start. 
    !
    ! At first we define...
    ! 1.) the minimum refinement level that is used to get the coarse mesh. 
    NLMIN = 2
    
    ! 2.) the maximum refinement level of the h-adaptive refinement routines.
    ! This is relative to NLMIN: Starting from level NLMIN, the h-adaptivity
    ! routines refine the mesh up to NLMAXhRefinement times.
    NLMAXhRefinement = 6

    ! 3.) the maximum number of refinement sweeps.
    MAXhRefinementSteps = 8
    
    ! +------------------------------------------------------------------------
    ! | BOUNDARY AND TRIANGULATION
    ! +------------------------------------------------------------------------

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    ! Set rboundary to NULL to create a new structure on the heap.
    call boundary_read_prm(rboundary, './pre/QUAD.prm')
        
    ! Now read in the basic triangulation.
    call tria_readTriFile2D (rtriangulation, './pre/QUAD.tri', rboundary)
    
    ! Refine it to get the coarse mesh.
    call tria_quickRefine2LevelOrdering (NLMIN-1,rtriangulation,rboundary)
    
    ! Create information about adjacencies and everything one needs from
    ! a triangulation.
    call tria_initStandardMeshFromRaw (rtriangulation,rboundary)

    ! +------------------------------------------------------------------------
    ! | SETUP OF THE H-ADAPTIVITY
    ! +------------------------------------------------------------------------
       
    rhadapt%iSpec               = 2**1
    rhadapt%nsubdividemax       = NLMAXhRefinement
    rhadapt%iadaptationStrategy = HADAPT_REDGREEN
    
    ! The following tolerance values are the bounds for the monitor function
    ! gethadaptMonitorFunction below. If this monitor function returns a value
    ! > drefinementTolerance to an element, that element is refined.
    ! A value < dcoarseningTolerance on the other hand results in coarsening.
    rhadapt%drefinementTolerance = 0.5
    rhadapt%dcoarseningTolerance = 0.1
    call hadapt_initFromTriangulation(rhadapt,rtriangulation)

    ! +------------------------------------------------------------------------
    ! | SETUP OF THE DISCRETISATION
    ! +------------------------------------------------------------------------

    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we only have one block.
    call spdiscr_initBlockDiscr (rdiscretisation,1,&
                                 rtriangulation, rboundary)

    ! Repeat the procedure until the maximum number of refinement
    ! steps has been reached. This will be checked below.
    do 
      
      ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
      ! structures for every component of the solution vector.
      ! Initialise the first element of the list to specify the element
      ! and cubature rule for this solution component:
      call spdiscr_initDiscr_triquad (rdiscretisation%RspatialDiscr(1), &
          EL_P1,EL_Q1,CUB_G3_T,CUB_G2X2,rtriangulation, &
          rboundary)
      
      ! Now as the discretisation is set up, we can start to generate
      ! the structure of the system matrix which is to solve.
      ! We create a scalar matrix, based on the discretisation structure
      ! for our one and only solution component.
      call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
          LSYSSC_MATRIX9,rmatrix)
      
      ! +------------------------------------------------------------------------
      ! | DISCRETISATION OF MATRICES/VECTORS
      ! +------------------------------------------------------------------------
      
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
      rform%Dcoefficients(1)  = 1.0 
      rform%Dcoefficients(2)  = 1.0 
      
      ! Now we can build the matrix entries.
      ! We specify the callback function coeff_Laplace for the coefficients.
      ! As long as we use constant coefficients, this routine is not used.
      ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
      ! the framework will call the callback routine to get analytical
      ! data.
      call bilf_buildMatrixScalar (rform,.true.,rmatrix,coeff_Laplace_2D)
      
      ! The same has to be done for the right hand side of the problem.
      ! At first set up the corresponding linear form (f,Phi_j):
      rlinform%itermCount = 1
      rlinform%Idescriptors(1) = DER_FUNC
      
      ! ... and then discretise the RHS to get a discrete version of it.
      ! Again we simply create a scalar vector based on the one and only
      ! discretisation structure.
      ! This scalar vector will later be used as the one and only first
      ! component in a block vector.
      call linf_buildVectorScalar (rdiscretisation%RspatialDiscr(1),&
          rlinform,.true.,rrhs,coeff_RHS_2D)
      
      ! The linear solver only works for block matrices/vectors - but above,
      ! we created scalar ones. So the next step is to make a 1x1 block
      ! system from the matrices/vectors above which the linear solver
      ! understands.
      call lsysbl_createMatFromScalar (rmatrix,rmatrixBlock,rdiscretisation)
      call lsysbl_createVecFromScalar (rrhs,rrhsBlock,rdiscretisation)
      
      ! +------------------------------------------------------------------------
      ! | DISCRETISATION AND IMPLEMENTATION OF BOUNDARY CONDITIONS
      ! +------------------------------------------------------------------------
      
      ! Now we have the raw problem. What is missing is the definition of the boundary
      ! conditions.
      ! For implementing boundary conditions, we use a 'filter technique with
      ! discretised boundary conditions'. This means, we first have to calculate
      ! a discrete version of the analytic BC, which we can implement into the
      ! solution/RHS vectors using the corresponding filter.
      !
      ! Create a t_discreteBC structure where we store all discretised boundary
      ! conditions.
      call bcasm_initDiscreteBC(rdiscreteBC)
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
      
      ! We use this boundary region and specify that we want to have Dirichlet
      ! boundary there. The following call does the following:
      ! - Create Dirichlet boundary conditions on the region rboundaryRegion.
      !   We specify icomponent='1' to indicate that we set up the
      !   Dirichlet BC's for the first (here: one and only) component in the 
      !   solution vector.
      ! - Discretise the boundary condition so that the BC's can be applied
      !   to matrices and vectors
      ! - Add the calculated discrete BC's to rdiscreteBC for later use.
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues_2D)
                               
      ! Now to the edge 2 of boundary component 1 the domain. 
      call boundary_createRegion(rboundary,1,2,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues_2D)
                               
      ! Edge 3 of boundary component 1.
      call boundary_createRegion(rboundary,1,3,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues_2D)
      
      ! Edge 4 of boundary component 1. That's it.
      call boundary_createRegion(rboundary,1,4,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues_2D)
      
      ! Hang the pointer into the vector and matrix. That way, these
      ! boundary conditions are always connected to that matrix and that
      ! vector.
      rmatrixBlock%p_rdiscreteBC => rdiscreteBC
      rrhsBlock%p_rdiscreteBC => rdiscreteBC
      
      ! Now we have block vectors for the RHS and the matrix. What we
      ! need additionally is a block vector for the solution and
      ! temporary data. Create them using the RHS as template.
      ! Fill the solution vector with 0:
      call lsysbl_createVecBlockIndirect (rrhsBlock, rvectorBlock, .true.)
      call lsysbl_createVecBlockIndirect (rrhsBlock, rtempBlock, .false.)
      
      ! Next step is to implement boundary conditions into the RHS,
      ! solution and matrix. This is done using a vector/matrix filter
      ! for discrete boundary conditions.
      ! The discrete boundary conditions are already attached to the
      ! vectors/matrix. Call the appropriate vector/matrix filter that
      ! modifies the vectors/matrix according to the boundary conditions.
      call vecfil_discreteBCrhs (rrhsBlock)
      call vecfil_discreteBCsol (rvectorBlock)
      call matfil_discreteBC (rmatrixBlock)
      
      ! +------------------------------------------------------------------------
      ! | INVOKE THE SOLVER
      ! +------------------------------------------------------------------------
      
      nullify(p_rpreconditioner)
      call linsol_initUMFPACK4 (p_rsolverNode)
      
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
      call linsol_setMatrices(p_RsolverNode,Rmatrices)
      
      ! Initialise structure/data of the solver. This allows the
      ! solver to allocate memory / perform some precalculation
      ! to the problem.
      call linsol_initStructure (p_rsolverNode, ierror)
      if (ierror .ne. LINSOL_ERR_NOERROR) stop
      call linsol_initData (p_rsolverNode, ierror)
      if (ierror .ne. LINSOL_ERR_NOERROR) stop
      
      ! Finally solve the system. As we want to solve Ax=b with
      ! b being the real RHS and x being the real solution vector,
      ! we use linsol_solveAdaptively. If b is a defect
      ! RHS and x a defect update to be added to a solution vector,
      ! we would have to use linsol_precondDefect instead.
      call linsol_solveAdaptively (p_rsolverNode,rvectorBlock,rrhsBlock,rtempBlock)
      
      ! Do we have to perform one step of h-adaptivity?
      if (rhadapt%nRefinementSteps .ge. MAXhRefinementSteps) exit
      
      ! +----------------------------------------------------------------------
      ! | COMPUTE INDICATOR FOR H-ADAPTIVITY
      ! +----------------------------------------------------------------------     

      ! Perform a posteriori error estimation
      call lsyssc_createVector(rindicator,rtriangulation%NEL,.true.)
      call gethadaptMonitorFunction_2D(rtriangulation,rvectorBlock%RvectorBlock(1),&
          EL_UNDEFINED,3,rindicator)
      
      ! Output error
      call lsyssc_getbase_double(rindicator,p_Ddata)
      !
      ! Get the path for writing postprocessing files from the environment variable
      ! $UCDDIR. If that does not exist, write to the directory "./gmv".
      if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = './gmv'

      ! Start UCD export to GMV file:
      call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,TRIM(sucddir)//'/error8.'//&
          trim(sys_siL(rhadapt%nRefinementSteps,3))//'.gmv')
      call ucd_addVariableElementBased (rexport,'error',UCD_VAR_STANDARD, p_Ddata)
      call ucd_write (rexport)
      call ucd_release (rexport)

      ! Perform one step h-adaptivity
      call hadapt_refreshAdaptation(rhadapt,rtriangulation)
      call hadapt_performAdaptation(rhadapt,rindicator)
      
      ! Release the indicator vector again
      call lsyssc_releaseVector(rindicator)
      
      ! Generate raw mesh from adaptivity structure
      call hadapt_generateRawMesh(rhadapt,rtriangulation)
      
      ! Create information about adjacencies and everything one needs from
      ! a triangulation.
      call tria_initStandardMeshFromRaw (rtriangulation,rboundary)
      
      ! +----------------------------------------------------------------------
      ! | TEMPORAL CLEANUP
      ! +----------------------------------------------------------------------

      ! Release solver data and structure
      call linsol_doneData (p_rsolverNode)
      call linsol_doneStructure (p_rsolverNode)
      
      ! Release the solver node and all subnodes attached to it (if at all):
      call linsol_releaseSolver (p_rsolverNode)
      
      ! Release the block matrix/vectors
      call lsysbl_releaseVector (rtempBlock)
      call lsysbl_releaseVector (rvectorBlock)
      call lsysbl_releaseVector (rrhsBlock)
      call lsysbl_releaseMatrix (rmatrixBlock)

      ! Release the scalar matrix/rhs vector which were used to create
      ! the block matrices/vectors. These must exist as long as the
      ! block matrices/vectors exist, as the block matrices/vectors are
      ! only 'copies' of the scalar ones, sharing the same handles!
      call lsyssc_releaseVector (rrhs)
      call lsyssc_releaseMatrix (rmatrix)

      ! Release our discrete version of the boundary conditions
      call bcasm_releaseDiscreteBC (rdiscreteBC)
      
    end do
    
    ! +------------------------------------------------------------------------
    ! | POSTPROCESSING
    ! +------------------------------------------------------------------------
    
    ! That's it, rvectorBlock now contains our solution. We can now
    ! start the postprocessing. 
    !
    ! Get the path for writing postprocessing files from the environment variable
    ! $UCDDIR. If that does not exist, write to the directory "./gmv".
    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = './gmv'

    ! Start UCD export to GMV file:
    call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       TRIM(sucddir)//'/u2d_1_hadapt.gmv')
    
    call lsyssc_getbase_double (rvectorBlock%RvectorBlock(1),p_Ddata)
    call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
      
    ! Write the file to disc, that's it.
    call ucd_write (rexport)
    call ucd_release (rexport)
    
    ! Calculate the error to the reference function.
    call pperr_scalar (rvectorBlock%RvectorBlock(1),PPERR_L2ERROR,derror,&
                       getReferenceFunction_2D)
    call output_line ('L2-error: ' // sys_sdEL(derror,10) )

    call pperr_scalar (rvectorBlock%RvectorBlock(1),PPERR_H1ERROR,derror,&
                       getReferenceFunction_2D)
    call output_line ('H1-error: ' // sys_sdEL(derror,10) )
    
    ! +------------------------------------------------------------------------
    ! | CLEANUP
    ! +------------------------------------------------------------------------

    ! We are finished - but not completely!
    ! Now, clean up so that all the memory is available again.
    !
    ! Release adaptivity structure
    call hadapt_releaseAdaptation(rhadapt)

    ! Release solver data and structure
    call linsol_doneData (p_rsolverNode)
    call linsol_doneStructure (p_rsolverNode)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    call linsol_releaseSolver (p_rsolverNode)
    
    ! Release the block matrix/vectors
    call lsysbl_releaseVector (rtempBlock)
    call lsysbl_releaseVector (rvectorBlock)
    call lsysbl_releaseVector (rrhsBlock)
    call lsysbl_releaseMatrix (rmatrixBlock)

    ! Release the scalar matrix/rhs vector which were used to create
    ! the block matrices/vectors. These must exist as long as the
    ! block matrices/vectors exist, as the block matrices/vectors are
    ! only 'copies' of the scalar ones, sharing the same handles!
    call lsyssc_releaseVector (rrhs)
    call lsyssc_releaseMatrix (rmatrix)
    
    ! Release our discrete version of the boundary conditions
    call bcasm_releaseDiscreteBC (rdiscreteBC)

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    call spdiscr_releaseBlockDiscr(rdiscretisation)
    
    ! Release the triangulation. 
    call tria_done (rtriangulation)
    
    ! Finally release the domain, that's it.
    call boundary_release (rboundary)
    
  end subroutine

end module

