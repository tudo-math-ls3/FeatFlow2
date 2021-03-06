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
  use derivatives
  use element
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use filtersupport
  use linearsystemscalar
  use linearsystemblock
  use matrixfilters
  use vectorfilters
  use discretebc
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use scalarpde
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
  ! 7.) Write solution to VTK file
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

    ! Path to the mesh
    character(len=SYS_STRLEN) :: spredir

    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    type(t_blockDiscretisation) :: rdiscretisation
    
    ! Cubature info structure which encapsules the cubature formula
    type(t_scalarCubatureInfo) :: rcubatureInfo
    
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform
    
    ! A matrix, a RHS vector, a solution vector and a temporary vector. 
    ! The RHS vector accepts the RHS of the problem, the solution vector
    ! accepts the solution. All are block vectors with only one block.
    type(t_matrixBlock) :: rmatSystem
    type(t_vectorBlock) :: rvecSol,rvecRhs,rvecTmp
    
    ! A set of variables describing the discrete boundary conditions.
    type(t_boundaryRegion) :: rboundaryRegion
    type(t_discreteBC), target :: rdiscreteBC

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode,p_rpreconditioner

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
    
    ! Output block for UCD output to VTK file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir
    real(DP), dimension(:), pointer :: p_Ddata

    ! Ok, let us start.
    !
    ! At first we define...
    ! 1.) the minimum refinement level that is used to get the coarse mesh.
    NLMIN = 2
    
    ! 2.) the maximum refinement level of the h-adaptive refinement routines.
    ! This is relative to NLMIN: Starting from level NLMIN, the h-adaptivity
    ! routines refine the mesh up to NLMAXhRefinement times.
    NLMAXhRefinement = 8

    ! 3.) the maximum number of refinement sweeps.
    MAXhRefinementSteps = 10
    
    ! +------------------------------------------------------------------------
    ! | BOUNDARY AND TRIANGULATION
    ! +------------------------------------------------------------------------

    ! Get the path $PREDIR from the environment, where to read .prm/.tri files
    ! from. If that does not exist, write to the directory "./pre".
    if (.not. sys_getenv_string("PREDIR", spredir)) spredir = "./pre"

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    ! Set rboundary to NULL to create a new structure on the heap.
    call boundary_read_prm(rboundary, trim(spredir)//"/QUAD.prm")
        
    ! Now read in the basic triangulation.
    call tria_readTriFile2D (rtriangulation, trim(spredir)//"/QUAD.tri", rboundary)
    
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
    rhadapt%drefinementTolerance = 0.2
    rhadapt%dcoarseningTolerance = 0.1

    ! Specify that parameters have been set correctly and perform initialisation
    rhadapt%iSpec = ior(rhadapt%iSpec, HADAPT_HAS_PARAMETERS)
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
      ! for this solution component:
      call spdiscr_initDiscr_triquad (rdiscretisation%RspatialDiscr(1), &
          EL_P1,EL_Q1,rtriangulation, rboundary)

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Set up an cubature info structure to tell the code which cubature
      ! formula to use
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
                   
      ! Create an assembly information structure which tells the code
      ! the cubature formula to use. Standard: Gauss 3x3.
      call spdiscr_createDefCubStructure(&  
          rdiscretisation%RspatialDiscr(1),rcubatureInfo,CUB_GEN_AUTO_G2)

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Create a 1x1 block matrix with the operator
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      
      ! Now as the discretisation is set up, we can start to generate
      ! the structure of the system matrix which is to solve.
      ! At first, create a basic 1x1 block matrix based on the discretisation.
      call lsysbl_createMatrix (rdiscretisation,rmatSystem)
      
      ! We create a scalar matrix, based on the discretisation structure
      ! for our one and only solution component.
      call bilf_createMatrixStructure (rmatSystem, 1, 1, LSYSSC_MATRIX9)
      
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
      rform%Dcoefficients(1)  = 1.0
      rform%Dcoefficients(2)  = 1.0
      
      ! Now we can build the matrix entries.
      ! We specify the callback function coeff_Laplace for the coefficients.
      ! As long as we use constant coefficients, this routine is not used.
      ! By specifying ballCoeffConstant = .FALSE. above,
      ! the framework will call the callback routine to get analytical
      ! data.
      call bilf_buildMatrixScalar (&
          rform,.true.,rmatSystem%RmatrixBlock(1,1),rcubatureInfo)

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Create RHS and solution vectors
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        
      ! Next step: Create a RHS vector, a solution vector and a temporary
      ! vector. All are filled with zero.
      call lsysbl_createVector (rdiscretisation,rvecRhs,.true.)
      call lsysbl_createVector (rdiscretisation,rvecSol,.true.)
      call lsysbl_createVector (rdiscretisation,rvecTmp,.true.)

      ! Set up a linear form structure for the assembly of the
      ! the right hand side.
      ! At first set up the corresponding linear form (f,Phi_j):
      rlinform%itermCount = 1
      rlinform%Idescriptors(1) = DER_FUNC
      
      ! ... and then discretise the RHS to get a discrete version of it.
      ! Again we simply create a scalar vector based on the one and only
      ! discretisation structure.
      ! This scalar vector will later be used as the one and only first
      ! component in a block vector.
      call linf_buildVectorScalar (&
          rlinform,.true.,rvecRhs%RvectorBlock(1),rcubatureInfo,coeff_RHS_2D)

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Assembly of matrices/vectors finished
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Discretise the boundary conditions and apply them to the matrix/RHS/sol.
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      
      ! Now we have the raw problem. What is missing is the definition of the boundary
      ! conditions.
      ! For implementing boundary conditions, we use a `filter technique with
      ! discretised boundary conditions`. This means, we first have to calculate
      ! a discrete version of the analytic BC, which we can implement into the
      ! solution/RHS vectors using the corresponding filter.
      !
      ! Create a t_discreteBC structure where we store all discretised boundary
      ! conditions.
      call bcasm_initDiscreteBC(rdiscreteBC)
      !
      ! We "know" already (from the problem definition) that we have four boundary
      ! segments in the domain. Each of these, we want to use for enforcing
      ! some kind of boundary condition.
      !
      ! We ask the boundary routines to create a "boundary region" - which is
      ! simply a part of the boundary corresponding to a boundary segment.
      ! A boundary region roughly contains the type, the min/max parameter value
      ! and whether the endpoints are inside the region or not.
      call boundary_createRegion(rboundary,1,1,rboundaryRegion)
      
      ! We use this boundary region and specify that we want to have Dirichlet
      ! boundary there. The following call does the following:
      ! - Create Dirichlet boundary conditions on the region rboundaryRegion.
      !   We specify icomponent="1" to indicate that we set up the
      !   Dirichlet BC`s for the first (here: one and only) component in the
      !   solution vector.
      ! - Discretise the boundary condition so that the BC`s can be applied
      !   to matrices and vectors
      ! - Add the calculated discrete BC`s to rdiscreteBC for later use.
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
      
      ! Edge 4 of boundary component 1. That is it.
      call boundary_createRegion(rboundary,1,4,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues_2D)

      ! Next step is to implement boundary conditions into the RHS,
      ! solution and matrix. This is done using a vector/matrix filter
      ! for discrete boundary conditions.
      call vecfil_discreteBCrhs (rvecRhs,rdiscreteBC)
      call vecfil_discreteBCsol (rvecSol,rdiscreteBC)
      call matfil_discreteBC (rmatSystem,rdiscreteBC)
      
      ! +------------------------------------------------------------------------
      ! | INVOKE THE SOLVER
      ! +------------------------------------------------------------------------
      
      nullify(p_rpreconditioner)
      call linsol_initUMFPACK4 (p_rsolverNode)
      
      ! Set the output level of the solver to 2 for some output
      p_rsolverNode%ioutputLevel = 2
      
      ! Attach the system matrix to the solver.
      call linsol_setMatrix(p_RsolverNode,rmatSystem)
      
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

      ! Finally solve the system. As we want to solve Ax=b with
      ! b being the real RHS and x being the real solution vector,
      ! we use linsol_solveAdaptively. If b is a defect
      ! RHS and x a defect update to be added to a solution vector,
      ! we would have to use linsol_precondDefect instead.
      call linsol_solveAdaptively (p_rsolverNode,rvecSol,rvecRhs,rvecTmp)
      
      ! Do we have to perform one step of h-adaptivity?
      if (rhadapt%nRefinementSteps .ge. MAXhRefinementSteps) exit
      
      ! +----------------------------------------------------------------------
      ! | COMPUTE INDICATOR FOR H-ADAPTIVITY
      ! +----------------------------------------------------------------------

      ! Perform a posteriori error estimation
      call lsyssc_createVector(rindicator,rtriangulation%NEL,.true.)
      call gethadaptMonitorFunction_2D(rtriangulation,rvecSol%RvectorBlock(1),&
          EL_UNDEFINED,3,rindicator)
      
      ! Output error
      call lsyssc_getbase_double(rindicator,p_Ddata)
      !
      ! Get the path for writing postprocessing files from the environment variable
      ! $UCDDIR. If that does not exist, write to the directory "./gmv".
      if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = "./gmv"

      ! Start UCD export to VTK file:
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,trim(sucddir)//"/error8."//&
          trim(sys_siL(rhadapt%nRefinementSteps,3))//".vtk")
      call ucd_addVariableElementBased (rexport,"error",UCD_VAR_STANDARD, p_Ddata)
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
      call lsysbl_releaseVector (rvecTmp)
      call lsysbl_releaseVector (rvecSol)
      call lsysbl_releaseVector (rvecRhs)
      call lsysbl_releaseMatrix (rmatSystem)

      ! Release the cubature info structure.
      call spdiscr_releaseCubStructure(rcubatureInfo)

      ! Release our discrete version of the boundary conditions
      call bcasm_releaseDiscreteBC (rdiscreteBC)
      
    end do
    
    ! +------------------------------------------------------------------------
    ! | POSTPROCESSING
    ! +------------------------------------------------------------------------
    
    ! That is it, rvecSol now contains our solution. We can now
    ! start the postprocessing.
    !
    ! Get the path for writing postprocessing files from the environment variable
    ! $UCDDIR. If that does not exist, write to the directory "./gmv".
    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = "./gmv"

    ! Start UCD export to VTK file:
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       trim(sucddir)//"/u2d_1_hadapt.vtk")
    
    ! Add the solution to the UCD exporter
    call ucd_addVectorByVertex (rexport, "sol", UCD_VAR_STANDARD, &
        rvecSol%RvectorBlock(1))
      
    ! Write the file to disc, that is it.
    call ucd_write (rexport)
    call ucd_release (rexport)
    
    ! Calculate the error to the reference function.
    call pperr_scalar (PPERR_L2ERROR,derror,rvecSol%RvectorBlock(1),&
                       getReferenceFunction_2D, rcubatureInfo=rcubatureInfo)
    call output_line ("L2-error: " // sys_sdEL(derror,10) )

    call pperr_scalar (PPERR_H1ERROR,derror,rvecSol%RvectorBlock(1),&
                       getReferenceFunction_2D, rcubatureInfo=rcubatureInfo)
    call output_line ("H1-error: " // sys_sdEL(derror,10) )
    
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
    call lsysbl_releaseVector (rvecTmp)
    call lsysbl_releaseVector (rvecSol)
    call lsysbl_releaseVector (rvecRhs)
    call lsysbl_releaseMatrix (rmatSystem)

    ! Release the cubature info structure.
    call spdiscr_releaseCubStructure(rcubatureInfo)

    ! Release our discrete version of the boundary conditions
    call bcasm_releaseDiscreteBC (rdiscreteBC)

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    call spdiscr_releaseBlockDiscr(rdiscretisation)
    
    ! Release the triangulation.
    call tria_done (rtriangulation)
    
    ! Finally release the domain, that is it.
    call boundary_release (rboundary)
    
  end subroutine

end module

