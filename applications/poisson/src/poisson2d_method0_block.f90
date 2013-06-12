!##############################################################################
!# ****************************************************************************
!# <name> poisson2d_method0_block </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple Poisson
!# problem with constant coefficients on a simple domain.
!# For the assembly, the block assembly routines are used.
!#
!# WARNING: This test does not use the routines from poisson_callback
!#    for the RHS; only boundary values are taken from there!
!# </purpose>
!##############################################################################

module poisson2d_method0_block

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
  use element
  use ucd
  use pprocerror
  
  use blockmatassemblybase
  use blockmatassembly
  use blockmatassemblystdop
  use feevaluation2
  use collection
    
  use poisson2d_callback
  
  implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine poisson2d_0_block
  
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

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    type(t_filterChain), dimension(1), target :: RfilterChain
    
    ! Number of filters in the filter chain
    integer :: nfilters
    
    ! NLMAX receives the level where we want to solve.
    integer :: NLMAX
    
    ! Error indicator during initialisation of the solver
    integer :: ierror
    
    ! Error of FE function to reference function
    real(DP) :: derror
    
    ! A collection for passing values
    type(t_collection) :: rcollection

    ! Vector evaluation structure for the calculation of errors.
    type(t_fev2Vectors) :: revalVectors
    
    ! Output block for UCD output to VTK file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir

    ! Ok, let us start.
    !
    ! We want to solve our Poisson problem on level...
    NLMAX = 7
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Read the domain, read the mesh, refine
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Get the path $PREDIR from the environment, where to read .prm/.tri files
    ! from. If that does not exist, write to the directory "./pre".
    if (.not. sys_getenv_string("PREDIR", spredir)) spredir = "./pre"

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call boundary_read_prm(rboundary, trim(spredir)//"/QUAD.prm")
        
    ! Now read in the basic triangulation.
    call tria_readTriFile2D (rtriangulation, trim(spredir)//"/QUAD.tri", rboundary)
     
    ! Refine it.
    call tria_quickRefine2LevelOrdering (NLMAX-1,rtriangulation,rboundary)
    
    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    call tria_initStandardMeshFromRaw (rtriangulation,rboundary)
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up a discretisation structure which tells the code which
    ! finite element to use
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we only have one block.
    call spdiscr_initBlockDiscr (rdiscretisation,1,&
                                 rtriangulation, rboundary)
    
    ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
    ! structures for every component of the solution vector.
    ! Initialise the first element of the list to specify the element
    ! for this solution component:
    call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
                                   EL_Q1,rtriangulation, rboundary)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up an cubature info structure to tell the code which cubature
    ! formula to use
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
                 
    ! Create an assembly information structure which tells the code
    ! the cubature formula to use. Standard: Gauss 3x3.
    call spdiscr_createDefCubStructure(&  
        rdiscretisation%RspatialDiscr(1),rcubatureInfo,CUB_GEN_AUTO_G3)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Create a 1x1 block matrix with the operator
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! At first, create a basic 1x1 block matrix based on the discretisation.
    call lsysbl_createMatBlockByDiscr (rdiscretisation,rmatSystem)
    
    ! We create a scalar matrix, based on the discretisation structure
    ! for our one and only solution component.
    call bilf_createMatrixStructure (rmatSystem, 1, 1, LSYSSC_MATRIX9)
        
    ! And now to the entries of the matrix.
    !
    ! Allocate memory for the entries, fill with zero.
    call lsysbl_allocEmptyMatrix (rmatSystem,LSYSSC_SETM_ZERO)
    
    ! Pass the constant dnu=1.0 in front of the Laplace via rcollection.
    rcollection%DquickAccess(1) = 1.0_DP
    
    ! Laplace is to be assembled at position (y,x) = (1,1) in the block matrix
    rcollection%IquickAccess(1) = 1
    rcollection%IquickAccess(2) = 1

    ! The routine "bma_fcalc_laplace" can be found in "blockmatassembly.f90" as
    ! an example for a matrix calculation routine!
    call bma_buildMatrix (rmatSystem,BMA_CALC_STANDARD,&
          bma_fcalc_laplace,rcubatureInfo=rcubatureInfo,rcollection=rcollection)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Create RHS and solution vectors
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        
    ! Next step: Create a RHS vector, a solution vector and a temporary
    ! vector. All are filled with zero.
    call lsysbl_createVectorBlock (rdiscretisation,rvecRhs,.true.)
    call lsysbl_createVectorBlock (rdiscretisation,rvecSol,.true.)
    call lsysbl_createVectorBlock (rdiscretisation,rvecTmp,.true.)

    ! Discretise the RHS to get a discrete version of it.
    !
    ! The routine "bma_fcalc_rhsBubble" can be found in "blockmatassembly.f90" as
    ! an example for a RHS routine!
    call bma_buildVector (rvecRhs,BMA_CALC_STANDARD,&
        bma_fcalc_rhsBubble,rcubatureInfo=rcubatureInfo,rcollection=rcollection)

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

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up a linear solver
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! During the linear solver, the boundary conditions are also
    ! frequently imposed to the vectors. But as the linear solver
    ! does not work with the actual solution vectors but with
    ! defect vectors instead.
    ! So, set up a filter chain that filters the defect vector
    ! during the solution process to implement discrete boundary conditions.
    call filter_clearFilterChain (RfilterChain,nfilters)
    call filter_newFilterDiscBCDef (RfilterChain,nfilters,rdiscreteBC)

    ! Create a BiCGStab-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    nullify(p_rpreconditioner)
    call linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,RfilterChain)
    
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
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Solve the system
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! Finally solve the system. As we want to solve Ax=b with
    ! b being the real RHS and x being the real solution vector,
    ! we use linsol_solveAdaptively. If b is a defect
    ! RHS and x a defect update to be added to a solution vector,
    ! we would have to use linsol_precondDefect instead.
    call linsol_solveAdaptively (p_rsolverNode,rvecSol,rvecRhs,rvecTmp)
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Postprocessing of the solution
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! That is it, rvecSol now contains our solution. We can now
    ! start the postprocessing.
    !
    ! Get the path for writing postprocessing files from the environment variable
    ! $UCDDIR. If that does not exist, write to the directory "./gmv".
    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = "./gmv"

    ! Start UCD export to VTK file:
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       trim(sucddir)//"/u2d_0_block.vtk")
    
    ! Add the solution to the UCD exporter
    call ucd_addVectorByVertex (rexport, "sol", UCD_VAR_STANDARD, &
        rvecSol%RvectorBlock(1))
    
    ! Write the file to disc, that is it.
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Calculate the error to the reference function
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Set up revalVectors with the solution vectors.
    call fev2_addVectorToEvalList(revalVectors,rvecSol%RvectorBlock(1),1)

    ! The routines "bma_fcalc_bubbleL2error" and "bma_fcalc_bubbleH1error"
    ! can be found in the module "blockmatassembly.f90" as examples
    ! for L2/H1 reference functions.
    !
    ! L2-error, squared
    call bma_buildIntegral(derror,BMA_CALC_STANDARD,bma_fcalc_bubbleL2error,&
        rcollection=rcollection,revalVectors=revalVectors,&
        rcubatureInfo=rcubatureInfo)
    
    derror = sqrt(derror)

    call output_line ("L2-error: " // sys_sdEL(derror,10) )

    ! H1-error, squared
    call bma_buildIntegral(derror,BMA_CALC_STANDARD,bma_fcalc_bubbleH1error,&
        rcollection=rcollection,revalVectors=revalVectors,&
        rcubatureInfo=rcubatureInfo)
    
    derror = sqrt(derror)
    
    call output_line ("H1-error: " // sys_sdEL(derror,10) )
    
    ! Cleanup
    call fev2_releaseVectorList(revalVectors)

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
