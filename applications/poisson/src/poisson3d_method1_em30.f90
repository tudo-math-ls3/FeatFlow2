!##############################################################################
!# ****************************************************************************
!# <name> poisson3d_method1_em30 </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple Poisson
!# problem with constant coefficients on a simple domain.
!# This module is a (provisional) equivalent to the 2D example
!# poisson2d_method7.
!# The element used here is EM30. For VTK output, it is interpolated
!# to Q1 and then written to the output file.
!# </purpose>
!##############################################################################

module poisson3d_method1_em30

  use fsystem
  use genoutput
  use storage
  use linearsolver
  use boundary
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
  use matrixio
  use meshregion
  use spdiscprojection
    
  use poisson3d_callback
  
  implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine poisson3d_1_em30
  
!<description>
  ! This is an all-in-one poisson solver for directly solving a Poisson
  ! problem without making use of special features like collections
  ! and so on. The routine performs the following tasks:
  !
  ! 1.) Read in triangulation
  ! 2.) Set up RHS
  ! 3.) Set up matrix
  ! 4.) Create solver structure
  ! 5.) Solve the problem
  ! 6.) Convert the solution vector into VTK style
  ! 7.) Write solution to VTK file
  ! 8.) Release all variables, finish
!</description>

!</subroutine>

    ! Definitions of variables.
    !
    ! We need a couple of variables for this problem. Let us see...
    !
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

    ! An object for saving the boundary mesh region
    type(t_meshregion) :: rmeshRegion

    ! A variable describing the discrete boundary conditions.
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
    
    ! Output block for UCD output to VTK file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir
    
    ! Declarations for projecting a vector to the Q1 space for VTK export
    type(t_vectorBlock) :: rprjVector
    type(t_blockDiscretisation) :: rprjDiscretisation
    type(t_discreteBC), target :: rdiscreteBC_Q1

    ! Ok, let us start.
    !
    ! We want to solve our Poisson problem on level...
    NLMAX = 4
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Read the domain, read the mesh, refine
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Get the path $PREDIR from the environment, where to read .prm/.tri files
    ! from. If that does not exist, write to the directory "./pre".
    if (.not. sys_getenv_string("PREDIR", spredir)) spredir = "./pre"

    ! At first, read in the basic triangulation.
    call tria_readTriFile3D (rtriangulation, trim(spredir)//"/CUBE.tri")
    
    ! Refine it.
    call tria_quickRefine2LevelOrdering (NLMAX-1,rtriangulation)
    
    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    call tria_initStandardMeshFromRaw (rtriangulation)
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up a discretisation structure which tells the code which
    ! finite element to use
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we only have one block.
    call spdiscr_initBlockDiscr (rdiscretisation,1,rtriangulation)
    
    ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
    ! structures for every component of the solution vector.
    ! Initialise the first element of the list to specify the element
    ! for this solution component:
    call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
                                   EL_EM30_3D,rtriangulation)

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
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
        LSYSSC_MATRIX9,rmatSystem%RmatrixBlock(1,1))
    
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
    rform%BconstantCoeff = .true.
    rform%Dcoefficients(1)  = 1.0
    rform%Dcoefficients(2)  = 1.0
    rform%Dcoefficients(3)  = 1.0

    ! Now we can build the matrix entries.
    ! We specify the callback function coeff_Laplace for the coefficients.
    ! As long as we use constant coefficients, this routine is not used.
    ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
    ! the framework will call the callback routine to get analytical
    ! data.
    call bilf_buildMatrixScalar (&
        rform,.true.,rmatSystem%RmatrixBlock(1,1),rcubatureInfo)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Create RHS and solution vectors
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        
    ! Next step: Create a RHS vector, a solution vector and a temporary
    ! vector. All are filled with zero.
    call lsysbl_createVectorBlock (rdiscretisation,rvecRhs,.true.)
    call lsysbl_createVectorBlock (rdiscretisation,rvecSol,.true.)
    call lsysbl_createVectorBlock (rdiscretisation,rvecTmp,.true.)

    ! Set up a linear form structure for the assembly of the
    ! the right hand side.
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC3D
    
    ! ... and then discretise the RHS to get a discrete version of it.
    ! Again we simply create a scalar vector based on the one and only
    ! discretisation structure.
    ! This scalar vector will later be used as the one and only first
    ! component in a block vector.
    call linf_buildVectorScalar (&
        rlinform,.true.,rvecRhs%RvectorBlock(1),rcubatureInfo,coeff_RHS_3D)

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
    
    ! In contrast to the 2D examples, we currently do not have an analytic
    ! description of the domain`s boundary, therefore we need a discrete
    ! (mesh-dependent) description of the mesh`s boundary. This can be done
    ! using mesh-regions.
    !
    ! Create a t_discreteBC structure where we store all discretised boundary
    ! conditions.
    call bcasm_initDiscreteBC(rdiscreteBC)
    
    ! Create a mesh region describing the mesh`s boundary based on the
    ! nodal-property-array of the current triangulation.
    call mshreg_createFromNodalProp(rmeshRegion, rtriangulation, MSHREG_IDX_ALL)
    
    ! Describe Dirichlet BCs on that mesh region
    call bcasm_newDirichletBConMR(rdiscretisation, 1, rdiscreteBC, rmeshRegion,&
                                  getBoundaryValuesMR_3D)
    
    ! We will not free the mesh region structure here, as we will need it
    ! later for calculating the Q1 solution.

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
    !CALL linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,RfilterChain)
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

    ! -------------------------------------------------------------------------
    ! Solution vector projection for projecting a solution vector of
    ! an arbitrary element to Q1 such that it can be written into a VTK file
    ! -------------------------------------------------------------------------

    ! Now, Our vector block is off the way VTK`s style. We need to convert the vector
    ! to Q1 as that is the format, Paraview understands. So the task is to
    ! create a Q1 solution from rvecSol.
    !
    ! Step 1: Create a discretisation structure for Q1, based on our
    ! previous discretisation structure:
    call spdiscr_duplicateBlockDiscr(rvecSol%p_rblockDiscr,rprjDiscretisation)
    call spdiscr_deriveSimpleDiscrSc (&
                 rvecSol%p_rblockDiscr%RspatialDiscr(1), &
                 EL_Q1_3D, CUB_G3_3D, rprjDiscretisation%RspatialDiscr(1))

    ! Step 2: Setup a new solution vector based on this discretisation,
    ! allocate memory.
    call lsysbl_createVecBlockByDiscr (rprjDiscretisation,rprjVector,.false.)
 
    ! Step 3: Take the original solution vector and convert it according to the
    ! new discretisation:
    call spdp_projectSolution(rvecSol,rprjVector)

    ! Step 4: Discretise the boundary condition according to the Q1
    ! discretisation:
    call bcasm_initDiscreteBC(rdiscreteBC_Q1)
    
    ! Describe Dirichlet BCs on the mesh region which was created before
    call bcasm_newDirichletBConMR(rprjDiscretisation, 1, rdiscreteBC_Q1, rmeshRegion,&
                                  getBoundaryValuesMR_3D)
    
    ! Free the mesh region structure as we will not need it anymore
    call mshreg_done(rmeshRegion)

    ! Implement the boundary conditions
    call vecfil_discreteBCsol (rprjVector,rdiscreteBC_Q1)

    ! Now we have a Q1 solution in rprjVector.
    !
    ! Step 6: Write the VTK file.
    ! That is it, rvecSol now contains our solution in VTK style.
    ! We can now start the postprocessing.
      
    ! Call the UCD library to write out a VTK file for our solution.
    !
    ! Get the path for writing postprocessing files from the environment variable
    ! $UCDDIR. If that does not exist, write to the directory "./gmv".
    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = "./gmv"

    ! Start UCD export to VTK file:
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       trim(sucddir)//"/u3d_1_em30.vtk")
    
    ! Add the solution to the UCD exporter
    call ucd_addVectorByVertex (rexport, "sol", UCD_VAR_STANDARD, &
        rvecSol%RvectorBlock(1))
    
    ! Write the file to disc, that is it.
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! Step 7: Release our discrete version of the boundary conditions of Q1.
    call bcasm_releaseDiscreteBC(rdiscreteBC_Q1)
    
    ! Release the block matrix/vectors rprjVector
    call lsysbl_releaseVector(rprjVector)
    
    ! Release the discretisation structure.
    call spdiscr_releaseBlockDiscr (rprjDiscretisation)
    
    ! -------------------------------------------------------------------------
    ! Projection and VTK export finished.
    ! -------------------------------------------------------------------------

    ! Calculate the error to the reference function.
    call pperr_scalar (PPERR_L2ERROR,derror,rvecSol%RvectorBlock(1),&
                       getReferenceFunction_3D, rcubatureInfo=rcubatureInfo)
    call output_line ("L2-error: " // sys_sdEL(derror,10) )

    call pperr_scalar (PPERR_H1ERROR,derror,rvecSol%RvectorBlock(1),&
                       getReferenceFunction_3D, rcubatureInfo=rcubatureInfo)
    call output_line ("H1-error: " // sys_sdEL(derror,10) )
    
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
    
  end subroutine

end module
