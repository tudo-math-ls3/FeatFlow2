!##############################################################################
!# ****************************************************************************
!# <name> heatcond_method1 </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstation program how to solve a simple heat conduction
!# problem with constant coefficients on a simple domain.
!#
!# The head conduction equation solved in this module is the simplest at all:
!#
!#    d/dt u  - Laplace u  =  f
!#
!# With the RHS f and the boundary conditions being constant in time.
!# The equation is discretised in time by implicit Euler, prividing the
!# following discrete equation:
!#
!#    (1/dt M  +  L) u_{n+1}  =  f  +  1/dt M u_n
!#
!# with M=Mass matrix, L=-Laplace, dt=time step size.
!#
!# The whole solution process is implemented into one routine, using a standard
!# linear solver (not multigrid) for solving the linear subproblems in every
!# timestep. The time domain is [0,T] and discretised by a simple explicit
!# Euler.
!#
!# The module bases on the standard poisson example and shows which
!# modifications are necessary to create a heat equation solver from a poisson
!# solver. Boundary conditions, matrices and vectors are all reassembled in
!# every timestep
!# </purpose>
!##############################################################################

module heatcond_method1

  use fsystem
  use genoutput
  use storage
  use boundary
  use cubature
  use derivatives
  use matrixfilters
  use vectorfilters
  use discretebc
  use bcassembly
  use sortstrategy
  use triangulation
  use element
  use spatialdiscretisation
  use coarsegridcorrection
  use filtersupport
  use linearsystemscalar
  use linearsystemblock
  use scalarpde
  use bilinearformevaluation
  use linearformevaluation
  use linearsolver
  use discretebc
  use ucd
    
  use heatcond_callback
  
  implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine heatcond1
  
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
    type(t_matrixBlock) :: rmatrixBlock
    type(t_vectorBlock) :: rvectorBlock,rrhsBlock,rtempBlock

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
    
    ! Output block for UCD output to VTK file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir
    real(DP), dimension(:), pointer :: p_Ddata
    
    ! Time step size, number of timesteps.
    real(DP) :: dtstep
    integer :: ntimesteps

    ! Time and time step counter
    real(DP) :: dtime
    integer :: itimestep
    
    ! Ok, let us start.
    !
    ! We want to solve our heat equation problem on level...
    NLMAX = 6
    
    ! Initialise time step size and number of timesteps
    dtstep = 0.01_DP
    ntimesteps = 10

    ! We start at time 0.0.
    dtime = 0.0_DP

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
    call lsysbl_createMatBlockByDiscr (rdiscretisation,rmatrixBlock)
    
    ! We create a scalar matrix, based on the discretisation structure
    ! for our one and only solution component.
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
        LSYSSC_MATRIX9,rmatrixBlock%RmatrixBlock(1,1))
                                     
    ! Use the discretisation to create a solution vector.
    ! Fill it with zero -- that is out initial condition!
    call lsysbl_createVectorBlock(rdiscretisation,rvectorBlock,.true.)
    
    ! Create a temporary vector
    call lsysbl_createVectorBlock(rdiscretisation,rtempBlock,.true.)
    
    ! Create an empty RHS vector filled with 0.
    call lsysbl_createVectorBlock(rdiscretisation,rrhsBlock,.true.)

    ! Start the timeloop
    do itimestep=1,ntimesteps
       
      ! Next time step.
      dtime = dtime + dtstep
      
      call output_separator(OU_SEP_MINUS)
      call output_line ("Time step "//trim(sys_siL(itimestep,6))// &
                        "     Time "//trim(sys_sdL(dtime,5)))
      call output_lbrk ()

      ! STEP 1: Form the right hand side:  dtimestep*f + M u_{old}

      ! To assemble the basic RHS f, set up the corresponding linear
      ! form (f,Phi_j):
      rlinform%itermCount = 1
      rlinform%Idescriptors(1) = DER_FUNC
      
      ! Discretise the RHS. We simply create a scalar vector
      ! based on the one and only discretisation structure.
      call linf_buildVectorScalar (&
         rlinform,.true.,rrhsBlock%RvectorBlock(1),rcubatureInfo,coeff_RHS)

      ! And now to the entries of the mass matrix.
      ! For assembling of the entries,
      ! we need a bilinear form, which first has to be set up manually.
      ! We specify the bilinear form (Psi_j, Phi_i) for the
      ! scalar system matrix in 2D.
      rform%itermCount = 1
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_FUNC
      
      rform%ballCoeffConstant = .true.
      rform%BconstantCoeff = .true.
      rform%Dcoefficients(1)  = 1.0
      call bilf_buildMatrixScalar (&
          rform,.true.,rmatrixBlock%RmatrixBlock(1,1), rcubatureInfo)
      
      ! Now form the actual RHS by matrix vector multiplication!
      ! dtimestep*f + M u_{old}
      call lsysbl_blockMatVec(rmatrixBlock,rvectorBlock,rrhsBlock,1.0_DP,dtstep)
      
      ! STEP 2: Assemble the system matrix (M + dtimestep*Laplace)
      !
      ! For this purpose, set up the corresponding bilinear form
      ! (grad Psi_j, grad Phi_i):
      rform%itermCount = 2
      rform%Idescriptors(1,1) = DER_DERIV_X
      rform%Idescriptors(2,1) = DER_DERIV_X
      rform%Idescriptors(1,2) = DER_DERIV_Y
      rform%Idescriptors(2,2) = DER_DERIV_Y

      ! In the standard case, we have constant coefficients.
      rform%ballCoeffConstant = .true.
      rform%BconstantCoeff = .true.
      rform%Dcoefficients(1)  = dtstep
      rform%Dcoefficients(2)  = dtstep

      ! Now we can build the matrix entries.
      ! Note that we set bclear=FALSE in this call, so the Laplace part
      ! is added to the existing (!) mass matrix!
      ! So this results in (M + dtstep*Laplace), as we set
      ! the coefficient rform%Dcoefficients to dtstep above!
      call bilf_buildMatrixScalar (&
          rform,.true.,rmatrixBlock%RmatrixBlock(1,1),rcubatureInfo)
      
      ! STEP 3: Create boundary conditions.
      !
      
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
                                        getBoundaryValues)
                               
      ! Now to the edge 2 of boundary component 1 the domain.
      call boundary_createRegion(rboundary,1,2,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues)
                               
      ! Edge 3 of boundary component 1.
      call boundary_createRegion(rboundary,1,3,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues)
      
      ! Edge 4 of boundary component 1. That is it.
      call boundary_createRegion(rboundary,1,4,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues)

      ! Next step is to implement boundary conditions into the RHS,
      ! solution and matrix. This is done using a vector/matrix filter
      ! for discrete boundary conditions.
      ! The discrete boundary conditions are already attached to the
      ! vectors/matrix. Call the appropriate vector/matrix filter that
      ! modifies the vectors/matrix according to the boundary conditions.
      call vecfil_discreteBCrhs (rrhsBlock,rdiscreteBC)
      call vecfil_discreteBCsol (rvectorBlock,rdiscreteBC)
      call matfil_discreteBC (rmatrixBlock,rdiscreteBC)
      
      ! STEP 6: Solve the system
      !
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
      call linsol_setMatrix(p_RsolverNode,rmatrixBlock)
      
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
      call linsol_solveAdaptively (p_rsolverNode,rvectorBlock,rrhsBlock,rtempBlock)
      
      ! STEP 7: Postprocessing
      !
      ! That is it, rvectorBlock now contains our solution. We can now
      ! start the postprocessing.
      !
      ! Get the path for writing postprocessing files from the environment variable
      ! $UCDDIR. If that does not exist, write to the directory "./gmv".
      if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = "./gmv"

      ! Start UCD export to VTK file:
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                        trim(sucddir)//"/u1.vtk."//trim(sys_si0L(itimestep,5)))
      
      call lsyssc_getbase_double (rvectorBlock%RvectorBlock(1),p_Ddata)
      call ucd_addVariableVertexBased (rexport,"sol",UCD_VAR_STANDARD, p_Ddata)
      
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
      
      ! Release our discrete version of the boundary conditions
      call bcasm_releaseDiscreteBC (rdiscreteBC)

    end do

    ! Release the block matrix/vectors
    call lsysbl_releaseVector (rrhsBlock)
    call lsysbl_releaseVector (rtempBlock)
    call lsysbl_releaseVector (rvectorBlock)
    call lsysbl_releaseMatrix (rmatrixBlock)

    ! Release the cubature info structure.
    call spdiscr_releaseCubStructure(rcubatureInfo)

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    call spdiscr_releaseBlockDiscr(rdiscretisation)
    
    ! Release the triangulation.
    call tria_done (rtriangulation)
    
    ! Finally release the domain, that is it.
    call boundary_release (rboundary)
    
  end subroutine

end module
