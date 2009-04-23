!##############################################################################
!# ****************************************************************************
!# <name> poisson1d_method1_mg </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple 1D Poisson
!# problem with constant coefficients on a simple domain.
!#
!# This module is based on poisson1d_method0_simple, but using a multi-grid 
!# solver.
!# </purpose>
!##############################################################################

module poisson1d_method1_mg

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
  use genoutput
  use matrixio
  use vectorio
  use meshregion
    
  use poisson1d_callback
  
  implicit none
  
!<types>

!<typeblock description="Type block defining all information about one level">

  type t_level
  
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation (structure of the
    ! solution, trial/test functions,...)
    type(t_blockDiscretisation) :: rdiscretisation
    
    ! A system matrix for that specific level. The matrix will receive the 
    ! discrete Laplace operator.
    type(t_matrixBlock) :: rmatrix

    ! A variable describing the discrete boundary conditions.    
    type(t_discreteBC) :: rdiscreteBC
  
  end type
  
!</typeblock>

!</types>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine poisson1d_1_mg
  
!<description>
  ! This is an all-in-one poisson solver for directly solving a Poisson
  ! problem without making use of special features like collections
  ! and so on. The routine performs the following tasks:
  !
  ! 1.) Create triangulation
  ! 2.) Set up RHS
  ! 3.) Set up matrix
  ! 4.) Create solver structure
  ! 5.) Solve the problem
  ! 6.) Write solution to GMV file
  ! 7.) Release all variables, finish
!</description>

!</subroutine>

    ! Definitions of variables.
    !
    ! We need a couple of variables for this problem. Let's see...
    !
    ! An array of problem levels for the multigrid solver
    type(t_level), dimension(:), pointer :: Rlevels
    
    ! An object for saving the boundary mesh region
    type(t_meshregion) :: rmeshRegion
    
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform
    
    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_vectorBlock) :: rvectorBlock,rrhsBlock,rtempBlock

    ! A solver node that accepts parameters for the linear solver    
    type(t_linsolNode), pointer :: p_rsolverNode,p_rsmoother

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(:), pointer :: Rmatrices

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    type(t_filterChain), dimension(1), target :: RfilterChain
    type(t_filterChain), dimension(:), pointer :: p_RfilterChain
    
    ! One level of multigrid
    type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo

    ! Error indicator during initialisation of the solver
    integer :: ierror    
    
    ! Error of FE function to reference function
    real(DP) :: derror

    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir
    real(DP), dimension(:), pointer :: p_Ddata

    ! A temporary variable for the Level-loops
    integer :: i
    
    ! The number of sub-intervals for the discretisation
    ! This defines the coarse mesh on level 1!
    integer :: nintervals = 16
    
    ! The number of levels. The coarse mesh on level 1 with nintervals intervals is
    ! refined nlevel-1 times to get the fine mesh
    integer :: nlevels = 4
    
    ! Ok, let's start. 
    !
    ! Allocate memory for all levels
    allocate(Rlevels(nlevels))
    
    ! At first, create the basic coarse grid triangulation.
    ! Our domain is [0, 1], divided into nintervals sub-intervals.
    call tria_createRawTria1D(Rlevels(1)%rtriangulation, 0.0_DP, 1.0_DP, nintervals)
    
    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    call tria_initStandardMeshFromRaw (Rlevels(1)%rtriangulation)

    ! Now refine the grid for the fine levels.
    do i = 2, nlevels
    
      ! Refine the grid using the 2-Level-Ordering algorithm
      call tria_refine2LevelOrdering(&
          Rlevels(i-1)%rtriangulation, Rlevels(i)%rtriangulation)
      
      ! Create a standard mesh
      call tria_initStandardMeshFromRaw(Rlevels(i)%rtriangulation)
    
    end do
    
    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we only have one block.
    ! Do this for all levels
    do i = 1, nlevels
      call spdiscr_initBlockDiscr (&
          Rlevels(i)%rdiscretisation, 1, Rlevels(i)%rtriangulation)
    end do
    
    ! In the next step, we will define the element type and the cubature
    ! formula that is to be used. For this 1D poisson-example we currently
    ! have 2 possible element types: linear and quadratic ones.
    ! For linear elements the trapezoidal formula satisfies all our
    ! needs, for quadratic elements we should choose a 3-point Gauss-formula.
    !
    ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
    ! structures for every component of the solution vector.
    ! Initialise the first element of the list to specify the element
    ! and cubature rule for this solution component:
    do i = 1, nlevels
      call spdiscr_initDiscr_simple (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1), &
      ! Setting up a linear element and trapezoidal rule would be...
          EL_P1_1D,CUB_TRZ_1D,Rlevels(i)%rtriangulation)
      ! Setting up a quadratic element and 3-point Gauss rule would be...
          !EL_P2_1D,CUB_G3_1D,Rlevels(i)%rtriangulation)
    end do

    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! We create a scalar matrix, based on the discretisation structure
    ! for our one and only solution component.
    do i = 1, nlevels
      
      ! Initialise the block matrix with default values based on
      ! the discretisation.
      call lsysbl_createMatBlockByDiscr (&
          Rlevels(i)%rdiscretisation,Rlevels(i)%rmatrix)    

      ! Now as the discretisation is set up, we can start to generate
      ! the structure of the system matrix which is to solve.
      ! We create that directly in the block (1,1) of the block matrix
      ! using the discretisation structure of the first block.
      call bilf_createMatrixStructure ( &
           Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
           LSYSSC_MATRIX9,Rlevels(i)%rmatrix%RmatrixBlock(1,1))

      ! Update the structural information of the block matrix, as we manually
      ! changed one of the submatrices:
      call lsysbl_updateMatStrucInfo (Rlevels(i)%rmatrix)
      
      ! And now to the entries of the matrix. For assembling of the entries,
      ! we need a bilinear form, which first has to be set up manually.
      ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
      ! scalar system matrix in 1D.
      rform%itermCount = 1
      rform%Idescriptors(1,1) = DER_DERIV_X
      rform%Idescriptors(2,1) = DER_DERIV_X

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
      call bilf_buildMatrixScalar (rform,.true.,&
           Rlevels(i)%rmatrix%RmatrixBlock(1,1),coeff_Laplace_1D)

    end do
    ! Keep in mind that now p_rdisc and p_rmatrix are pointers to the
    ! discretisation and matrix structures on the finest level!
    
    ! Although we could manually create the solution/RHS vector,
    ! the easiest way to set up the vector structure is
    ! to create it by using our matrix as template:
    call lsysbl_createVecBlockIndMat (Rlevels(nlevels)%rmatrix,rrhsBlock, .false.)

    ! The vector structure is ready but the entries are missing. 
    ! So the next thing is to calculate the content of that vector.
    !
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    
    ! ... and then discretise the RHS to get a discrete version of it.
    call linf_buildVectorScalar (&
        Rlevels(nlevels)%rdiscretisation%RspatialDiscr(1),&
        rlinform,.true.,rrhsBlock%RvectorBlock(1),coeff_RHS_1D)
    
    ! Now we have the raw problem. What is missing is the definition of the boundary
    ! conditions.
    ! For implementing boundary conditions, we use a 'filter technique with
    ! discretised boundary conditions'. This means, we first have to calculate
    ! a discrete version of the analytic BC, which we can implement into the
    ! solution/RHS vectors using the corresponding filter.
    !
    ! In contrast to the 2D poisson examples, we will directly set the
    ! dirichlet boundary conditions by hand instead of discretising an analytic
    ! boundary condition function using a boundary structure.
    !
    do i = 1, nlevels
    
      ! Initialise the discrete BC structure
      call bcasm_initDiscreteBC(Rlevels(i)%rdiscreteBC)

      ! Create a mesh region describing the mesh's boundary based on the
      ! nodal-property-array of the current triangulation.
      call mshreg_createFromNodalProp(rmeshRegion, Rlevels(i)%rtriangulation, &
                                      MSHREG_IDX_ALL)
      
      ! Describe Dirichlet BCs on that mesh region
      call bcasm_newDirichletBConMR(Rlevels(i)%rdiscretisation, 1,&
        Rlevels(i)%rdiscreteBC,rmeshRegion,getBoundaryValuesMR_1D)
      
      ! Free the mesh region structure as we won't need it anymore
      call mshreg_done(rmeshRegion)

      ! Hang the pointer into the matrix. That way, these
      ! boundary conditions are always connected to that matrix.
      Rlevels(i)%rmatrix%p_rdiscreteBC => Rlevels(i)%rdiscreteBC
      
      ! Also implement the boundary conditions into the matrix.
      call matfil_discreteBC (Rlevels(i)%rmatrix)
  
    end do
    
    ! Our one and only right-hand-side also needs to know the boundary conditions.
    rrhsBlock%p_rdiscreteBC => Rlevels(nlevels)%rdiscreteBC
                             
    ! Now we have block vectors for the RHS and the matrix. What we
    ! need additionally is a block vector for the solution and
    ! temporary data. Create them using the RHS as template.
    ! Fill the solution vector with 0:
    call lsysbl_createVecBlockIndirect (rrhsBlock, rvectorBlock, .true.)
    call lsysbl_createVecBlockIndirect (rrhsBlock, rtempBlock, .false.)
    
    ! Next step is to implement boundary conditions into the RHS and
    ! solution vectors. This is done using a vector filter
    ! for discrete boundary conditions.
    ! The discrete boundary conditions are already attached to the
    ! vectors. Call the appropriate vector filter that
    ! modifies the vectors according to the boundary conditions.
    call vecfil_discreteBCrhs (rrhsBlock)
    call vecfil_discreteBCsol (rvectorBlock)

    ! During the linear solver, the boundary conditions are also
    ! frequently imposed to the vectors. But as the linear solver
    ! does not work with the actual solution vectors but with
    ! defect vectors instead.
    ! So, set up a filter chain that filters the defect vector
    ! during the solution process to implement discrete boundary conditions.
    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

    ! Now we have to build up the level information for multigrid.
    !
    ! Create a Multigrid-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    p_RfilterChain => RfilterChain
    call linsol_initMultigrid2 (p_rsolverNode,nlevels,p_RfilterChain)
    
    ! Set up a coarse grid solver.
    ! The coarse grid in multigrid is always grid 1!
    call linsol_getMultigrid2Level (p_rsolverNode,1,p_rlevelInfo)
    call linsol_initUMFPACK4 (p_rlevelInfo%p_rcoarseGridSolver)
    
    ! Now set up the other levels...
    do i = 2, nlevels
    
      ! Create a Jacobi smoother
      call linsol_initJacobi(p_rsmoother)
      
      ! We will use 4 smoothing steps with damping parameter 0.7
      call linsol_convertToSmoother(p_rsmoother, 4, 0.7_DP)
      
      ! And add this multi-grid level. We will use the same smoother
      ! for pre- and post-smoothing.
      call linsol_getMultigrid2Level (p_rsolverNode,i,p_rlevelInfo)
      p_rlevelInfo%p_rpresmoother => p_rsmoother
      p_rlevelInfo%p_rpostsmoother => p_rsmoother
      
    end do

    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2
    
    ! Attach the system matrices to the solver.
    !
    ! We copy our matrices to a big matrix array and transfer that
    ! to the setMatrices routines. This intitialises then the matrices
    ! on all levels according to that array. Note that this does not
    ! allocate new memory, we create only 'links' to existing matrices
    ! into Rmatrices(:)!
    allocate(Rmatrices(nlevels))
    do i=1,nlevels
      call lsysbl_duplicateMatrix (Rlevels(i)%rmatrix,&
          Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    end do
    
    call linsol_setMatrices(p_RsolverNode,Rmatrices(1:nlevels))

    ! We can release Rmatrices immediately -- as long as we don't
    ! release rproblem%RlevelInfo(i)%rmatrix!
    do i=1,nlevels
      call lsysbl_releaseMatrix (Rmatrices(i))
    end do
    deallocate(Rmatrices)
    
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
    
    ! That's it, rvectorBlock now contains our solution. We can now
    ! start the postprocessing. 
    !
    ! Get the path for writing postprocessing files from the environment variable
    ! $UCDDIR. If that does not exist, write to the directory "./gmv".
    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = './gmv'

    ! Start UCD export to GMV file:
    call ucd_startGMV (rexport,UCD_FLAG_STANDARD,&
        Rlevels(nlevels)%rtriangulation,TRIM(sucddir)//'/u1d_1_mg.gmv')
    
    call lsyssc_getbase_double (rvectorBlock%RvectorBlock(1),p_Ddata)
    call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that's it.
    call ucd_write (rexport)
    call ucd_release (rexport)
    
    ! Calculate the error to the reference function.
    call pperr_scalar (rvectorBlock%RvectorBlock(1),PPERR_L2ERROR,derror,&
                       getReferenceFunction_1D)
    call output_line ('L2-error: ' // sys_sdEL(derror,10) )
    call pperr_scalar (rvectorBlock%RvectorBlock(1),PPERR_H1ERROR,derror,&
                       getReferenceFunction_1D)
    call output_line ('H1-error: ' // sys_sdEL(derror,10) )

    ! We are finished - but not completely!
    ! Now, clean up so that all the memory is available again.
    !
    ! Release solver data and structure
    call linsol_doneData (p_rsolverNode)
    call linsol_doneStructure (p_rsolverNode)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    call linsol_releaseSolver (p_rsolverNode)
    
    ! Release the block matrices/vectors
    call lsysbl_releaseVector (rtempBlock)
    call lsysbl_releaseVector (rvectorBlock)
    call lsysbl_releaseVector (rrhsBlock)
    do i = nlevels, 1, -1
      call lsysbl_releaseMatrix (Rlevels(i)%rmatrix)
    end do

    ! Release our discrete version of the boundary conditions
    do i = nlevels, 1, -1
      call bcasm_releaseDiscreteBC (Rlevels(i)%rdiscreteBC)
    end do

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    do i = nlevels, 1, -1
      call spdiscr_releaseBlockDiscr(Rlevels(i)%rdiscretisation)
    end do
    
    ! Release the triangulation. 
    do i = nlevels, 1, -1
      call tria_done (Rlevels(i)%rtriangulation)
    end do

    ! Release level information
    deallocate(Rlevels)
    
    ! That's it!
    
  end subroutine

end module
