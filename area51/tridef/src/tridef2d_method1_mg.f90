!##############################################################################
!# ****************************************************************************
!# <name> poisson2d_method1_mg </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple Poisson
!# problem with constant coefficients on a simple domain.
!#
!# This module is based on poisson2d_method0_simple, but using a multi-grid
!# solver.
!# </purpose>
!##############################################################################

module tridef2d_method1_mg

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
  use ucd
  use pprocerror
  use genoutput
  use collection
  use griddeform
  use gridsmooth
  use tridef2d_callback
  use statistics
  
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

  subroutine poisson2d_1_mg
  
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
    ! An array of problem levels for the multigrid solver
    type(t_level), dimension(:), pointer :: Rlevels

    ! An object for saving the domain:
    type(t_boundary) :: rboundary
    
    ! Path to the mesh
    character(len=SYS_STRLEN) :: spredir

    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform
    
    ! A couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_vectorBlock) :: rvectorBlock,rrhsBlock,rtempBlock

    ! A variable that is used to specify a region on the boundary.
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
    type(t_filterChain), dimension(:), pointer :: p_RfilterChain
    
    ! One level of multigrid
    type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo
    type(t_triangulation),target :: rtriangulation2,rtriangulation3
    ! NLMIN receives the level of the coarse grid.
    integer :: NLMIN

    ! NLMAX receives the level where we want to solve.
    integer :: NLMAX
    
    ! Error indicator during initialisation of the solver
    integer :: ierror
    
    ! Error of FE function to reference function
    real(DP) :: derror
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir
    real(DP), dimension(:), pointer :: p_Ddata
    
    ! Some temporary variables
    integer :: i
    
    type(t_griddefWork) :: rgriddefWork
    
    ! dummy
    integer :: idummy,iloop,idupFlag,iiter

    ! grid deform structure
    type(t_griddefInfo) :: rgriddefInfo

    type(t_timer) :: rtimerGridDeformation

    ! Ok, let's start.
    !
    ! We want to solve our Poisson problem on level...
    NLMIN = 3
    NLMAX = 7
    
    ! Allocate memory for all levels
    allocate(Rlevels(NLMIN:NLMAX))
    
    ! Get the path $PREDIR from the environment, where to read .prm/.tri files
    ! from. If that does not exist, write to the directory "./pre".
    if (.not. sys_getenv_string("PREDIR", spredir)) spredir = './pre'

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call boundary_read_prm(rboundary, trim(spredir)//'/QUAD.prm')
        
    ! Now read in the basic triangulation into our coarse level.
    call tria_readTriFile2D (Rlevels(NLMIN)%rtriangulation, &
                             trim(spredir)//'/QUAD.tri', rboundary)
    
    ! Refine it.
    call tria_quickRefine2LevelOrdering (NLMIN-1,&
        Rlevels(NLMIN)%rtriangulation,rboundary)
    
    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    call tria_initStandardMeshFromRaw (Rlevels(NLMIN)%rtriangulation,&
        rboundary)
    
    ! Now refine the grid for the fine levels.
    do i = NLMIN, NLMAX

      if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = './gmv'
      ! Start UCD export to GMV file:
      ! write out the mesh before refinement
      call ucd_startGMV (rexport,UCD_FLAG_STANDARD,&
          Rlevels(i)%rtriangulation,TRIM(sucddir)//'/u2d_gridef_pre'//&
          trim(sys_siL(i,3))//'.gmv')


      ! griddeformation setup
      call griddef_deformationInit(rgriddefInfo,Rlevels(i)%rtriangulation,&
                                   NLMIN,&
                                   i,rboundary,2)

      do iloop=NLMIN,i
        call griddef_buildHGrid(rgriddefInfo,Rlevels(iloop)%rtriangulation,iloop)
      end do
      
      ! pre-smooth the mesh
!      do iiter=1,3
!        call gsmth_umbrella(rgriddefInfo%p_rhLevels(i)%rtriangulation)
!      end do
      
      
      ! Write the file to disc, that's it.
      call ucd_write (rexport)
      call ucd_release (rexport)
      

      ! Deform
      call griddef_performDeformation(rgriddefInfo, rgriddefWork,idummy,&
                                      .TRUE., .FALSE., .FALSE., &
                                      .FALSE., i, 0, 0,&
                                      tridef2d_monitorfct)
                                      
      ! Start UCD export to GMV file:
      ! write out the mesh after refinement
      call ucd_startGMV (rexport,UCD_FLAG_STANDARD,&
          Rlevels(i)%rtriangulation,TRIM(sucddir)//'/u2d_gridef_post'//&
          trim(sys_siL(i,3))//'.gmv')
      
      ! Write the file to disc, that's it.
      call ucd_write (rexport)
      call ucd_release (rexport)

      
      ! Release Deformation structures
      call griddef_DeformationDone(rgriddefInfo,rgriddefWork)

      if(i .lt. NLMAX)then
        ! Refine the grid that was just deformed to get the next level mesh
        call tria_refine2LevelOrdering(Rlevels(i)%rtriangulation,&
            Rlevels(i+1)%rtriangulation,rboundary)
        
        ! Create a standard mesh out of the deformed mesh and save it
        call tria_initStandardMeshFromRaw(Rlevels(i+1)%rtriangulation,&
          rboundary)

      end if


      ! postsmoothing
!      do iiter=1,1
!        call gsmth_umbrella(Rlevels(i)%rtriangulation)
!      end do
    
    end do ! end for levels
    
    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    

    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we only have one block.
    ! Do this for all levels
    do i = NLMIN, NLMAX
      call spdiscr_initBlockDiscr (Rlevels(i)%rdiscretisation, 1, &
                                   Rlevels(i)%rtriangulation, rboundary)
    end do
    
    ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
    ! structures for every component of the solution vector.
    ! Initialise the first element of the list to specify the element
    ! and cubature rule for this solution component:
    do i = NLMIN, NLMAX
      call spdiscr_initDiscr_simple (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1), &
          EL_E011,CUB_G2X2,Rlevels(i)%rtriangulation, rboundary)
    end do
                 
    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! We create a scalar matrix, based on the discretisation structure
    ! for our one and only solution component.
    do i = NLMIN, NLMAX

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
      call bilf_buildMatrixScalar (rform,.true.,&
           Rlevels(i)%rmatrix%RmatrixBlock(1,1),coeff_Laplace_2D)
    
    end do
      
    ! Although we could manually create the solution/RHS vector,
    ! the easiest way to set up the vector structure is
    ! to create it by using our matrix as template:
    call lsysbl_createVecBlockIndMat (Rlevels(NLMAX)%rmatrix,rrhsBlock, .false.)

    ! The vector structure is ready but the entries are missing.
    ! So the next thing is to calculate the content of that vector.
    !
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    
    ! ... and then discretise the RHS to get a discrete version of it.
    ! Again we simply create a scalar vector based on the one and only
    ! discretisation structure.
    ! This scalar vector will later be used as the one and only first
    ! component in a block vector.
    call linf_buildVectorScalar (&
        Rlevels(NLMAX)%rdiscretisation%RspatialDiscr(1),&
        rlinform,.true.,rrhsBlock%RvectorBlock(1),coeff_RHS_2D)
    
    do i = NLMIN, NLMAX
    
      ! Initialise the discrete BC structure
      call bcasm_initDiscreteBC(Rlevels(i)%rdiscreteBC)

      ! On edge 1 of boundary component 1 add Dirichlet boundary conditions.
      call boundary_createRegion(rboundary,1,1,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                        rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                        getBoundaryValues_2D)
                               
      ! Now to the edge 2 of boundary component 1 the domain.
      call boundary_createRegion(rboundary,1,2,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                        rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                        getBoundaryValues_2D)
                               
      ! Edge 3 of boundary component 1.
      call boundary_createRegion(rboundary,1,3,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                        rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                        getBoundaryValues_2D)
      
      ! Edge 4 of boundary component 1. That's it.
      call boundary_createRegion(rboundary,1,4,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                        rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                        getBoundaryValues_2D)
      
      ! Hang the pointer into the matrix. That way, these
      ! boundary conditions are always connected to that matrix.
      Rlevels(i)%rmatrix%p_rdiscreteBC => Rlevels(i)%rdiscreteBC
  
      ! Also implement the boundary conditions into the matrix.
      call matfil_discreteBC (Rlevels(i)%rmatrix)
      
    end do

    ! Our right-hand-side also needs to know the boundary conditions.
    rrhsBlock%p_rdiscreteBC => Rlevels(NLMAX)%rdiscreteBC
    
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
    
    ! During the linear solver, the boundary conditions are also
    ! frequently imposed to the vectors. But as the linear solver
    ! does not work with the actual solution vectors but with
    ! defect vectors instead.
    ! So, set up a filter chain that filters the defect vector
    ! during the solution process to implement discrete boundary conditions.
    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

    ! Now we have to build up the level information for multigrid.

    ! Create a Multigrid-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    p_RfilterChain => RfilterChain
    call linsol_initMultigrid2 (p_rsolverNode,NLMAX-NLMIN+1,p_RfilterChain)
    
    ! Set up a coarse grid solver.
    ! The coarse grid in multigrid is always grid 1!
    call linsol_getMultigrid2Level (p_rsolverNode,1,p_rlevelInfo)
    call linsol_initUMFPACK4 (p_rlevelInfo%p_rcoarseGridSolver)
    
    ! Now set up the other levels...
    do i = NLMIN+1, NLMAX
    
      ! Create a Jacobi smoother
      !CALL linsol_initJacobi(p_rsmoother)
      
      ! Create an ILU(0) smoother
      call linsol_initMILUs1x1 (p_rsmoother,0,0.0_DP)
      
      ! We will use 4 smoothing steps with damping parameter 0.7
      call linsol_convertToSmoother(p_rsmoother, 4, 0.7_DP)
      
      ! And add this multi-grid level. We will use the same smoother
      ! for pre- and post-smoothing.
      call linsol_getMultigrid2Level (p_rsolverNode,i-NLMIN+1,p_rlevelInfo)
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
    allocate(Rmatrices(NLMIN:NLMAX))
    do i = NLMIN, NLMAX
      call lsysbl_duplicateMatrix (Rlevels(i)%rmatrix,&
          Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    end do
    
    call linsol_setMatrices(p_RsolverNode,Rmatrices(NLMIN:NLMAX))

    ! We can release Rmatrices immediately -- as long as we don't
    ! release Rlevels(i)%rmatrix!
    do i=NLMIN,NLMAX
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
        Rlevels(NLMAX)%rtriangulation,TRIM(sucddir)//'/u2d_1_mg.gmv')
    
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
    
    ! We are finished - but not completely!
    ! Now, clean up so that all the memory is available again.
    !
    ! Release solver data and structure
    call linsol_doneData (p_rsolverNode)
    call linsol_doneStructure (p_rsolverNode)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    call linsol_releaseSolver (p_rsolverNode)
    
    ! Release the block matrix/vectors
    call lsysbl_releaseVector (rtempBlock)
    call lsysbl_releaseVector (rvectorBlock)
    call lsysbl_releaseVector (rrhsBlock)
    do i = NLMAX, NLMIN, -1
      call lsysbl_releaseMatrix (Rlevels(i)%rmatrix)
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
    
    ! Finally release the domain, that's it.
    call boundary_release (rboundary)

  end subroutine
 
!******************************************************************************
  
  !<subroutine>
    subroutine tridef2d_monitorfct(DvertexCoords,Dentries)
  
  
    !<description>
      ! In this function we build the nodewise area distribution out
      ! of an elementwise distribution
    !</description>

    !<inputoutput>
     real(DP), dimension(:,:) :: DvertexCoords
     real(DP), dimension(:) :: Dentries
    !</inputoutput>

    !</subroutine>
    ! local variables
     real(dp),dimension(:,:),allocatable :: Dpoints
     integer :: ive,i1,ipoints
     integer :: iMethod
     real(DP) :: Dist,t,dt,dmin
     iMethod = 1
      
     ipoints = ubound(Dentries,1)
      
      
     select case(iMethod)
       ! case(0) 0.5 + x
       case(0)
         ! loop over all vertices and compute the monitor function
        do ive=1, ipoints
          Dentries(ive) = 0.5_dp + DvertexCoords(1,ive)
         end do
         
       ! Case(1) circle
       case(1)
         ! loop over all vertices and compute the monitor function
         do ive=1,ipoints
           Dist = sqrt((0.5_dp - DvertexCoords(1,ive))**2 + (0.5_dp - DvertexCoords(2,ive))**2)
           ! Good now define the monitor function
           Dist = abs(Dist - 0.2_dp)/0.2_dp
           Dist=max(dist,0.1_dp)
           Dist=min(1.0_dp,dist)
           Dentries(ive)=Dist
         end do
       ! Case(2) do nothing
       case(2)
        do ive=1,ipoints
          Dentries(ive) = 1.0_dp
        end do
       ! Case(3) Elipse
       case(3)
       
        allocate(Dpoints(2,10000))
        dt = 6.28_dp/real(10000)
        t  = 0.0_dp
        do i1=1,10000
         Dpoints(1,i1) = 0.5_dp + 0.1_dp * cos(t)
         Dpoints(2,i1) = 0.5_dp + 0.2_dp * sin(t)
         t = t + dt
        end do
        ! loop over all points on the elipse
        do ive=1,ipoints
          dmin = 10000.0_dp
          do i1=1,10000
            Dist = sqrt((Dpoints(1,i1)-DvertexCoords(1,ive))**2 + (Dpoints(2,i1)-DvertexCoords(2,ive))**2)
            dmin =min(Dist,dmin)
          end do
          Dentries(ive) = dmin
        end do
       !Case 4
       case(4)
         do ive=1,ipoints
           Dist = sqrt((0.5_dp - DvertexCoords(1,ive))**2 + (0.2_dp - DvertexCoords(2,ive))**2)
           ! Good now define the monitor function
           Dist = abs(Dist - 0.05_dp)/0.05_dp
           Dist=max(dist,0.1_dp)
           Dist=min(1.0_dp,dist)
           Dentries(ive)=Dist
         end do
       case default
     end select
      
  end subroutine
  
  ! ***************************************************************************
  
  subroutine gd_MultigridDeformation(rgriddefInfo, rgriddefWork, &
                                        h_Dcontrib,&
                                        bstartNew, blevelHasChanged, bterminate, &
                                        bdegenerated, imgLevelCalc, iiteradapt, ibcIdx,&
                                        def_monitorfct)
  !<description>
    ! This subroutine is the main routine for the grid deformation process, as all
    ! necessary steps are included here. For performing grid deformation, it is sufficient
    ! to define a monitor function and call this routine or to have an error distribution
    ! at hand and call this subroutine.
  !</description>

  !<input>
    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.

    ! structure containing all parameter settings for grid deformation
    type(t_griddefInfo), intent(inout) :: rgriddefInfo
    
    type(t_griddefWork), intent(inout) :: rgriddefWork
    
    ! if true, start from scratch: new vectors, new boundary conditions structures
    logical, intent(in) :: bstartNew

    ! if true, adjust the vectors after level change since a previous deformation call
    logical, intent(in) :: blevelHasChanged

    logical, intent(in) :: bterminate
    
    ! number of adaptive iteration
    integer, intent(in) :: iiterAdapt

    ! multigrid level on which the simulation was computed
    integer, intent(in) :: imgLevelCalc

    ! index of boundary condition related to the deformation PDE
    integer, intent(in):: ibcIdx
  !</input>

    ! flag for grid checking: if true, the deformation process would lead to
    ! a grid with tangled elements
    logical , intent(in):: bdegenerated

  !<inoutput>
    ! handle of vector with elementwise error contributions
    integer, intent(inout):: h_Dcontrib
  !</inoutput>
  
    ! A callback routine for the monitor function
    include 'intf_monitorfct.inc'
    optional :: def_monitorfct
    
  !</subroutine>
  
  ! Deform
  call griddef_performDeformation(rgriddefInfo, rgriddefWork,h_Dcontrib,&
                                  .TRUE., .FALSE., .FALSE., &
                                  .FALSE., imgLevelCalc, 0, 0,&
                                  def_monitorfct)
  
    
  end subroutine ! griddef_MultigridDeformation
 
  
    
end module
