!##############################################################################
!# ****************************************************************************
!# <name> poisson2d_method1_deflgmres </name>
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

module poisson2d_method1_deflgmres

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
  use convection
  use spdiscprojection
!  use sortstrategy
    
  use poisson2d_callback
  
  implicit none

!<types>

!<typeblock description="Type block defining all information about one level">

  type t_level
  
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation (structure of the
    ! solution, trial/test functions,...)
    type(t_blockDiscretisation) :: rdiscretisation
    
    ! Cubature info structure which encapsules the cubature formula
    type(t_scalarCubatureInfo) :: rcubatureInfo
    
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

  subroutine poisson2d_1_deflgmres_quad
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
    type(t_linsolNode), pointer :: p_rsolverNode,p_rcoarseGridSolver,p_rsmoother,p_rprecond
    type(t_linsolNode), pointer :: p_rspecDefl

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(:), pointer :: Rmatrices

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    type(t_filterChain), dimension(1), target :: RfilterChain

    ! One level of multigrid
    type(t_linsolDeflGMRESLevelInfo), pointer :: p_rlevelInfo

    ! NLMIN receives the level of the coarse grid.
    integer :: NLMIN

    ! NLMAX receives the level where we want to solve.
    integer :: NLMAX
    
    ! Error indicator during initialisation of the solver
    integer :: ierror
    
    ! Error of FE function to reference function
    real(DP) :: derror,dtemp,dtemp2
    
    type(t_matrixScalar) :: rtempMatrix
    
    type(t_convStreamDiff2) :: rsdconfig
    type(t_vectorBlock) :: rsdconvection
    
    type(t_blockDiscretisation) :: rvelDiscr
    
    ! Output block for UCD output to VTK file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir
    real(DP), dimension(:), pointer :: p_Ddata, p_Ddata2
    integer, dimension(:), pointer :: p_Kld
    
    ! Some temporary variables
    integer :: i,j,k,l
    integer :: h_Ipermutation
    
    integer :: ipeclet
    real(DP), dimension(6), parameter :: Dpeclet = (/1.0_DP,10.0_DP,20.0_DP,50.0_DP,100.0_DP,200.0_DP/)

    ! Ok, let us start.
    !
    ! We want to solve our Poisson problem on level...
    NLMIN = 3

    do ipeclet = 1,size(Dpeclet)

      do NLMAX = 5,9
      
        ! Allocate memory for all levels
        allocate(Rlevels(NLMIN:NLMAX))

        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Read the domain, read the mesh, refine
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        
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
        
        ! And create information about adjacencies and everything one needs
        ! from a triangulation.
        call tria_initStandardMeshFromRaw (Rlevels(NLMIN)%rtriangulation,&
            rboundary)
        
        ! Now refine the grid for the fine levels.
        do i = NLMIN+1, NLMAX

          ! Refine the grid using the 2-Level-Ordering algorithm
          call tria_refine2LevelOrdering(Rlevels(i-1)%rtriangulation,&
              Rlevels(i)%rtriangulation,rboundary)
          
          ! Create a standard mesh
          call tria_initStandardMeshFromRaw(Rlevels(i)%rtriangulation,&
              rboundary)
        
        end do

        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Set up discretisation structures which tells the code which
        ! finite element to use
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

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
        ! for this solution component:
        do i = NLMIN, NLMAX
          call spdiscr_initDiscr_simple (&
              Rlevels(i)%rdiscretisation%RspatialDiscr(1), &
              EL_Q1,Rlevels(i)%rtriangulation, rboundary)
        end do

        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Set up an cubature info structure to tell the code which cubature
        ! formula to use
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        ! Create an assembly information structure on each level which tells the code
        ! the cubature formula to use. Standard: Gauss 3x3.
        do i = NLMIN, NLMAX
          call spdiscr_createDefCubStructure(&  
              Rlevels(i)%rdiscretisation%RspatialDiscr(1),Rlevels(i)%rcubatureInfo,&
              CUB_GEN_AUTO_G3)
        end do

        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Create a 1x1 block matrix with the operator on every level
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

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
          
          ! And now to the entries of the matrix. For assembling of the entries,
          ! we need a bilinear form, which first has to be set up manually.
          ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
          ! scalar system matrix in 2D.
          rform%itermCount = 3
          rform%Idescriptors(1,1) = DER_DERIV_X
          rform%Idescriptors(2,1) = DER_DERIV_X
          rform%Idescriptors(1,2) = DER_DERIV_Y
          rform%Idescriptors(2,2) = DER_DERIV_Y
          rform%Idescriptors(1,3) = DER_DERIV_Y
          rform%Idescriptors(2,3) = DER_FUNC

          ! In the standard case, we have constant coefficients:
          rform%ballCoeffConstant = .true.
          rform%BconstantCoeff = .true.
          rform%Dcoefficients(1)  = 1.0_DP/Dpeclet(ipeclet)
          rform%Dcoefficients(2)  = 1.0_DP/Dpeclet(ipeclet)
          rform%Dcoefficients(3)  = 1.0_DP

          ! Now we can build the matrix entries.
          ! We specify the callback function coeff_Laplace for the coefficients.
          ! As long as we use constant coefficients, this routine is not used.
          ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
          ! the framework will call the callback routine to get analytical
          ! data.
          call bilf_buildMatrixScalar (rform,.true.,&
              Rlevels(i)%rmatrix%RmatrixBlock(1,1),Rlevels(i)%rcubatureInfo)
               
!          call sstrat_calcCuthillMcKee (Rlevels(i)%rmatrix%RmatrixBlock(1,1), h_Ipermutation)
!          call lsyssc_sortMatrix (Rlevels(i)%rmatrix%RmatrixBlock(1,1),.true.,&
!              SSTRAT_CM,h_Ipermutation)
          
    !      call spdiscr_initBlockDiscr (rvelDiscr, 2, &
    !                                   Rlevels(i)%rtriangulation, rboundary)
    !      call spdiscr_duplicateDiscrSc (Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
    !          rvelDiscr%RspatialDiscr(1), .true.)
    !      call spdiscr_duplicateDiscrSc (Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
    !          rvelDiscr%RspatialDiscr(2), .true.)
    !
    !      call lsysbl_createVectorBlock(rvelDiscr,rsdconvection,.true.)
    !      call lsyssc_clearVector(rsdconvection%RvectorBlock(2),1.0_DP)
    !      rsdconfig%dupsam = 1.0_DP
    !      rsdconfig%dnu = rform%Dcoefficients(1)
    !      rsdconfig%ddelta = 1.0_DP
    !      call conv_streamDiff2Blk2dMat (rsdconfig,Rlevels(i)%rmatrix,rsdconvection,&
    !          rcubatureInfo=Rlevels(i)%rcubatureInfo)
    !      call lsysbl_releaseVector(rsdconvection)
    !      call spdiscr_releaseBlockDiscr(rvelDiscr)
               
        end do
          
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Create RHS and solution vectors
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        ! Next step: Create a RHS vector, a solution vector and a temporary
        ! vector. All are filled with zero.
        call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,rrhsBlock,.true.)
        call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,rvectorBlock,.true.)
        call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,rtempBlock,.true.)
          
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
            rlinform,.true.,rrhsBlock%RvectorBlock(1),Rlevels(NLMAX)%rcubatureInfo,coeff_RHS_2D)
        
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Discretise the boundary conditions and apply them to the matrix/RHS/sol.
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        
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
          
          ! Edge 4 of boundary component 1. That is it.
          call boundary_createRegion(rboundary,1,4,rboundaryRegion)
          call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                            rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                            getBoundaryValues_2D)
          
          ! Assign the BC`s to the matrix. That way, these
          ! boundary conditions are always connected to that matrix.
          call lsysbl_assignDiscreteBC(Rlevels(i)%rmatrix,Rlevels(i)%rdiscreteBC)
      
          ! Also implement the boundary conditions into the matrix.
          call matfil_discreteBC (Rlevels(i)%rmatrix)
          
        end do

        ! Our right-hand-side/solution/temp vectors also needs to 
        ! know the boundary conditions.
        call lsysbl_assignDiscreteBC(rrhsBlock,Rlevels(NLMAX)%rdiscreteBC)
        call lsysbl_assignDiscreteBC(rvectorBlock,Rlevels(NLMAX)%rdiscreteBC)
        call lsysbl_assignDiscreteBC(rtempBlock,Rlevels(NLMAX)%rdiscreteBC)

        ! Next step is to implement boundary conditions into the RHS,
        ! solution and matrix. This is done using a vector/matrix filter
        ! for discrete boundary conditions.
        ! The discrete boundary conditions are already attached to the
        ! vectors. Call the appropriate vector filter that
        ! modifies the vectors according to the boundary conditions.
        call vecfil_discreteBCrhs (rrhsBlock)
        call vecfil_discreteBCsol (rvectorBlock)
        
!        call lsyssc_sortVectorInSitu (rrhsBlock%RvectorBlock(1),&
!            rtempBlock%RvectorBlock(1),-SSTRAT_CM,h_Ipermutation)
!        call lsyssc_sortVectorInSitu (rvectorBlock%RvectorBlock(1),&
!            rtempBlock%RvectorBlock(1),-SSTRAT_CM,h_Ipermutation)
        
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Set up a linear solver
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

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
        call linsol_initDeflGMRES (p_rsolverNode,NLMAX-NLMIN+1,Rfilter=RfilterChain)

        ! Now set up the other levels...
        do i = NLMIN, NLMAX
        
          ! And add this multi-grid level. We will use the same smoother
          ! for pre- and post-smoothing.
          call linsol_getDeflGMRESLevel (p_rsolverNode,i-NLMIN+1,p_rlevelInfo)
          
          call linsol_initJacobi (p_rlevelInfo%p_rpreconditioner)
          !call linsol_initSOR (p_rlevelInfo%p_rpreconditioner,1.0_DP)
          !call linsol_initMILUs1x1 (p_rlevelInfo%p_rpreconditioner,0,0.0_DP)
          
    !      call lsyssc_calcGerschgorin (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
    !          dtemp,p_rlevelInfo%dmaxEigenvalue)
              
    !      call lsyssc_duplicateMatrix (&
    !          Rlevels(i)%rmatrix%RmatrixBlock(1,1),rtempMatrix,LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)
    !      call lsyssc_lumpMatrixScalar (rtempMatrix,LSYSSC_LUMP_STD,.false.)
    !      call lsyssc_calcGerschgorin (rtempMatrix,dtemp,dtemp2)
    !      p_rlevelInfo%dmaxEigenvalue = p_rlevelInfo%dmaxEigenvalue / dtemp
              
          p_rlevelInfo%dmaxEigenvalue = 1.0_DP
          if (i .eq. NLMIN) then
            p_rlevelInfo%ikrylowDim = 4
          else
            p_rlevelInfo%ikrylowDim = 2
          end if
          
          ! Use 5 iterations. This should actually be testet with =1 !
          ! Use =5 for getting the reference results from the first test!
          p_rlevelInfo%niterations = 5
          
        end do
        p_rsolverNode%p_rsubnodeDeflGMRES%drelax = 1.0_DP
        !p_rsolverNode%p_rsubnodeDeflGMRES%brightPrec = .true.
        
        ! Set the output level of the solver to 2 for some output
        p_rsolvernode%ioutputLevel = 0
        p_rsolvernode%nmaxIterations = 1000
        p_rsolvernode%depsRel = 1E-6
        
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

        ! We can release Rmatrices immediately -- as long as we do not
        ! release Rlevels(i)%rmatrix!
        do i=NLMIN,NLMAX
          call lsysbl_releaseMatrix (Rmatrices(i))
        end do
        deallocate(Rmatrices)

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
        call linsol_solveAdaptively (p_rsolverNode,rvectorBlock,rrhsBlock,rtempBlock)
        
        print *,"Level=",NLMAX,", Pe=",1.0_DP/Dpeclet(ipeclet),", ite=",p_rsolverNode%iiterations
        
!        call lsyssc_sortVectorInSitu (rvectorBlock%RvectorBlock(1),rtempBlock%RvectorBlock(1),-SSTRAT_CM)
        
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Postprocessing of the solution
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        
        ! That is it, rvectorBlock now contains our solution. We can now
        ! start the postprocessing.
        !
        ! Get the path for writing postprocessing files from the environment variable
        ! $UCDDIR. If that does not exist, write to the directory "./gmv".
        if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = './gmv'
    
        ! Start UCD export to VTK file:
        call ucd_startVTK (rexport,UCD_FLAG_STANDARD,&
            Rlevels(NLMAX)%rtriangulation,trim(sucddir)//'/u2d_1_deflgmres.vtk')
        
        call lsyssc_getbase_double (rvectorBlock%RvectorBlock(1),p_Ddata)
        call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
        
        ! Write the file to disc, that is it.
        call ucd_write (rexport)
        call ucd_release (rexport)
        
    !    ! Calculate the error to the reference function.
    !    call pperr_scalar (PPERR_L2ERROR,derror,rvectorBlock%RvectorBlock(1),&
    !                       getReferenceFunction_2D)
    !    call output_line ('L2-error: ' // sys_sdEL(derror,10) )
    !
    !    call pperr_scalar (PPERR_H1ERROR,derror,rvectorBlock%RvectorBlock(1),&
    !                       getReferenceFunction_2D)
    !    call output_line ('H1-error: ' // sys_sdEL(derror,10) )
        
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
        
        ! Release the block matrices/vectors
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

        ! Release the cubature info structures
        do i=NLMAX,NLMIN,-1
          call spdiscr_releaseCubStructure(Rlevels(i)%rcubatureInfo)
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

        ! Finally release the domain, that is it.
        call boundary_release (rboundary)
        
      end do
      
      print *
    
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getBoundaryValues_b1_2D (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
                                      cinfoNeeded,iwhere,dwhere, Dvalues, rcollection)
  
  use collection
  use spatialdiscretisation
  use discretebc
  
!<description>
  ! This subroutine is called during the discretisation of boundary
  ! conditions. It calculates a special quantity on the boundary, which is
  ! then used by the discretisation routines to generate a discrete
  ! 'snapshot' of the (actually analytic) boundary conditions.
!</description>
  
!<input>
  ! Component specifier.
  ! For Dirichlet boundary:
  !   Icomponents(1) defines the number of the boundary component, the value
  !   should be calculated for (e.g. 1=1st solution component, e.g. X-velocitry,
  !   2=2nd solution component, e.g. Y-velocity,...)
  integer, dimension(:), intent(in)                           :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  
  ! Boundary region that is currently being processed.
  type(t_boundaryRegion), intent(in)                          :: rboundaryRegion
  
  ! The element number on the boundary which is currently being processed
  integer, intent(in)                                         :: ielement
  
  ! The type of information, the routine should calculate. One of the
  ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
  ! to return one or multiple information value in the result array.
  integer, intent(in)                                         :: cinfoNeeded
  
  ! A reference to a geometric object where information should be computed.
  ! cinfoNeeded=DISCBC_NEEDFUNC :
  !   iwhere = number of the point in the triangulation or
  !          = 0, if only the parameter value of the point is known; this
  !               can be found in dwhere,
  ! cinfoNeeded=DISCBC_NEEDDERIV :
  !   iwhere = number of the point in the triangulation or
  !          = 0, if only the parameter value of the point is known; this
  !               can be found in dwhere,
  ! cinfoNeeded=DISCBC_NEEDINTMEAN :
  !   iwhere = number of the edge where the value integral mean value
  !            should be computed
  integer, intent(in)                                          :: iwhere

  ! A reference to a geometric object where information should be computed.
  ! cinfoNeeded=DISCBC_NEEDFUNC :
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDDERIV :
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDINTMEAN :
  !   dwhere = 0 (not used)
  real(DP), intent(in)                                        :: dwhere
    
  ! Optional: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional                 :: rcollection

!</input>

!<output>
  ! This array receives the calculated information. If the caller
  ! only needs one value, the computed quantity is put into Dvalues(1).
  ! If multiple values are needed, they are collected here (e.g. for
  ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
  real(DP), dimension(:), intent(out)                         :: Dvalues
!</output>
  
!</subroutine>

    ! To get the X/Y-coordinates of the boundary point, use:
    !
    REAL(DP) :: dx,dy
    
    CALL boundary_getCoords(rdiscretisation%p_rboundary, &
        rboundaryRegion%iboundCompIdx, dwhere, dx, dy)

    ! Return zero Dirichlet boundary values for all situations.
    Dvalues(1) = 0.0_DP !dx*dy
    
    if ((dx .gt. 0.0_DP) .and. (dx .lt. 2.2_DP) .and. &
        (dy .gt. 0.0_DP) .and. (dy .lt. 0.41_DP)) then
      Dvalues(1) = 1.0_DP
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine poisson2d_1_deflgmres_bench1
  
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
    type(t_linsolNode), pointer :: p_rsolverNode,p_rcoarseGridSolver,p_rsmoother,p_rprecond
    type(t_linsolNode), pointer :: p_rspecDefl

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(:), pointer :: Rmatrices

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    type(t_filterChain), dimension(1), target :: RfilterChain

    ! One level of multigrid
    type(t_linsolDeflGMRESLevelInfo), pointer :: p_rlevelInfo

    ! NLMIN receives the level of the coarse grid.
    integer :: NLMIN

    ! NLMAX receives the level where we want to solve.
    integer :: NLMAX
    
    ! Error indicator during initialisation of the solver
    integer :: ierror
    
    ! Error of FE function to reference function
    real(DP) :: derror,dtemp,dtemp2
    
    type(t_matrixScalar) :: rtempMatrix
    
    type(t_convStreamDiff2) :: rsdconfig
    type(t_vectorBlock) :: rsdconvection
    
    type(t_blockDiscretisation) :: rvelDiscr
    
    ! Output block for UCD output to VTK file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir
    real(DP), dimension(:), pointer :: p_Ddata, p_Ddata2
    integer, dimension(:), pointer :: p_Kld
    
    ! Some temporary variables
    integer :: i,j,k,l
    integer :: h_Ipermutation
    
    integer :: ipeclet
    real(DP), dimension(8), parameter :: Dpeclet = (/1.0_DP,10.0_DP,20.0_DP,50.0_DP,100.0_DP,200.0_DP,500.0_DP,1000.0_DP/)

    ! Ok, let us start.
    !
    ! We want to solve our Poisson problem on level...
    NLMIN = 1

    do ipeclet = 1,size(Dpeclet)

      do NLMAX = 2,5
      
        ! Allocate memory for all levels
        allocate(Rlevels(NLMIN:NLMAX))

        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Read the domain, read the mesh, refine
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        
        ! Get the path $PREDIR from the environment, where to read .prm/.tri files
        ! from. If that does not exist, write to the directory "./pre".
        if (.not. sys_getenv_string("PREDIR", spredir)) spredir = './pre'

        ! At first, read in the parametrisation of the boundary and save
        ! it to rboundary.
        call boundary_read_prm(rboundary, trim(spredir)//'/bench1.prm')
            
        ! Now read in the basic triangulation into our coarse level.
        call tria_readTriFile2D (Rlevels(NLMIN)%rtriangulation, &
                                trim(spredir)//'/bench1.tri', rboundary)
        
        ! Refine it.
        call tria_quickRefine2LevelOrdering (NLMIN-1,&
            Rlevels(NLMIN)%rtriangulation,rboundary)
        
        ! And create information about adjacencies and everything one needs
        ! from a triangulation.
        call tria_initStandardMeshFromRaw (Rlevels(NLMIN)%rtriangulation,&
            rboundary)
        
        ! Now refine the grid for the fine levels.
        do i = NLMIN+1, NLMAX

          ! Refine the grid using the 2-Level-Ordering algorithm
          call tria_refine2LevelOrdering(Rlevels(i-1)%rtriangulation,&
              Rlevels(i)%rtriangulation,rboundary)
          
          ! Create a standard mesh
          call tria_initStandardMeshFromRaw(Rlevels(i)%rtriangulation,&
              rboundary)
        
        end do

        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Set up discretisation structures which tells the code which
        ! finite element to use
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

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
        ! for this solution component:
        do i = NLMIN, NLMAX
          call spdiscr_initDiscr_simple (&
              Rlevels(i)%rdiscretisation%RspatialDiscr(1), &
              EL_Q1,Rlevels(i)%rtriangulation, rboundary)
        end do

        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Set up an cubature info structure to tell the code which cubature
        ! formula to use
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        ! Create an assembly information structure on each level which tells the code
        ! the cubature formula to use. Standard: Gauss 3x3.
        do i = NLMIN, NLMAX
          call spdiscr_createDefCubStructure(&  
              Rlevels(i)%rdiscretisation%RspatialDiscr(1),Rlevels(i)%rcubatureInfo,&
              CUB_GEN_AUTO_G4)
        end do

        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Create a 1x1 block matrix with the operator on every level
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

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
          
          ! And now to the entries of the matrix. For assembling of the entries,
          ! we need a bilinear form, which first has to be set up manually.
          ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
          ! scalar system matrix in 2D.
          rform%itermCount = 3
          rform%Idescriptors(1,1) = DER_DERIV_X
          rform%Idescriptors(2,1) = DER_DERIV_X
          rform%Idescriptors(1,2) = DER_DERIV_Y
          rform%Idescriptors(2,2) = DER_DERIV_Y
          rform%Idescriptors(1,3) = DER_DERIV_X
          rform%Idescriptors(2,3) = DER_FUNC

          ! In the standard case, we have constant coefficients:
          rform%ballCoeffConstant = .true.
          rform%BconstantCoeff = .true.
          rform%Dcoefficients(1)  = 1.0_DP/Dpeclet(ipeclet)
          rform%Dcoefficients(2)  = 1.0_DP/Dpeclet(ipeclet)
          rform%Dcoefficients(3)  = 1.0_DP

          ! Now we can build the matrix entries.
          ! We specify the callback function coeff_Laplace for the coefficients.
          ! As long as we use constant coefficients, this routine is not used.
          ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
          ! the framework will call the callback routine to get analytical
          ! data.
          call bilf_buildMatrixScalar (rform,.true.,&
              Rlevels(i)%rmatrix%RmatrixBlock(1,1),Rlevels(i)%rcubatureInfo)
               
!          call sstrat_calcCuthillMcKee (Rlevels(i)%rmatrix%RmatrixBlock(1,1), h_Ipermutation)
!          call lsyssc_sortMatrix (Rlevels(i)%rmatrix%RmatrixBlock(1,1),.true.,&
!              SSTRAT_CM,h_Ipermutation)
          
    !      call spdiscr_initBlockDiscr (rvelDiscr, 2, &
    !                                   Rlevels(i)%rtriangulation, rboundary)
    !      call spdiscr_duplicateDiscrSc (Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
    !          rvelDiscr%RspatialDiscr(1), .true.)
    !      call spdiscr_duplicateDiscrSc (Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
    !          rvelDiscr%RspatialDiscr(2), .true.)
    !
    !      call lsysbl_createVectorBlock(rvelDiscr,rsdconvection,.true.)
    !      call lsyssc_clearVector(rsdconvection%RvectorBlock(2),1.0_DP)
    !      rsdconfig%dupsam = 1.0_DP
    !      rsdconfig%dnu = rform%Dcoefficients(1)
    !      rsdconfig%ddelta = 1.0_DP
    !      call conv_streamDiff2Blk2dMat (rsdconfig,Rlevels(i)%rmatrix,rsdconvection,&
    !          rcubatureInfo=Rlevels(i)%rcubatureInfo)
    !      call lsysbl_releaseVector(rsdconvection)
    !      call spdiscr_releaseBlockDiscr(rvelDiscr)
               
        end do
          
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Create RHS and solution vectors
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        ! Next step: Create a RHS vector, a solution vector and a temporary
        ! vector. All are filled with zero.
        call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,rrhsBlock,.true.)
        call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,rvectorBlock,.true.)
        call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,rtempBlock,.true.)
          
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
            rlinform,.true.,rrhsBlock%RvectorBlock(1),Rlevels(NLMAX)%rcubatureInfo,coeff_RHS_2D)
        
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Discretise the boundary conditions and apply them to the matrix/RHS/sol.
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        
        do i = NLMIN, NLMAX
        
          ! Initialise the discrete BC structure
          call bcasm_initDiscreteBC(Rlevels(i)%rdiscreteBC)

          ! On edge 1 of boundary component 1 add Dirichlet boundary conditions.
          call boundary_createRegion(rboundary,1,1,rboundaryRegion)
          call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                            rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                            getBoundaryValues_b1_2D)
                                   
          ! Now to the edge 2 of boundary component 1 the domain.
          !call boundary_createRegion(rboundary,1,2,rboundaryRegion)
          !call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
          !                                  rboundaryRegion,Rlevels(i)%rdiscreteBC,&
          !                                  getBoundaryValues_b1_2D)
                                   
          ! Edge 3 of boundary component 1.
          call boundary_createRegion(rboundary,1,3,rboundaryRegion)
          call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                            rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                            getBoundaryValues_b1_2D)
          
          ! Edge 4 of boundary component 1. That is it.
          call boundary_createRegion(rboundary,1,4,rboundaryRegion)
          call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                            rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                            getBoundaryValues_b1_2D)
          
          call boundary_createRegion(rboundary,2,1,rboundaryRegion)
          call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                            rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                            getBoundaryValues_b1_2D)
          
          ! Assign the BC`s to the matrix. That way, these
          ! boundary conditions are always connected to that matrix.
          call lsysbl_assignDiscreteBC(Rlevels(i)%rmatrix,Rlevels(i)%rdiscreteBC)
      
          ! Also implement the boundary conditions into the matrix.
          call matfil_discreteBC (Rlevels(i)%rmatrix)
          
        end do

        ! Our right-hand-side/solution/temp vectors also needs to 
        ! know the boundary conditions.
        call lsysbl_assignDiscreteBC(rrhsBlock,Rlevels(NLMAX)%rdiscreteBC)
        call lsysbl_assignDiscreteBC(rvectorBlock,Rlevels(NLMAX)%rdiscreteBC)
        call lsysbl_assignDiscreteBC(rtempBlock,Rlevels(NLMAX)%rdiscreteBC)

        ! Next step is to implement boundary conditions into the RHS,
        ! solution and matrix. This is done using a vector/matrix filter
        ! for discrete boundary conditions.
        ! The discrete boundary conditions are already attached to the
        ! vectors. Call the appropriate vector filter that
        ! modifies the vectors according to the boundary conditions.
        call vecfil_discreteBCrhs (rrhsBlock)
        call vecfil_discreteBCsol (rvectorBlock)
        
!        call lsyssc_sortVectorInSitu (rrhsBlock%RvectorBlock(1),&
!            rtempBlock%RvectorBlock(1),-SSTRAT_CM,h_Ipermutation)
!        call lsyssc_sortVectorInSitu (rvectorBlock%RvectorBlock(1),&
!            rtempBlock%RvectorBlock(1),-SSTRAT_CM,h_Ipermutation)
        
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Set up a linear solver
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

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
        call linsol_initDeflGMRES (p_rsolverNode,NLMAX-NLMIN+1,Rfilter=RfilterChain)

        ! Now set up the other levels...
        do i = NLMIN, NLMAX
        
          ! And add this multi-grid level. We will use the same smoother
          ! for pre- and post-smoothing.
          call linsol_getDeflGMRESLevel (p_rsolverNode,i-NLMIN+1,p_rlevelInfo)
          
          call linsol_initJacobi (p_rlevelInfo%p_rpreconditioner)
          !call linsol_initSOR (p_rlevelInfo%p_rpreconditioner,1.0_DP)
          !call linsol_initMILUs1x1 (p_rlevelInfo%p_rpreconditioner,0,0.0_DP)
          
    !      call lsyssc_calcGerschgorin (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
    !          dtemp,p_rlevelInfo%dmaxEigenvalue)
              
    !      call lsyssc_duplicateMatrix (&
    !          Rlevels(i)%rmatrix%RmatrixBlock(1,1),rtempMatrix,LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)
    !      call lsyssc_lumpMatrixScalar (rtempMatrix,LSYSSC_LUMP_STD,.false.)
    !      call lsyssc_calcGerschgorin (rtempMatrix,dtemp,dtemp2)
    !      p_rlevelInfo%dmaxEigenvalue = p_rlevelInfo%dmaxEigenvalue / dtemp
              
          p_rlevelInfo%dmaxEigenvalue = 1.0_DP
          if (i .eq. NLMIN) then
            p_rlevelInfo%ikrylowDim = 4
          else
            p_rlevelInfo%ikrylowDim = 2
          end if

          ! Use 5 iterations. This should actually be testet with =1 !!!
          ! Use =5 for getting the reference results from the first test!
          p_rlevelInfo%niterations = 5
          
        end do
        p_rsolverNode%p_rsubnodeDeflGMRES%drelax = 1.0_DP
        !p_rsolverNode%p_rsubnodeDeflGMRES%brightPrec = .true.
        
        ! Set the output level of the solver to 2 for some output
        p_rsolvernode%ioutputLevel = 0
        p_rsolvernode%nmaxIterations = 1000
        p_rsolvernode%depsRel = 1E-6
        
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

        ! We can release Rmatrices immediately -- as long as we do not
        ! release Rlevels(i)%rmatrix!
        do i=NLMIN,NLMAX
          call lsysbl_releaseMatrix (Rmatrices(i))
        end do
        deallocate(Rmatrices)

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
        call linsol_solveAdaptively (p_rsolverNode,rvectorBlock,rrhsBlock,rtempBlock)
        
        print *,"Level=",NLMAX,", Pe=",1.0_DP/Dpeclet(ipeclet),", ite=",p_rsolverNode%iiterations
        
!        call lsyssc_sortVectorInSitu (rvectorBlock%RvectorBlock(1),rtempBlock%RvectorBlock(1),-SSTRAT_CM)
        
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Postprocessing of the solution
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        
        ! That is it, rvectorBlock now contains our solution. We can now
        ! start the postprocessing.
        !
        ! Get the path for writing postprocessing files from the environment variable
        ! $UCDDIR. If that does not exist, write to the directory "./gmv".
        if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = './gmv'
    
        ! Start UCD export to VTK file:
        call ucd_startVTK (rexport,UCD_FLAG_STANDARD,&
            Rlevels(NLMAX)%rtriangulation,trim(sucddir)//'/u2d_1_deflgmres.vtk')
        
        nullify(p_Ddata)
        call spdp_projectToVertices (rvectorBlock%RvectorBlock(1), p_Ddata, DER_FUNC)
        !call lsyssc_getbase_double (rvectorBlock%RvectorBlock(1),p_Ddata)
        call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
        deallocate(p_Ddata)
        
        ! Write the file to disc, that is it.
        call ucd_write (rexport)
        call ucd_release (rexport)
        
!         print *,"ok."
!         read (*,*)
        
    !    ! Calculate the error to the reference function.
    !    call pperr_scalar (PPERR_L2ERROR,derror,rvectorBlock%RvectorBlock(1),&
    !                       getReferenceFunction_2D)
    !    call output_line ('L2-error: ' // sys_sdEL(derror,10) )
    !
    !    call pperr_scalar (PPERR_H1ERROR,derror,rvectorBlock%RvectorBlock(1),&
    !                       getReferenceFunction_2D)
    !    call output_line ('H1-error: ' // sys_sdEL(derror,10) )
        
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
        
        ! Release the block matrices/vectors
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

        ! Release the cubature info structures
        do i=NLMAX,NLMIN,-1
          call spdiscr_releaseCubStructure(Rlevels(i)%rcubatureInfo)
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

        ! Finally release the domain, that is it.
        call boundary_release (rboundary)
        
      end do
      
      print *
    
    end do

  end subroutine

end module
