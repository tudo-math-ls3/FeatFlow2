!##############################################################################
!# ****************************************************************************
!# <name> dg2d_method0_simple </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple Poisson
!# problem with constant coefficients on a simple domain.
!# </purpose>
!##############################################################################

module dg2d_multigridscalar

  use fsystem
  use stdoperators
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
  use genoutput
  use dg2d_routines
  use collection
  use linearalgebra
  use paramlist
  use matrixio

  implicit none
  
  
  !<types>

!<typeblock description="Type block defining all information about one level">

  type t_problem_lvl
  
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation
    
    ! An object for additional triangulation data for DG
    type(t_additionalTriaData) :: raddTriaData
    
    ! An object specifying the discretisation (structure of the
    ! solution, trial/test functions,...)
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
    ! A system matrix for that specific level. The matrix will receive the
    ! discrete Laplace operator.
    type(t_matrixBlock) :: rmatrix

    ! A variable describing the discrete boundary conditions.
    type(t_discreteBC) :: rdiscreteBC
  
  end type
  
!</typeblock>


!<typeblock description="Application-specific type block for poisson problem">

  type t_problem
  
    ! Minimum refinement level; = Level i in RlevelInfo
    integer :: ilvmin
    
    ! Maximum refinement level
    integer :: ilvmax

    ! An object for saving the domain:
    type(t_boundary) :: rboundary

    ! A solution vector and a RHS vector on the finest level.
    type(t_vectorBlock) :: rvector,rrhs

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode

    ! An array of t_problem_lvl structures, each corresponding
    ! to one level of the discretisation.
    type(t_problem_lvl), dimension(:), pointer :: RlevelInfo
    
    ! A collection object that saves structural data and some
    ! problem-dependent information which is e.g. passed to
    ! callback routines.
    type(t_collection) :: rcollection

  end type

!</typeblock>

!</types>
  

contains


  ! ***************************************************************************

!<subroutine>

  subroutine dgmgs_initParamTriang (rproblem,rparlist)
  
!<description>
  ! This routine initialises the parametrisation and triangulation of the
  ! domain. The corresponding .prm/.tri files are read from disc and
  ! the triangulation is refined as described by the parameter ilv.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

  ! local variables
  integer :: i
    ! Minimum refinement level of the mesh; = coarse grid = level 1
  integer :: ilvmin
  
  ! Maximum refinement level
  integer :: ilvmax
  
  ! Filename of parametrisation
  character(LEN=SYS_STRLEN) :: sstring
  
  ! We want to solve our problem on level... Default=1
    call parlst_getvalue_int(rparlist, 'TRIANGULATION', 'NLMIN', ilvmin)
  
  ! We want to solve our problem on level... Default=1
    call parlst_getvalue_int(rparlist, 'TRIANGULATION', 'NLMAX', ilvmax)
  
    ! Initialise the level in the problem structure
    rproblem%ilvmin = ilvmin
    rproblem%ilvmax = ilvmax
    allocate(rproblem%RlevelInfo(ilvmin:ilvmax))
    
    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call parlst_getvalue_string (rparlist, 'TRIANGULATION', &
         'prmname', sstring)
    call boundary_read_prm(rproblem%rboundary, trim(sstring))
        
    ! Now read in the basic triangulation.
    call parlst_getvalue_string (rparlist, 'TRIANGULATION', &
         'triname', sstring)
    call tria_readTriFile2D (rproblem%RlevelInfo(rproblem%ilvmin)%rtriangulation, &
         trim(sstring),rproblem%rboundary)
        
    ! Refine the mesh up to the minimum level
    call tria_quickRefine2LevelOrdering(rproblem%ilvmin-1,&
        rproblem%RlevelInfo(rproblem%ilvmin)%rtriangulation,rproblem%rboundary)

    ! Create information about adjacencies and everything one needs from
    ! a triangulation. Afterwards, we have the coarse mesh.
    call tria_initStandardMeshFromRaw (&
        rproblem%RlevelInfo(rproblem%ilvmin)%rtriangulation,rproblem%rboundary)
        
    ! Calculate additional triangulation data, as the normal vectors and local edge numbers
    call addTriaData(rproblem%RlevelInfo(rproblem%ilvmin)%rtriangulation,&
                     rproblem%RlevelInfo(rproblem%ilvmin)%raddTriaData)
    
    ! Now, refine to level up to nlmax.
    do i=rproblem%ilvmin+1,rproblem%ilvmax
      call tria_refine2LevelOrdering (rproblem%RlevelInfo(i-1)%rtriangulation,&
          rproblem%RlevelInfo(i)%rtriangulation, rproblem%rboundary)
      call tria_initStandardMeshFromRaw (rproblem%RlevelInfo(i)%rtriangulation,&
          rproblem%rboundary)
      call addTriaData(rproblem%RlevelInfo(i)%rtriangulation,&
                       rproblem%RlevelInfo(i)%raddTriaData)
    end do

  end subroutine
  
  
  
  ! ***************************************************************************

!<subroutine>

  subroutine dgmgs_initDiscretisation (rproblem,rparlist)
  
!<description>
  ! This routine initialises the discretisation structure of the underlying
  ! problem and saves it to the problem structure.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

  ! local variables
  integer :: I, ielementType
  
    ! An object for saving the domain:
    type(t_boundary), pointer :: rboundary
    
    ! An object for saving the triangulation on the domain
    type(t_triangulation), pointer :: p_rtriangulation

    ! An object for the block discretisation on one level
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
    ! Type of finite element to use
    call parlst_getvalue_int(rparlist, 'TRIANGULATION', 'FEkind', ielementType, 0)
    
    select case (ielementType)
    case (0)
       ielementType = EL_DG_T0_2D
    case (1)
       ielementType = EL_DG_T1_2D
    case (2)
       ielementType = EL_DG_T2_2D
    case (3)
       ielementType = EL_DG_Q1_2D
    case (4)
       ielementType = EL_DG_Q2_2D
    end select
    
    do i=rproblem%ilvmin,rproblem%ilvmax
      ! Ask the problem structure to give us the boundary and triangulation.
      ! We need it for the discretisation.
      rboundary => rproblem%rboundary
      p_rtriangulation => rproblem%RlevelInfo(i)%rtriangulation
      
      ! Now we can start to initialise the discretisation. At first, set up
      ! a block discretisation structure that specifies the blocks in the
      ! solution vector. In this simple problem, we only have one block.
      allocate(p_rdiscretisation)
      call spdiscr_initBlockDiscr (p_rdiscretisation,1,&
                                   p_rtriangulation, rboundary)

      ! Save the discretisation structure to our local LevelInfo structure
      ! for later use.
      rproblem%RlevelInfo(i)%p_rdiscretisation => p_rdiscretisation

      ! p_rdiscretisation%Rdiscretisations is a list of scalar
      ! discretisation structures for every component of the solution vector.
      ! Initialise the first element of the list to specify the element
      ! and cubature rule for this solution component:
      call spdiscr_initDiscr_simple ( &
                  p_rdiscretisation%RspatialDiscr(1), &
                  ielementType,CUB_G6_2D, &
                  p_rtriangulation, rboundary)

    end do
                                   
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine dgmgs_initMatVec (rproblem,rparlist)
  
!<description>
  ! Calculates the system matrix and RHS vector of the linear system
  ! by discretising the problem with the default discretisation structure
  ! in the problem structure.
  ! Sets up a solution vector for the linear system.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

  ! local variables
  integer :: i, ielementtype
  
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform, rlinformedge, rlinformconv
    
    ! A pointer to the system matrix and the RHS/solution vectors.
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_vectorBlock), pointer :: p_rrhs,p_rvector
    type(t_vectorScalar) :: rvectorSolTemp

    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
  
      
    do i=rproblem%ilvmin,rproblem%ilvmax
      ! Ask the problem structure to give us the discretisation structure
      p_rdiscretisation => rproblem%RlevelInfo(i)%p_rdiscretisation
      
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      
      ! Initialise the block matrix with default values based on
      ! the discretisation.
      call lsysbl_createMatBlockByDiscr (p_rdiscretisation,p_rmatrix)

      ! Now as the discretisation is set up, we can start to generate
      ! the structure of the system matrix which is to solve.
      ! We create that directly in the block (1,1) of the block matrix
      ! using the discretisation structure of the first block.
      call bilf_createMatrixStructure (&
                p_rdiscretisation%RspatialDiscr(1),LSYSSC_MATRIX9,&
                p_rmatrix%RmatrixBlock(1,1),&
                p_rdiscretisation%RspatialDiscr(1),&
                BILF_MATC_EDGEBASED)
                
       ! Build mass matrix
       rform%itermCount = 1
       rform%Idescriptors(1,1) = DER_FUNC
       rform%Idescriptors(2,1) = DER_FUNC
       rform%ballCoeffConstant = .true.
       rform%BconstantCoeff = .true.
       rform%Dcoefficients(1)  = 1.0
       call bilf_buildMatrixScalar (rform,.true.,p_rmatrix%RmatrixBlock(1,1))
                
      ! Create temporary empty solution vector which is needed to build the matrices
      call lsyssc_createVecIndMat (p_rmatrix%RmatrixBlock(1,1),rvectorSolTemp,.true.)
                
       
       ! First calculate the matrix for the cell terms
       rform%itermCount = 2
       rform%Idescriptors(1,1) = DER_FUNC
       rform%Idescriptors(2,1) = DER_DERIV_X
       rform%Idescriptors(1,2) = DER_FUNC
       rform%Idescriptors(2,2) = DER_DERIV_Y
       rform%ballCoeffConstant = .false.
       rform%BconstantCoeff = .false.
       !rcollection%p_rvectorQuickAccess1 => rsolBlock
       
       ! Type of finite element to use
       call parlst_getvalue_int(rparlist, 'TRIANGULATION', 'FEkind', ielementType)
       
       if (ielementtype.ne.0) then
          call bilf_buildMatrixScalar2 (rform, .true., p_rmatrix%RmatrixBlock(1,1),&
               fcoeff_MatrixScalarMgCell)!,rcollection)!,rscalarAssemblyInfo)
       else
          call lsyssc_clearMatrix(p_rmatrix%RmatrixBlock(1,1))
       end if

       call lsyssc_scaleMatrix (p_rmatrix%RmatrixBlock(1,1),-1.0_DP)


       ! Next calculate the edge terms
       rform%itermCount = 2
       rform%Idescriptors(1,1) = DER_FUNC
       rform%Idescriptors(2,1) = DER_FUNC
       rform%Idescriptors(1,2) = DER_FUNC
       rform%Idescriptors(2,2) = DER_FUNC
       rform%ballCoeffConstant = .false.
       rform%BconstantCoeff = .false.
       !rcollection%p_rvectorQuickAccess1 => rsolBlock
       !rcollection%Dquickaccess(1) = dt
       call bilf_dg_buildMatrixScEdge2D_ss (rform, CUB_G5_1D, .false., p_rmatrix%RmatrixBlock(1,1),&
            rvectorSolTemp, rproblem%RlevelInfo(i)%raddTriaData,&
            flux_dg_MatrixScalarMgEdge)!,&
            !rcollection)!, cconstrType)
       
       ! Deallocate temporary solution vector
       call lsyssc_releaseVector (rvectorSolTemp)

    end do

    ! (Only) on the finest level, we need to calculate a RHS vector
    ! and to allocate a solution vector.
    
    p_rrhs    => rproblem%rrhs
    p_rvector => rproblem%rvector

    ! Although we could manually create the solution/RHS vector,
    ! the easiest way to set up the vector structure is
    ! to create it by using our matrix as template:
    call lsysbl_createVecBlockIndMat (p_rmatrix,p_rrhs, .false.)
    call lsysbl_createVecBlockIndMat (p_rmatrix,p_rvector, .false.)

    ! Clear the solution vector on the finest level.
    call lsysbl_clearVector(rproblem%rvector)


    ! Create RHS-Vector

    ! At first set up the corresponding linear form (f,Phi_j):
    rlinformedge%itermCount = 1
    rlinformedge%Idescriptors(1) = DER_FUNC2D
    rproblem%rcollection%p_rvectorQuickAccess1 => rproblem%rvector
    
    ! Now use the dg-function for the edge terms
    call linf_dg_buildVectorScalarEdge2d (rlinformedge, CUB_G5_1D, .true.,&
         p_rrhs%RvectorBlock(1),&
         p_rvector%RvectorBlock(1),&
         rproblem%RlevelInfo(rproblem%ilvmax)%raddTriaData,&
         flux_dg_VectorScalarMgEdge,&
         rproblem%rcollection)

    call lsysbl_scaleVector (p_rrhs,-1.0_DP)

    ! Type of finite element to use
    call parlst_getvalue_int(rparlist, 'TRIANGULATION', 'FEkind', ielementType)

    ! Then add the cell terms
    if (ielementtype.ne.0) then
      ! Set up linear form
       rlinformconv%itermCount = 2
       rlinformconv%Idescriptors(1) = DER_DERIV_X
       rlinformconv%Idescriptors(2) = DER_DERIV_Y
       rproblem%rcollection%p_rvectorQuickAccess1 => rproblem%rvector
       ! Call the linearformroutine
       call linf_buildVectorScalar2 (rlinformconv, .false., p_rrhs%RvectorBlock(1),&
                                     fcoeff_VectorScalarMgEdge, rproblem%rcollection)
    end if


    
  end subroutine
  
  
  
  ! ***************************************************************************

!<subroutine>

  subroutine dgmgs_solve (rproblem)
  
!<description>
  ! Solves the given problem by applying a linear solver.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

  ! local variables
    integer :: ilvmin,ilvmax
    integer :: i

    ! Error indicator during initialisation of the solver
    integer :: ierror
  
    ! A filter chain to filter the vectors and the matrix during the
    ! solution process.
    type(t_filterChain), dimension(1), target :: RfilterChain

    ! A pointer to the system matrix and the RHS vector as well as
    ! the discretisation
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_vectorBlock), pointer :: p_rrhs,p_rvector
    type(t_vectorBlock), target :: rtempBlock

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode,p_rsmoother
    type(t_linsolNode), pointer :: p_rcoarseGridSolver,p_rpreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(:), pointer :: Rmatrices
    
    ! One level of multigrid
    type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo
    
    ilvmin = rproblem%ilvmin
    ilvmax = rproblem%ilvmax
    
    ! Get our right hand side / solution / matrix on the finest
    ! level from the problem structure.
    p_rrhs    => rproblem%rrhs
    p_rvector => rproblem%rvector
    p_rmatrix => rproblem%RlevelInfo(ilvmax)%rmatrix
    
    ! Create a temporary vector we need that for some preparation.
    call lsysbl_createVecBlockIndirect (p_rrhs, rtempBlock, .false.)

    ! Now we have to build up the level information for multigrid.
    !
    ! Create a Multigrid-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    call linsol_initMultigrid2 (p_rsolverNode,ilvmax-ilvmin+1,RfilterChain)
    
    ! Then set up smoothers / coarse grid solver:
    do i=ilvmin,ilvmax
      
      ! On the coarsest grid, set up a coarse grid solver and no smoother
      ! On finer grids, set up a smoother but no coarse grid solver.
      nullify(p_rpreconditioner)
      nullify(p_rsmoother)
      nullify(p_rcoarseGridSolver)
      
      ! Get the level
      call linsol_getMultigrid2Level (p_rsolverNode,1,p_rlevelInfo)
      
      if (i .eq. ilvmin) then
        ! Set up a BiCGStab solver with ILU preconditioning as coarse grid solver
        ! would be:
        ! CALL linsol_initMILUs1x1 (p_rpreconditioner,0,0.0_DP)
        ! CALL linsol_initBiCGStab (p_rcoarseGridSolver,p_rpreconditioner,RfilterChain)
        
        ! Set up UMFPACK coarse grid solver.
        call linsol_initUMFPACK4 (p_rcoarseGridSolver)
        
      else
        ! Setting up Jacobi smoother for multigrid would be:
        ! CALL linsol_initJacobi (p_rsmoother)
        
        ! Setting up VANCA smoother for multigrid would be:
        ! CALL linsol_initVANCA (p_rsmoother)

        ! Set up an ILU smoother for multigrid with damping parameter 0.7,
        ! 4 smoothing steps:
        call linsol_initMILUs1x1 (p_rsmoother,0,0.0_DP)
        call linsol_convertToSmoother (p_rsmoother,4,0.7_DP)
      end if
    
      ! And add this multi-grid level. We will use the same smoother
      ! for pre- and post-smoothing.
      call linsol_getMultigrid2Level (p_rsolverNode,i-ilvmin+1,p_rlevelInfo)
      p_rlevelInfo%p_rcoarseGridSolver => p_rcoarseGridSolver
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
    allocate(Rmatrices(ilvmin:ilvmax))
    do i=ilvmin,ilvmax
      call lsysbl_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrix,&
          Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    end do
    
    call linsol_setMatrices(p_RsolverNode,Rmatrices(ilvmin:ilvmax))
    
    ! We can release Rmatrices immediately -- as long as we do not
    ! release rproblem%RlevelInfo(i)%rmatrix!
    do i=ilvmin,ilvmax
      call lsysbl_releaseMatrix (Rmatrices(i))
    end do
    deallocate(Rmatrices)
    
    ! Initialise structure/data of the solver. This allows the
    ! solver to allocate memory / perform some precalculation
    ! to the problem.
    call linsol_initStructure (p_rsolverNode,ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) stop
    call linsol_initData (p_rsolverNode,ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) stop
    
    ! Finally solve the system. As we want to solve Ax=b with
    ! b being the real RHS and x being the real solution vector,
    ! we use linsol_solveAdaptively. If b is a defect
    ! RHS and x a defect update to be added to a solution vector,
    ! we would have to use linsol_precondDefect instead.
    call linsol_solveAdaptively (p_rsolverNode,p_rvector,p_rrhs,rtempBlock)
    
    ! Release solver data and structure
    call linsol_doneData (p_rsolverNode)
    call linsol_doneStructure (p_rsolverNode)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    call linsol_releaseSolver (p_rsolverNode)
    
    ! Release the temporary vector
    call lsysbl_releaseVector (rtempBlock)

  end subroutine
  
  
    ! ***************************************************************************

!<subroutine>

  subroutine dgmgs_postprocessing (rproblem)
  
!<description>
  ! Writes the solution into a GMV file.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

  ! local variables
  
    ! We need some more variables for postprocessing.
    real(DP), dimension(:), pointer :: p_Ddata
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir, sofile

    ! A pointer to the solution vector and to the triangulation.
    type(t_vectorBlock), pointer :: p_rvector
    type(t_triangulation), pointer :: p_rtriangulation
    
    real(DP) :: derror

    ! Get the solution vector from the problem structure
    p_rvector => rproblem%rvector
    
    write(*,*) ''
    write(*,*) 'Writing solution to file'
    
    sofile = './gmv/u2d'

    ! Output solution to gmv file
    call dg2gmv(p_rvector%Rvectorblock(1),3,sofile,-1)

    ! Output solution to vtk file
    call dg2vtk(p_rvector%Rvectorblock(1),3,sofile,-1)

  end subroutine
  
  
    ! ***************************************************************************

!<subroutine>

  subroutine dgmgs_doneMatVec (rproblem)
  
!<description>
  ! Releases system matrix and vectors.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    integer :: ihandle,i

    ! Release matrices and vectors on all levels
    do i=rproblem%ilvmax,rproblem%ilvmin,-1
      ! Delete the matrix
      call lsysbl_releaseMatrix (rproblem%RlevelInfo(i)%rmatrix)
      ! And the additional triangulation data
      call releaseTriaData(rproblem%RlevelInfo(i)%raddTriaData)
    end do

    ! Delete solution/RHS vector
    call lsysbl_releaseVector (rproblem%rvector)
    call lsysbl_releaseVector (rproblem%rrhs)

  end subroutine
  
    ! ***************************************************************************

!<subroutine>

  subroutine dgmgs_doneDiscretisation (rproblem)
  
!<description>
  ! Releases the discretisation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i

    do i=rproblem%ilvmax,rproblem%ilvmin,-1
      ! Delete the block discretisation together with the associated
      ! scalar spatial discretisations....
      call spdiscr_releaseBlockDiscr(rproblem%RlevelInfo(i)%p_rdiscretisation)
      
      ! and remove the allocated block discretisation structure from the heap.
      deallocate(rproblem%RlevelInfo(i)%p_rdiscretisation)
    end do
    
  end subroutine
    
  ! ***************************************************************************

!<subroutine>

  subroutine dgmgs_doneParamTriang (rproblem)
  
!<description>
  ! Releases the triangulation and parametrisation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i

    do i=rproblem%ilvmax,rproblem%ilvmin,-1
      ! Release the triangulation
      call tria_done (rproblem%RlevelInfo(i)%rtriangulation)
    end do
    
    deallocate(rproblem%RlevelInfo)
    
    ! Finally release the domain.
    call boundary_release (rproblem%rboundary)
    
  end subroutine
  
  



  ! ***************************************************************************

  !<subroutine>

  subroutine dg2d_mgsc

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
    ! We need a couple of variables for this problem. Let us see...
    !
    ! An object for saving the domain:
    type(t_boundary) :: rboundary

    ! Number of Variables
    ! Shallow water : 3 (h, hu, hv)
    ! (h=Waterheights, u/v=speed in x/y-direction)
    ! Euler: 4 (rho, rho u, rho v, rho E)
    integer, parameter :: nvar2d = 1

    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! Path to the mesh
    character(len=SYS_STRLEN) :: spredir

    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    type(t_blockDiscretisation) :: rdiscretisation

    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinformconv, rlinformedge, rlinformIC, rlinformSource, rlinform

    ! A scalar matrix and vector. The vector accepts the RHS of the problem
    ! in scalar form.
    type(t_matrixScalar) :: rmatrixMC, rmatrixiMC, rmatrixA
    type(t_vectorScalar) :: rrhs,rsol,redge,rconv,rsoltemp,rrhstemp,rsolUp,rsolold,rdef

    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_matrixBlock) :: rmatrixBlock, rmatrixiBlock, rmatrixABlock
    type(t_vectorBlock), target :: rvectorBlock,rrhsBlock,rtempBlock,rsolBlock,redgeBlock,rconvBlock,rsolTempBlock,rsolUpBlock,rsolOldBlock,rsolLimiterBlock,rsolSaveBlock,rsourceTermBlock
    type(t_vectorBlock), target :: rk0, rk1, rk2, rk3, rdefBlock, rimf1, rimf2


    ! A set of variables describing the discrete boundary conditions.
    type(t_boundaryRegion) :: rboundaryRegion
    type(t_discreteBC), target :: rdiscreteBC

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode, p_rpreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(1) :: Rmatrices

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    type(t_filterChain), dimension(1), target :: RfilterChain

!    ! NLMAX receives the level where we want to solve.
!    integer :: NLMAX

    ! Error indicator during initialisation of the solver
    integer :: ierror

    ! Error of FE function to reference function
    real(DP) :: derror

    integer :: ilimiting

    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir
    real(DP), dimension(:), pointer :: p_Ddata

    ! Command line and name of the paramater file
    character(LEN=SYS_STRLEN) :: cbuffer
    character(LEN=SYS_STRLEN) :: sparameterfileName

    ! Strings describing initial condition and inlet condition for the function parser
    character(len=SYS_STRLEN) :: sstring, sic, sinlet

    real(DP) :: ttime, dt, ttfinal

    real(DP), dimension(2) :: vel

    real(dp) :: dL2updnorm

    integer :: ielementType

    type(t_collection) :: rcollection

    ! Parameter list
    type(t_parlist) :: rparlist

    character(LEN=*), dimension(2), parameter ::&
         cvariables = (/ (/'x'/), (/'y'/) /)

    ! How many extra points for output
    integer :: iextraPoints

    integer :: ivar

    ! For time measurement
    real(dp) :: dtime1, dtime2

    type(t_additionalTriaData) :: raddTriaData

    integer :: ilimiter

    real(dp) :: ddefnorm

    integer :: idef

    ! The profiler to measure the time
    type(t_profiler) :: rprofiler

    integer :: itimestepping

    integer :: iwithoutlimiting

    real(dp) :: dCFL, dgravconst

    ! Name of output file
    character (LEN=SYS_STRLEN) :: sofile

    integer :: imakeVideo, ifilenumber, ioutputtype

    real(dp) :: dvideotimestep, dvideotime

    integer, dimension(6) :: IdofGlob

    integer :: iel, NEL

    integer :: irhstype

    real(dp) :: dtheta

    integer :: iinsertSourceTerm = 0

    integer :: ilimitEveryStep = 1

    real(dp) , dimension(:), pointer :: p_DiMCdata, p_DMCdata

    integer :: i, j
    
    integer :: inumSolverCalls = 0, iiterations =0
    
        ! NLMIN receives the minimal level where to discretise for supporting
    ! the solution process.
    ! NLMAX receives the level where we want to solve.
    integer :: NLMIN,NLMAX
    
    ! A problem structure for our problem
    type(t_problem), pointer :: p_rproblem
    
!    integer :: i
    
    ! Start time measurement
    call cpu_time(dtime1)
    
    ! Get command line arguments and extract name of parameter file
    if (command_argument_count() .eq. 0) then
       call output_lbrk()
       call output_line('Using standart parameterfile: ./dat/1.dat')
       call output_lbrk()
       sparameterfileName = './dat/1.dat'
    else
       call get_command_argument(command_argument_count(), cbuffer)
       sparameterfileName = adjustl(cbuffer)
    end if
    
    ! Read parameter file
    call parlst_init(rparlist)
    call parlst_readFromFile(rparlist,sparameterfileName)

    ! Allocate the problem structure
    allocate(p_rproblem)
   
    ! Initialise triangulation and parametrisation
    call dgmgs_initParamTriang (p_rproblem,rparlist)
    
    ! Initialise discretisation
    call dgmgs_initDiscretisation (p_rproblem,rparlist)
    
    ! Calculate matrix and vectors
    call dgmgs_initMatVec (p_rproblem,rparlist)
    
    ! Solve
    call dgmgs_solve (p_rproblem)
    
    ! Write solution to file
    call dgmgs_postprocessing (p_rproblem)
    
    ! Release memory
    call dgmgs_doneMatVec (p_rproblem)
    call dgmgs_doneDiscretisation (p_rproblem)
    call dgmgs_doneParamTriang (p_rproblem)

  end subroutine

end module
