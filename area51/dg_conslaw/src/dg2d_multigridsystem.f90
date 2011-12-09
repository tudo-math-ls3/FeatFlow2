!##############################################################################
!# ****************************************************************************
!# <name> dg2d_multigridsystem </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module solves the transport equation
!# 
!#   Q_t + \nabla \cdot F(Q) = 0
!#    
!# using an implicit dg discretisation with
!# multigrid solver.
!#
!# </purpose>
!##############################################################################

module dg2d_multigridsystem

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
  use multilevelprojection

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
    
    ! A vector on the level
    type(t_vectorBlock) :: rvector

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
    
    ! Number of variables (components of the system)
    integer :: nvar2d
    
    ! Variables for the time stepping (current time, final time and timestep)
    real(dp) :: dtcurrent, dtfinal, dtstep
    
    ! The type of element used for the discretisation
    integer :: ielementType

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
    
    ! Parameter list
    type(t_parlist) :: rparlist

  end type

!</typeblock>

!</types>
  

contains

  ! ***************************************************************************

!<subroutine>

  subroutine dgmgs_initParlist (rproblem)
  
!<description>
  ! This routine reads in the parameter file.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

  ! Command line and name of the paramater file
  character(LEN=SYS_STRLEN) :: cbuffer
  character(LEN=SYS_STRLEN) :: sparameterfileName


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
    call parlst_init(rproblem%rparlist)
    call parlst_readFromFile(rproblem%rparlist,sparameterfileName)
    
    ! Set timesteppingparameters
    rproblem%dtcurrent = 0.0_dp
    call parlst_getvalue_double(rproblem%rparlist, 'TIMESTEPPING', 'dt', rproblem%dtstep)
    call parlst_getvalue_double(rproblem%rparlist, 'TIMESTEPPING', 'ttfinal', rproblem%dtfinal)

  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine dgmgs_initParamTriang (rproblem)
  
!<description>
  ! This routine initialises the parametrisation and triangulation of the
  ! domain. The corresponding .prm/.tri files are read from disc and
  ! the triangulation is refined as described by the parameter ilv.
!</description>

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
  
  ! We want to solve our problem on level...
  call parlst_getvalue_int(rproblem%rparlist, 'TRIANGULATION', 'NLMIN', ilvmin)
  
  ! We want to solve our problem on level...
  call parlst_getvalue_int(rproblem%rparlist, 'TRIANGULATION', 'NLMAX', ilvmax)
  
  ! Initialise the level in the problem structure
  rproblem%ilvmin = ilvmin
  rproblem%ilvmax = ilvmax
  allocate(rproblem%RlevelInfo(ilvmin:ilvmax))
    
  ! At first, read in the parametrisation of the boundary and save
  ! it to rboundary.
  call parlst_getvalue_string (rproblem%rparlist, 'TRIANGULATION', &
                               'prmname', sstring)
  call boundary_read_prm(rproblem%rboundary, trim(sstring))
        
  ! Now read in the basic triangulation.
  call parlst_getvalue_string (rproblem%rparlist, 'TRIANGULATION', &
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

  subroutine dgmgs_initDiscretisation (rproblem)
  
!<description>
  ! This routine initialises the discretisation structure of the underlying
  ! problem and saves it to the problem structure.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

  ! local variables
  integer :: i, ielementType, nvar2d, ivar
  
    ! An object for saving the domain:
    type(t_boundary), pointer :: rboundary
    
    ! An object for saving the triangulation on the domain
    type(t_triangulation), pointer :: p_rtriangulation

    ! An object for the block discretisation on one level
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
    ! Type of finite element to use
    call parlst_getvalue_int(rproblem%rparlist, 'TRIANGULATION', 'FEkind', ielementType, 0)
    
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
    
    rproblem%ielementType = ielementType
    
    ! Get nuber of components of the system
    call parlst_getvalue_int(rproblem%rparlist, 'PROBLEM', 'nvar2d', nvar2d)
    rproblem%nvar2d = nvar2d
    
    do i=rproblem%ilvmin,rproblem%ilvmax
      ! Ask the problem structure to give us the boundary and triangulation.
      ! We need it for the discretisation.
      rboundary => rproblem%rboundary
      p_rtriangulation => rproblem%RlevelInfo(i)%rtriangulation
      
      ! Now we can start to initialise the discretisation. At first, set up
      ! a block discretisation structure that specifies the blocks in the
      ! solution vector.
      allocate(p_rdiscretisation)
      call spdiscr_initBlockDiscr (p_rdiscretisation, nvar2d,&
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
      
      ! Now copy this initialised block into the other ones
      ! But only create a derived copy, which shares the handles of the original one
      do ivar = 2, nvar2d
         call spdiscr_duplicateDiscrSc (p_rdiscretisation%Rspatialdiscr(1), &
                                        p_rdiscretisation%Rspatialdiscr(ivar), &
                                        .true.)
      end do

    end do
                                   
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine dgmgs_projectInitialConditions (rproblem)
  
!<description>
  ! Sets initial conditions for solution vector.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

  ! local variables
  ! A solver node that accepts parameters for the linear solver
  type(t_linsolNode), pointer :: p_rsolverNode, p_rpreconditioner
  
  ! A bilinear and linear form describing the analytic problem to solve
  type(t_bilinearForm) :: rform
  type(t_linearForm)   :: rlinformIC
  
  ! Pointer to a block matrix
  type(t_matrixBlock), pointer :: p_rmatrixBlock
  
  ! RHS and temp vector
  type(t_vectorBlock) :: rvectorRHS, rtempBlock
  
  ! Pointer to the discretisation
  type(t_blockDiscretisation), pointer :: p_rdiscretisation
  
  ! For the solver
  type(t_matrixBlock), dimension(1) :: Rmatrices
  
  integer :: ivar, ierror
  
  
  ! Now set the initial conditions via L2 projection
  write(*,*) 'Projecting initial condition'
  
  p_rmatrixBlock => rproblem%RlevelInfo(rproblem%ilvmax)%rmatrix
  p_rdiscretisation => rproblem%RlevelInfo(rproblem%ilvmax)%p_rdiscretisation  
  
  !!! First build the mass matrix !!!
  
  ! First create an empty block matrix structur with nvar2d x nvar2d blocks
  call lsysbl_createEmptyMatrix (p_rmatrixBlock, rproblem%nvar2d)
  
  ! Now create the mass matrix in the first component of the block matrix.
  ! First the matrix structure.
  call bilf_createMatrixStructure (p_rdiscretisation%RspatialDiscr(1),&
                                   LSYSSC_MATRIX9,p_rmatrixBlock%RmatrixBlock(1,1))

  ! And now to the entries of the matrix. For assembling of the entries,
  ! we need a bilinear form, which first has to be set up manually.
  ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
  ! scalar system matrix in 2D.

  rform%itermCount = 1
  rform%Idescriptors(1,1) = DER_FUNC
  rform%Idescriptors(2,1) = DER_FUNC

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
  call bilf_buildMatrixScalar (rform,.true.,p_rmatrixBlock%RmatrixBlock(1,1),&
                               coeff_Laplace_2D)  

  ! Now copy the mass matrix into the other diagonal blocks
  do ivar = 2, rproblem%nvar2d
     call lsyssc_duplicateMatrix (p_rmatrixBlock%RmatrixBlock(1,1), &
                                  p_rmatrixBlock%Rmatrixblock(ivar,ivar), &
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)
  end do
  
  ! Finalize matrix
  call lsysbl_updateMatStrucInfo (p_rmatrixBlock)
  
  
  
  !!! Now build RHS !!!
  
  ! Create structure of RHS, temp and solution vector
  call lsysbl_createVecBlockIndMat (p_rmatrixBlock,rproblem%rvector,.true.)
  call lsysbl_createVecBlockIndMat (p_rmatrixBlock,rvectorRHS,.true.)
  call lsysbl_createVecBlockIndMat (p_rmatrixBlock,rtempBlock,.true.)
  

  ! Set up linear form
  rlinformIC%itermCount = 1
  rlinformIC%Idescriptors(1) = DER_FUNC2D

  ! Build RHS for each component
  do ivar = 1, rproblem%nvar2d
    rproblem%rcollection%IquickAccess(1) = ivar
    !rrhsBlock%p_rblockDiscr%RspatialDiscr(ivar)%RelementDistr(1)%ccubTypeLinForm=CUB_G6_2D
    call linf_buildVectorScalar2 (rlinformIC, .true., rvectorRHS%RvectorBlock(ivar),&
                                  Euler_coeff_RHS_IC, rproblem%rcollection)
    !rrhsBlock%p_rblockDiscr%RspatialDiscr(ivar)%RelementDistr(1)%ccubTypeLinForm=CUB_G3x3
  end do
  
  
  
  !!! Create a solver !!!
  nullify(p_rpreconditioner)
  call linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner)
  !call linsol_initUMFPACK4 (p_rsolverNode)

  ! Set the output level of the solver to 2 for some output
  p_rsolverNode%ioutputLevel = 2

  ! The linear solver stops, when this relative or absolut norm of
  ! the residual is reached.
  p_rsolverNode%depsRel = 1.0e-14
  p_rsolverNode%depsAbs = 1.0e-14

  ! Attach the system matrix to the solver.
  Rmatrices = (/p_rmatrixBlock/)
  call linsol_setMatrices(p_rsolverNode,Rmatrices)

  ! Initialise structure/data of the solver. This allows the
  ! solver to allocate memory / perform some precalculation
  ! to the problem.
  call linsol_initStructure (p_rsolverNode, ierror)
  if (ierror .ne. LINSOL_ERR_NOERROR) stop
  call linsol_initData (p_rsolverNode, ierror)
  if (ierror .ne. LINSOL_ERR_NOERROR) stop
  
  
  !!! Solve for initial conditions !!!
  call linsol_solveAdaptively (p_rsolverNode,rproblem%rvector,rvectorRHS,rtempBlock)


  !!! Release solver !!!
  ! Release solver data and structure
  call linsol_doneData (p_rsolverNode)
  call linsol_doneStructure (p_rsolverNode)

  ! Release the solver node and all subnodes attached to it (if at all):
  call linsol_releaseSolver (p_rsolverNode)
  

  !!! Free memory !!!
  call lsysbl_releaseVector (rvectorRHS)
  call lsysbl_releaseVector (rtempBlock)
  call lsysbl_releaseMatrix (p_rmatrixBlock)
                                     
  end subroutine
  
  
  
  
  
  
  
  
  
  
  
  
  
  ! ***************************************************************************

!<subroutine>

  subroutine dgmgs_doTimeStepping (rproblem)
  
!<description>
  ! Calculates the system matrix and RHS vector for the multigrid solver
  ! and performs timesteps until the final time is reached.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

  ! local variables
  integer :: ilvmin, ilvmax
  integer :: ilevel, nlmax
  integer :: nvar2d
  integer :: i,j
  integer :: ierror
  
  ! A solver node that accepts parameters for the linear solver
  type(t_linsolNode), pointer :: p_rsolverNode
  type(t_linsolNode), pointer :: p_rcoarseGridSolver, p_rpreconditioner, p_rsmoother
  
  ! Pointer to the block discretisations
  type(t_blockDiscretisation), pointer :: p_rdiscrCoarse,p_rdiscrFine
  
  ! Temp vector
  type(t_vectorScalar) :: rvectorTemp1
  
  ! Mass matrices for all levels
  type(t_matrixScalar), dimension(:), allocatable :: rmatrixMC
  
  ! Bilinear form
  type(t_bilinearForm) :: rform
  
  ! Linear forms
  type(t_linearForm) :: rlinformedge, rlinformcell
  
  ! A pointer to the system matrix and the RHS vector as well as
  ! the discretisation
  type(t_matrixBlock), pointer :: p_rmatrix
  type(t_vectorBlock), pointer :: p_rrhs,p_rvector
  type(t_vectorBlock), target :: rtempBlock
  
  ! An array for the system matrix(matrices) during the initialisation of
  ! the linear solver.
  type(t_matrixBlock), dimension(:), pointer :: Rmatrices
  
  ! One level of multigrid
  type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo
  


  ilvmin = rproblem%ilvmin
  ilvmax = rproblem%ilvmax
  nvar2d = rproblem%nvar2d
  
  
    
  !!! Create a Multigrid-solver !!!
  call linsol_initMultigrid2 (p_rsolverNode,ilvmax-ilvmin+1)
  
  
  
  !!! Now initialise the projections !!!
  nlmax = p_rsolverNode%p_rsubnodeMultigrid2%nlevels
  
  ! Loop through all levels. Whereever the projection structure
  ! is missing, create it
  do ilevel = 2,nlmax
    ! Initialise the projection structure
    call linsol_initProjMG2LvByDiscr (p_rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel),&
                                      rproblem%RlevelInfo(ilvmin+ilevel-2)%p_rdiscretisation,&
                                      rproblem%RlevelInfo(ilvmin+ilevel-1)%p_rdiscretisation)
  end do
  
  
  
  !!! Project the initial solution to all levels !!!
  do ilevel = ilvmax, ilvmin, -1
    
    ! First copy initial solution or project the solution to lower level
    if (ilevel == ilvmax) then
      ! We're on the highest level - just copy the initial data
      call lsysbl_copyVector(rproblem%rvector,rproblem%RlevelInfo(ilevel)%rvector)
    else
      ! We're on a lower level and have to project the solution from the higher level
      call lsysbl_createVecBlockByDiscr (rproblem%RlevelInfo(ilevel)%p_rdiscretisation, rproblem%RlevelInfo(ilevel)%rvector)
      call lsyssc_createVecByDiscr (rproblem%RlevelInfo(ilevel)%p_rdiscretisation%RspatialDiscr(1),rvectorTemp1)
      call mlprj_performInterpolation (p_rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel-ilvmin+2)%p_rprojection,&
                                       rproblem%RlevelInfo(ilevel)%rvector, &
                                       rproblem%RlevelInfo(ilevel+1)%rvector, &
                                       rvectorTemp1)
      call lsyssc_releaseVector(rvectorTemp1)
    end if
  end do
  
  
  
  !!! Build mass matrices on all levels !!!
  
  ! Allocate the matrices
  allocate(rmatrixMC(ilvmin:ilvmax))
  
  ! Set up the bilinearform
  rform%itermCount = 1
  rform%Idescriptors(1,1) = DER_FUNC
  rform%Idescriptors(2,1) = DER_FUNC
  rform%ballCoeffConstant = .true.
  rform%BconstantCoeff = .true.
  
  do ilevel = ilvmin, ilvmax
    ! Create matrix structure
    call bilf_createMatrixStructure (rproblem%RlevelInfo(ilevel)%p_rdiscretisation%RspatialDiscr(1),&
                                     LSYSSC_MATRIX9,&
                                     rmatrixMC(ilevel),&
                                     rproblem%RlevelInfo(ilevel)%p_rdiscretisation%RspatialDiscr(1),&
                                     BILF_MATC_EDGEBASED)
    ! Calculate entries
    call bilf_buildMatrixScalar (rform,.true.,rmatrixMC(ilevel))
  end do
  
  
  
  !!! Build the RHS on the finest level !!!
  ! Create structure of RHS vector (by copying the solution vector into it)
  call lsysbl_createVecBlockByDiscr (rproblem%RlevelInfo(ilvmax)%p_rdiscretisation, rproblem%rrhs)

  ! Set up the corresponding linear form (f,Phi_j):
  rlinformedge%itermCount = 1
  rlinformedge%Idescriptors(1) = DER_FUNC2D
   
  ! Now use the dg-function for the edge terms
  call linf_dg_buildVectorBlockEdge2d (rlinformedge, CUB_G5_1D, .true.,&
       rproblem%rrhs,rproblem%rvector,&
       rproblem%RlevelInfo(ilvmax)%raddTriaData,&
       Euler_flux_dg_buildVectorBlEdge2D_sim_Newton,&
       rproblem%rcollection)

  call lsysbl_scaleVector (rproblem%rrhs,-1.0_DP)


  ! Then add the cell term
  if(rproblem%ielementType .ne. EL_DG_T0_2D) then
     ! Set up linear form
     rlinformcell%itermCount = 2
     rlinformcell%Idescriptors(1) = DER_DERIV_X
     rlinformcell%Idescriptors(2) = DER_DERIV_Y
     rproblem%rcollection%p_rvectorQuickAccess1 => rproblem%rvector
     ! Call the linearformroutine
     call linf_buildVectorBlock2 (rlinformcell, .false., rproblem%rrhs,&
                                  Euler_flux_sys_block,rproblem%rcollection)
  end if
  
  

  
  !!! Now build the nonlinear matrices on all levels !!!
  do ilevel = ilvmin, ilvmax
    
    ! Create structure of the block matrix by copying MC to all blocks
    call lsysbl_createEmptyMatrix (rproblem%RlevelInfo(ilevel)%rmatrix,nvar2d,nvar2d)
    do i = 1, nvar2d
       do j = 1, nvar2d
          call lsyssc_copyMatrix (rmatrixMC(ilevel),rproblem%RlevelInfo(ilevel)%rmatrix%RmatrixBlock(i,j))
       end do
    end do
    call lsysbl_updateMatStrucInfo (rproblem%RlevelInfo(ilevel)%rmatrix)
    call lsysbl_clearMatrix(rproblem%RlevelInfo(ilevel)%rmatrix)
    
    
    ! Calculate the matrix for the cell terms
    ! Set up the bilinear form
    rform%itermCount = 2
    rform%Idescriptors(1,1) = DER_FUNC
    rform%Idescriptors(2,1) = DER_DERIV_X
    rform%Idescriptors(1,2) = DER_FUNC
    rform%Idescriptors(2,2) = DER_DERIV_Y
    rform%ballCoeffConstant = .false.
    rform%BconstantCoeff = .false.
    rproblem%rcollection%p_rvectorQuickAccess1 => rproblem%RlevelInfo(ilevel)%rvector
    if (rproblem%ielementtype.ne.EL_DG_T0_2D) then
       call bilf_buildMatrixBlock2 (rform, .true., rproblem%RlevelInfo(ilevel)%rmatrix,&
                                    fcoeff_buildMatrixBl_sim_iEuler,rproblem%rcollection)
    else
       call lsysbl_clearMatrix(rproblem%RlevelInfo(ilevel)%rmatrix)
    end if

    call lsysbl_scaleMatrix (rproblem%RlevelInfo(ilevel)%rmatrix,-1.0_DP)


    ! Calculate the matrix for the edge terms
    ! Set up the bilinear form
    rform%itermCount = 2
    rform%Idescriptors(1,1) = DER_FUNC
    rform%Idescriptors(2,1) = DER_FUNC
    rform%Idescriptors(1,2) = DER_FUNC
    rform%Idescriptors(2,2) = DER_FUNC
    rform%ballCoeffConstant = .false.
    rform%BconstantCoeff = .false.
    rproblem%rcollection%p_rvectorQuickAccess1 => rproblem%RlevelInfo(ilevel)%rvector
    rproblem%rcollection%Dquickaccess(1) = rproblem%dtstep
    call bilf_dg_buildMatrixBlEdge2D_ss (rform, CUB_G5_1D, .false.,&
                                         rproblem%RlevelInfo(ilevel)%rmatrix,&
                                         rproblem%RlevelInfo(ilevel)%rvector,&
                                         rproblem%RlevelInfo(ilevel)%raddTriaData,&
                                         flux_dg_buildMatrixBlEdge2D_sim_iEuler_Newton_ss,&
                                         rproblem%rcollection)
    
    ! Add the mass matrix on the main diagonal
    do i = 1, rproblem%nvar2d
      call lsyssc_matrixLinearComb (rmatrixMC(ilevel),rproblem%RlevelInfo(ilevel)%rmatrix%RmatrixBlock(i,i),&
                                    1.0_dp/rproblem%dtstep,1.0_dp,&
                                    .false.,.false.,.true.,.false.)
    end do
    
  end do
  
  
  
  !!! Create (rest of) multigrid solver and solve !!!
  ! Get our right hand side / solution / matrix on the finest
  ! level from the problem structure.
  p_rrhs    => rproblem%rrhs
  p_rvector => rproblem%rvector
  p_rmatrix => rproblem%RlevelInfo(ilvmax)%rmatrix
  
  ! Create a temporary vector we need that for some preparation.
  call lsysbl_createVecBlockIndirect (p_rrhs, rtempBlock, .false.)
  
  ! Then set up smoothers / coarse grid solver:
  do i=ilvmin,ilvmax
    
    ! On the coarsest grid, set up a coarse grid solver and no smoother
    ! On finer grids, set up a smoother but no coarse grid solver.
    nullify(p_rpreconditioner)
    nullify(p_rsmoother)
    nullify(p_rcoarseGridSolver)
    
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
      !call linsol_initMILUs1x1 (p_rsmoother,0,0.0_DP)
      !call linsol_convertToSmoother (p_rsmoother,4,1.0_DP)
      
      
        nullify(p_rpreconditioner)
!        call linsol_initJacobi (p_rpreconditioner)
!        call linsol_initMILUs1x1 (p_rpreconditioner,0,0.0_DP)
      
!        call linsol_initJacobi (p_rsmoother)
!        call linsol_initGMRES (p_rsmoother,4,p_rpreconditioner)
        call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner)!,RfilterChain)
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
                                 Rmatrices(i),LSYSSC_DUP_SHARE,&
                                 LSYSSC_DUP_SHARE)
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
  
  
  
  ! Calculate defect and solve for solution update NOW
  
  
  
  
  !!! Free memory !!!

  ! Release the (solution) vectors on all levels
  do ilevel = ilvmin, ilvmax
    call lsysbl_releaseVector(rproblem%RlevelInfo(ilevel)%rvector)
  end do
  
  ! Release RHS vector
  call lsysbl_releaseVector (rproblem%rrhs)
  
  ! Release the solver node and all subnodes attached to it (if at all):
  call linsol_releaseSolver (p_rsolverNode)
  
  ! Deallocate mass matrices
  do ilevel = ilvmin, ilvmax
    call lsyssc_releaseMatrix(rmatrixMC(ilevel))
  end do
  deallocate(rmatrixMC)
  
  
  end subroutine
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ! ***************************************************************************

!<subroutine>

  subroutine dgmgs_initMatVec (rproblem)
  
!<description>
  ! Calculates the system matrix and RHS vector of the linear system
  ! by discretising the problem with the default discretisation structure
  ! in the problem structure.
  ! Sets up a solution vector for the linear system.
!</description>

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
  
      
    do i=rproblem%ilvmin,rproblem%ilvmax,-1
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
       call parlst_getvalue_int(rproblem%rparlist, 'TRIANGULATION', 'FEkind', ielementType)
       
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
    call parlst_getvalue_int(rproblem%rparlist, 'TRIANGULATION', 'FEkind', ielementType)

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
    
    ! Später wieder löschen
    type(t_vectorBlock) :: rvectorCoarse, rrhsCoarse, rtempCoarse
    character(LEN=SYS_STRLEN) :: sofile
    type(t_linearForm) :: rlinformIC
    type(t_bilinearForm) :: rform
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    integer :: ivar
    type(t_linsolNode), pointer :: p_rsolverNode1
    
    

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
        call linsol_convertToSmoother (p_rsmoother,4,1.0_DP)
        
        
!        nullify(p_rpreconditioner)
!        call linsol_initJacobi (p_rpreconditioner)
!        call linsol_initMILUs1x1 (p_rpreconditioner,0,0.0_DP)
        
!        call linsol_initJacobi (p_rsmoother)
!        call linsol_initGMRES (p_rsmoother,4,p_rpreconditioner)
!        call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner)!,RfilterChain)
!        call linsol_convertToSmoother (p_rsmoother,4,0.7_DP)
        
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
    
    
    
    
   
!    ! Solve with MILU
!    
!          call linsol_doneData (p_rsolverNode)
!          call linsol_doneStructure (p_rsolverNode)
!
!          ! Release the solver node and all subnodes attached to it (if at all):
!          call linsol_releaseSolver (p_rsolverNode)
!    
!    
!!       nullify(p_rpreconditioner)
!       call linsol_initMILUs1x1 (p_rpreconditioner,0,0.0_DP)
!       CALL linsol_initDefCorr (p_rsolverNode,p_rpreconditioner)
!       p_rsolverNode%domega = 1.0
!       p_rsolverNode%nmaxIterations = 5000
!       
!       
!       allocate(Rmatrices(1))
!    
!    i=rproblem%ilvmax
!    
!      call lsysbl_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrix,&
!          Rmatrices(1),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!    p_rsolverNode%ioutputLevel = 2
!       call linsol_setMatrices(p_rsolverNode,Rmatrices)
!    
!       call linsol_initStructure (p_rsolverNode,ierror)
!       if (ierror .ne. LINSOL_ERR_NOERROR) stop
!       call linsol_initData (p_rsolverNode,ierror)
!       if (ierror .ne. LINSOL_ERR_NOERROR) stop
!       
!       
!       call linsol_solveAdaptively (p_rsolverNode,p_rvector,p_rrhs,rtempBlock)
    
    
    
    
    








!    ! Test Projection of RHS (are the bounday conditions projected rightly)
!    
!
!      i=rproblem%ilvmax-1
!      ! Ask the problem structure to give us the discretisation structure
!      p_rdiscretisation => rproblem%RlevelInfo(i)%p_rdiscretisation
!         
!      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
!
!       nullify(p_rpreconditioner)
!       call linsol_initUMFPACK4 (p_rsolverNode1)
!       
!       p_rsolverNode1%ioutputLevel = 2
!       p_rsolverNode1%nmaxIterations = 5000
!       
!       allocate(Rmatrices(1))
!    
!      call lsysbl_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrix,&
!          Rmatrices(1),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!    
!       call linsol_setMatrices(p_rsolverNode1,Rmatrices)
!    
!    call linsol_initStructure (p_rsolverNode1,ierror)
!    if (ierror .ne. LINSOL_ERR_NOERROR) stop
!    call linsol_initData (p_rsolverNode1,ierror)
!    if (ierror .ne. LINSOL_ERR_NOERROR) stop
!
!    ! Create temporary empty solution vector which is needed to build the matrices
!    call lsysbl_createVecBlockIndMat (rproblem%RlevelInfo(i)%rmatrix,rvectorCoarse,.true.)
!    call lsysbl_createVecBlockIndMat (rproblem%RlevelInfo(i)%rmatrix,rrhsCoarse,.true.)
!    call lsysbl_createVecBlockIndMat (rproblem%RlevelInfo(i)%rmatrix,rtempCoarse,.true.)
!    
!!    call Test_mlprj_performInterpolation (rvectorCoarse, &
!!                                         p_rvector)
!
!!    call mlprj_initProjectionDiscr (rlevelInfo%p_rprojection,rdiscrCoarse)
!    
!    call mlprj_performRestriction (p_rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(2)%p_rprojection,rrhsCoarse, &
!                                       rproblem%rrhs,rtempBlock%rvectorblock(1),1)
!
!    call linsol_solveAdaptively (p_rsolverNode1,rvectorCoarse,rrhsCoarse,rtempCoarse)
!    
!    sofile = './gmv/testrestr' 
!    
!    ! Output solution to vtk file
!    call dg2vtk(rvectorCoarse%Rvectorblock(1),3,sofile,-1)









    
    
    
    
    
    
    
!    ! Test Projections
!    
!    call linsol_doneData (p_rsolverNode)
!          call linsol_doneStructure (p_rsolverNode)
!
!          ! Release the solver node and all subnodes attached to it (if at all):
!          call linsol_releaseSolver (p_rsolverNode)
!
!      i=rproblem%ilvmax
!      ! Ask the problem structure to give us the discretisation structure
!      p_rdiscretisation => rproblem%RlevelInfo(i)%p_rdiscretisation
!         
!      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
!      
!      ! Initialise the block matrix with default values based on
!      ! the discretisation.
!      call lsysbl_createMatBlockByDiscr (p_rdiscretisation,p_rmatrix)    
!
!      ! Now as the discretisation is set up, we can start to generate
!      ! the structure of the system matrix which is to solve.
!      ! We create that directly in the block (1,1) of the block matrix
!      ! using the discretisation structure of the first block.
!      call bilf_createMatrixStructure (&
!                p_rdiscretisation%RspatialDiscr(1),LSYSSC_MATRIX9,&
!                p_rmatrix%RmatrixBlock(1,1),&
!                p_rdiscretisation%RspatialDiscr(1),&
!                BILF_MATC_EDGEBASED)
!                
!       ! Build mass matrix
!       rform%itermCount = 1
!       rform%Idescriptors(1,1) = DER_FUNC
!       rform%Idescriptors(2,1) = DER_FUNC
!       rform%ballCoeffConstant = .true.
!       rform%BconstantCoeff = .true.
!       rform%Dcoefficients(1)  = 1.0
!       call bilf_buildMatrixScalar (rform,.true.,p_rmatrix%RmatrixBlock(1,1))
!    
!    
!    
!    
!    
!       ! Now set the initial conditions via L2 projection
!       rlinformIC%itermCount = 1
!       rlinformIC%Idescriptors(1) = DER_FUNC2D
!       !rcollection%SquickAccess(2) = cvariables
!       !rcollection%SquickAccess(1) = sic
!
!     write(*,*) 'Projecting initial condition'
!
!       do ivar = 1, 1
!
!!          rcollection%IquickAccess(1) = ivar
!
!          !rrhsBlock%p_rblockDiscr%RspatialDiscr(ivar)%RelementDistr(1)%ccubTypeLinForm=CUB_G6_2D
!          call linf_buildVectorScalar2 (rlinformIC, .true., p_rrhs%RvectorBlock(ivar),&
!               Euler_coeff_RHS_IC)!, rcollection)
!          !rrhsBlock%p_rblockDiscr%RspatialDiscr(ivar)%RelementDistr(1)%ccubTypeLinForm=CUB_G3x3
!
!       end do
!       
!       
!       nullify(p_rpreconditioner)
!       call linsol_initBiCGStab (p_rsolvernode,p_rpreconditioner)
!       
!       allocate(Rmatrices(1))
!    
!      call lsysbl_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrix,&
!          Rmatrices(1),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!    
!       call linsol_setMatrices(p_rsolverNode,Rmatrices)
!    
!    call linsol_initStructure (p_rsolverNode,ierror)
!    if (ierror .ne. LINSOL_ERR_NOERROR) stop
!    call linsol_initData (p_rsolverNode,ierror)
!    if (ierror .ne. LINSOL_ERR_NOERROR) stop
!       
!       
!       call linsol_solveAdaptively (p_rsolverNode,p_rvector,p_rrhs,rtempBlock)
!       !call linsol_solveAdaptively (p_rsolverNode,rsolBlock,rrhsBlock,rtempBlock)  
!    
!    
!    
!    
!    
!    
!    
!    
!    
!    ! Create temporary empty solution vector which is needed to build the matrices
!    call lsysbl_createVecBlockIndMat (rproblem%RlevelInfo(ilvmax-1)%rmatrix,rvectorCoarse,.true.)
!    
!!    call Test_mlprj_performInterpolation (rvectorCoarse, &
!!                                         p_rvector)
!    
!    sofile = './gmv/fine' 
!    
!    ! Output solution to vtk file
!    call dg2vtk(p_rvector%Rvectorblock(1),3,sofile,-1)
!    
!    
!    call Test_mlprj_performRestriction (rvectorCoarse, &
!                                         p_rvector)
!    
!    sofile = './gmv/coarse' 
!    
!    ! Output solution to vtk file
!    call dg2vtk(rvectorCoarse%Rvectorblock(1),3,sofile,-1)
!    
!    
!!    
!!    call Test_mlprj_performProlongation (rvectorCoarse, &
!!                                        p_rvector)
!!                                        
!!    sofile = './gmv/againfine' 
!!    
!!    ! Output solution to vtk file
!!    call dg2vtk(p_rvector%Rvectorblock(1),3,sofile,-1)
!    
!    
!    pause
    
    
    
    
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

  subroutine dgmgs_doneParlist (rproblem)
  
!<description>
  ! Releases the parameter list.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

  call parlst_done (rproblem%rparlist)

!</subroutine>

      
  end subroutine



  ! ***************************************************************************

  !<subroutine>

  subroutine dg2d_mgsys

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

    ! Number of Variables
    ! Shallow water : 3 (h, hu, hv)
    ! (h=Waterheights, u/v=speed in x/y-direction)
    ! Euler: 4 (rho, rho u, rho v, rho E)
    integer, parameter :: nvar2d = 1

    ! For time measurement
    real(dp) :: dtime1, dtime2
    
    ! A problem structure for our problem
    type(t_problem), pointer :: p_rproblem
    
!    integer :: i
    
    ! Start time measurement
    call cpu_time(dtime1)

    ! Allocate the problem structure
    allocate(p_rproblem)
    
    ! Initialise triangulation and parametrisation
    call dgmgs_initParlist (p_rproblem)
   
    ! Initialise triangulation and parametrisation
    call dgmgs_initParamTriang (p_rproblem)
    
    ! Initialise discretisation
    call dgmgs_initDiscretisation (p_rproblem)
    
    ! Set initial conditions via L2 projection
    call dgmgs_projectInitialConditions (p_rproblem)
    
    ! Perform the timestepping
    call dgmgs_doTimeStepping (p_rproblem)
    
!    ! Calculate matrix and vectors
!    call dgmgs_initMatVec (p_rproblem)
!    
!    ! Solve
!    call cpu_time(dtime1)
!    call dgmgs_solve (p_rproblem)
!    call cpu_time(dtime2)
!    write(*,*) 'Solver took: ', dtime2-dtime1
    
    ! Write solution to file
    call dgmgs_postprocessing (p_rproblem)
    
    ! Release memory
    call dgmgs_doneMatVec (p_rproblem)
    call dgmgs_doneDiscretisation (p_rproblem)
    call dgmgs_doneParamTriang (p_rproblem)
    call dgmgs_doneParList (p_rproblem)
    
    ! Deallocate the problem structure
    deallocate(p_rproblem)

  end subroutine

end module
