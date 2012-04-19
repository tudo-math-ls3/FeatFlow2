!##############################################################################
!# ****************************************************************************
!# <name> stokes2d_vortex_noslip </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains a multigrid solver for the 2D steady Stokes equation
!# with no-slip boundary conditions and pressure integral-mean filtering
!# on a simple square domain.
!#
!# The analytical solution of the problem solved in this module is given by
!#   u_1(x,y) =  (1 - cos(2*pi*x)) * sin(2*pi*y)
!#   u_2(x,y) = -(1 - cos(2*pi*y)) * sin(2*pi*x)
!#     p(x,y) = c * (4/pi - sin(pi*x) - sin(pi*y))
!#
!# </purpose>
!##############################################################################

module stokes2d_vortex_noslip

  use fsystem
  use storage
  use genoutput
  use linearsolver
  use boundary
  use cubature
  use derivatives
  use matrixfilters
  use vectorfilters
  use discretebc
  use bcassembly
  use triangulation
  use domainintegration
  use element
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use coarsegridcorrection
  use multileveloperators
  use multilevelprojection
  use spdiscprojection
  use filtersupport
  use scalarpde
  use linearformevaluation
  use discretebc
  use ucd
  use collection, only: t_collection
  use pprocerror
  use bilinearformevaluation
  use stdoperators
  
  use stokes2d_aux

  implicit none

  private
  public :: stokes2d_vtx_noslip

!<types>

!<typeblock description="Type block defining all information about one level">

  type t_level
  
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtria

    ! An object specifying the discretisation (structure of the
    ! solution, trial/test functions,...)
    type(t_blockDiscretisation) :: rdiscr
    
    ! Cubature info structure which encapsules the cubature formula
    type(t_scalarCubatureInfo) :: rcubInfo
    
    ! A system matrix for that specific level.
    type(t_matrixBlock) :: rmatSys

    ! Prolongation matrix for velocity
    type(t_matrixScalar) :: rmatProlVelocity

    ! Prolongation matrix for pressure
    type(t_matrixScalar) :: rmatProlPressure

    ! Multilevel projection structure
    type(t_interlevelProjectionBlock) :: rprojection

    ! A variable describing the discrete boundary conditions.
    type(t_discreteBC) :: rdiscreteBC
  
  end type
  
!</typeblock>

!</types>

contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine stokes2d_vtx_noslip
  
!<description>
  ! This is the main routine of the 2D steady Stokes solver with no-slip boundary
  ! conditions and pressure integral-mean filtering.
!</description>

!</subroutine>

  ! Definitions of variables.
  type(t_boundary) :: rboundary
  type(t_level), dimension(:), pointer :: Rlevels
  type(t_vectorBlock) :: rvecSol, rvecRhs, rvecTmp
  type(t_boundaryRegion) :: rrgn
  type(t_linearForm) :: rform
  type(t_linsolNode), pointer :: p_rsolverNode, p_rcoarseGridSolver,p_rsmoother
  type(t_matrixBlock), dimension(:), pointer :: Rmatrices
  type(t_filterChain), dimension(2), target :: RfilterChain
  type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo
  integer :: NLMIN, NLMAX, ierror, i
  real(DP) :: dnu, dc, dL2, dH1, ddiv
  type(t_ucdExport) :: rexport
  character(len=SYS_STRLEN) :: spredir, sucddir
  real(DP), dimension(:), pointer :: p_Du, p_Dv, p_Dp
  integer(I32) :: ccubature, celemVelocity, celemPressure
  type(t_collection) :: rcollection
  type(t_errorScVec) :: rerrorU, rerrorP
  real(DP), dimension(2), target :: DerrorUL2, DerrorUH1
  real(DP), dimension(1), target :: DerrorPL2

    ! We want to solve our Stokes problem on level...
    NLMIN = 2
    NLMAX = 5
    
    ! Viscosity parameter:
    dnu = 1.0_DP

    ! Pressure multiplier
    dc = 1.0_DP

    ! FE spaces
    celemVelocity = EL_Q2_2D
    celemPressure = EL_QP1_2D

    ! cubature rule
    ccubature = CUB_G5_2D

    ! Allocate memory for all levels
    allocate(Rlevels(NLMIN:NLMAX))

    ! Read the parametrisation
    call output_line('Reading boundary parametrisation...')
    if (.not. sys_getenv_string("PREDIR", spredir)) spredir = './pre'
    call boundary_read_prm(rboundary, trim(spredir)//'/QUAD.prm')

    ! Now read in the basic triangulation.
    call output_line('Reading coarse mesh...')
    call tria_readTriFile2D (Rlevels(NLMIN)%rtria, trim(spredir)//'/QUAD.tri', rboundary)

    ! Refine the mesh up to the minimum level
    call output_line('Creating mesh hierarchy...')
    call tria_quickRefine2LevelOrdering (NLMIN-1, Rlevels(NLMIN)%rtria, rboundary)
    call tria_initStandardMeshFromRaw (Rlevels(NLMIN)%rtria, rboundary)
    
    ! Now refine the grid for the fine levels.
    do i = NLMIN+1, NLMAX
      call tria_refine2LevelOrdering(Rlevels(i-1)%rtria, Rlevels(i)%rtria, rboundary)
      call tria_initStandardMeshFromRaw(Rlevels(i)%rtria, rboundary)
    end do

    ! Loop over all levels
    call output_line('Discretising Stokes operators...')
    do i = NLMIN, NLMAX

      call output_line('Discretising level ' // trim(sys_sil(i,4)) // '...')

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Discretisation Setup
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Create the block discretisation
      call spdiscr_initBlockDiscr (Rlevels(i)%rdiscr, 3, Rlevels(i)%rtria, rboundary)

      ! Set up velocity spaces
      call spdiscr_initDiscr_simple (Rlevels(i)%rdiscr%RspatialDiscr(1), &
          celemVelocity, ccubature, Rlevels(i)%rtria, rboundary)

      call spdiscr_duplicateDiscrSc(&
          Rlevels(i)%rdiscr%RspatialDiscr(1),&
          Rlevels(i)%rdiscr%RspatialDiscr(2))

      ! Set up pressure space
      call spdiscr_deriveSimpleDiscrSc (Rlevels(i)%rdiscr%RspatialDiscr(1), &
          celemPressure, ccubature, Rlevels(i)%rdiscr%RspatialDiscr(3))
    
      ! Set up cubature info structure
      call spdiscr_createDefCubStructure(Rlevels(i)%rdiscr%RspatialDiscr(1), &
          Rlevels(i)%rcubInfo, ccubature)

      ! Create block matrix
      call lsysbl_createMatBlockByDiscr (Rlevels(i)%rdiscr, Rlevels(i)%rmatSys)

      ! Houston, we have a saddle point problem...
      Rlevels(i)%rmatSys%imatrixSpec = LSYSBS_MSPEC_SADDLEPOINT

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Structural Matrix Assembly
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Assemble A-matrix structure
      call bilf_createMatrixStructure(Rlevels(i)%rdiscr%RspatialDiscr(1),&
          LSYSSC_MATRIX9, Rlevels(i)%rmatSys%RmatrixBlock(1,1))

      call lsyssc_duplicateMatrix (Rlevels(i)%rmatSys%RmatrixBlock(1,1),&
          Rlevels(i)%rmatSys%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)

      call lsyssc_assignDiscrDirectMat (Rlevels(i)%rmatSys%RmatrixBlock(2,2),&
          Rlevels(i)%rdiscr%RspatialDiscr(2))

      ! Assemble B-matrix structure
      call bilf_createMatrixStructure (Rlevels(i)%rdiscr%RspatialDiscr(3),&
          LSYSSC_MATRIX9, Rlevels(i)%rmatSys%RmatrixBlock(1,3),&
          Rlevels(i)%rdiscr%RspatialDiscr(1))

      call lsyssc_duplicateMatrix (Rlevels(i)%rmatSys%RmatrixBlock(1,3),&
          Rlevels(i)%rmatSys%RmatrixBlock(2,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)

      ! Transpose B-matrix to obtain D-matrix structure
      call lsyssc_transposeMatrix (Rlevels(i)%rmatSys%RmatrixBlock(1,3), &
          Rlevels(i)%rmatSys%RmatrixBlock(3,1), LSYSSC_TR_STRUCTURE)

      call lsyssc_duplicateMatrix (Rlevels(i)%rmatSys%RmatrixBlock(3,1),&
          Rlevels(i)%rmatSys%RmatrixBlock(3,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Matrix Content Assembly
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Assemble laplace matrix for X-velocity
      call stdop_assembleLaplaceMatrix (Rlevels(i)%rmatSys%RmatrixBlock(1,1),&
          .true., dnu, Rlevels(i)%rcubInfo)

      ! Copy the matrix into the Y-velocity block
      call lsyssc_duplicateMatrix (Rlevels(i)%rmatSys%RmatrixBlock(1,1),&
          Rlevels(i)%rmatSys%RmatrixBlock(2,2),LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPY)

      ! Assemble the B-matrices
      call stdop_assembleSimpleMatrix (Rlevels(i)%rmatSys%RmatrixBlock(1,3),&
        DER_FUNC, DER_DERIV_X, -1.0_DP, .true., Rlevels(i)%rcubInfo)
      call stdop_assembleSimpleMatrix (Rlevels(i)%rmatSys%RmatrixBlock(2,3),&
        DER_FUNC, DER_DERIV_Y, -1.0_DP, .true., Rlevels(i)%rcubInfo)

      ! And transpose the B-matrices to obtain the D-matrices
      call lsyssc_transposeMatrix (Rlevels(i)%rmatSys%RmatrixBlock(1,3), &
          Rlevels(i)%rmatSys%RmatrixBlock(3,1), LSYSSC_TR_CONTENT)
      call lsyssc_transposeMatrix (Rlevels(i)%rmatSys%RmatrixBlock(2,3), &
          Rlevels(i)%rmatSys%RmatrixBlock(3,2), LSYSSC_TR_CONTENT)

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Boundary Condition Assembly
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Initialise the discrete BC structure
      call bcasm_initDiscreteBC(Rlevels(i)%rdiscreteBC)

      ! Bottom edge: (u1,u2) = 0
      call boundary_createRegion(rboundary,1,1,rrgn)
      rrgn%iproperties = BDR_PROP_WITHSTART
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscr,1,rrgn,&
          Rlevels(i)%rdiscreteBC, funcZeroBC2D)
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscr,2,rrgn,&
          Rlevels(i)%rdiscreteBC, funcZeroBC2D)

      ! Right edge: (u1,u2) = 0
      call boundary_createRegion(rboundary,1,2,rrgn)
      rrgn%iproperties = BDR_PROP_WITHSTART
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscr,1,rrgn,&
          Rlevels(i)%rdiscreteBC, funcZeroBC2D)
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscr,2,rrgn,&
          Rlevels(i)%rdiscreteBC, funcZeroBC2D)

      ! Top edge: (u1,u2) = 0
      call boundary_createRegion(rboundary,1,3,rrgn)
      rrgn%iproperties = BDR_PROP_WITHSTART
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscr,1,rrgn,&
          Rlevels(i)%rdiscreteBC, funcZeroBC2D)
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscr,2,rrgn,&
          Rlevels(i)%rdiscreteBC, funcZeroBC2D)

      ! Left edge: (u1,u2) = 0
      call boundary_createRegion(rboundary,1,4,rrgn)
      rrgn%iproperties = BDR_PROP_WITHSTART
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscr,1,rrgn,&
          Rlevels(i)%rdiscreteBC, funcZeroBC2D)
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscr,2,rrgn,&
          Rlevels(i)%rdiscreteBC, funcZeroBC2D)

      ! Assign BCs to system matrix
      Rlevels(i)%rmatSys%p_rdiscreteBC => Rlevels(i)%rdiscreteBC

      ! Filter system matrix
      call matfil_discreteBC(Rlevels(i)%rmatSys)

    end do

    call output_line('Discretising right-hand-side...')

    ! Create three vectors
    call lsysbl_createVecBlockIndMat (Rlevels(NLMAX)%rmatSys, rvecSol, .true.)
    call lsysbl_createVecBlockIndMat (Rlevels(NLMAX)%rmatSys, rvecRhs, .true.)
    call lsysbl_createVecBlockIndMat (Rlevels(NLMAX)%rmatSys, rvecTmp, .true.)

    ! Assemble the RHS vector
    rform%itermCount = 1
    rform%Idescriptors(1) = DER_FUNC
    rcollection%DquickAccess(1) = dnu
    rcollection%DquickAccess(2) = dc
    call linf_buildVectorScalar (rform,.true., rvecRhs%RvectorBlock(1), &
        Rlevels(NLMAX)%rcubInfo, funcRhsX2D, rcollection)
    call linf_buildVectorScalar (rform,.true., rvecRhs%RvectorBlock(2), &
        Rlevels(NLMAX)%rcubInfo, funcRhsY2D, rcollection)

    ! Assign BCs to vectors
    call lsysbl_assignDiscreteBC(rvecSol, Rlevels(NLMAX)%rdiscreteBC)
    call lsysbl_assignDiscreteBC(rvecRhs, Rlevels(NLMAX)%rdiscreteBC)
    call lsysbl_assignDiscreteBC(rvecTmp, Rlevels(NLMAX)%rdiscreteBC)

    ! Filter solution and rhs vectors
    call vecfil_discreteBCsol(rvecSol)
    call vecfil_discreteBCrhs(rvecRhs)

    ! Initialise multi-level projection for coarse level
    call output_line('Assembling multilevel projections...')
    call mlprj_initProjectionMat (Rlevels(NLMIN)%rprojection, Rlevels(NLMIN)%rmatSys)

    do i = NLMIN + 1, NLMAX
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Prolongation Matrix Assembly
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Create prolongation matrix structures
      call mlop_create2LvlMatrixStruct(Rlevels(i-1)%rdiscr%RspatialDiscr(1), &
          Rlevels(i)%rdiscr%RspatialDiscr(1), LSYSSC_MATRIX9, Rlevels(i)%rmatProlVelocity)
      call mlop_create2LvlMatrixStruct(Rlevels(i-1)%rdiscr%RspatialDiscr(3), &
          Rlevels(i)%rdiscr%RspatialDiscr(3), LSYSSC_MATRIX9, Rlevels(i)%rmatProlPressure)

      ! Assemble prolongation matrices
      call mlop_build2LvlProlMatrix (Rlevels(i-1)%rdiscr%RspatialDiscr(1),&
          Rlevels(i)%rdiscr%RspatialDiscr(1), .true., Rlevels(i)%rmatProlVelocity)
      call mlop_build2LvlProlMatrix (Rlevels(i-1)%rdiscr%RspatialDiscr(3),&
          Rlevels(i)%rdiscr%RspatialDiscr(3), .true., Rlevels(i)%rmatProlPressure)

      ! Initialise multi-level projection
      call mlprj_initProjectionMat (Rlevels(i)%rprojection, Rlevels(i)%rmatSys)
      call mlprj_initMatrixProjection(&
          Rlevels(i)%rprojection%RscalarProjection(1,1), Rlevels(i)%rmatProlVelocity)
      call mlprj_initMatrixProjection(&
          Rlevels(i)%rprojection%RscalarProjection(1,2), Rlevels(i)%rmatProlVelocity)
      call mlprj_initMatrixProjection(&
          Rlevels(i)%rprojection%RscalarProjection(1,3), Rlevels(i)%rmatProlPressure)

    end do

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Multigrid solver setup
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    call output_line('Setting up multigrid solver...')

    ! Define a filter chain; the second entry is an integral-mean filter for the pressure
    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL
    RfilterChain(2)%ifilterType = FILTER_TOL20
    RfilterChain(2)%itoL20component = 3

    ! Initialise multigrid
    call linsol_initMultigrid2 (p_rsolverNode,NLMAX-NLMIN+1,RfilterChain)
    
    ! Use BiCGStab-Vanka as a coarse grid solver
    call linsol_initVANKA (p_rsmoother)
    call linsol_initBiCGStab (p_rcoarseGridSolver,p_rsmoother,RfilterChain)
    
    ! The coarse grid in multigrid is always grid 1!
    call linsol_getMultigrid2Level (p_rsolverNode,1,p_rlevelInfo)
    p_rlevelInfo%p_rcoarseGridSolver => p_rcoarseGridSolver

    ! Now set up the other levels...
    do i = NLMIN+1, NLMAX
    
      ! Set up the VANKA smoother.
      call linsol_initVANKA (p_rsmoother)
      
      ! We will use 4 smoothing steps with damping parameter 0.7
      call linsol_convertToSmoother(p_rsmoother, 4, 0.7_DP)
      
      ! And add this multi-grid level. We will use the same smoother
      ! for pre- and post-smoothing.
      call linsol_getMultigrid2Level (p_rsolverNode,i-NLMIN+1,p_rlevelInfo)
      p_rlevelInfo%p_rpresmoother => p_rsmoother
      p_rlevelInfo%p_rpostsmoother => p_rsmoother

      ! Attach our user-defined projection to the level.
      call linsol_initProjMultigrid2Level(p_rlevelInfo,Rlevels(i)%rprojection)

    end do
    
    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2

    ! Attach the system matrix to the solver.
    allocate(Rmatrices(NLMIN:NLMAX))
    do i = NLMIN, NLMAX
      call lsysbl_duplicateMatrix (Rlevels(i)%rmatSys,&
          Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    end do
    call linsol_setMatrices(p_RsolverNode,Rmatrices(NLMIN:NLMAX))
    do i=NLMIN,NLMAX
      call lsysbl_releaseMatrix (Rmatrices(i))
    end do
    deallocate(Rmatrices)
    
    ! Initialise solver
    call output_line('Initialising multigrid solver...')
    call linsol_initStructure (p_rsolverNode, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) then
      call output_line("Failed to initialise solver structure!", OU_CLASS_ERROR)
      call sys_halt()
    end if
    call linsol_initData (p_rsolverNode, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) then
      call output_line("Failed to initialise solver data!", OU_CLASS_ERROR)
      call sys_halt()
    end if

    ! Solve...
    call output_lbrk()
    call output_line('Solving...')
    call linsol_solveAdaptively (p_rsolverNode,rvecSol,rvecRhs,rvecTmp)

    ! filter pressure
    call vecfil_normaliseToL20Sca(rvecSol%RvectorBlock(3))

    ! Release solver data and structure
    call linsol_doneData (p_rsolverNode)
    call linsol_doneStructure (p_rsolverNode)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Error Analysis
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    call output_lbrk()
    call output_line('Performing error analysis...')

    ! Store the viscosity parameter nu in the collection's quick access array
    rcollection%DquickAccess(1) = dnu
    rcollection%DquickAccess(2) = dc

    ! Set up the error structure for velocity
    rerrorU%p_RvecCoeff => rvecSol%RvectorBlock(1:2)
    rerrorU%p_DerrorL2 => DerrorUL2
    rerrorU%p_DerrorH1 => DerrorUH1

    ! Set up the error structure for pressure
    rerrorP%p_RvecCoeff => rvecSol%RvectorBlock(3:3)
    rerrorP%p_DerrorL2 => DerrorPL2

    ! Calculate errors of velocity and pressure against analytic solutions.
    call pperr_scalarVec(rerrorU, funcVelocity2D, rcollection);
    call pperr_scalarVec(rerrorP, funcPressure2D, rcollection);

    ! Calculate errors of velocity field
    dL2 = sqrt(DerrorUL2(1)**2 + DerrorUL2(2)**2)
    dH1 = sqrt(DerrorUH1(1)**2 + DerrorUH1(2)**2)

    ! Calculate divergence of velocity field
    call aux_calcDiv2D(ddiv, rvecSol%RvectorBlock(1:2))

    ! Print the errors.
    call output_lbrk()
    call output_line('Error Analysis')
    call output_line('--------------')
    call output_line('|u - u_h|_L2 = ' // trim(sys_sdEL(dL2, 10)) // ' ( ' &
                                       // trim(sys_sdEL(DerrorUL2(1), 10)) // ' , ' &
                                       // trim(sys_sdEL(DerrorUL2(2), 10)) // ' )')
    call output_line('|u - u_h|_H1 = ' // trim(sys_sdEL(dH1, 10)) // ' ( ' &
                                       // trim(sys_sdEL(DerrorUH1(1), 10)) // ' , ' &
                                       // trim(sys_sdEL(DerrorUH1(2), 10)) // ' )')
    call output_line('|div u_h|_L2 = ' // trim(sys_sdEL(ddiv, 10)))
    call output_line('|p - p_h|_L2 = ' // trim(sys_sdEL(derrorPL2(1), 10)))

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! UCD Output
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    call output_lbrk()
    call output_line('Writing VTK output...')

    ! Start UCD export to VTK file:
    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = './ucd'
    call ucd_startVTK (rexport, UCD_FLAG_STANDARD, Rlevels(NLMAX)%rtria, &
        trim(sucddir)//'/u2d_vortex_noslip.vtk')

    ! Allocate temporary memory for projection
    allocate(p_Du(Rlevels(NLMAX)%rtria%NVT))
    allocate(p_Dv(Rlevels(NLMAX)%rtria%NVT))
    allocate(p_Dp(Rlevels(NLMAX)%rtria%NEL))

    ! Project velocity and pressure
    call spdp_projectToVertices (rvecSol%RvectorBlock(1), p_Du)
    call spdp_projectToVertices (rvecSol%RvectorBlock(2), p_Dv)
    call spdp_projectToCells(rvecSol%RvectorBlock(3), p_Dp)

    ! Write velocity and pressure
    call ucd_addVarVertBasedVec(rexport, 'velocity', p_Du, p_Dv)
    call ucd_addVariableElementBased (rexport,'pressure',UCD_VAR_STANDARD, p_Dp)
    
    ! Write the file to disc, that is it.
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! deallocate temporary arrays
    deallocate(p_Dp)
    deallocate(p_Dv)
    deallocate(p_Du)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Cleanup
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    call output_line('Cleaning up...')

    ! Release the solver node and all subnodes attached to it (if at all):
    call linsol_releaseSolver (p_rsolverNode)
    
    ! Release projections
    do i = NLMAX, NLMIN+1, -1
      call mlprj_doneProjection(Rlevels(i)%rprojection)
      call lsyssc_releaseMatrix (Rlevels(i)%rmatProlVelocity)
      call lsyssc_releaseMatrix (Rlevels(i)%rmatProlPressure)
    end do
    call mlprj_doneProjection(Rlevels(NLMIN)%rprojection)

    ! Release the block matrix/vectors
    call lsysbl_releaseVector (rvecTmp)
    call lsysbl_releaseVector (rvecRhs)
    call lsysbl_releaseVector (rvecSol)
    do i = NLMAX, NLMIN, -1
      call lsysbl_releaseMatrix (Rlevels(i)%rmatSys)
      call bcasm_releaseDiscreteBC (Rlevels(i)%rdiscreteBC)
      call spdiscr_releaseCubStructure(Rlevels(i)%rcubInfo)
      call spdiscr_releaseBlockDiscr(Rlevels(i)%rdiscr)
      call tria_done (Rlevels(i)%rtria)
    end do
    
    deallocate(Rlevels)
    
    ! Finally release the domain, that is it.
    call boundary_release (rboundary)

  end subroutine

  ! ***********************************************************************************************

  subroutine funcZeroBC2D (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
                           cinfoNeeded,iwhere,dwhere, Dvalues, rcollection)
  integer, dimension(:), intent(in)                           :: Icomponents
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  type(t_boundaryRegion), intent(in)                          :: rboundaryRegion
  integer, intent(in)                                         :: ielement
  integer, intent(in)                                         :: cinfoNeeded
  integer, intent(in)                                         :: iwhere
  real(DP), intent(in)                                        :: dwhere
  type(t_collection), intent(inout), optional                 :: rcollection
  real(DP), dimension(:), intent(out)                         :: Dvalues

    Dvalues(1) = 0.0_DP

  end subroutine

  ! ***********************************************************************************************

  subroutine funcRhsX2D (rdiscretisation,rform,nelements,npointsPerElement,Dpoints, &
                         IdofsTest,rdomainIntSubset,Dcoefficients,rcollection)
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  type(t_linearForm), intent(in)                              :: rform
  integer, intent(in)                                         :: nelements
  integer, intent(in)                                         :: npointsPerElement
  real(DP), dimension(:,:,:), intent(in)  :: Dpoints
  integer, dimension(:,:), intent(in) :: IdofsTest
  type(t_domainIntSubset), intent(in)              :: rdomainIntSubset
  type(t_collection), intent(inout), optional      :: rcollection
  real(DP), dimension(:,:,:), intent(out) :: Dcoefficients

  real(DP) :: dnu,dc

    ! fetch coefficients from collection
    dnu = rcollection%DquickAccess(1)
    dc = rcollection%DquickAccess(2)

    ! f_1(x,y) = -pi * (c*cos(pi*x) - 4*pi*nu*sin(2*pi*y) * (1 - 2*cos(2*pi*x)))
    Dcoefficients(1,:,:) = -SYS_PI * (dc*cos(SYS_PI * Dpoints(1,:,:)) - &
      4.0_DP * SYS_PI * dnu * sin(2.0_DP * SYS_PI * Dpoints(2,:,:)) * &
      (1.0_DP - 2.0_DP * cos(2.0_DP * SYS_PI * Dpoints(1,:,:))))

  end subroutine

  ! ***********************************************************************************************

  subroutine funcRhsY2D (rdiscretisation,rform,nelements,npointsPerElement,Dpoints, &
                         IdofsTest,rdomainIntSubset,Dcoefficients,rcollection)
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  type(t_linearForm), intent(in)                              :: rform
  integer, intent(in)                                         :: nelements
  integer, intent(in)                                         :: npointsPerElement
  real(DP), dimension(:,:,:), intent(in)  :: Dpoints
  integer, dimension(:,:), intent(in) :: IdofsTest
  type(t_domainIntSubset), intent(in)              :: rdomainIntSubset
  type(t_collection), intent(inout), optional      :: rcollection
  real(DP), dimension(:,:,:), intent(out) :: Dcoefficients

  real(DP) :: dnu,dc

    ! fetch coefficients from collection
    dnu = rcollection%DquickAccess(1)
    dc = rcollection%DquickAccess(2)

    ! f_2(x,y) = -pi * (c*cos(pi*y) + 4*pi*nu*sin(2*pi*x) * (1 - 2*cos(2*pi*y)))
    Dcoefficients(1,:,:) = -SYS_PI * (dc * cos(SYS_PI * Dpoints(2,:,:)) + &
      4.0_DP * SYS_PI * dnu * sin(2.0_DP * SYS_PI * Dpoints(1,:,:)) * &
      (1.0_DP - 2.0_DP * cos(2.0_DP * SYS_PI * Dpoints(2,:,:))))

  end subroutine

  ! ***********************************************************************************************

  subroutine funcVelocity2D (icomponent, cderivative, rdiscretisation, &
                nelements, npointsPerElement, Dpoints, rdomainIntSubset, &
                Dvalues, rcollection)
  integer, intent(in)                          :: icomponent
  integer, intent(in)                          :: cderivative
  type(t_spatialDiscretisation), intent(in)    :: rdiscretisation
  integer, intent(in)                          :: nelements
  integer, intent(in)                          :: npointsPerElement
  real(DP), dimension(:,:,:), intent(in)       :: Dpoints
  type(t_domainIntSubset), intent(in)          :: rdomainIntSubset
  type(t_collection), intent(inout), optional  :: rcollection
  real(DP), dimension(:,:), intent(out)        :: Dvalues

  ! define 2*pi; this cannot be a parameter as SYS_PI is not a parameter...
  real(DP) :: pi2

    pi2 = 2.0_DP * SYS_PI

    select case(icomponent)
    case (1)
      ! X-velocity
      select case(cderivative)
      case (DER_FUNC2D)
        ! u_1(x,y) = (1 - cos(2*pi*x)) * sin(2*pi*y)
        Dvalues(:,:) = (1.0_DP - cos(pi2 * Dpoints(1,:,:))) * sin(pi2 * Dpoints(2,:,:))

      case (DER_DERIV2D_X)
        ! d_x u_1(x,y) = 2 * pi * sin(2*pi*x) * sin(2*pi*y)
        Dvalues(:,:) = pi2 * sin(pi2 * Dpoints(1,:,:)) * sin(pi2 * Dpoints(2,:,:))

      case (DER_DERIV2D_Y)
        ! d_y u_1(x,y) = 2 * pi *(1- cos(2*pi*x)) * cos(2*pi*y)
        Dvalues(:,:) = pi2 * (1.0_DP - cos(pi2 * Dpoints(1,:,:))) * cos(pi2 * Dpoints(2,:,:))

      end select

    case (2)
      ! Y-velocity
      select case(cderivative)
      case (DER_FUNC2D)
        ! u_2(x,y) = -(1 - cos(2*pi*y)) * sin(2*pi*x)
        Dvalues(:,:) = - (1.0_DP - cos(pi2 * Dpoints(2,:,:))) * sin(pi2 * Dpoints(1,:,:))

      case (DER_DERIV2D_X)
        ! d_x u_2(x,y) = -2 * pi * (1 - cos(pi*y)) * sin(pi*y)
        Dvalues(:,:) = -pi2 * (1.0_DP - cos(pi2 * Dpoints(2,:,:))) * cos(pi2 * Dpoints(1,:,:))

      case (DER_DERIV2D_Y)
        ! d_y u_2(x,y) = -2 * pi * sin(pi*x) * sin(pi*y)
        Dvalues(:,:) = -pi2 * sin(pi2 * Dpoints(1,:,:)) * sin(pi2 * Dpoints(2,:,:))

      end select

    end select

  end subroutine

  ! ***********************************************************************************************

  subroutine funcPressure2D (icomponent, cderivative, rdiscretisation, &
                nelements, npointsPerElement, Dpoints, rdomainIntSubset, &
                Dvalues, rcollection)
  integer, intent(in)                          :: icomponent
  integer, intent(in)                          :: cderivative
  type(t_spatialDiscretisation), intent(in)    :: rdiscretisation
  integer, intent(in)                          :: nelements
  integer, intent(in)                          :: npointsPerElement
  real(DP), dimension(:,:,:), intent(in)       :: Dpoints
  type(t_domainIntSubset), intent(in)          :: rdomainIntSubset
  type(t_collection), intent(inout), optional  :: rcollection
  real(DP), dimension(:,:), intent(out)        :: Dvalues

  real(DP) :: dc

    dc = rcollection%DquickAccess(2)

    ! p(x,y) = c*(4/pi - sin(pi*x) - sin(pi*y))
    Dvalues(:,:) = dc*(4.0_DP / SYS_PI - sin(SYS_PI * Dpoints(1,:,:)) - sin(SYS_PI * Dpoints(2,:,:)))

  end subroutine

end module
