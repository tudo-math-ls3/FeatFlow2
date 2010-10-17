!##############################################################################
!# ****************************************************************************
!# <name> spacetimelinsol </name>
!# ****************************************************************************
!#
!# <purpose>
!# Lightweight linear solver library for space-time systems.
!# </purpose>
!##############################################################################

module spacetimelinsol

  use fsystem
  use storage
  use genoutput
  use linearalgebra
  use linearsystemscalar
  use linearsystemblock
  use linearsolver
  use matrixfilters
  use matrixio
  
  use physics
  
  use spacetimematvec
  use spacetimehierarchy
  use fespacehierarchy
  use spacetimeinterlevelprojection
  
  use postprocessing
  
  use statistics
  
  implicit none

!<constantblock description="Linear space-time solver tyes">
  ! None
  integer, parameter :: STLS_TYPE_NONE       = 0
  
  ! Jacobi
  integer, parameter :: STLS_TYPE_JACOBI     = 1
  
  ! FB-GS
  integer, parameter :: STLS_TYPE_FBSIM      = 2

  ! FB-GS2
  integer, parameter :: STLS_TYPE_FBGS       = 3
  
  ! Defect correction
  integer, parameter :: STLS_TYPE_DEFCORR    = 4
  
  ! BiCGStab 
  integer, parameter :: STLS_TYPE_BICGSTAB   = 5

  ! Two-grid
  integer, parameter :: STLS_TYPE_MULTIGRID  = 6
!</constantblock>

!<constantblock description="Type of preconditioner in space">

  ! No preconditioning; Used for UMFPACK solvers e.g. 
  ! which do not use preconditioning.
  integer, parameter, public :: STLS_PC_NONE = -1

  ! Jacobi preconditioning
  integer, parameter, public :: STLS_PC_JACOBI = 0

  ! Diagonal VANKA preconditioner
  integer, parameter, public :: STLS_PC_VANKA = 0
  
  ! ILU-0 preconditioner
  integer, parameter, public :: STLS_PC_ILU0 = 1

  ! SSOR preconditioner
  integer, parameter, public :: STLS_PC_SSOR = 2

  ! Jacobi preconditioner + BiCGStab
  integer, parameter, public :: STLS_PC_BICGSTABJACOBI = 3

  ! General VANKA preconditioner + BiCGStab
  integer, parameter, public :: STLS_PC_BICGSTABVANKA = 3

  ! UMFPACK
  integer, parameter, public :: STLS_PC_UMFPACK = 4

  ! Full VANKA preconditioner + BiCGStab
  integer, parameter, public :: STLS_PC_BICGSTABFULLVANKA = 5
  
!</constantblock>

!<constantblock description="Problem type">

  ! Standard problem. All diagonal blocks are invertible.
  integer, parameter, public :: STLS_PR_STANDARD = 0

  ! 2D Saddle point, 2 equations.
  ! (Optimal control of Stokes, Navier--Stokes,...
  integer, parameter, public :: STLS_PC_2DSADDLEPT2EQ = 1

!</constantblock>
  
  ! Matrix pointer
  type t_spacetimeMatrixPointer
    ! Pointer to a matrix
    type(t_spaceTimeMatrix), pointer :: p_rmatrix
  end type
  
  ! Parameters for space-solver
  type t_spaceSolverParams
  
    ! Type of the linear solver in space if used.
    ! =LINSOL_ALG_UNDEFINED: not used.
    ! =LINSOL_ALG_UMFPACK4 : UMFPACK
    integer :: cspaceSolverType = LINSOL_ALG_UNDEFINED
    
    ! Relative stopping criterion solver
    real(DP) :: depsrel = 1E-5

    ! Absolute stopping criterion solver
    real(DP) :: depsabs = 1E-14
    
    ! Max. number of iterations
    integer :: nmaxiterations = 100
    
    ! Output level
    integer :: ioutputLevel = 0
    
    ! Problem type. Decides upon which smoothers are allowed.
    integer :: cproblemtype = STLS_PR_STANDARD

    ! ONLY MG: Type of smoother.
    integer :: csmoother = 0

    ! ONLY MG: Number of (post-)smoothing steps
    integer :: nsmpost = 4
    
    ! ONLY MG: Damping parameter for the smoother
    real(DP) :: domegaSmoother = 1.0_DP

    ! ONLY MG: Type identifier for the space smoother.
    ! A STLS_PC_xxxx constant.
    integer :: cspacesmoother = STLS_PC_JACOBI

    ! ONLY MG: Type identifier for the space coarse grid solver.
    ! A STLS_PC_xxxx constant.
    integer :: cspacecoarsegridsolver = STLS_PC_JACOBI

    ! ONLY MG: Output level on the coarse grid
    integer :: ioutputLevelCoarse = 0
    
    ! ONLY MG: Relative stopping criterion solver, coarse grid solver
    real(DP) :: depsrelCoarse = 1E-10

    ! ONLY MG: Absolute stopping criterion solver, coarse grid solver
    real(DP) :: depsabsCoarse = 1E-14
    
    ! ONLY MG: Max. number of iterations, coarse grid solver
    integer :: nmaxiterationsCoarse = 1000
    
    ! Minimum difference in the solutions before the solution is treated as
    ! converged.
    real(DP) :: depsdiff = 1E-6_DP

    ! Minimum difference in the solutions before the solution is treated as
    ! converged. Coarse grid.
    real(DP) :: depsdiffCoarse = 1E-6_DP

  end type

  ! Linear solver structure.
  type t_spacetimelinsol
  
    ! Output: Status of the solver.
    ! =0: solver converged.
    ! =1: solver did not converge
    ! =2: solver diverged.
    integer :: csolverStatus = 0
    
    ! Type of the solver.
    integer :: csolverType = STLS_TYPE_NONE
    
    ! Relaxation parameter; algorithm dependent.
    real(DP) :: drelax = 1.0_DP
    
    ! Damping parameter for the return vector.
    real(DP) :: domega = 1.0_DP

    ! Relative stopping criterion
    real(DP) :: depsRel = 1E-12_DP

    ! Absolute stopping criterion
    real(DP) :: depsAbs = 1E-14_DP

    ! Relative difference stopping criterion
    real(DP) :: depsRelDiff = 1E-4_DP

    ! Absolute difference stopping criterion
    real(DP) :: depsAbsDiff = 0.0_DP
    
    ! Absolute divergence criterion
    real(DP) :: ddivAbs = 1E20_DP

    ! Relative divergence criterion
    real(DP) :: ddivRel = 1E10_DP
    
    ! Whether or not to check the residual for convergence.
    ! =0: do not check.
    ! =1: check: relative or absolute stopping criterion must be fulfilled.
    ! =2: check: relative and absolute stopping criterion must be fulfilled.
    integer :: iresCheck = 1
    
    ! Number of levels
    integer :: nlevels = 1
    
    ! Type of cycle.
    ! =0: F-cycle, =1: V-cycle, =2: W-cycle
    integer :: icycle = 1
    
    ! Minimum/Maximum number of steps
    integer :: nminIterations = 1
    integer :: nmaxIterations = 1000
    
    ! Output level
    ! =0: no output
    ! =1: basic output
    ! =2: standard output
    ! =3: extended output
    integer :: ioutputLevel = 2
    
    ! Algorithm-specific option field.
    integer :: ialgOptions = 0
    
    ! ONLY MULTIGRID: Specifies the type of the adaptive coarse grid correction used.
    ! =0: deactivated.
    ! =1: adaptive coarse grid correction by energy minimisation.
    integer :: iadcgcorr = 0
    
    ! Min/max correction factor
    real(DP) :: dadcgcorrMin = 0.5_DP
    real(DP) :: dadcgcorrMax = 2.0_DP
    
    ! Determins the number of iterations after which the
    ! algorithm is reinitialised (if it supports reinitialisation).
    ! =0: deactivate.
    integer :: niteReinit = 0
    
    ! ONLY BICGSTAB: If set to TRUE, BiCGStab measures the real residuum for
    ! determining the stopping criterion.
    logical :: bstopOnRealResiduum = .false.
    
    ! System matrix if used as 1-level solver.
    type(t_spaceTimeMatrix), pointer :: p_rmatrix => null()
    
    ! Preconditioner (if present), or null
    type(t_spacetimelinsol), pointer :: p_rpreconditioner => null()
    
    ! Underlying space-time hierarchy.
    type(t_spacetimeHierarchy), pointer :: p_rspaceTimeHierarchy => null()
    
    ! Level, this solver represents in the space-time hierarchy.
    integer :: ilevel = 0
    
    ! ONLY MULTIGRID: Smoother on each level
    type(t_spacetimelinsol), dimension(:), pointer :: p_Rpresmoothers => null()
    type(t_spacetimelinsol), dimension(:), pointer :: p_Rpostsmoothers => null()

    ! DEBUG!!! Preconditioners
    type(t_spacetimelinsol), dimension(:), pointer :: p_Rpreconditioners => null()

    ! ONLY MULTIGRID: Coarse grid solver
    type(t_spacetimelinsol), pointer :: p_rcoarsegridsolver => null()
    
    ! ONLY MULTIGRID: System matrices if used as multigrid solver.
    type(t_spacetimeMatrixPointer), dimension(:), pointer :: p_rmatrices => null()
    
    ! ----------
    ! Statistics
    ! ----------
    
    ! Sums up the total time of the solver. Is reset in initData!
    type(t_timer) :: rtotalTime
    
    ! ------------------------------------
    ! Parameters maintaines by each solver
    ! ------------------------------------
    
    ! General purpose temporary space-time vectors
    type(t_spaceTimeVector) :: rspaceTimeTemp1
    type(t_spaceTimeVector) :: rspaceTimeTemp2
    
    ! General purpose temporary space vectors
    type(t_vectorBlock) :: rspaceTemp1
    type(t_vectorBlock) :: rspaceTemp2
    type(t_vectorBlock) :: rspaceTemp3
    
    ! General purpose temporary space matrices/vectors
    type(t_matrixBlock), dimension(:), pointer :: p_RspaceMatrices => null()
    type(t_vectorBlock), dimension(:), pointer :: p_RspaceVectors => null()
    type(t_spacetimeVector), dimension(:), pointer :: p_Rvectors1 => null()
    type(t_spacetimeVector), dimension(:), pointer :: p_Rvectors2 => null()
    type(t_spacetimeVector), dimension(:), pointer :: p_Rvectors3 => null()
    type(t_spacetimeVector), dimension(:), pointer :: p_Rvectors4 => null()
    type(t_spacetimeVector), dimension(:), pointer :: p_Rvectors5 => null()
    
    ! Temporary vectors to use during the solution process
    type(t_spaceTimeVector), dimension(:), pointer :: p_RtempVectors => null()
    
    ! Parameters defining the linear solver in space
    type(t_spaceSolverParams) :: rspaceSolverParams
    
    ! Linear solver in space.
    type(t_linsolNode), pointer :: p_rspaceSolver => null()
    
    ! Boundary conditions in space
    type(t_discreteBC), dimension(:), pointer :: p_RdiscreteBC => null()
    
    ! Spatial template matrices on all space levels.
    type(t_matvecTemplates), dimension(:), pointer :: p_RmatVecTempl => null()
    
    ! ONLY MULTIGRID: Interlevel projection hierarchy for switching between levels.
    type(t_sptiProjHierarchy), pointer :: p_rspaceTimeProjection => null()
    
  end type
  

contains

  ! ***************************************************************************

  subroutine stls_initspacesolver (p_rsolver, cproblemtype, csolvertype, &
      depsrel, depsabs, depsdiff, nmaxiterations, ioutputlevel)
  
  ! Initialises a space solver.
  
  ! Type of the problem
  integer, intent(in) :: cproblemtype
  
  ! Type of the solver
  integer, intent(in) :: csolvertype
  
  ! Stopping criteria
  real(DP), intent(in), optional :: depsrel
  real(DP), intent(in), optional :: depsabs
  real(DP), intent(in), optional :: depsdiff
  integer, intent(in), optional :: nmaxiterations
  
  ! Output level
  integer, intent(in), optional :: ioutputlevel
  
  ! OUTPUT: Pointer to the solver
  type(t_linsolNode), pointer :: p_rsolver
  
    ! local variables
    type(t_linsolNode), pointer :: p_rpreconditioner
  
    ! Create the solver depending on the problem
    select case (cproblemtype)
    case (STLS_PR_STANDARD)
    
      select case (csolvertype)
      case (STLS_PC_JACOBI)  
        call linsol_initJacobi (p_rpreconditioner)
        call linsol_initDefCorr (p_rsolver,p_rpreconditioner)
        
      case (STLS_PC_BICGSTABJACOBI)
        call linsol_initJacobi (p_rpreconditioner)
        call linsol_initBiCGStab (p_rsolver,p_rpreconditioner)
        
      case (STLS_PC_ILU0)
        call linsol_initILU0 (p_rsolver)
        
      case (STLS_PC_SSOR)
        call linsol_initSSOR (p_rsolver)

      case (STLS_PC_UMFPACK)
        call linsol_initUMFPACK4 (p_rsolver)
      end select
      
      p_rsolver%istoppingCriterion = LINSOL_STOP_ONEOF
      
      if (present(depsRel)) &
        p_rsolver%depsRel = depsRel
      if (present(depsAbs)) &
        p_rsolver%depsAbs = depsAbs
      if (present(depsdiff)) &
        p_rsolver%depsdiff = depsdiff
      if (present(nmaxIterations)) &
        p_rsolver%nmaxIterations = nmaxIterations
      if (present(ioutputlevel)) &
        p_rsolver%ioutputlevel = ioutputlevel
      
    case (STLS_PC_2DSADDLEPT2EQ)

      select case (csolvertype)
      case (STLS_PC_VANKA)  
        call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG2)
        call linsol_initDefCorr (p_rsolver,p_rpreconditioner)
        
      case (STLS_PC_BICGSTABVANKA)
        call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG2)
        call linsol_initBiCGStab (p_rsolver,p_rpreconditioner)

      case (STLS_PC_BICGSTABFULLVANKA)
        call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_GENERAL)
        call linsol_initBiCGStab (p_rsolver,p_rpreconditioner)
        
      case (STLS_PC_ILU0)
        call output_line ("Invalid solver.")
        call sys_halt()
        
      case (STLS_PC_SSOR)
        call output_line ("Invalid solver.")
        call sys_halt()

      case (STLS_PC_UMFPACK)
        call linsol_initUMFPACK4 (p_rsolver)
      end select

      p_rsolver%istoppingCriterion = LINSOL_STOP_ONEOF
      if (present(depsRel)) &
        p_rsolver%depsRel = depsRel
      if (present(depsAbs)) &
        p_rsolver%depsAbs = depsAbs
      if (present(depsdiff)) &
        p_rsolver%depsdiff = depsdiff
      if (present(nmaxIterations)) &
        p_rsolver%nmaxIterations = nmaxIterations
      if (present(ioutputlevel)) &
        p_rsolver%ioutputlevel = ioutputlevel
      
    end select
  
  end subroutine

  ! ***************************************************************************

  subroutine stls_init (rsolver,csolverType,rspaceTimeHierarchy,ilevel,&
      ballocTempMatrices,ballocTempVectors,rpreconditioner,&
      rspaceSolverParams,RmatVecTempl)
  
  ! Basic initialisation of a solver node.
  
  ! Solver structure to be initialised
  type(t_spacetimelinsol), intent(out) :: rsolver
  
  ! Type identifier of the solver.
  integer, intent(in) :: csolverType
  
  ! Underlying space-time hierarchy
  type(t_spacetimeHierarchy), intent(in), target :: rspaceTimeHierarchy
  
  ! Level of the solver  in the space-time hierarchy
  integer, intent(in) :: ilevel
  
  ! TRUE = the solver needs temp matrices in space.
  ! FALSE = automatically decide depending on the solver in space.
  logical, intent(in) :: ballocTempMatrices
  
  ! TRUE = the solver needs temp vectors in space.
  ! FALSE = temp vectors in space not needed.
  logical, intent(in) :: ballocTempVectors
  
  ! OPTIONAL: Solver structure of a preconditioner
  type(t_spacetimelinsol), intent(in), target, optional :: rpreconditioner
  
  ! OPTIONAL: Structure defining the linear solver.
  type(t_spaceSolverParams), intent(in), optional :: rspaceSolverParams
  
  ! OPTINAL: Spatial matrix templates. Necessary if a matrix-based
  ! preconditioner in space is used.
  type(t_matvecTemplates), dimension(:), target, optional :: RmatVecTempl
  
    ! local variables
    integer :: ispacelevel,ilev
    type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo
  
    rsolver%csolverType = csolverType
    rsolver%ilevel = ilevel
    rsolver%p_rspaceTimeHierarchy => rspaceTimeHierarchy
    
    ! Allocate memory for the system matrices
    allocate(rsolver%p_rmatrices(ilevel))
    
    ! Attach the preconditioner if present
    if (present(rpreconditioner)) then
      rsolver%p_rpreconditioner => rpreconditioner
    end if

    if (present(RmatVecTempl)) then    
      ! Remember matrix templates
      rsolver%p_RmatVecTempl => RmatVecTempl
    end if
    
    if (present(rspaceSolverParams) .or. ballocTempMatrices) then

      ! Get the space level associated to the space-time level    
      call sth_getLevel (rspaceTimeHierarchy,ilevel,ispaceLevel=ispaceLevel)
      
      ! Allocate matrices for the space solver
      allocate(rsolver%p_RspaceMatrices(ispaceLevel))

      ! Initialise the solver in space.
      rsolver%rspaceSolverParams = rspaceSolverParams
      select case (rspaceSolverParams%cspaceSolverType)
      case (LINSOL_ALG_UMFPACK4) 
        call linsol_initUMFPACK4 (rsolver%p_rspaceSolver)
        !rsolver%p_rspaceSolver%p_rsubnodeUmfpack4%imatrixDebugOutput = 1
        
      case (LINSOL_ALG_MULTIGRID2)
      
        ! Initialise the MG solver
        call linsol_initMultigrid2 (rsolver%p_rspaceSolver,ispaceLevel)
        
        ! Initialise the levels
        do ilev = 1,ispaceLevel
        
          call linsol_getMultigrid2Level (rsolver%p_rspaceSolver,ilev,p_rlevelInfo)
          if (ilev .eq. 1) then
            ! Coarse grid solver
            if (ispaceLevel .eq. 1) then
              ! Only one level
              call stls_initspacesolver (p_rlevelInfo%p_rcoarseGridSolver, &
                  rspaceSolverParams%cproblemtype, rspaceSolverParams%cspacecoarsegridsolver, &
                  rspaceSolverParams%depsRel, rspaceSolverParams%depsAbs, &
                  rspaceSolverParams%depsdiff,&
                  rspaceSolverParams%nmaxIterations, rspaceSolverParams%ioutputLevel)
            else            
              call stls_initspacesolver (p_rlevelInfo%p_rcoarseGridSolver, &
                  rspaceSolverParams%cproblemtype, rspaceSolverParams%cspacecoarsegridsolver, &
                  rspaceSolverParams%depsRelCoarse, rspaceSolverParams%depsAbsCoarse, &
                  rspaceSolverParams%depsdiffCoarse,&
                  rspaceSolverParams%nmaxIterationsCoarse, rspaceSolverParams%ioutputLevelCoarse)
            end if
          else
            ! Smoother
            call stls_initspacesolver (p_rlevelInfo%p_rpostsmoother, &
                rspaceSolverParams%cproblemtype, rspaceSolverParams%cspacesmoother)            

            call linsol_convertToSmoother (p_rlevelInfo%p_rpostsmoother,&
                rspaceSolverParams%nsmPost,rspaceSolverParams%domegaSmoother)
          end if
      
          ! Stopping criteria
          rsolver%p_rspaceSolver%depsRel = rspaceSolverParams%depsRel
          rsolver%p_rspaceSolver%depsAbs = rspaceSolverParams%depsAbs
          rsolver%p_rspaceSolver%depsdiff = rspaceSolverParams%depsdiff
          rsolver%p_rspaceSolver%nmaxIterations = rspaceSolverParams%nmaxIterations
          rsolver%p_rspaceSolver%ioutputlevel = rspaceSolverParams%ioutputlevel
          rsolver%p_rspaceSolver%istoppingCriterion = LINSOL_STOP_ONEOF
      
        end do
      
      end select
        
      ! Prepare boundary conditions
      allocate(rsolver%p_RdiscreteBC(ispaceLevel))
      do ilev = 1,ispaceLevel
        call bcasm_initDiscreteBC(rsolver%p_RdiscreteBC(ilev))
      end do
      
    end if
  
    if (ballocTempVectors) then
    
      ! Get the space level of the space-time level
      call sth_getLevel (rspaceTimeHierarchy,ilevel,ispaceLevel=ispaceLevel)
    
      ! Allocate space vectors for the space solver
      allocate(rsolver%p_RspaceVectors(ispaceLevel))

    end if
  
  end subroutine

  ! ***************************************************************************

  subroutine stls_done (rsolver)
  
  ! Releases a solver node
  
  ! Solver structure to be released
  type(t_spacetimelinsol), intent(inout) :: rsolver
  
    integer :: ilev
  
    rsolver%rspaceSolverParams%cspaceSolverType = LINSOL_ALG_UNDEFINED
    if (associated(rsolver%p_RdiscreteBC)) then
      ! Release boundary conditions
      do ilev = 1,size(rsolver%p_RdiscreteBC)
        call bcasm_releaseDiscreteBC(rsolver%p_RdiscreteBC(ilev))
      end do
      deallocate(rsolver%p_RdiscreteBC)
    end if

    if (associated(rsolver%p_rspaceSolver)) then
      ! Release solver
      call linsol_releaseSolver (rsolver%p_rspaceSolver)
    end if

    if (associated(rsolver%p_RspaceMatrices)) then
      deallocate(rsolver%p_RspaceMatrices)
    end if

    if (associated(rsolver%p_RspaceVectors)) then
      deallocate(rsolver%p_RspaceVectors)
    end if
  
  end subroutine

  ! ***************************************************************************

  recursive subroutine stls_setMatrix (rsolver,ilevel,rmatrix)
  
  ! Assigns the system matrix of level ilevel to the solver
  
  ! Solver structure to be initialised
  type(t_spacetimelinsol), intent(inout) :: rsolver
  
  ! Level, the matrix corresponds to.
  integer, intent(in) :: ilevel
  
  ! Space-time system matrix.
  type(t_spaceTimeMatrix), intent(in), target :: rmatrix
  
    integer :: i
  
    if (ilevel .eq. rsolver%ilevel) then
      rsolver%p_rmatrix => rmatrix
    end if
    
    if (ilevel .le. rsolver%ilevel) then
      
      ! Pass to subsolvers
      if (associated (rsolver%p_rpreconditioner)) then
        call stls_setMatrix(rsolver%p_rpreconditioner,ilevel,rmatrix)
      end if

      if (associated (rsolver%p_rcoarsegridsolver)) then
        call stls_setMatrix(rsolver%p_rcoarsegridsolver,ilevel,rmatrix)
      end if
      
      if (associated (rsolver%p_rpresmoothers)) then
        do i=1,rsolver%ilevel
          if (rsolver%p_rpresmoothers(i)%csolverType .ne. STLS_TYPE_NONE) then
            call stls_setMatrix(rsolver%p_rpresmoothers(i),ilevel,rmatrix)
          end if
        end do
      end if

      if (associated (rsolver%p_Rpostsmoothers)) then
        do i=1,rsolver%ilevel
          if (rsolver%p_Rpostsmoothers(i)%csolverType .ne. STLS_TYPE_NONE) then
            call stls_setMatrix(rsolver%p_rpostsmoothers(i),ilevel,rmatrix)
          end if
        end do
      end if

      if (associated (rsolver%p_Rpreconditioners)) then
        do i=1,rsolver%ilevel
          if (rsolver%p_Rpreconditioners(i)%csolverType .ne. STLS_TYPE_NONE) then
            call stls_setMatrix(rsolver%p_Rpreconditioners(i),ilevel,rmatrix)
          end if
        end do
      end if

      ! Remember the matrix
      rsolver%p_Rmatrices(ilevel)%p_rmatrix => rmatrix

    end if
  
  end subroutine

  ! ***************************************************************************

  recursive subroutine stls_initData (rsolver)
  
  ! Solver preparation right before solving.
  
  ! Solver structure to be initialised
  type(t_spacetimelinsol), intent(inout) :: rsolver
  
    ! Local variables
    type(t_feSpaceLevel), pointer :: p_rfeSpaceLevel
    type(t_blockDiscretisation), pointer :: p_rspaceDiscr
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr
    integer :: ilev,i,ispaceLevel
  
    if (rsolver%csolverType .eq. STLS_TYPE_NONE) then
      call sys_halt()
    end if
  
    call sth_getLevel (rsolver%p_rspaceTimeHierarchy,&
        rsolver%ilevel,p_rfeSpaceLevel=p_rfeSpaceLevel,p_rtimeDiscr=p_rtimeDiscr,&
        ispaceLevel=ispaceLevel)
    p_rspaceDiscr => p_rfeSpaceLevel%p_rdiscretisation
  
    ! Some initialisation, based on the solver type.
    select case (rsolver%csolverType)
      
    case (STLS_TYPE_JACOBI)
      
      ! Allocate temp vectors.
      call lsysbl_createVectorBlock (p_rspaceDiscr,rsolver%rspaceTemp1)
    
    case (STLS_TYPE_FBSIM)
      
      ! Allocate temp vectors.
      call lsysbl_createVectorBlock (p_rspaceDiscr,rsolver%rspaceTemp1)
      call lsysbl_createVectorBlock (p_rspaceDiscr,rsolver%rspaceTemp2)
      call lsysbl_createVectorBlock (p_rspaceDiscr,rsolver%rspaceTemp3)
      call sptivec_initVector (rsolver%rspaceTimeTemp1,p_rtimeDiscr,p_rspaceDiscr)

    case (STLS_TYPE_FBGS)
      
      ! Allocate temp vectors.
      call lsysbl_createVectorBlock (p_rspaceDiscr,rsolver%rspaceTemp1)
      call lsysbl_createVectorBlock (p_rspaceDiscr,rsolver%rspaceTemp2)
      call lsysbl_createVectorBlock (p_rspaceDiscr,rsolver%rspaceTemp3)
      call sptivec_initVector (rsolver%rspaceTimeTemp1,p_rtimeDiscr,p_rspaceDiscr)
    
    case (STLS_TYPE_DEFCORR)

      ! Allocate temp vectors.
      call sptivec_initVector (rsolver%rspaceTimeTemp1,p_rtimeDiscr,p_rspaceDiscr)
      call sptivec_initVector (rsolver%rspaceTimeTemp2,p_rtimeDiscr,p_rspaceDiscr)
      
      if (associated(rsolver%p_rpreconditioner)) then
        call stls_initData(rsolver%p_rpreconditioner)
      end if
    
    case (STLS_TYPE_BICGSTAB)

      ! Allocate temp vectors
      call lsysbl_createVectorBlock (p_rspaceDiscr,rsolver%rspaceTemp1)

      allocate (rsolver%p_rtempVectors(7))
      do i=1,7
        call sptivec_initVector (rsolver%p_rtempVectors(i),p_rtimeDiscr,p_rspaceDiscr)
      end do

      if (associated(rsolver%p_rpreconditioner)) then
        call stls_initData(rsolver%p_rpreconditioner)
      end if
    
    case (STLS_TYPE_MULTIGRID)
    
      ! Allocate temp vectors on all levels.
      allocate(rsolver%p_Rvectors1(rsolver%ilevel))
      allocate(rsolver%p_Rvectors2(rsolver%ilevel))
      allocate(rsolver%p_Rvectors3(rsolver%ilevel))
      allocate(rsolver%p_Rvectors4(rsolver%ilevel))
      
      if (rsolver%iadcgcorr .eq. 2) then
        allocate(rsolver%p_Rvectors5(rsolver%ilevel))
      end if
      
      do ilev = 1,rsolver%ilevel

        call sth_getLevel (rsolver%p_rspaceTimeHierarchy,&
            ilev,p_rfeSpaceLevel=p_rfeSpaceLevel,p_rtimeDiscr=p_rtimeDiscr)
        p_rspaceDiscr => p_rfeSpaceLevel%p_rdiscretisation

        call sptivec_initVector (rsolver%p_Rvectors1(ilev),p_rtimeDiscr,p_rspaceDiscr)
        call sptivec_initVector (rsolver%p_Rvectors2(ilev),p_rtimeDiscr,p_rspaceDiscr)
        call sptivec_initVector (rsolver%p_Rvectors3(ilev),p_rtimeDiscr,p_rspaceDiscr)
        call sptivec_initVector (rsolver%p_Rvectors4(ilev),p_rtimeDiscr,p_rspaceDiscr)

        if (associated(rsolver%p_Rvectors5)) then
          call sptivec_initVector (rsolver%p_Rvectors5(ilev),p_rtimeDiscr,p_rspaceDiscr)
        end if
      end do
    
      if (associated(rsolver%p_rcoarseGridSolver)) then
        call stls_initData(rsolver%p_rcoarseGridSolver)
      end if

      if (associated(rsolver%p_RpreSmoothers)) then
        do ilev = 2,rsolver%ilevel
          call stls_initData(rsolver%p_RpreSmoothers(ilev))
        end do
      end if

      if (associated(rsolver%p_RpostSmoothers)) then
        if (.not. associated (rsolver%p_RpreSmoothers,rsolver%p_RpostSmoothers)) then
          do ilev = 2,rsolver%ilevel
            call stls_initData(rsolver%p_RpostSmoothers(ilev))
          end do
        end if
      end if

      if (associated(rsolver%p_Rpreconditioners)) then
        do ilev = 2,rsolver%ilevel
          call stls_initData(rsolver%p_Rpreconditioners(ilev))
        end do
      end if
    
    end select
    
    ! Temp matrices.
    if (associated(rsolver%p_RspaceMatrices)) then
      do ilev = 1,ispaceLevel
        
        p_rfeSpaceLevel => rsolver%p_rspaceTimeHierarchy%p_rfeHierarchy%p_rfeSpaces(ilev)
        p_rspaceDiscr => p_rfeSpaceLevel%p_rdiscretisation
        
        call stmat_allocSubmatrix (rsolver%p_rmatrix%cmatrixType,&
            rsolver%p_rmatrix%p_rphysics,rsolver%p_RmatVecTempl(ilev),&
            rsolver%p_RspaceMatrices(ilev),rsolver%p_RdiscreteBC(ilev))
            
      end do
      
    end if
    
    ! Temp vectors
    if (associated(rsolver%p_RspaceVectors)) then
      do ilev = 1,ispaceLevel
        
        p_rfeSpaceLevel => rsolver%p_rspaceTimeHierarchy%p_rfeHierarchy%p_rfeSpaces(ilev)
        p_rspaceDiscr => p_rfeSpaceLevel%p_rdiscretisation
        
        call lsysbl_createVecBlockByDiscr (p_rspaceDiscr,rsolver%p_RspaceVectors(ilev),.false.)
        
      end do
      
    end if
  
    ! Spatial solver
    if (associated(rsolver%p_rspaceSolver)) then
      
      select case (rsolver%p_rspaceSolver%calgorithm)
      case (LINSOL_ALG_UMFPACK4,&
            LINSOL_ALG_MULTIGRID2)
        ! Attach the temp matrices to the linear solver.
        call linsol_setMatrices (rsolver%p_rspaceSolver,rsolver%p_RspaceMatrices)
      end select
      
    end if
    
    ! Reset the total solver time
    call stat_clearTimer (rsolver%rtotalTime)

  end subroutine

  ! ***************************************************************************

  recursive subroutine stls_doneData (rsolver)
  
  ! Some cleanup after solving
  
  ! Solver structure 
  type(t_spacetimelinsol), intent(inout) :: rsolver
  
    ! local variables
    integer :: ilev,i
  
    ! Solver dependent cleanup
    select case (rsolver%csolverType)
      
    case (STLS_TYPE_JACOBI)
      
      ! Deallocate temp vectors.
      call lsysbl_releaseVector (rsolver%rspaceTemp1)
    
    case (STLS_TYPE_FBSIM)
      
      ! Deallocate temp vectors.
      call sptivec_releaseVector (rsolver%rspaceTimeTemp1)
      call lsysbl_releaseVector (rsolver%rspaceTemp3)
      call lsysbl_releaseVector (rsolver%rspaceTemp2)
      call lsysbl_releaseVector (rsolver%rspaceTemp1)

    case (STLS_TYPE_FBGS)
      
      ! Deallocate temp vectors.
      call sptivec_releaseVector (rsolver%rspaceTimeTemp1)
      call lsysbl_releaseVector (rsolver%rspaceTemp3)
      call lsysbl_releaseVector (rsolver%rspaceTemp2)
      call lsysbl_releaseVector (rsolver%rspaceTemp1)
    
    case (STLS_TYPE_DEFCORR)
      call sptivec_releaseVector (rsolver%rspaceTimeTemp2)
      call sptivec_releaseVector (rsolver%rspaceTimeTemp1)
    
      if (associated(rsolver%p_rpreconditioner)) then
        call stls_doneData(rsolver%p_rpreconditioner)
      end if
    
    case (STLS_TYPE_BICGSTAB)

      do i=1,7
        call sptivec_releaseVector (rsolver%p_rtempVectors(i))
      end do
      deallocate (rsolver%p_rtempVectors)
    
      call lsysbl_releaseVector (rsolver%rspaceTemp1)
    
      if (associated(rsolver%p_rpreconditioner)) then
        call stls_doneData(rsolver%p_rpreconditioner)
      end if

    case (STLS_TYPE_MULTIGRID)
    
      if (associated(rsolver%p_rcoarseGridSolver)) then
        call stls_doneData(rsolver%p_rcoarseGridSolver)
      end if

      if (associated(rsolver%p_RpreSmoothers)) then
        do ilev = 2,rsolver%ilevel
          call stls_doneData(rsolver%p_RpreSmoothers(ilev))
        end do
      end if

      if (associated(rsolver%p_RpostSmoothers)) then
        do ilev = 2,rsolver%ilevel
          call stls_doneData(rsolver%p_RpostSmoothers(ilev))
        end do
      end if

      if (associated(rsolver%p_Rpreconditioners)) then
        do ilev = 2,rsolver%ilevel
          call stls_doneData(rsolver%p_Rpreconditioners(ilev))
        end do
      end if
    
      ! Deallocate temp vectors on all levels.
      do ilev = 1,rsolver%ilevel
        call sptivec_releaseVector (rsolver%p_Rvectors4(ilev))
        call sptivec_releaseVector (rsolver%p_Rvectors3(ilev))
        call sptivec_releaseVector (rsolver%p_Rvectors2(ilev))
        call sptivec_releaseVector (rsolver%p_Rvectors1(ilev))
      end do
      deallocate(rsolver%p_Rvectors4)
      deallocate(rsolver%p_Rvectors3)
      deallocate(rsolver%p_Rvectors2)
      deallocate(rsolver%p_Rvectors1)
      
      if (associated (rsolver%p_Rvectors5)) then
        do ilev = 1,rsolver%ilevel
          call sptivec_releaseVector (rsolver%p_Rvectors5(ilev))
        end do
        deallocate(rsolver%p_Rvectors5)
      end if
    
    end select
    
    if (associated(rsolver%p_RspaceMatrices)) then
      ! Release matrices
      do ilev = 1,size(rsolver%p_RspaceMatrices)
        call lsysbl_releaseMatrix(rsolver%p_RspaceMatrices(ilev))
      end do
    end if

    if (associated(rsolver%p_RspaceVectors)) then
      ! Release matrices
      do ilev = 1,size(rsolver%p_RspaceVectors)
        call lsysbl_releaseVector(rsolver%p_RspaceVectors(ilev))
      end do
    end if

  end subroutine

  ! ***************************************************************************

  recursive subroutine stls_precondDefect (rsolver, rd)
  
  ! General preconditioning to a defect vector rd
  
  ! Solver structure
  type(t_spacetimelinsol), intent(inout) :: rsolver
  
  ! Defect vector to apply preconditioning to.
  type(t_spaceTimeVector), intent(inout) :: rd
  
    ! By default, the solver worked.
    rsolver%csolverStatus = 0
    
    ! Compute the time for solving
    call stat_startTimer(rsolver%rtotalTime)

    select case (rsolver%csolverType)
    case (STLS_TYPE_DEFCORR)
      call stls_precondDefCorr (rsolver, rd)
    case (STLS_TYPE_BICGSTAB)
      call stls_precondBiCGStab (rsolver, rd)
    case (STLS_TYPE_JACOBI)
      call stls_precondBlockJacobi (rsolver, rd)
    case (STLS_TYPE_FBSIM)
      call stls_precondBlockFBSIM (rsolver, rd)
    case (STLS_TYPE_FBGS)
      call stls_precondBlockFBGS (rsolver, rd)
    case (STLS_TYPE_MULTIGRID)
      call stls_precondMultigrid (rsolver, rd)
    end select
    
    call stat_stopTimer(rsolver%rtotalTime)
  
  end subroutine

  ! ***************************************************************************

  recursive subroutine stls_solveAdaptively (rsolver, rx,rb,rd)
  
  ! Solves the system Ax=b
  
  ! Solver structure
  type(t_spacetimelinsol), intent(inout) :: rsolver

  ! Initial solution vector. Is overwritten by the new solution
  type(t_spaceTimeVector), intent(inout) :: rx
  
  ! RHS vector.
  type(t_spaceTimeVector), intent(inout) :: rb

  ! Temporary vector.
  type(t_spaceTimeVector), intent(inout) :: rd
  
    ! Create the defect
    call sptivec_copyVector (rb,rd)
    call stmv_matvec (rsolver%p_rmatrix,rx, rd, -1.0_DP, 1.0_DP)
    call spop_applyBC (rsolver%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, rd)

    ! Apply preconditioning
    call stls_precondDefect(rsolver,rd)
    
    ! Add to the new solution, that's it.
    if (rsolver%csolverStatus .ne. 2) then
      call sptivec_vectorLinearComb (rd,rx,1.0_DP,1.0_DP)
    end if
    ! else: Divergence.

  end subroutine

  ! ***************************************************************************

  subroutine stls_initDefCorr (rsolver,rspaceTimeHierarchy,ilevel,rpreconditioner)
  
  ! Initialise a defect correction solver.
  
  ! Solver structure to be initialised
  type(t_spacetimelinsol), intent(out) :: rsolver
  
  ! Underlying space-time hierarchy
  type(t_spacetimeHierarchy), intent(in), target :: rspaceTimeHierarchy
  
  ! Level of the solver  in the space-time hierarchy
  integer, intent(in) :: ilevel
  
  ! OPTIONAL: Solver structure of a preconditioner
  type(t_spacetimelinsol), intent(in), target, optional :: rpreconditioner

    ! Basic initialisation
    call stls_init(rsolver,STLS_TYPE_DEFCORR,rspaceTimeHierarchy,&
        ilevel,.false.,.false.,rpreconditioner)
  
  end subroutine

  ! ***************************************************************************

  recursive subroutine stls_precondDefCorr (rsolver, rd)
  
  ! General preconditioning to a defect vector rd
  
  ! Solver structure
  type(t_spacetimelinsol), intent(inout) :: rsolver
  
  ! Defect vector to apply preconditioning to.
  type(t_spaceTimeVector), intent(inout) :: rd
  
    ! local variables
    integer :: ite
    real(DP) :: dresInit, dresCurrent, drho, dresMeanDev, dresMeanLast
    real(DP), dimension(3) :: DresLast
    
    ! One temp vector is used for the residuum, one for the solution.
    call sptivec_clearVector (rsolver%rspaceTimeTemp1)
    
    ! Initial residuum.
    call sptivec_copyVector (rd,rsolver%rspaceTimeTemp2)
    call spop_applyBC (rsolver%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, rsolver%rspaceTimeTemp2)
    dresInit = sptivec_vectorNorm (rsolver%rspaceTimeTemp2,LINALG_NORML2)
    dresCurrent = dresInit
    DresLast(:) = 0.0_DP
    
    do ite = 1,rsolver%nmaxiterations+1
      if (rsolver%ioutputlevel .ge. 2) then
        call output_line("Space-Time DefCorr: Step "//trim(sys_siL(ite-1,10))//", Residuum = "//&
            trim(sys_sdEP(dresCurrent,15,7)))
      end if
      
      if (ite .gt. rsolver%nmaxiterations) then
        ! Not completed.
        rsolver%csolverStatus = 1
        exit
      end if

      if (ite .ge. rsolver%nminiterations+1) then
        ! Check the stopping criterion.
        select case (rsolver%iresCheck)
        case (1)
          if ((dresCurrent .lt. rsolver%depsAbs) .or. &
              (dresCurrent .lt. rsolver%depsRel*dresInit)) then
            exit
          end if
        case (2)
          if ((dresCurrent .lt. rsolver%depsAbs) .and. &
              (dresCurrent .lt. rsolver%depsRel*dresInit)) then
            exit
          end if
        end select
        
        ! Calculate mean deviation in the last iterations.
        dresMeanDev = sum(abs(dresCurrent-DresLast(:))) / min(ite,size(DresLast))
        dresMeanLast = sum(DresLast(:)) / min(ite,size(DresLast))
        
        if ( (dresMeanDev .lt. rsolver%depsRelDiff*dresMeanLast) .or. &
            (dresMeanDev .lt. rsolver%depsAbsDiff) ) then
          exit
        end if
      end if
      
      if ((.not. (dresCurrent .lt. rsolver%ddivAbs)) .or. &
          (.not. (dresCurrent .lt. rsolver%ddivRel*dresInit))) then
        if (rsolver%ioutputlevel .ge. 1) then
          call output_lbrk()
          call output_line("Space-Time DefCorr: Divergence detected! Iteration stopped.");
        end if
        rsolver%csolverStatus = 2
        exit
      end if
      
      ! Defect preconditioning
      if (associated (rsolver%p_rpreconditioner)) then
        call stls_precondDefect(rsolver%p_rpreconditioner,rsolver%rspaceTimeTemp2)
        
        if (rsolver%p_rpreconditioner%csolverStatus .eq. 2) then
          ! Divergence. Stop immediately.
          rsolver%csolverStatus = 2
          exit
        end if
      end if
      
      ! Add to the new solution.
      call sptivec_vectorLinearComb (rsolver%rspaceTimeTemp2,rsolver%rspaceTimeTemp1,1.0_DP,1.0_DP)
      
      ! New defect
      call sptivec_copyVector (rd,rsolver%rspaceTimeTemp2)
      call stmv_matvec (rsolver%p_rmatrix, &
          rsolver%rspaceTimeTemp1, rsolver%rspaceTimeTemp2, -1.0_DP, 1.0_DP)
      call spop_applyBC (rsolver%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, rsolver%rspaceTimeTemp2)
      
      DresLast = eoshift(DresLast,shift=-1,boundary=dresCurrent)
      dresCurrent = sptivec_vectorNorm (rsolver%rspaceTimeTemp2,LINALG_NORML2)
      
    end do
    
    ! Convergence rate
    if (ite .gt. 0) then
      if (dresInit .eq. 0.0_DP) then
        drho = dresCurrent ** (1.0_DP/real(ite-1,DP))
      else
        drho = (dresCurrent / dresInit) ** (1.0_DP/real(ite-1,DP))
      end if
    else
      drho = 0.0_DP
    end if
    
    if (rsolver%ioutputlevel .ge. 2) then
      call output_lbrk()
      call output_line("Initial Residual: "//trim(sys_sdEP(dresInit,15,7)))
      call output_line("Final Residual  : "//trim(sys_sdEP(dresCurrent,15,7)))
      call output_line("Iterations      : "//trim(sys_siL(ite-1,10)))
      call output_line("Convergence rate: "//trim(sys_sdEP(drho,15,7)))
    else if (rsolver%ioutputlevel .ge. 1) then
      call output_line("Space-Time DefCorr: Iterations/Conv. rate: "//&
          trim(sys_siL(ite-1,10))//"/"//trim(sys_sdEP(drho,15,7)))
    end if
    
    ! Put the (weighted) solution into rd as returm value
    call sptivec_vectorLinearComb (rsolver%rspaceTimeTemp1,rd,rsolver%domega,0.0_DP)
  
  end subroutine

  ! ***************************************************************************

  subroutine stls_initBiCGStab (rsolver,rspaceTimeHierarchy,ilevel,rpreconditioner)
  
  ! Initialise a defect correction solver.
  
  ! Solver structure to be initialised
  type(t_spacetimelinsol), intent(out) :: rsolver
  
  ! Underlying space-time hierarchy
  type(t_spacetimeHierarchy), intent(in), target :: rspaceTimeHierarchy
  
  ! Level of the solver  in the space-time hierarchy
  integer, intent(in) :: ilevel
  
  ! OPTIONAL: Solver structure of a preconditioner
  type(t_spacetimelinsol), intent(in), target, optional :: rpreconditioner

    ! Basic initialisation
    call stls_init(rsolver,STLS_TYPE_BICGSTAB,rspaceTimeHierarchy,&
        ilevel,.false.,.false.,rpreconditioner)
  
  end subroutine

  ! ***************************************************************************

  recursive subroutine stls_precondBiCGStab (rsolver, rd)
  
  ! General preconditioning to a defect vector rd
  
  ! Solver structure
  type(t_spacetimelinsol), intent(inout) :: rsolver
  
  ! Defect vector to apply preconditioning to.
  type(t_spaceTimeVector), intent(inout) :: rd
  
  ! local variables
  real(DP) :: dalpha,dbeta,domega0,domega1,domega2,dresreal
  real(DP) :: drho1,drho0,dresscale,dresunprec
  integer :: ite
  logical :: bstopOnRealResiduum
  real(DP) :: dresInit, dresCurrent, drho, dresFinal,dresMeanDev,dresMeanLast
  real(DP), dimension(3) :: DresLast

  ! The system matrix
  type(t_spaceTimeMatrix), pointer :: p_rmatrix
  
  ! Minimum number of iterations, print-sequence for residuals
  integer :: nminIterations
  
  ! Whether to filter/prcondition
  logical bprec
  
  ! Pointers to temporary vectors - named for easier access
  type(t_spaceTimeVector), pointer :: p_DR,p_DR0,p_DP,p_DPA,p_DSA,p_rx,p_rres
  type(t_spacetimelinsol), pointer :: p_rpreconditioner
  
    ! Solve the system!
  
    ! Status reset
    rsolver%csolverStatus = 0
    
    ! Getch some information
    p_rmatrix => rsolver%p_rmatrix
    bstopOnRealResiduum = rsolver%bstopOnRealResiduum

    ! Check the parameters
    if (rd%NEQtime .eq. 0) then
    
      ! Parameters wrong
      rsolver%csolverStatus = 2
      return
    end if

    ! Minimum number of iterations
 
    nminIterations = max(rsolver%nminIterations,0)
      
    ! Use preconditioning? Filtering?

    bprec = associated(rsolver%p_rpreconditioner)
    
    ! Set pointers to the temporary vectors
    p_DR   => rsolver%p_RtempVectors(1)
    p_DR0  => rsolver%p_RtempVectors(2)
    p_DP   => rsolver%p_RtempVectors(3)
    p_DPA  => rsolver%p_RtempVectors(4)
    p_DSA  => rsolver%p_RtempVectors(5)
    p_rx   => rsolver%p_RtempVectors(6)
    
    p_rres => rsolver%p_RtempVectors(7)
    
    if (bprec) then
      p_rpreconditioner => rsolver%p_rpreconditioner
    end if
    
    ! rd is our RHS. p_rx points to a new vector which will be our
    ! iteration vector. At the end of this routine, we replace
    ! rd by p_rx.
    ! Clear our iteration vector p_rx.
    call sptivec_clearVector (p_rx)
      
    ! Initialize used vectors with zero
      
    call sptivec_clearVector(p_DP)
    call sptivec_clearVector(p_DPA)
    
    ! Initialise the iteration vector with zero.

    ! Initialization

    drho0  = 1.0_DP
    dalpha = 1.0_DP
    domega0 = 1.0_DP

    ! Copy our RHS rd to p_DR. As the iteration vector is 0, this
    ! is also our initial defect.

    call sptivec_copyVector(rd,p_DR)
    
    ! Filter the defect for boundary conditions in space and time.
    call spop_applyBC (rsolver%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, p_DR)

    ! Get the norm of the residuum.
    !
    ! If we measure the real residuum, remember the current residuum
    ! in p_rres.
    if (bstopOnRealResiduum) then
      call sptivec_copyVector(rd,p_rres)
      dresCurrent = sptivec_vectorNorm (p_rres,LINALG_NORML2)
    else
      dresCurrent = sptivec_vectorNorm (p_DR,LINALG_NORML2)
    end if
    
    if (.not.((dresCurrent .ge. 1E-99_DP) .and. &
              (dresCurrent .le. 1E99_DP))) dresCurrent = 0.0_DP

    if (bprec) then
      ! Perform preconditioning with the assigned preconditioning
      ! solver structure.
      call stls_precondDefect (p_rpreconditioner,p_DR)
      
      if (p_rpreconditioner%csolverStatus .eq. 2) then
        ! Divergence. Stop immediately.
        rsolver%csolverStatus = 2
        return
      end if
      
      if (.not. bstopOnRealResiduum) then
        ! We scale the absolute stopping criterion by the difference
        ! between the preconditioned and unpreconditioned defect --
        ! to encounter the difference in the residuals.
        ! This is of course an approximation to 
        dresunprec = dresCurrent
        dresCurrent = sptivec_vectorNorm (p_DR,LINALG_NORML2)
        
        if (rsolver%ioutputLevel .ge. 2) then
          call output_line ('Space-Time-BiCGStab: Iteration '// &
              trim(sys_siL(0,10))//',  !!RES(unscaled)!! = '//&
              trim(sys_sdEL(dresCurrent,15)) )
        end if
        
        if (.not.((dresCurrent .ge. 1E-99_DP) .and. &
                  (dresCurrent .le. 1E99_DP))) dresCurrent = 1.0_DP
        dresscale = dresunprec / dresCurrent
      else
      
        if (rsolver%ioutputLevel .ge. 3) then
          dresreal = sptivec_vectorNorm (p_DR,LINALG_NORML2)
          call output_line ('Space-Time-BiCGStab: Iteration '// &
              trim(sys_siL(0,10))//',  !!RES(precond)!! = '//&
              trim(sys_sdEL(dresreal,15)) )
        end if

        dresscale = 1.0_DP
        
      end if
    else
      dresscale = 1.0_DP
    end if
    
    ! Initialize starting residuum
      
    dresInit = dresCurrent
    dresLast(:) = 0.0_DP
    dresFinal = dresCurrent

    ! Check if out initial defect is zero. This may happen if the filtering
    ! routine filters "everything out"!
    ! In that case we can directly stop our computation.

    if ( dresInit .lt. SYS_EPSREAL ) then
     
      ! final defect is 0, as initialised in the output variable above

      call sptivec_clearVector(p_rx)
      ite = 0
      dresFinal = dresCurrent
          
    else

      if (rsolver%ioutputLevel .ge. 2) then
        if (bprec .and. .not. bstopOnRealResiduum) then
          call output_line ('Space-Time-BiCGStab: Iteration '// &
              trim(sys_siL(0,10))//',  !!RES(scaled)!! = '//&
              trim(sys_sdEL(dresInit*dresscale,15)) )
              
          if (rsolver%ioutputLevel .ge. 3) then
          
            ! Compute the real residual.
            call sptivec_copyVector (rd,p_rres)
            call stmv_matvec (rsolver%p_rmatrix, &
                p_rx, p_rres, -1.0_DP, 1.0_DP)
            call spop_applyBC (rsolver%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, p_rres)
            dresreal = sptivec_vectorNorm (p_rres,LINALG_NORML2)
            
            call output_line ('Space-Time-BiCGStab: Iteration '// &
                trim(sys_siL(ITE,10))//',  !!RES(real)!! = '//&
                trim(sys_sdEL(dresreal,15)) )
            
          end if
        else
          call output_line ('Space-Time-BiCGStab: Iteration '// &
              trim(sys_siL(0,10))//',  !!RES!! = '//&
              trim(sys_sdEL(dresInit,15)) )
              
        end if
      end if

      call sptivec_copyVector(p_DR,p_DR0)

      ! Perform at most nmaxIterations loops to get a new vector

      do ite = 1,rsolver%nmaxIterations
      
        if (rsolver%niteReinit .gt. 0) then
          if ((ite .gt. 1) .and. (mod(ite,rsolver%niteReinit) .eq. 1)) then
            if (rsolver%ioutputLevel .ge. 2) then
              call output_line ('Space-Time-BiCGStab: Reinitialisation.')
            end if

            ! Reinitialisation. Reompute the residual and reset dr/dp.
            call sptivec_copyVector (rd,p_DR)
            
            call stmv_matvec (rsolver%p_rmatrix, &
                p_rx, p_DR, -1.0_DP, 1.0_DP)
            
            ! Filter the defect for boundary conditions in space and time.
            call spop_applyBC (rsolver%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, p_DR)
            
            if (bprec) then
              ! Perform preconditioning with the assigned preconditioning
              ! solver structure.
              call stls_precondDefect (p_rpreconditioner,p_DR)
              
              if (p_rpreconditioner%csolverStatus .eq. 2) then
                ! Divergence. Stop immediately.
                rsolver%csolverStatus = 2
                exit
              end if
              
            end if
            
            call sptivec_copyVector(p_DR,p_DR0)

            call sptivec_clearVector(p_DP)
            call sptivec_clearVector(p_DPA)
            
            drho0  = 1.0_DP
            dalpha = 1.0_DP
            domega0 = 1.0_DP
            
          end if
        end if

        drho1 = sptivec_scalarProduct (p_DR0,p_DR) 

        if (drho0*domega0 .eq. 0.0_DP) then
          ! Should not happen
          if (rsolver%ioutputLevel .ge. 2) then
            call output_line ('Space-Time-BiCGStab: Iteration prematurely stopped! '//&
                 'Correction vector is zero!')
          end if

          ! Some tuning for the output, then cancel.

          rsolver%csolverStatus = -1
          exit
          
        end if

        dbeta=(drho1*dalpha)/(drho0*domega0)
        drho0 = drho1

        call sptivec_vectorLinearComb (p_DR ,p_DP,1.0_DP,dbeta)
        call sptivec_vectorLinearComb (p_DPA ,p_DP,-dbeta*domega0,1.0_DP)

        ! Filter the defect for boundary conditions in space and time.
        call spop_applyBC (rsolver%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, p_DP)

        call stmv_matvec (rsolver%p_rmatrix, &
            p_DP, p_DPA, 1.0_DP, 0.0_DP)
    
        call spop_applyBC (rsolver%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, p_DPA)
    
        if (bprec) then
          ! Perform preconditioning with the assigned preconditioning
          ! solver structure.
          call stls_precondDefect (p_rpreconditioner,p_DPA)
          
          if (p_rpreconditioner%csolverStatus .eq. 2) then
            ! Divergence. Stop immediately.
            rsolver%csolverStatus = 2
            exit
          end if
          
        end if

        dalpha = sptivec_scalarProduct (p_DR0,p_DPA)
        
        if (dalpha .eq. 0.0_DP) then
          ! We are below machine exactness - we can't do anything more...
          ! May happen with very small problems with very few unknowns!
          if (rsolver%ioutputLevel .ge. 2) then
            call output_line ('Space-Time-BiCGStab: Convergence failed, ALPHA=0!')
            rsolver%csolverStatus = -2
            exit
          end if
        end if
        
        dalpha = drho1/dalpha

        call sptivec_vectorLinearComb (p_DPA,p_DR,-dalpha,1.0_DP)

        call stmv_matvec (rsolver%p_rmatrix, &
            p_DR,p_DSA, 1.0_DP, 0.0_DP)
                
        ! Filter the defect for boundary conditions in space and time.
        call spop_applyBC (rsolver%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, p_DSA)

        if (bprec) then
          ! Perform preconditioning with the assigned preconditioning
          ! solver structure.
          call stls_precondDefect (p_rpreconditioner,p_DSA)
          
          if (p_rpreconditioner%csolverStatus .eq. 2) then
            ! Divergence. Stop immediately.
            rsolver%csolverStatus = 2
            exit
          end if
          
        end if
        
        domega1 = sptivec_scalarProduct (p_DSA,p_DR)
        domega2 = sptivec_scalarProduct (p_DSA,p_DSA)
        
        if (domega1 .eq. 0.0_DP) then
          domega0 = 0.0_DP
        else
          if (domega2 .eq. 0.0_DP) then
            if (rsolver%ioutputLevel .ge. 2) then
              call output_line ('Space-Time-BiCGStab: Convergence failed: omega=0!')
              rsolver%csolverStatus = -2
              exit
            end if
          end if
          domega0 = domega1/domega2
        end if

        call sptivec_vectorLinearComb (p_DP ,p_rx,dalpha,1.0_DP)
        call sptivec_vectorLinearComb (p_DR ,p_rx,domega0,1.0_DP)
        
        call sptivec_vectorLinearComb (p_DSA,p_DR,-domega0,1.0_DP)

        call spop_applyBC (rsolver%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, p_DSA)

        ! Get the norm of the new (final?) residuum
        if (bstopOnRealResiduum) then
          ! Calculate the real residuum.
          call sptivec_copyVector (rd,p_rres)
          call stmv_matvec (rsolver%p_rmatrix, &
              p_rx, p_rres, -1.0_DP, 1.0_DP)
          call spop_applyBC (rsolver%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, p_rres)
          dresCurrent = sptivec_vectorNorm (p_rres,LINALG_NORML2)
        else
          ! Take the preconditioned residuum
          dresCurrent = sptivec_vectorNorm (p_DR,LINALG_NORML2)
        end if
     
        DresLast = eoshift(DresLast,shift=-1,boundary=dresFinal)
        dresFinal = dresCurrent

        if (ite .ge. rsolver%nminiterations+1) then
          ! Check the stopping criterion.
          select case (rsolver%iresCheck)
          case (1)
            if ((dresCurrent .lt. rsolver%depsAbs) .or. &
                (dresCurrent .lt. rsolver%depsRel*dresInit)) then
              exit
            end if
          case (2)
            if ((dresCurrent .lt. rsolver%depsAbs) .and. &
                (dresCurrent .lt. rsolver%depsRel*dresInit)) then
              exit
            end if
          end select
          
          ! Calculate mean deviation in the last iterations.
          dresMeanDev = sum(abs(dresCurrent-DresLast(:))) / min(ite,size(DresLast))
          dresMeanLast = sum(DresLast(:)) / min(ite,size(DresLast))
          
          if ( (dresMeanDev .lt. rsolver%depsRelDiff*dresMeanLast) .or. &
              (dresMeanDev .lt. rsolver%depsAbsDiff) ) then
            exit
          end if
        end if
        
        if ((.not. (dresCurrent .lt. rsolver%ddivAbs)) .or. &
            (.not. (dresCurrent .lt. rsolver%ddivRel*dresInit))) then
          if (rsolver%ioutputlevel .ge. 1) then
            call output_lbrk()
            call output_line("Space-Time-BiCGStab: Divergence detected! Iteration stopped.");
          end if
          rsolver%csolverStatus = 2
          exit
        end if

        ! print out the current residuum

        if (rsolver%ioutputLevel .ge. 2) then
          if (bprec .and. .not. bstopOnRealResiduum) then
            call output_line ('Space-Time-BiCGStab: Iteration '// &
                trim(sys_siL(ITE,10))//',  !!RES(scaled)!! = '//&
                trim(sys_sdEL(dresFinal*dresscale,15)) )

            if (rsolver%ioutputLevel .ge. 3) then
            
              ! Compute the real residual.
              call sptivec_copyVector (rd,p_rres)
              call stmv_matvec (rsolver%p_rmatrix, &
                  p_rx, p_rres, -1.0_DP, 1.0_DP)
              call spop_applyBC (rsolver%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, p_rres)
              dresreal = sptivec_vectorNorm (p_rres,LINALG_NORML2)
              
              call output_line ('Space-Time-BiCGStab: Iteration '// &
                  trim(sys_siL(ITE,10))//',  !!RES(real)!! = '//&
                  trim(sys_sdEL(dresreal,15)) )
              
            end if

          else
            if (bstopOnRealResiduum .and. rsolver%ioutputLevel .ge. 3) then
              dresreal = sptivec_vectorNorm (p_DR,LINALG_NORML2)
              call output_line ('Space-Time-BiCGStab: Iteration '// &
                  trim(sys_siL(ITE,10))//',  !!RES(precond)!! = '//&
                  trim(sys_sdEL(dresreal,15)) )
            end if                

            call output_line ('Space-Time-BiCGStab: Iteration '// &
                trim(sys_siL(ITE,10))//',  !!RES!! = '//&
                trim(sys_sdEL(dresFinal,15)) )
          end if
        end if

      end do

      ! Set ITE to NIT to prevent printing of "NIT+1" of the loop was
      ! completed

      if (ite .gt. rsolver%nmaxIterations) &
        ite = rsolver%nmaxIterations

      ! Finish - either with an error or if converged.
      ! Print the last residuum.

      if ((rsolver%ioutputLevel .ge. 2) .and. &
          (ite .ge. 1) .and. (ITE .lt. rsolver%nmaxIterations) .and. &
          (rsolver%csolverStatus .ge. 0)) then
          
        if (bprec .and. .not. bstopOnRealResiduum) then
          call output_line ('Space-Time-BiCGStab: Iteration '// &
              trim(sys_siL(ITE,10))//',  !!RES(scaled)!! = '//&
              trim(sys_sdEL(dresFinal*dresscale,15)) )

          if (rsolver%ioutputLevel .ge. 3) then
          
            ! Compute the real residual.
            call sptivec_copyVector (rd,p_rres)
            call stmv_matvec (rsolver%p_rmatrix, &
                p_rx, p_rres, -1.0_DP, 1.0_DP)
            call spop_applyBC (rsolver%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, p_rres)
            dresreal = sptivec_vectorNorm (p_rres,LINALG_NORML2)
            
            call output_line ('Space-Time-BiCGStab: Iteration '// &
                trim(sys_siL(ITE,10))//',  !!RES(real)!! = '//&
                trim(sys_sdEL(dresreal,15)) )
            
          end if

        else
          if (bstopOnRealResiduum .and. rsolver%ioutputLevel .ge. 3) then
            dresreal = sptivec_vectorNorm (p_DR,LINALG_NORML2)
            call output_line ('Space-Time-BiCGStab: Iteration '// &
                trim(sys_siL(ITE,10))//',  !!RES(precond)!! = '//&
                trim(sys_sdEL(dresreal,15)) )
          end if                

          call output_line ('Space-Time-BiCGStab: Iteration '// &
              trim(sys_siL(ITE,10))//',  !!RES!! = '//&
              trim(sys_sdEL(dresFinal,15)) )
        end if
      end if

    end if

    ! Overwrite our previous RHS by the new correction vector p_rx.
    ! This completes the preconditioning.
    call sptivec_copyVector (p_rx,rd)
    call sptivec_scaleVector (rd,rsolver%domega)
      
    ! Don't calculate anything if the final residuum is out of bounds -
    ! would result in NaN's,...
      
    if (dresFinal .lt. 1E99_DP) then
    
      ! If the initial defect was zero, the solver immediately
      ! exits - and so the final residuum is zero and we performed
      ! no steps; so the resulting convergence rate stays zero.
      ! In the other case the convergence rate computes as
      ! (final defect/initial defect) ** 1/nit :

      drho = 0.0_DP
      if ((dresFinal .gt. SYS_EPSREAL) .and. (ite .gt. 0)) then
        drho = (dresFinal / dresInit) ** &
               (1.0_DP/real(ite,DP))
      end if

      if (rsolver%ioutputLevel .ge. 2) then
        call output_lbrk()
        call output_line ('Space-Time-BiCGStab statistics:')
        call output_lbrk()
        call output_line ('Iterations              : '//&
             trim(sys_siL(ite,10)) )
        call output_line ('!!INITIAL RES!!         : '//&
             trim(sys_sdEL(dresInit,15)) )
        call output_line ('!!RES!!                 : '//&
             trim(sys_sdEL(dresFinal,15)) )
        if (dresInit .gt. SYS_EPSREAL) then     
          call output_line ('!!RES!!/!!INITIAL RES!! : '//&
            trim(sys_sdEL(dresFinal / dresInit,15)) )
        else
          call output_line ('!!RES!!/!!INITIAL RES!! : '//&
               trim(sys_sdEL(0.0_DP,15)) )
        end if
        call output_lbrk ()
        call output_line ('Rate of convergence     : '//&
             trim(sys_sdEL(drho,15)) )

      end if

      if (rsolver%ioutputLevel .eq. 1) then
        call output_line (&
              'Space-Time-BiCGStab: Iterations/Rate of convergence: '//&
              trim(sys_siL(ite,10))//' /'//&
              trim(sys_sdEL(drho,15)) )
      end if
      
    else
      ! DEF=Infinity; RHO=Infinity, set to 1
      drho = 1.0_DP
    end if  

  end subroutine

  ! ***************************************************************************

  subroutine stls_initBlockJacobi (rsolver,rspaceTimeHierarchy,ilevel,domega,&
      rspaceSolverParams,RmatVecTempl)
  
  ! Initialise a block Jacobi correction solver.
  
  ! Solver structure to be initialised
  type(t_spacetimelinsol), intent(out) :: rsolver
  
  ! Underlying space-time hierarchy
  type(t_spacetimeHierarchy), intent(in), target :: rspaceTimeHierarchy
  
  ! Level of the solver  in the space-time hierarchy
  integer, intent(in) :: ilevel

  ! Damping parameter
  real(DP), intent(in) :: domega

  ! Structure defining the linear solver.
  type(t_spaceSolverParams), intent(in) :: rspaceSolverParams

  ! OPTINAL: Spatial matrix templates. Necessary if a matrix-based
  ! preconditioner in space is used.
  type(t_matvecTemplates), dimension(:), intent(in), target :: RmatVecTempl

    ! Basic initialisation
    call stls_init(rsolver,STLS_TYPE_JACOBI,rspaceTimeHierarchy,ilevel,&
        .false.,.false.,rspaceSolverParams=rspaceSolverParams,RmatVecTempl=RmatVecTempl)
        
    rsolver%domega = domega
  
  end subroutine

  ! ***************************************************************************

  recursive subroutine stls_precondBlockJacobi (rsolver, rd)
  
  ! General preconditioning to a defect vector rd
  
  ! Solver structure
  type(t_spacetimelinsol), intent(inout) :: rsolver
  
  ! Defect vector to apply preconditioning to.
  type(t_spaceTimeVector), intent(inout) :: rd

    ! local variables
    integer :: istep, ierror, ispaceLevel, ilev
    real(DP), dimension(:), pointer :: p_Ddata
    
    ! Initialise the linear solver.
    call linsol_initStructure(rsolver%p_rspaceSolver,ierror)
    if (ierror .ne. 0) then
      call output_line ("Spatial system cannot be symbolically factorised in timestep "&
          //trim(sys_siL(istep,10)))
      call sys_halt()
    end if

    call lsysbl_getbase_double (rsolver%rspaceTemp1,p_Ddata)
    
    ! Get the space level
    call sth_getLevel(rsolver%p_rspaceTimeHierarchy,rsolver%ilevel,ispaceLevel=ispaceLevel)

    ! We loop through the timesteps and apply the spatial solver in each timestep.
    do istep = 1,rd%NEQtime
      
      ! Load the timestep.
      call sptivec_getTimestepData(rd,istep,rsolver%rspaceTemp1)
      
      do ilev = 1,ispaceLevel
        ! Assemble diagonal submatrix of that timestep on all levels
        call stmat_getSubmatrix (rsolver%p_rmatrix, ilev, istep, istep, &
            rsolver%p_RspaceMatrices(ilev))
            
        ! Apply the boundary conditions to the matrix
        call stmat_implementDefBCSubmatrix (rsolver%p_rmatrix%p_rboundaryCond, &
            ilev, rsolver%p_RspaceMatrices(ilev), istep, istep, &
            rsolver%p_RdiscreteBC(ilev))
      end do
      
      !call matio_writeBlockMatrixHR (rsolver%p_RspaceMatrices(ispaceLevel), "matrix",&
      !    .true., 0, "matrix.txt", "(E10.3)")

      ! Apply the space solver
      call linsol_initData(rsolver%p_rspaceSolver,ierror)
      if (ierror .ne. 0) then
        call output_line ("Error "//trim(sys_siL(ierror,10))&
            //". Spatial system cannot be factorised in timestep "//&
            trim(sys_siL(istep,10)),OU_CLASS_ERROR,OU_MODE_STD,'stls_precondBlockJacobi')
        call sys_halt()
      end if
      
      call linsol_precondDefect (rsolver%p_rspaceSolver,rsolver%rspaceTemp1)
      call linsol_doneData(rsolver%p_rspaceSolver,ierror)
      
      ! Scale and store back.
      call lsysbl_scaleVector (rsolver%rspaceTemp1,rsolver%domega)
      call sptivec_setTimestepData(rd,istep,rsolver%rspaceTemp1)
      
    end do
    
    call linsol_doneStructure(rsolver%p_rspaceSolver,ierror)

  end subroutine

  ! ***************************************************************************

  subroutine stls_initBlockFBSIM (rsolver,rspaceTimeHierarchy,ilevel,drelax,&
      rspaceSolverParams,RmatVecTempl)
  
  ! Initialise a block GS correction solver, working on the decoupled solution.
  
  ! Solver structure to be initialised
  type(t_spacetimelinsol), intent(out) :: rsolver
  
  ! Underlying space-time hierarchy
  type(t_spacetimeHierarchy), intent(in), target :: rspaceTimeHierarchy
  
  ! Level of the solver  in the space-time hierarchy
  integer, intent(in) :: ilevel

  ! Relaxation parameter
  real(DP), intent(in) :: drelax

  ! Structure defining the linear solver.
  type(t_spaceSolverParams), intent(in) :: rspaceSolverParams

  ! OPTINAL: Spatial matrix templates. Necessary if a matrix-based
  ! preconditioner in space is used.
  type(t_matvecTemplates), dimension(:), intent(in), target :: RmatVecTempl

    ! Basic initialisation
    call stls_init(rsolver,STLS_TYPE_FBSIM,rspaceTimeHierarchy,ilevel,&
        .false.,.false.,rspaceSolverParams=rspaceSolverParams,RmatVecTempl=RmatVecTempl)
        
    rsolver%drelax = drelax
  
  end subroutine

  ! ***************************************************************************

  recursive subroutine stls_precondBlockFBSIM (rsolver, rd)
  
  ! General preconditioning to a defect vector rd
  
  ! Solver structure
  type(t_spacetimelinsol), intent(inout) :: rsolver
  
  ! Defect vector to apply preconditioning to.
  type(t_spaceTimeVector), intent(inout) :: rd

    ! local variables
    integer :: istep, ierror,jstep,i,ispaceLevel,ilev
    real(DP), dimension(:), pointer :: p_Dx,p_Dd,p_Ddata
    
    ! Get the space level
    call sth_getLevel(rsolver%p_rspaceTimeHierarchy,rsolver%ilevel,ispaceLevel=ispaceLevel)

    ! Initialise the linear solver.
    call linsol_initStructure(rsolver%p_rspaceSolver,ierror)
    if (ierror .ne. 0) then
      call output_line ("Spatial system cannot be symbolically factorised in timestep "&
          //trim(sys_siL(istep,10)))
      call sys_halt()
    end if
    
    call lsysbl_getbase_double (rsolver%rspaceTemp1,p_Dd)
    call lsysbl_getbase_double (rsolver%rspaceTemp2,p_Dx)
    call lsysbl_getbase_double (rsolver%rspaceTemp3,p_Ddata)
    
    ! Clear the output vector.
    call sptivec_clearVector (rsolver%rspaceTimeTemp1)
    
    do jstep = 1,1
    
      ! Forward sweep.
        
      ! We loop through the timesteps and apply the spatial solver in each timestep.
      do istep = 1,rd%NEQtime
        
        ! Load the timestep.
        call sptivec_getTimestepData(rd,istep,rsolver%rspaceTemp1)
        
        if (istep .gt. 1) then
          ! Load the previous timestep.
          call sptivec_getTimestepData(rsolver%rspaceTimeTemp1,istep-1,rsolver%rspaceTemp2)
          
          ! Subtract the primal.
          call stmat_getSubmatrix (rsolver%p_rmatrix, ispaceLevel, istep, istep-1, &
              rsolver%p_RspaceMatrices(ispaceLevel))

          call stmat_implementDefBCSubmatrix (rsolver%p_rmatrix%p_rboundaryCond, &
              ispaceLevel,rsolver%p_RspaceMatrices(ispaceLevel), istep, istep-1, &
              rsolver%p_RdiscreteBC(ispaceLevel))

          call lsysbl_blockMatVec (rsolver%p_RspaceMatrices(ispaceLevel), &
              rsolver%rspaceTemp2, rsolver%rspaceTemp1, -1.0_DP, 1.0_DP)
        end if
        
        ! Subtract the diagonal, set up the preconditioner matrix on all space levels
        
        do ilev = 1,ispaceLevel
          ! Assemble diagonal submatrix of that timestep on all levels
          call stmat_getSubmatrix (rsolver%p_rmatrix, ilev, istep, istep, &
              rsolver%p_RspaceMatrices(ilev))
              
          ! Apply the boundary conditions to the matrix
          call stmat_implementDefBCSubmatrix (rsolver%p_rmatrix%p_rboundaryCond, &
              ilev, rsolver%p_RspaceMatrices(ilev), istep, istep, &
              rsolver%p_RdiscreteBC(ilev))
        end do

        ! Subtract the diagonal.
        call sptivec_getTimestepData(rsolver%rspaceTimeTemp1,istep,rsolver%rspaceTemp2)
        
        call lsysbl_blockMatVec (rsolver%p_RspaceMatrices(ispaceLevel), &
            rsolver%rspaceTemp2, rsolver%rspaceTemp1, -1.0_DP, 1.0_DP)
        
        ! Remove the primal part from the matrix
        do ilev = 1,ispaceLevel
          call stmat_reduceDiagToPrimal (rsolver%p_RspaceMatrices(ilev))
        end do
        
        ! Apply the space solver
        call linsol_initData(rsolver%p_rspaceSolver,ierror)
        if (ierror .ne. 0) then
          call output_line ("Error "//trim(sys_siL(ierror,10))&
              //". Spatial system cannot be factorised in timestep "//&
              trim(sys_siL(istep,10)),OU_CLASS_ERROR,OU_MODE_STD,'stls_precondBlockFBSIM')
              
          !call matio_writeBlockMatrixHR (rsolver%p_RspaceMatrices(ispaceLevel), &
          !    "matrix", .true., 0, "matrix.txt", "(E10.3)") 
              
          call sys_halt()
        end if
        
        call linsol_precondDefect (rsolver%p_rspaceSolver,rsolver%rspaceTemp1)
        call linsol_doneData(rsolver%p_rspaceSolver,ierror)
        
        do i=1,rsolver%rspaceTemp1%nblocks/2
          call lsyssc_vectorLinearComb (rsolver%rspaceTemp1%RvectorBlock(i),&
              rsolver%rspaceTemp2%RvectorBlock(i),rsolver%drelax,1.0_DP)
        end do
        
        call sptivec_setTimestepData(rsolver%rspaceTimeTemp1,istep,rsolver%rspaceTemp2)
        
      end do
      
      ! Backward sweep.
        
      ! We loop through the timesteps and apply the spatial solver in each timestep.
      do istep = rd%NEQtime,1,-1
        
        ! Load the timestep.
        call sptivec_getTimestepData(rd,istep,rsolver%rspaceTemp1)
        
        if (istep .lt. rd%NEQtime) then
          ! Load the previous timestep.
          call sptivec_getTimestepData(rsolver%rspaceTimeTemp1,istep+1,rsolver%rspaceTemp2)
          
          ! Subtract the primal.
          call stmat_getSubmatrix (rsolver%p_rmatrix, ispaceLevel, istep, istep+1, &
              rsolver%p_RspaceMatrices(ispaceLevel))

          call stmat_implementDefBCSubmatrix (rsolver%p_rmatrix%p_rboundaryCond, &
              ispaceLevel,rsolver%p_RspaceMatrices(ispaceLevel), istep, istep+1, &
              rsolver%p_RdiscreteBC(ispaceLevel))

          call lsysbl_blockMatVec (rsolver%p_RspaceMatrices(ispaceLevel), &
              rsolver%rspaceTemp2, rsolver%rspaceTemp1, -1.0_DP, 1.0_DP)
        end if
        
        ! Subtract the diagonal, set up the preconditioner matrix on all space levels
        
        do ilev = 1,ispaceLevel
          ! Assemble diagonal submatrix of that timestep on all levels
          call stmat_getSubmatrix (rsolver%p_rmatrix, ilev, istep, istep, &
              rsolver%p_RspaceMatrices(ilev))
              
          ! Apply the boundary conditions to the matrix
          call stmat_implementDefBCSubmatrix (rsolver%p_rmatrix%p_rboundaryCond, &
              ilev, rsolver%p_RspaceMatrices(ilev), istep, istep, &
              rsolver%p_RdiscreteBC(ilev))
        end do
        
        ! Subtract the diagonal.
        call sptivec_getTimestepData(rsolver%rspaceTimeTemp1,istep,rsolver%rspaceTemp2)
        
        call lsysbl_blockMatVec (rsolver%p_RspaceMatrices(ispaceLevel), &
            rsolver%rspaceTemp2, rsolver%rspaceTemp1, -1.0_DP, 1.0_DP)
        
        ! Remove the dual part from the matrix
        do ilev = 1,ispaceLevel
          call stmat_reduceDiagToDual (rsolver%p_RspaceMatrices(ilev))
        end do
        
        ! Apply the space solver
        call linsol_initData(rsolver%p_rspaceSolver,ierror)
        if (ierror .ne. 0) then
          call output_line ("Error "//trim(sys_siL(ierror,10))&
              //". Spatial system cannot be factorised in timestep "//&
              trim(sys_siL(istep,10)),OU_CLASS_ERROR,OU_MODE_STD,'stls_precondBlockFBSIM')
          call sys_halt()
        end if
        
        call linsol_precondDefect (rsolver%p_rspaceSolver,rsolver%rspaceTemp1)
        call linsol_doneData(rsolver%p_rspaceSolver,ierror)
        
        do i=rsolver%rspaceTemp1%nblocks/2+1,rsolver%rspaceTemp1%nblocks
          call lsyssc_vectorLinearComb (rsolver%rspaceTemp1%RvectorBlock(i),&
              rsolver%rspaceTemp2%RvectorBlock(i),rsolver%drelax,1.0_DP)
        end do
        
        call sptivec_setTimestepData(rsolver%rspaceTimeTemp1,istep,rsolver%rspaceTemp2)
      end do
      
    end do
    
    call linsol_doneStructure(rsolver%p_rspaceSolver,ierror)

    ! Write back the output vector.
    call sptivec_vectorLinearComb (rsolver%rspaceTimeTemp1,rd,rsolver%domega,0.0_DP)

  end subroutine

  ! ***************************************************************************

  subroutine stls_initBlockFBGS (rsolver,rspaceTimeHierarchy,ilevel,drelax,&
      rspaceSolverParams,icoupling,RmatVecTempl)
  
  ! Initialise a block GS correction solver, working on the coupled solution.
  
  ! Solver structure to be initialised
  type(t_spacetimelinsol), intent(out) :: rsolver
  
  ! Underlying space-time hierarchy
  type(t_spacetimeHierarchy), intent(in), target :: rspaceTimeHierarchy
  
  ! Level of the solver  in the space-time hierarchy
  integer, intent(in) :: ilevel

  ! Relaxation parameter
  real(DP), intent(in) :: drelax

  ! Structure defining the linear solver.
  type(t_spaceSolverParams), intent(in) :: rspaceSolverParams

  ! Specifies the coupling in the algorithm.
  ! =0: standard GS coupling.
  ! =1: extended GS coupling, also using the lower diagonal in the backward sweep.
  integer :: icoupling

  ! OPTINAL: Spatial matrix templates. Necessary if a matrix-based
  ! preconditioner in space is used.
  type(t_matvecTemplates), dimension(:), intent(in), target :: RmatVecTempl

    ! Basic initialisation
    call stls_init(rsolver,STLS_TYPE_FBGS,rspaceTimeHierarchy,ilevel,&
        .false.,.false.,rspaceSolverParams=rspaceSolverParams,RmatVecTempl=RmatVecTempl)
        
    rsolver%drelax = drelax
    
    rsolver%ialgoptions = icoupling
  
  end subroutine

  ! ***************************************************************************

  recursive subroutine stls_precondBlockFBGS (rsolver, rd)
  
  ! General preconditioning to a defect vector rd
  
  ! Solver structure
  type(t_spacetimelinsol), intent(inout) :: rsolver
  
  ! Defect vector to apply preconditioning to.
  type(t_spaceTimeVector), intent(inout) :: rd

    ! local variables
    integer :: istep, ierror,jstep,ispaceLevel,ilev
    real(DP), dimension(:), pointer :: p_Dx,p_Dd,p_Ddata
    
    ! Get the space level
    call sth_getLevel(rsolver%p_rspaceTimeHierarchy,rsolver%ilevel,ispaceLevel=ispaceLevel)

    ! Initialise the linear solver.
    call linsol_initStructure(rsolver%p_rspaceSolver,ierror)
    if (ierror .ne. 0) then
      call output_line ("Spatial system cannot be symbolically factorised in timestep "&
          //trim(sys_siL(istep,10)))
      call sys_halt()
    end if
    
    call lsysbl_getbase_double (rsolver%rspaceTemp1,p_Dd)
    call lsysbl_getbase_double (rsolver%rspaceTemp2,p_Dx)
    call lsysbl_getbase_double (rsolver%rspaceTemp3,p_Ddata)
    
    ! Clear the output vector.
    call sptivec_clearVector (rsolver%rspaceTimeTemp1)
    
    do jstep = 1,1
    
      ! Forward sweep.
        
      ! We loop through the timesteps and apply the spatial solver in each timestep.
      do istep = 1,rd%NEQtime
        
        ! Load the timestep.
        call sptivec_getTimestepData(rd,istep,rsolver%rspaceTemp1)
        
        if (istep .gt. 1) then
          ! Load the previous timestep.
          call sptivec_getTimestepData(rsolver%rspaceTimeTemp1,istep-1,rsolver%rspaceTemp2)
          
          ! Subtract the primal.
          call stmat_getSubmatrix (rsolver%p_rmatrix, ispaceLevel, istep, istep-1, &
              rsolver%p_RspaceMatrices(ispaceLevel))

          call stmat_implementDefBCSubmatrix (rsolver%p_rmatrix%p_rboundaryCond, &
              ispaceLevel,rsolver%p_RspaceMatrices(ispaceLevel), istep, istep-1, &
              rsolver%p_RdiscreteBC(ispaceLevel))

          call lsysbl_blockMatVec (rsolver%p_RspaceMatrices(ispaceLevel), &
              rsolver%rspaceTemp2, rsolver%rspaceTemp1, -1.0_DP, 1.0_DP)
        end if

!        if (istep .lt. rd%NEQtime) then
!          ! Load the previous timestep.
!          call sptivec_getTimestepData(rsolver%rspaceTimeTemp1,istep+1,rsolver%rspaceTemp2)
!          
!          ! Subtract the primal.
!          call stmat_getSubmatrix (rsolver%p_rmatrix, ispaceLevel, istep, istep+1, &
!              rsolver%p_RspaceMatrices(ispaceLevel))
!
!          call stmat_implementDefBCSubmatrix (rsolver%p_rmatrix%p_rboundaryCond, &
!              ispaceLevel,rsolver%p_RspaceMatrices(ispaceLevel), istep, istep+1, &
!              rsolver%p_RdiscreteBC(ispaceLevel))
!
!          call lsysbl_blockMatVec (rsolver%p_RspaceMatrices(ispaceLevel), &
!              rsolver%rspaceTemp2, rsolver%rspaceTemp1, -1.0_DP, 1.0_DP)
!        end if
        
        ! Subtract the diagonal, set up the preconditioner matrix on all levels
        
        do ilev = 1,ispaceLevel
          ! Assemble diagonal submatrix of that timestep on all levels
          call stmat_getSubmatrix (rsolver%p_rmatrix, ilev, istep, istep, &
              rsolver%p_RspaceMatrices(ilev))
              
          ! Apply the boundary conditions to the matrix
          call stmat_implementDefBCSubmatrix (rsolver%p_rmatrix%p_rboundaryCond, &
              ilev, rsolver%p_RspaceMatrices(ilev), istep, istep, &
              rsolver%p_RdiscreteBC(ilev))
        end do

        ! Subtract the diagonal.
        call sptivec_getTimestepData(rsolver%rspaceTimeTemp1,istep,rsolver%rspaceTemp2)
        
        call lsysbl_blockMatVec (rsolver%p_RspaceMatrices(ispaceLevel), &
            rsolver%rspaceTemp2, rsolver%rspaceTemp1, -1.0_DP, 1.0_DP)
        
        ! Apply the space solver
        call linsol_initData(rsolver%p_rspaceSolver,ierror)
        if (ierror .ne. 0) then
          call output_line ("Error "//trim(sys_siL(ierror,10))&
              //". Spatial system cannot be factorised in timestep "//&
              trim(sys_siL(istep,10)),OU_CLASS_ERROR,OU_MODE_STD,'stls_precondBlockFBGS')
              
          !call matio_writeBlockMatrixHR (rsolver%p_RspaceMatrices(ispaceLevel), &
          !    "matrix", .true., 0, "matrix.txt", "(E10.3)") 
              
          call sys_halt()
        end if
        
        call linsol_precondDefect (rsolver%p_rspaceSolver,rsolver%rspaceTemp1)
        call linsol_doneData(rsolver%p_rspaceSolver,ierror)
        
        call lsysbl_vectorLinearComb (rsolver%rspaceTemp1,rsolver%rspaceTemp2,rsolver%drelax,1.0_DP)
        
        call sptivec_setTimestepData(rsolver%rspaceTimeTemp1,istep,rsolver%rspaceTemp2)
        
      end do
      
      ! Backward sweep.
        
      ! We loop through the timesteps and apply the spatial solver in each timestep.
      do istep = rd%NEQtime,1,-1
        
        ! Load the timestep.
        call sptivec_getTimestepData(rd,istep,rsolver%rspaceTemp1)
        
        if (istep .gt. 1) then
          ! Load the previous timestep.
          call sptivec_getTimestepData(rsolver%rspaceTimeTemp1,istep-1,rsolver%rspaceTemp2)
          
          ! Subtract the primal.
          call stmat_getSubmatrix (rsolver%p_rmatrix, ispaceLevel, istep, istep-1, &
              rsolver%p_RspaceMatrices(ispaceLevel))

          call stmat_implementDefBCSubmatrix (rsolver%p_rmatrix%p_rboundaryCond, &
              ispaceLevel,rsolver%p_RspaceMatrices(ispaceLevel), istep, istep-1, &
              rsolver%p_RdiscreteBC(ispaceLevel))

          call lsysbl_blockMatVec (rsolver%p_RspaceMatrices(ispaceLevel), &
              rsolver%rspaceTemp2, rsolver%rspaceTemp1, -1.0_DP, 1.0_DP)
        end if

        if (rsolver%ialgoptions .eq. 1) then
          if (istep .lt. rd%NEQtime) then
            ! Load the previous timestep.
            call sptivec_getTimestepData(rsolver%rspaceTimeTemp1,istep+1,rsolver%rspaceTemp2)
            
            ! Subtract the primal.
            call stmat_getSubmatrix (rsolver%p_rmatrix, ispaceLevel, istep, istep+1, &
                rsolver%p_RspaceMatrices(ispaceLevel))

            call stmat_implementDefBCSubmatrix (rsolver%p_rmatrix%p_rboundaryCond, &
                ispaceLevel,rsolver%p_RspaceMatrices(ispaceLevel), istep, istep+1, &
                rsolver%p_RdiscreteBC(ispaceLevel))

            call lsysbl_blockMatVec (rsolver%p_RspaceMatrices(ispaceLevel), &
                rsolver%rspaceTemp2, rsolver%rspaceTemp1, -1.0_DP, 1.0_DP)
          end if
        end if
        
        do ilev = 1,ispaceLevel
          ! Assemble diagonal submatrix of that timestep on all levels
          call stmat_getSubmatrix (rsolver%p_rmatrix, ilev, istep, istep, &
              rsolver%p_RspaceMatrices(ilev))
              
          ! Apply the boundary conditions to the matrix
          call stmat_implementDefBCSubmatrix (rsolver%p_rmatrix%p_rboundaryCond, &
              ilev, rsolver%p_RspaceMatrices(ilev), istep, istep, &
              rsolver%p_RdiscreteBC(ilev))
        end do

        ! Subtract the diagonal.
        call sptivec_getTimestepData(rsolver%rspaceTimeTemp1,istep,rsolver%rspaceTemp2)
        
        call lsysbl_blockMatVec (rsolver%p_RspaceMatrices(ispaceLevel), &
            rsolver%rspaceTemp2, rsolver%rspaceTemp1, -1.0_DP, 1.0_DP)
        
        ! Apply the space solver
        call linsol_initData(rsolver%p_rspaceSolver,ierror)
        if (ierror .ne. 0) then
          call output_line ("Error "//trim(sys_siL(ierror,10))&
              //". Spatial system cannot be factorised in timestep "//&
              trim(sys_siL(istep,10)),OU_CLASS_ERROR,OU_MODE_STD,'stls_precondBlockFBGS')
          call sys_halt()
        end if
        
        call linsol_precondDefect (rsolver%p_rspaceSolver,rsolver%rspaceTemp1)
        call linsol_doneData(rsolver%p_rspaceSolver,ierror)
        
        call lsysbl_vectorLinearComb (rsolver%rspaceTemp1,rsolver%rspaceTemp2,rsolver%drelax,1.0_DP)
        
        call sptivec_setTimestepData(rsolver%rspaceTimeTemp1,istep,rsolver%rspaceTemp2)
      end do
      
    end do
    
    call linsol_doneStructure(rsolver%p_rspaceSolver,ierror)

    ! Write back the output vector.
    call sptivec_vectorLinearComb (rsolver%rspaceTimeTemp1,rd,rsolver%domega,0.0_DP)

  end subroutine

  ! ***************************************************************************

  subroutine stls_initMultigrid (rsolver,rspaceTimeHierarchy,ilevel,rspaceTimeProjection,&
      rcoarseGridSolver, RpreSmoothers, RpostSmoothers, Rpreconditioners)
  
  ! Initialise a multigrid solver.
  
  ! Solver structure to be initialised
  type(t_spacetimelinsol), intent(out) :: rsolver
  
  ! A space-time hierarchy that describes the discretisation in space and time
  ! for all levels.
  type(t_spaceTimeHierarchy), intent(in), target :: rspaceTimeHierarchy

  ! Projection hierarchy that describes how to do prolongation/restriction
  type(t_sptiProjHierarchy), intent(in), target :: rspaceTimeProjection

  ! Solver structure of the coarse grid solver.
  type(t_spacetimelinsol), intent(in), target :: rcoarseGridSolver
  
  ! OPTIONAL: Array of solver structures for the presmoother.
  ! if not present, no presmoothing is applied.
  type(t_spacetimelinsol), dimension(:), intent(in), optional, target :: RpreSmoothers

  ! OPTIONAL: Array of solver structures for the postsmoother.
  ! if not present, no presmoothing is applied.
  type(t_spacetimelinsol), dimension(:), intent(in), optional, target :: RpostSmoothers
  
  ! DEBUG!!!
  ! OPTIONAL: Array of solver structures for the preconditioners.
  type(t_spacetimelinsol), dimension(:), intent(in), optional, target :: Rpreconditioners

  ! Level of the solver in the space-time hierarchy
  integer, intent(in) :: ilevel
  
    ! Basic initialisation
    call stls_init(rsolver,STLS_TYPE_MULTIGRID,rspaceTimeHierarchy,ilevel,.false.,.true.)
    
    ! Remember the smoothers and the coarse grid solver.
    rsolver%p_rcoarseGridSolver => rcoarseGridSolver
    if (present(RpreSmoothers)) then
      rsolver%p_RpreSmoothers => RpreSmoothers
    else
      nullify(rsolver%p_RpreSmoothers)
    end if
    
    if (present(RpostSmoothers)) then
      rsolver%p_RpostSmoothers => RpostSmoothers
    else
      nullify(rsolver%p_RpostSmoothers)
    end if

    if (present(Rpreconditioners)) then
      rsolver%p_Rpreconditioners => Rpreconditioners
    else
      nullify(rsolver%p_Rpreconditioners)
    end if
    
    rsolver%p_rspaceTimeProjection => rspaceTimeProjection
        
  end subroutine

  ! ***************************************************************************
  
!<subroutine>
  
  subroutine stls_convertToSmoother (rsolver,nsmoothingSteps)
  
!<description>
  ! Converts a solver node to a smoother node. A smoother is a solver
  ! that performs a fixed number of iterations without respecting any
  ! residuum. nsmoothingSteps is the number of steps the smoother should
  ! perform.
!</description>
  
!<input>
  ! Number of steps the smoother should perform
  integer, intent(in)          :: nsmoothingSteps
!</input>
  
!<inputoutput>
  ! Solver node which should be configured as smoother.
  type(t_spacetimelinsol), intent(inout) :: rsolver
!</inputoutput>
  
!</subroutine>

    rsolver%depsRel = 0.0_DP
    rsolver%depsAbs = 0.0_DP
    rsolver%depsRelDiff = 0.0_DP
    rsolver%nminIterations = nsmoothingSteps
    rsolver%nmaxIterations = nsmoothingSteps
    rsolver%iresCheck  = 0
    
  end subroutine

  ! ***************************************************************************

  recursive subroutine stls_precondMultigrid (rsolver, rd)
  
  ! General preconditioning to a defect vector rd
  
  ! Solver structure
  type(t_spacetimelinsol), intent(inout) :: rsolver
  
  ! Defect vector to apply preconditioning to.
  type(t_spaceTimeVector), intent(inout) :: rd

    ! Set up the RHS by copying rd to p_Rvectors2.
    call sptivec_copyVector (rd,rsolver%p_Rvectors2(rsolver%ilevel))
    call spop_applyBC (rsolver%p_rmatrix%p_rboundaryCond, SPOP_DEFECT,&
        rsolver%p_Rvectors2(rsolver%ilevel))

    ! We forward the call to the actual multigrid solver which works
    ! recursively...
    call stls_precondMultigridInternal (rsolver, rsolver%ilevel, 0)

    ! Put the (weighted) solution into rd as returm value
    call sptivec_vectorLinearComb (rsolver%p_Rvectors1(rsolver%ilevel),rd,&
        rsolver%domega,0.0_DP)

  end subroutine

  ! ***************************************************************************

  recursive subroutine stls_precondMultigridInternal (rsolver, ilevel, niteFixed)
  
  ! General preconditioning to a defect vector in rsolver%p_Rvectors2.
  ! Result written into rsolver%p_Rvectors1.
  
  ! Solver structure
  type(t_spacetimelinsol), intent(inout) :: rsolver
  
  ! Current level
  integer, intent(in) :: ilevel
  
  ! Fixed number of iterations.
  ! =0: use stopping criteria from the solver structure.
  integer, intent(in) :: niteFixed
  
    ! local variables
    real(DP) :: dfactor
    integer :: ite,nite,ispacelevel,itimeLevel,ispacelevelcoarse
    integer :: nitecoarse
    real(DP) :: dresInit, dresCurrent, drho, drhoAsymp, dresMeanDev,dresMeanLast
    real(DP), dimension(3) :: DlastResiduals
    real(DP), dimension(3) :: DresLast
  
    ! Are we on the level of the coarse grid solver?
    if (ilevel .eq. rsolver%p_rcoarsegridsolver%ilevel) then
      if (rsolver%ioutputlevel .ge. 3) then
        call output_line("Space-Time Multigrid: Level "//trim(sys_siL(ilevel,10)))
      end if

      ! Call the coarse grid solver.
      call stls_precondDefect (rsolver%p_rcoarsegridsolver,rsolver%p_Rvectors2(ilevel))
      
      if (rsolver%p_rcoarsegridsolver%csolverStatus .eq. 2) then
        ! Coarse grid solver diverged. Stop here!
        rsolver%csolverStatus = 2
        return
      end if 
      
      ! Put the result to the solution.
      call sptivec_copyVector (rsolver%p_Rvectors2(ilevel),rsolver%p_Rvectors1(ilevel))
      
    else
      ! Do the iteration.
      
      ! p_Rvectors1 = solution
      ! p_Rvectors2 = rhs
      ! p_Rvectors3 = temp vector
      
      if (rsolver%ioutputlevel .ge. 3) then
        call output_line("Space-Time Multigrid: Level "//trim(sys_siL(ilevel,10)))
      end if
      
      ! Get the space and time level corresponding to the space-etime level.
      call sth_getLevel(rsolver%p_rspaceTimeHierarchy,ilevel,&
          ispaceLevel=ispaceLevel,itimeLevel=itimeLevel)
    
      ! and of the coarse grid level
      call sth_getLevel(rsolver%p_rspaceTimeHierarchy,ilevel-1,&
          ispaceLevel=ispaceLevelcoarse)
      
      ! Clear the solution vector p_Rvectors1.
      call sptivec_clearVector (rsolver%p_Rvectors1(ilevel))

      if (niteFixed .eq. 0) then
        ! Max. number of iterations.
        nite = rsolver%nmaxiterations+1
        
        ! If we are on a lower level, do as many iterations as prescribed by the cycle
        ! Otherwise, calculate the initial residuum for residuum checks.
        if (ilevel .eq. rsolver%ilevel) then
          ! Initial residuum.
          dresInit = sptivec_vectorNorm (rsolver%p_Rvectors2(ilevel),LINALG_NORML2)
          dresCurrent = dresInit
          DlastResiduals (:) = dresInit
        end if
      else
        ! Number of iterations fixed.
        nite = niteFixed
      end if
      
      DresLast(:) = 0.0_DP
      
      do ite = 1,nite
      
        if (niteFixed .eq. 0) then
          
          ! Stopping criterion check on the max. level.
        
          if (rsolver%ioutputlevel .ge. 2) then
            call output_line("Space-Time Multigrid: Step "//trim(sys_siL(ite-1,10))//", Residuum = "//&
                trim(sys_sdEP(dresCurrent,15,7)))
          end if

          if (ite .gt. rsolver%nmaxiterations) then
            ! Not completed.
            rsolver%csolverStatus = 1
            exit
          end if

          if (ite .ge. rsolver%nminiterations+1) then
            ! Check the stopping criterion.
            select case (rsolver%iresCheck)
            case (1)
              if ((dresCurrent .lt. rsolver%depsAbs) .or. &
                  (dresCurrent .lt. rsolver%depsRel*dresInit)) then
                exit
              end if
            case (2)
              if ((dresCurrent .lt. rsolver%depsAbs) .and. &
                  (dresCurrent .lt. rsolver%depsRel*dresInit)) then
                exit
              end if
            end select

            ! Calculate mean deviation in the last iterations.
            dresMeanDev = sum(abs(dresCurrent-DresLast(:))) / min(ite,size(DresLast))
            dresMeanLast = sum(DresLast(:)) / min(ite,size(DresLast))
            
            if ( (dresMeanDev .lt. rsolver%depsRelDiff*dresMeanLast) .or. &
                (dresMeanDev .lt. rsolver%depsAbsDiff) ) then
              exit
            end if
            
          end if
          
          if ((.not. (dresCurrent .lt. rsolver%ddivAbs)) .or. &
              (.not. (dresCurrent .lt. rsolver%ddivRel*dresInit))) then
            if (rsolver%ioutputlevel .ge. 1) then
              call output_lbrk()
              call output_line("Space-Time Multigrid: Divergence detected! Iteration stopped.");
            end if
            rsolver%csolverStatus = 2
            exit
          end if
          
        end if

!        ! RHS      
!        call stpp_postproc (rsolver%p_rmatrices(ilevel)%p_rmatrix%p_rphysics,rsolver%p_Rvectors2(ilevel))
!        print *,"RHS"
!        read (*,*)

!        ! Temp. defect
!        call sptivec_copyVector (rsolver%p_Rvectors2(ilevel),rsolver%p_Rvectors3(ilevel))
!        call stmv_matvec (rsolver%p_Rmatrices(ilevel)%p_rmatrix, &
!            rsolver%p_Rvectors1(ilevel),rsolver%p_Rvectors3(ilevel), -1.0_DP, 1.0_DP)
!        call spop_applyBC (rsolver%p_rmatrices(ilevel)%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, rsolver%p_Rvectors3(ilevel))
!        call stpp_postproc (rsolver%p_rmatrices(ilevel)%p_rmatrix%p_rphysics,rsolver%p_Rvectors3(ilevel))
!        call stpp_printDefectSubnormsDirect(rsolver%p_Rvectors3(ilevel))
!        print *,"current defect"
!        read (*,*)

        ! Presmoothing
        if (associated(rsolver%p_Rpresmoothers)) then
          if (rsolver%p_Rpresmoothers(ilevel)%csolverType .ne. STLS_TYPE_NONE) then
            call stls_solveAdaptively (rsolver%p_rpresmoothers(ilevel), &
                rsolver%p_Rvectors1(ilevel),rsolver%p_Rvectors2(ilevel),rsolver%p_Rvectors3(ilevel))

            if (rsolver%p_rpresmoothers(ilevel)%csolverStatus .eq. 2) then
              ! Smoother diverged. Stop here!
              rsolver%csolverStatus = 2
              exit
            end if
          end if
        end if

!        ! Temp. defect
!        call sptivec_copyVector (rsolver%p_Rvectors2(ilevel),rsolver%p_Rvectors3(ilevel))
!        call stmv_matvec (rsolver%p_Rmatrices(ilevel)%p_rmatrix, &
!            rsolver%p_Rvectors1(ilevel),rsolver%p_Rvectors3(ilevel), -1.0_DP, 1.0_DP)
!        call spop_applyBC (rsolver%p_rmatrices(ilevel)%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, rsolver%p_Rvectors3(ilevel))
!        call stpp_postproc (rsolver%p_rmatrices(ilevel)%p_rmatrix%p_rphysics,rsolver%p_Rvectors3(ilevel))
!        call stpp_printDefectSubnormsDirect(rsolver%p_Rvectors3(ilevel))
!        print *,"current defect after sm."
!        read (*,*)
            
!        ! DEBUG: Real solution.
!        call sptivec_copyVector (rsolver%p_Rvectors2(ilevel),rsolver%p_Rvectors3(ilevel))
!        call stmv_matvec (rsolver%p_Rmatrices(ilevel)%p_rmatrix, &
!            rsolver%p_Rvectors1(ilevel),rsolver%p_Rvectors3(ilevel), -1.0_DP, 1.0_DP)
!        call spop_applyBC (rsolver%p_rmatrices(ilevel)%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, rsolver%p_Rvectors3(ilevel))
!        call stls_precondDefect (rsolver%p_Rpreconditioners(ilevel),rsolver%p_Rvectors3(ilevel))
!        print *,"real sol."
!        call stpp_postproc (rsolver%p_rmatrices(ilevel)%p_rmatrix%p_rphysics,rsolver%p_Rvectors3(ilevel))
!        read (*,*)

        ! Build the defect into the temp vector.
        call sptivec_copyVector (rsolver%p_Rvectors2(ilevel),rsolver%p_Rvectors3(ilevel))
        call stmv_matvec (rsolver%p_Rmatrices(ilevel)%p_rmatrix, &
            rsolver%p_Rvectors1(ilevel),rsolver%p_Rvectors3(ilevel), -1.0_DP, 1.0_DP)
        call spop_applyBC (rsolver%p_rmatrices(ilevel)%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, rsolver%p_Rvectors3(ilevel))
        
        !call stpp_postproc (rsolver%p_rmatrices(ilevel)%p_rmatrix%p_rphysics,rsolver%p_Rvectors3(ilevel))

        ! Restriction.
        call sptipr_performRestriction (rsolver%p_rspaceTimeProjection,ilevel,&
            rsolver%p_Rvectors2(ilevel-1),rsolver%p_Rvectors3(ilevel),&
            rsolver%p_RspaceVectors(ispaceLevelcoarse),rsolver%p_RspaceVectors(ispaceLevel),&
            rsolver%p_Rvectors4(ilevel-1),rsolver%p_Rvectors4(ilevel))
        
        ! Boundary conditions.
        call spop_applyBC (rsolver%p_rmatrices(ilevel)%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, rsolver%p_Rvectors2(ilevel))

        ! Solve on the lower level.
        ! The number of coarse grid iterations depends on the cycle...
        select case (rsolver%icycle)
          case (0)
            ! F-cycle. The number of coarse grid iterations has to be calculated
            ! from the current iteration. 
            if (niteFixed .eq. 0) then
              niteCoarse = 2
            else
              ! 1st iteration: W-cycle. 2nc iteration: V-cycle.
              niteCoarse = niteFixed - ite + 1
            end if
            
          case (1,2)
            ! V/W cycle.
            ! The icycle identifier specifies the number of coarse 
            ! grid iterations.
            niteCoarse = rsolver%icycle
            
        end select
        call stls_precondMultigridInternal (rsolver, ilevel-1, niteCoarse)
        
        if (rsolver%ioutputlevel .ge. 3) then
          call output_line("Space-Time Multigrid: Level "//trim(sys_siL(ilevel,10)))
        end if
        
        if (rsolver%csolverStatus .eq. 2) then
          ! Coarse grid solver diverged. Stop here!
          exit
        end if
        
        ! Prolongation
        call sptipr_performProlongation (rsolver%p_rspaceTimeProjection,ilevel,&
            rsolver%p_Rvectors1(ilevel-1),rsolver%p_Rvectors3(ilevel),&
            rsolver%p_RspaceVectors(ispaceLevelcoarse),rsolver%p_RspaceVectors(ispaceLevel))
        
        ! Boundary conditions.
        call spop_applyBC (rsolver%p_rmatrices(ilevel)%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, rsolver%p_Rvectors3(ilevel))
        
!        print *,"cgr sol."
!        call stpp_postproc (rsolver%p_rmatrices(ilevel)%p_rmatrix%p_rphysics,rsolver%p_Rvectors3(ilevel))
!        read (*,*) 
        dfactor = 1.0_DP
!        if (sstring .eq. "y") then
!          print *,"ok"
!          dfactor = 2.0_DP
!        end if
        select case (rsolver%iadcgcorr)
        case (1)
          ! Get the factor by energy minimisatino.
          call stls_mgCalcCgCorrFactorEnergy(rsolver%p_Rmatrices(ilevel)%p_rmatrix,&
              rsolver%p_Rvectors1(ilevel),rsolver%p_Rvectors3(ilevel),&
              rsolver%p_Rvectors2(ilevel),rsolver%p_Rvectors4(ilevel),dfactor)

          if (rsolver%ioutputLevel .ge. 2) then
            call output_line ("Space-Time Multigrid: Coarse grid correction factor = "//&
                trim(sys_sdEL(dfactor,10)))
          end if
          
          dfactor = max(min(dfactor,rsolver%dadcgcorrMax),rsolver%dadcgcorrMin)

        case (2)
          ! Get the factor by energy minimisatino.
          call stls_mgCalcCgCorrFactorDef(rsolver%p_Rmatrices(ilevel)%p_rmatrix,&
              rsolver%p_Rvectors1(ilevel),rsolver%p_Rvectors3(ilevel),&
              rsolver%p_Rvectors2(ilevel),rsolver%p_Rvectors4(ilevel),&
              rsolver%p_Rvectors5(ilevel),dfactor)

          if (rsolver%ioutputLevel .ge. 2) then
            call output_line ("Space-Time Multigrid: Coarse grid correction factor = "//&
                trim(sys_sdEL(dfactor,10)))
          end if
          
          dfactor = max(min(dfactor,rsolver%dadcgcorrMax),rsolver%dadcgcorrMin)

        end select
        
        ! Coarse grid correction
        call sptivec_vectorLinearComb (rsolver%p_Rvectors3(ilevel),&
            rsolver%p_Rvectors1(ilevel),dfactor,1.0_DP)

        !call stpp_postproc (rsolver%p_rmatrices(ilevel)%p_rmatrix%p_rphysics,rsolver%p_Rvectors1(ilevel))

!          call sptivec_copyVector (rsolver%p_Rvectors2(ilevel),rsolver%p_Rvectors3(ilevel))
!          call stmv_matvec (rsolver%p_Rmatrices(ilevel)%p_rmatrix, &
!              rsolver%p_Rvectors1(ilevel),rsolver%p_Rvectors3(ilevel), -1.0_DP, 1.0_DP)
!          call spop_applyBC (rsolver%p_rmatrices(ilevel)%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, rsolver%p_Rvectors3(ilevel))
!
!        call stpp_postproc (rsolver%p_rmatrices(ilevel)%p_rmatrix%p_rphysics,rsolver%p_Rvectors3(ilevel))
        
!        if (ilevel .eq. rsolver%ilevel) then
!          call output_line ("Defect before smoothing:")
!          call sptivec_copyVector (rsolver%p_Rvectors2(ilevel),rsolver%p_Rvectors3(ilevel))
!          call stmv_matvec (rsolver%p_Rmatrices(ilevel)%p_rmatrix, &
!              rsolver%p_Rvectors1(ilevel),rsolver%p_Rvectors3(ilevel), -1.0_DP, 1.0_DP)
!          call spop_applyBC (rsolver%p_rmatrices(ilevel)%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, rsolver%p_Rvectors3(ilevel))
!          
!          call stpp_printDefectSubnormsDirect (rsolver%p_Rvectors3(ilevel))
!          
!          dresCurrent = sptivec_vectorNorm (rsolver%p_Rvectors3(ilevel),LINALG_NORML2)
!        end if
        
        ! Postsmoothing
        if (associated(rsolver%p_Rpostsmoothers)) then
          if (rsolver%p_Rpostsmoothers(ilevel)%csolverType .ne. STLS_TYPE_NONE) then
            call stls_solveAdaptively (rsolver%p_rpostsmoothers(ilevel), &
                rsolver%p_Rvectors1(ilevel),rsolver%p_Rvectors2(ilevel),rsolver%p_Rvectors3(ilevel))

            if (rsolver%p_rpostsmoothers(ilevel)%csolverStatus .eq. 2) then
              ! Smoother diverged. Stop here!
              rsolver%csolverStatus = 2
              exit
            end if
          end if
        end if
        
        !call stpp_postproc (rsolver%p_rmatrices(ilevel)%p_rmatrix%p_rphysics,rsolver%p_Rvectors1(ilevel))

        ! New defect on the max. level for residuum checks.
        if (ilevel .eq. rsolver%ilevel) then
          call sptivec_copyVector (rsolver%p_Rvectors2(ilevel),rsolver%p_Rvectors3(ilevel))
          call stmv_matvec (rsolver%p_Rmatrices(ilevel)%p_rmatrix, &
              rsolver%p_Rvectors1(ilevel),rsolver%p_Rvectors3(ilevel), -1.0_DP, 1.0_DP)
          call spop_applyBC (rsolver%p_rmatrices(ilevel)%p_rmatrix%p_rboundaryCond, &
              SPOP_DEFECT, rsolver%p_Rvectors3(ilevel))
          
!          call output_line ("Defect after smoothing:")
!          call stpp_printDefectSubnormsDirect (rsolver%p_Rvectors3(ilevel))
          
          DresLast(:) = eoshift(DresLast(:),shift=-1,boundary=dresCurrent)
          dresCurrent = sptivec_vectorNorm (rsolver%p_Rvectors3(ilevel),LINALG_NORML2)
          
          ! Shift the residuals
          DlastResiduals = eoshift(DlastResiduals,shift=1,boundary=dresCurrent)
        end if
        
      end do
      
      if (niteFixed .eq. 0) then
        ! Convergence rate
        if (ite .gt. 0) then
          if (dresInit .eq. 0.0_DP) then
            drho = dresCurrent ** (1.0_DP/real(ite-1,DP))
          else
            ! Convergence rate
            drho = (dresCurrent / dresInit) ** (1.0_DP/real(ite-1,DP))
            
            ! Asymptotic convergence rate
            drhoAsymp = (DlastResiduals(size(DlastResiduals)) / DlastResiduals(1)) ** &
                (1.0_DP/real(min(ite-1,size(DlastResiduals)-1),DP))
          end if
        else
          drho = 0.0_DP
        end if
        
        if (rsolver%ioutputlevel .ge. 2) then
          call output_lbrk()
          call output_line("Initial Residual        : "//trim(sys_sdEP(dresInit,15,7)))
          call output_line("Final Residual          : "//trim(sys_sdEP(dresCurrent,15,7)))
          call output_line("Iterations              : "//trim(sys_siL(ite-1,10)))
          call output_line("Convergence rate        : "//trim(sys_sdEP(drho,15,7)))
          call output_line("asympt. convergence rate: "//trim(sys_sdEP(drhoAsymp,15,7)))
        else if (rsolver%ioutputlevel .ge. 1) then
          call output_line("Space-Time Multigrid: Iterations/Conv. rate: "//&
              trim(sys_siL(ite-1,10))//"/"//trim(sys_sdEP(drho,15,7)))
        end if
      end if
      
    end if

  end subroutine
  
  ! ***************************************************************************

  subroutine stls_mgCalcCgCorrFactorEnergy(rmatrix,rsolution,rcorrection,rrhs,rtemp,dalpha)
  
  ! Calculates a factor for the adaptive coarse grid correction by 
  ! energy minimisation.
  
  ! Space-time matrix
  type(t_spaceTimeMatrix), intent(in) :: rmatrix
  
  ! Solution that should be corrected.
  type(t_spaceTimeVector), intent(in) :: rsolution

  ! prolongated correction vector
  type(t_spaceTimeVector), intent(in) :: rcorrection

  ! Temp vector
  type(t_spaceTimeVector), intent(inout) :: rtemp
  
  ! RHS vector
  type(t_spaceTimeVector), intent(in) :: rrhs
  
  ! OUT: Correction factor.
  real(DP), intent(out) :: dalpha
  
    ! local variables
    real(DP) :: d1,d2
    
    dalpha = 1.0_DP
    
    ! Calculate the nominator.
    call sptivec_copyVector (rrhs,rtemp)
    call stmv_matvec (rmatrix, rsolution, rtemp, -1.0_DP, 1.0_DP)
    call spop_applyBC (rmatrix%p_rboundaryCond, SPOP_DEFECT, rtemp)
    d1 = sptivec_scalarProduct(rtemp,rcorrection)
    if (d1 .eq. 0.0_DP) then
      ! Cancel. No correction or solution reached.
      return
    end if
    
    ! Calculate the denominator.
    call stmv_matvec (rmatrix, rcorrection, rtemp, 1.0_DP, 0.0_DP)
    call spop_applyBC (rmatrix%p_rboundaryCond, SPOP_DEFECT, rtemp)
    d2 = sptivec_scalarProduct(rtemp,rcorrection)
    if (d2 .eq. 0.0_DP) then
      ! Cancel. No correction.
      return
    end if
    
    ! Get the correction.
    dalpha = d1 / d2
    
  end subroutine

  ! ***************************************************************************

  subroutine stls_mgCalcCgCorrFactorDef(rmatrix,rsolution,rcorrection,rrhs,&
      rtemp1,rtemp2,dalpha)
  
  ! Calculates a factor for the adaptive coarse grid correction by 
  ! defect minimisation.
  
  ! Space-time matrix
  type(t_spaceTimeMatrix), intent(in) :: rmatrix
  
  ! Solution that should be corrected.
  type(t_spaceTimeVector), intent(in) :: rsolution

  ! prolongated correction vector
  type(t_spaceTimeVector), intent(in) :: rcorrection

  ! Temp vectors
  type(t_spaceTimeVector), intent(inout) :: rtemp1,rtemp2
  
  ! RHS vector
  type(t_spaceTimeVector), intent(in) :: rrhs
  
  ! OUT: Correction factor.
  real(DP), intent(out) :: dalpha
  
    ! local variables
    real(DP) :: d1,d2
    
    dalpha = 1.0_DP
    
    ! Calculate the nominator.
    call sptivec_copyVector (rrhs,rtemp1)
    call stmv_matvec (rmatrix, rsolution, rtemp1, -1.0_DP, 1.0_DP)
    call spop_applyBC (rmatrix%p_rboundaryCond, SPOP_DEFECT, rtemp1)
    
    call stmv_matvec (rmatrix, rcorrection, rtemp2, 1.0_DP, 0.0_DP)
    call spop_applyBC (rmatrix%p_rboundaryCond, SPOP_DEFECT, rtemp2)
    
    d1 = sptivec_scalarProduct(rtemp1,rtemp2)
    if (d1 .eq. 0.0_DP) then
      ! Cancel. No correction or solution reached.
      return
    end if
    
    ! Calculate the denominator.
    d2 = sptivec_scalarProduct(rtemp2,rtemp2)
    if (d2 .eq. 0.0_DP) then
      ! Cancel. No correction.
      return
    end if
    
    ! Get the correction.
    dalpha = d1 / d2
    
  end subroutine

end module
