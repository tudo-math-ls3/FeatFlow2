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
  
  implicit none

  ! Linear space-time solver tyes

  ! None
  integer, parameter :: STLS_TYPE_NONE       = 0
  
  ! Jacobi
  integer, parameter :: STLS_TYPE_JACOBI     = 1
  
  ! FB-GS
  integer, parameter :: STLS_TYPE_FBGS      = 2

  ! FB-GS2
  integer, parameter :: STLS_TYPE_FBGS2     = 3
  
  ! Defect correction
  integer, parameter :: STLS_TYPE_DEFCORR    = 4
  
  ! Two-grid
  integer, parameter :: STLS_TYPE_MULTIGRID  = 5

  ! Linear solver structure.
  type t_spacetimelinsol
    
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
    
    ! Whether or not to check the residual for convergence.
    ! =0: do not check.
    ! =1: check: relative or absolute stopping criterion must be fulfilled.
    ! =2: check: relative and absolute stopping criterion must be fulfilled.
    integer :: iresCheck = 1
    
    ! Number of levels
    integer :: nlevels = 1
    
    ! Minimum/Maximum number of steps
    integer :: nminIterations = 1
    integer :: nmaxIterations = 1000
    
    ! Output level
    ! =0: no output
    ! =1: basic output
    ! =2: standard output
    ! =3: extended output
    integer :: ioutputLevel = 2
    
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
    type(t_spaceTimeMatrix), dimension(:), pointer :: p_rmatrices => null()
    
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
    
    ! type of the linear solver in space if used.
    ! =LINSOL_ALG_UNDEFINED: not used.
    ! =LINSOL_ALG_UMFPACK4 : UMFPACK
    integer :: cspaceSolverType = LINSOL_ALG_UNDEFINED
    
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

  subroutine stls_init (rsolver,csolverType,rspaceTimeHierarchy,ilevel,&
      ballocTempMatrices,ballocTempVectors,rpreconditioner,cspaceSolverType,RmatVecTempl)
  
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
  
  ! OPTIONAL: Identifier for a linear solver in space that might be used
  ! for preconditioning.
  ! =LINSOL_ALG_UMFPACK4: UMFPACK-4
  integer, intent(in), optional :: cspaceSolverType
  
  ! OPTINAL: Spatial matrix templates. Necessary if a matrix-based
  ! preconditioner in space is used.
  type(t_matvecTemplates), dimension(:), target, optional :: RmatVecTempl
  
    ! local variables
    integer :: ispacelevel,ilev
  
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
    
    if (present(cspaceSolverType) .or. ballocTempMatrices) then
    
      ! Allocate matrices for the space solver
      allocate(rsolver%p_RspaceMatrices(ilevel))

      call sth_getLevel (rspaceTimeHierarchy,ilevel,ispaceLevel=ispaceLevel)
      
      ! Initialise the solver in space.
      rsolver%cspaceSolverType = cspaceSolverType
      select case (cspaceSolverType)
      case (LINSOL_ALG_UMFPACK4) 
        call linsol_initUMFPACK4 (rsolver%p_rspaceSolver)
      end select
        
      ! Prepare boundary conditions
      allocate(rsolver%p_RdiscreteBC(ilevel))
      do ilev = 1,ilevel
        call bcasm_initDiscreteBC(rsolver%p_RdiscreteBC(ilev))
      end do
      
    end if
  
    if (ballocTempVectors) then
    
      ! Allocate matrices for the space solver
      allocate(rsolver%p_RspaceVectors(ilevel))

    end if
  
  end subroutine

  ! ***************************************************************************

  subroutine stls_done (rsolver)
  
  ! Releases a solver node
  
  ! Solver structure to be released
  type(t_spacetimelinsol), intent(inout) :: rsolver
  
    integer :: ilev
  
    rsolver%cspaceSolverType = -1
    if (associated(rsolver%p_RdiscreteBC)) then
      ! Release boundary conditions
      do ilev = 1,rsolver%ilevel
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
    integer :: ispaceLevel,ilev
  
    if (rsolver%csolverType .eq. STLS_TYPE_NONE) then
      call sys_halt()
    end if
  
    call sth_getLevel (rsolver%p_rspaceTimeHierarchy,&
        rsolver%ilevel,p_rfeSpaceLevel=p_rfeSpaceLevel,p_rtimeDiscr=p_rtimeDiscr)
    p_rspaceDiscr => p_rfeSpaceLevel%p_rdiscretisation
  
    ! Some initialisation, based on the solver type.
    select case (rsolver%csolverType)
      
    case (STLS_TYPE_JACOBI)
      
      ! Allocate temp vectors.
      call lsysbl_createVectorBlock (p_rspaceDiscr,rsolver%rspaceTemp1)
    
    case (STLS_TYPE_FBGS)
      
      ! Allocate temp vectors.
      call lsysbl_createVectorBlock (p_rspaceDiscr,rsolver%rspaceTemp1)
      call lsysbl_createVectorBlock (p_rspaceDiscr,rsolver%rspaceTemp2)
      call lsysbl_createVectorBlock (p_rspaceDiscr,rsolver%rspaceTemp3)
      call sptivec_initVector (rsolver%rspaceTimeTemp1,p_rtimeDiscr,p_rspaceDiscr)

    case (STLS_TYPE_FBGS2)
      
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
    
    case (STLS_TYPE_MULTIGRID)
    
      ! Allocate temp vectors on all levels.
      allocate(rsolver%p_Rvectors1(rsolver%ilevel))
      allocate(rsolver%p_Rvectors2(rsolver%ilevel))
      allocate(rsolver%p_Rvectors3(rsolver%ilevel))
      do ilev = 1,rsolver%ilevel

        call sth_getLevel (rsolver%p_rspaceTimeHierarchy,&
            ilev,p_rfeSpaceLevel=p_rfeSpaceLevel,p_rtimeDiscr=p_rtimeDiscr)
        p_rspaceDiscr => p_rfeSpaceLevel%p_rdiscretisation

        call sptivec_initVector (rsolver%p_Rvectors1(ilev),p_rtimeDiscr,p_rspaceDiscr)
        call sptivec_initVector (rsolver%p_Rvectors2(ilev),p_rtimeDiscr,p_rspaceDiscr)
        call sptivec_initVector (rsolver%p_Rvectors3(ilev),p_rtimeDiscr,p_rspaceDiscr)
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
      do ilev = 1,rsolver%ilevel
        
        call sth_getLevel (rsolver%p_rspaceTimeHierarchy,ilev,&
            p_rfeSpaceLevel,p_rtimeDiscr,ispaceLevel)
        p_rspaceDiscr => p_rfeSpaceLevel%p_rdiscretisation
        
        call stmv_allocSubmatrix (rsolver%p_rmatrix%cmatrixType,&
            rsolver%p_rmatrix%p_rphysics,rsolver%p_RmatVecTempl(ispaceLevel),&
            rsolver%p_RspaceMatrices(ilev),rsolver%p_RdiscreteBC(ilev))
            
      end do
      
    end if
    
    ! Temp vectors
    if (associated(rsolver%p_RspaceVectors)) then
      do ilev = 1,rsolver%ilevel
        
        call sth_getLevel (rsolver%p_rspaceTimeHierarchy,ilev,&
            p_rfeSpaceLevel,p_rtimeDiscr,ispaceLevel)
        p_rspaceDiscr => p_rfeSpaceLevel%p_rdiscretisation
        
        call lsysbl_createVecBlockByDiscr (p_rspaceDiscr,rsolver%p_RspaceVectors(ilev),.false.)
        
      end do
      
    end if
  
    ! Spatial solver
    if (associated(rsolver%p_rspaceSolver)) then
      
      select case (rsolver%p_rspaceSolver%calgorithm)
      case (LINSOL_ALG_UMFPACK4)
        ! Attach the temp matrices to the linear solver.
        call linsol_setMatrices (rsolver%p_rspaceSolver,rsolver%p_RspaceMatrices)
      end select
      
    end if

  end subroutine

  ! ***************************************************************************

  recursive subroutine stls_doneData (rsolver)
  
  ! Some cleanup after solving
  
  ! Solver structure 
  type(t_spacetimelinsol), intent(inout) :: rsolver
  
    ! local variables
    integer :: ilev
  
    ! Solver dependent cleanup
    select case (rsolver%csolverType)
      
    case (STLS_TYPE_JACOBI)
      
      ! Deallocate temp vectors.
      call lsysbl_releaseVector (rsolver%rspaceTemp1)
    
    case (STLS_TYPE_FBGS)
      
      ! Deallocate temp vectors.
      call sptivec_releaseVector (rsolver%rspaceTimeTemp1)
      call lsysbl_releaseVector (rsolver%rspaceTemp3)
      call lsysbl_releaseVector (rsolver%rspaceTemp2)
      call lsysbl_releaseVector (rsolver%rspaceTemp1)

    case (STLS_TYPE_FBGS2)
      
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
        call sptivec_releaseVector (rsolver%p_Rvectors3(ilev))
        call sptivec_releaseVector (rsolver%p_Rvectors2(ilev))
        call sptivec_releaseVector (rsolver%p_Rvectors1(ilev))
      end do
      deallocate(rsolver%p_Rvectors3)
      deallocate(rsolver%p_Rvectors2)
      deallocate(rsolver%p_Rvectors1)
    
    end select
    
    if (associated(rsolver%p_RspaceMatrices)) then
      ! Release matrices
      do ilev = 1,rsolver%ilevel
        call lsysbl_releaseMatrix(rsolver%p_RspaceMatrices(ilev))
      end do
    end if

    if (associated(rsolver%p_RspaceVectors)) then
      ! Release matrices
      do ilev = 1,rsolver%ilevel
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
  
    select case (rsolver%csolverType)
    case (STLS_TYPE_DEFCORR)
      call stls_precondDefCorr (rsolver, rd)
    case (STLS_TYPE_JACOBI)
      call stls_precondBlockJacobi (rsolver, rd)
    case (STLS_TYPE_FBGS)
      call stls_precondBlockFBGS (rsolver, rd)
    case (STLS_TYPE_FBGS2)
      call stls_precondBlockFBGS2 (rsolver, rd)
    case (STLS_TYPE_MULTIGRID)
      call stls_precondMultigrid (rsolver, rd)
    end select
  
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
    call sptivec_vectorLinearComb (rd,rx,1.0_DP,1.0_DP)

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
    real(DP) :: dresInit, dresCurrent, drho
    
    ! One temp vector is used for the residuum, one for the solution.
    call sptivec_clearVector (rsolver%rspaceTimeTemp1)
    
    ! Initial residuum.
    call sptivec_copyVector (rd,rsolver%rspaceTimeTemp2)
    call spop_applyBC (rsolver%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, rsolver%rspaceTimeTemp2)
    dresInit = sptivec_vectorNorm (rsolver%rspaceTimeTemp2,LINALG_NORML2)
    dresCurrent = dresInit
    
    do ite = 1,rsolver%nmaxiterations+1
      if (rsolver%ioutputlevel .ge. 2) then
        call output_line("Space-Time DefCorr: Step "//trim(sys_siL(ite-1,10))//", Residuum = "//&
            trim(sys_sdEP(dresCurrent,15,7)))
      end if
      
      if (ite .gt. rsolver%nmaxiterations) &
        exit

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
      end if
      
      ! Defect preconditioning
      if (associated (rsolver%p_rpreconditioner)) then
        call stls_precondDefect(rsolver%p_rpreconditioner,rsolver%rspaceTimeTemp2)
      end if
      
      ! Add to the new solution.
      call sptivec_vectorLinearComb (rsolver%rspaceTimeTemp2,rsolver%rspaceTimeTemp1,1.0_DP,1.0_DP)
      
      ! New defect
      call sptivec_copyVector (rd,rsolver%rspaceTimeTemp2)
      call stmv_matvec (rsolver%p_rmatrix, &
          rsolver%rspaceTimeTemp1, rsolver%rspaceTimeTemp2, -1.0_DP, 1.0_DP)
      call spop_applyBC (rsolver%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, rsolver%rspaceTimeTemp2)
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
      call output_lbrk()
      call output_line("Iterations/Conv. rate: "//&
          trim(sys_siL(ite-1,10))//"/"//trim(sys_sdEP(drho,15,7)))
    end if
    
    ! Put the (weighted) solution into rd as returm value
    call sptivec_vectorLinearComb (rsolver%rspaceTimeTemp1,rd,rsolver%domega,0.0_DP)
  
  end subroutine

  ! ***************************************************************************

  subroutine stls_initBlockJacobi (rsolver,rspaceTimeHierarchy,ilevel,domega,&
      cspaceSolverType,RmatVecTempl)
  
  ! Initialise a block Jacobi correction solver.
  
  ! Solver structure to be initialised
  type(t_spacetimelinsol), intent(out) :: rsolver
  
  ! Underlying space-time hierarchy
  type(t_spacetimeHierarchy), intent(in), target :: rspaceTimeHierarchy
  
  ! Level of the solver  in the space-time hierarchy
  integer, intent(in) :: ilevel

  ! Damping parameter
  real(DP), intent(in) :: domega

  ! OPTIONAL: Identifier for a linear solver in space that might be used
  ! for preconditioning.
  ! =0: UMFPACK
  integer, intent(in) :: cspaceSolverType

  ! OPTINAL: Spatial matrix templates. Necessary if a matrix-based
  ! preconditioner in space is used.
  type(t_matvecTemplates), dimension(:), intent(in), target :: RmatVecTempl

    ! Basic initialisation
    call stls_init(rsolver,STLS_TYPE_JACOBI,rspaceTimeHierarchy,ilevel,&
        .false.,.false.,cspaceSolverType=cspaceSolverType,RmatVecTempl=RmatVecTempl)
        
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
    integer :: istep, ierror
    real(DP), dimension(:), pointer :: p_Ddata
    
    ! Initialise the linear solver.
    call linsol_initStructure(rsolver%p_rspaceSolver,ierror)
    if (ierror .ne. 0) then
      call output_line ("Spatial system cannot be symbolically factorised in timestep "&
          //trim(sys_siL(istep,10)))
      call sys_halt()
    end if

    call lsysbl_getbase_double (rsolver%rspaceTemp1,p_Ddata)

    ! We loop through the timesteps and apply the spatial solver in each timestep.
    do istep = 1,rd%NEQtime
      
      ! Load the timestep.
      call sptivec_getTimestepData(rd,istep,rsolver%rspaceTemp1)
      
      ! Assemble diagonal submatrix of that timestep.
      call stmv_getSubmatrix (rsolver%p_rmatrix, istep, istep, &
          rsolver%p_RspaceMatrices(rsolver%ilevel))
          
      ! Apply the boundary conditions to the matrix
      call stmv_implementDefBCSubmatrix (rsolver%p_rmatrix%p_rboundaryCond, &
          rsolver%p_RspaceMatrices(rsolver%ilevel), istep, istep, &
          rsolver%p_RdiscreteBC(rsolver%ilevel))
      
      !call matio_writeBlockMatrixHR (rsolver%p_RspaceMatrices(rsolver%ilevel), "matrix",&
      !    .true., 0, "matrix.txt", "(E10.3)")

      ! Apply the space solver
      call linsol_initData(rsolver%p_rspaceSolver,ierror)
      if (ierror .ne. 0) then
        call output_line ("Spatial system cannot be factorised in timestep "//&
            trim(sys_siL(istep,10)))
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

  subroutine stls_initBlockFBGS (rsolver,rspaceTimeHierarchy,ilevel,drelax,&
      cspaceSolverType,RmatVecTempl)
  
  ! Initialise a block GS correction solver, working on the decoupled solution.
  
  ! Solver structure to be initialised
  type(t_spacetimelinsol), intent(out) :: rsolver
  
  ! Underlying space-time hierarchy
  type(t_spacetimeHierarchy), intent(in), target :: rspaceTimeHierarchy
  
  ! Level of the solver  in the space-time hierarchy
  integer, intent(in) :: ilevel

  ! Relaxation parameter
  real(DP), intent(in) :: drelax

  ! OPTIONAL: Identifier for a linear solver in space that might be used
  ! for preconditioning.
  ! =0: UMFPACK
  integer, intent(in) :: cspaceSolverType

  ! OPTINAL: Spatial matrix templates. Necessary if a matrix-based
  ! preconditioner in space is used.
  type(t_matvecTemplates), dimension(:), intent(in), target :: RmatVecTempl

    ! Basic initialisation
    call stls_init(rsolver,STLS_TYPE_FBGS,rspaceTimeHierarchy,ilevel,&
        .false.,.false.,cspaceSolverType=cspaceSolverType,RmatVecTempl=RmatVecTempl)
        
    rsolver%drelax = drelax
  
  end subroutine

  ! ***************************************************************************

  recursive subroutine stls_precondBlockFBGS (rsolver, rd)
  
  ! General preconditioning to a defect vector rd
  
  ! Solver structure
  type(t_spacetimelinsol), intent(inout) :: rsolver
  
  ! Defect vector to apply preconditioning to.
  type(t_spaceTimeVector), intent(inout) :: rd

    ! local variables
    integer :: istep, ierror,jstep,i
    real(DP), dimension(:), pointer :: p_Dx,p_Dd,p_Ddata
    
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
    
    do jstep = 1,2
    
      ! Forward sweep.
        
      ! We loop through the timesteps and apply the spatial solver in each timestep.
      do istep = 1,rd%NEQtime
        
        ! Load the timestep.
        call sptivec_getTimestepData(rd,istep,rsolver%rspaceTemp1)
        
        if (istep .gt. 1) then
          ! Load the previous timestep.
          call sptivec_getTimestepData(rsolver%rspaceTimeTemp1,istep-1,rsolver%rspaceTemp2)
          
          ! Subtract the primal.
          call stmv_getSubmatrix (rsolver%p_rmatrix, istep, istep-1, &
              rsolver%p_RspaceMatrices(rsolver%ilevel))

          call stmv_implementDefBCSubmatrix (rsolver%p_rmatrix%p_rboundaryCond, &
              rsolver%p_RspaceMatrices(rsolver%ilevel), istep, istep-1, &
              rsolver%p_RdiscreteBC(rsolver%ilevel))

          call lsysbl_blockMatVec (rsolver%p_RspaceMatrices(rsolver%ilevel), &
              rsolver%rspaceTemp2, rsolver%rspaceTemp1, -1.0_DP, 1.0_DP)
        end if
        
        ! Subtract the diagonal, set up the preconditioner matrix
        
        ! Assemble diagonal submatrix of that timestep.
        call stmv_getSubmatrix (rsolver%p_rmatrix, istep, istep, &
            rsolver%p_RspaceMatrices(rsolver%ilevel))
            
        ! Apply the boundary conditions to the matrix
        call stmv_implementDefBCSubmatrix (rsolver%p_rmatrix%p_rboundaryCond, &
            rsolver%p_RspaceMatrices(rsolver%ilevel), istep, istep, &
            rsolver%p_RdiscreteBC(rsolver%ilevel))

        ! Subtract the diagonal.
        call sptivec_getTimestepData(rsolver%rspaceTimeTemp1,istep,rsolver%rspaceTemp2)
        
        call lsysbl_blockMatVec (rsolver%p_RspaceMatrices(rsolver%ilevel), &
            rsolver%rspaceTemp2, rsolver%rspaceTemp1, -1.0_DP, 1.0_DP)
        
        ! Remove the primal part from the matrix
        !call stmv_reduceDiagToDual (rsolver%p_RspaceMatrices(rsolver%ilevel))
        
        ! Apply the space solver
        call linsol_initData(rsolver%p_rspaceSolver,ierror)
        if (ierror .ne. 0) then
          call output_line ("Spatial system cannot be factorised in timestep "//&
              trim(sys_siL(istep,10)))
              
          call matio_writeBlockMatrixHR (rsolver%p_RspaceMatrices(rsolver%ilevel), &
              "matrix", .true., 0, "matrix.txt", "(E10.3)") 
              
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
      
      ! Backward sweep.
        
      ! We loop through the timesteps and apply the spatial solver in each timestep.
      do istep = rd%NEQtime,1,-1
        
        ! Load the timestep.
        call sptivec_getTimestepData(rd,istep,rsolver%rspaceTemp1)
        
        if (istep .lt. rd%NEQtime) then
          ! Load the previous timestep.
          call sptivec_getTimestepData(rsolver%rspaceTimeTemp1,istep+1,rsolver%rspaceTemp2)
          
          ! Subtract the primal.
          call stmv_getSubmatrix (rsolver%p_rmatrix, istep, istep+1, &
              rsolver%p_RspaceMatrices(rsolver%ilevel))

          call stmv_implementDefBCSubmatrix (rsolver%p_rmatrix%p_rboundaryCond, &
              rsolver%p_RspaceMatrices(rsolver%ilevel), istep, istep+1, &
              rsolver%p_RdiscreteBC(rsolver%ilevel))

          call lsysbl_blockMatVec (rsolver%p_RspaceMatrices(rsolver%ilevel), &
              rsolver%rspaceTemp2, rsolver%rspaceTemp1, -1.0_DP, 1.0_DP)
        end if
        
        ! Assemble diagonal submatrix of that timestep.
        call stmv_getSubmatrix (rsolver%p_rmatrix, istep, istep, &
            rsolver%p_RspaceMatrices(rsolver%ilevel))
            
        ! Apply the boundary conditions to the matrix
        call stmv_implementDefBCSubmatrix (rsolver%p_rmatrix%p_rboundaryCond, &
            rsolver%p_RspaceMatrices(rsolver%ilevel), istep, istep, &
            rsolver%p_RdiscreteBC(rsolver%ilevel))
        
        ! Subtract the diagonal.
        call sptivec_getTimestepData(rsolver%rspaceTimeTemp1,istep,rsolver%rspaceTemp2)
        
        call lsysbl_blockMatVec (rsolver%p_RspaceMatrices(rsolver%ilevel), &
            rsolver%rspaceTemp2, rsolver%rspaceTemp1, -1.0_DP, 1.0_DP)
        
        ! Remove the dual part from the matrix
        !call stmv_reduceDiagToPrimal (rsolver%p_RspaceMatrices(rsolver%ilevel))
        
        ! Apply the space solver
        call linsol_initData(rsolver%p_rspaceSolver,ierror)
        if (ierror .ne. 0) then
          call output_line ("Spatial system cannot be factorised in timestep "//&
              trim(sys_siL(istep,10)))
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
      
    end do
    
    call linsol_doneStructure(rsolver%p_rspaceSolver,ierror)

    ! Write back the output vector.
    call sptivec_vectorLinearComb (rsolver%rspaceTimeTemp1,rd,rsolver%domega,0.0_DP)

  end subroutine

  ! ***************************************************************************

  subroutine stls_initBlockFBGS2 (rsolver,rspaceTimeHierarchy,ilevel,drelax,&
      cspaceSolverType,RmatVecTempl)
  
  ! Initialise a block GS correction solver, working on the coupled solution.
  
  ! Solver structure to be initialised
  type(t_spacetimelinsol), intent(out) :: rsolver
  
  ! Underlying space-time hierarchy
  type(t_spacetimeHierarchy), intent(in), target :: rspaceTimeHierarchy
  
  ! Level of the solver  in the space-time hierarchy
  integer, intent(in) :: ilevel

  ! Relaxation parameter
  real(DP), intent(in) :: drelax

  ! OPTIONAL: Identifier for a linear solver in space that might be used
  ! for preconditioning.
  ! =0: UMFPACK
  integer, intent(in) :: cspaceSolverType

  ! OPTINAL: Spatial matrix templates. Necessary if a matrix-based
  ! preconditioner in space is used.
  type(t_matvecTemplates), dimension(:), intent(in), target :: RmatVecTempl

    ! Basic initialisation
    call stls_init(rsolver,STLS_TYPE_FBGS2,rspaceTimeHierarchy,ilevel,&
        .false.,.false.,cspaceSolverType=cspaceSolverType,RmatVecTempl=RmatVecTempl)
        
    rsolver%drelax = drelax
  
  end subroutine

  ! ***************************************************************************

  recursive subroutine stls_precondBlockFBGS2 (rsolver, rd)
  
  ! General preconditioning to a defect vector rd
  
  ! Solver structure
  type(t_spacetimelinsol), intent(inout) :: rsolver
  
  ! Defect vector to apply preconditioning to.
  type(t_spaceTimeVector), intent(inout) :: rd

    ! local variables
    integer :: istep, ierror,jstep
    real(DP), dimension(:), pointer :: p_Dx,p_Dd,p_Ddata
    
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
          call stmv_getSubmatrix (rsolver%p_rmatrix, istep, istep-1, &
              rsolver%p_RspaceMatrices(rsolver%ilevel))

          call stmv_implementDefBCSubmatrix (rsolver%p_rmatrix%p_rboundaryCond, &
              rsolver%p_RspaceMatrices(rsolver%ilevel), istep, istep-1, &
              rsolver%p_RdiscreteBC(rsolver%ilevel))

          call lsysbl_blockMatVec (rsolver%p_RspaceMatrices(rsolver%ilevel), &
              rsolver%rspaceTemp2, rsolver%rspaceTemp1, -1.0_DP, 1.0_DP)
        end if
        
        ! Subtract the diagonal, set up the preconditioner matrix
        
        ! Assemble diagonal submatrix of that timestep.
        call stmv_getSubmatrix (rsolver%p_rmatrix, istep, istep, &
            rsolver%p_RspaceMatrices(rsolver%ilevel))
            
        ! Apply the boundary conditions to the matrix
        call stmv_implementDefBCSubmatrix (rsolver%p_rmatrix%p_rboundaryCond, &
            rsolver%p_RspaceMatrices(rsolver%ilevel), istep, istep, &
            rsolver%p_RdiscreteBC(rsolver%ilevel))

        ! Subtract the diagonal.
        call sptivec_getTimestepData(rsolver%rspaceTimeTemp1,istep,rsolver%rspaceTemp2)
        
        call lsysbl_blockMatVec (rsolver%p_RspaceMatrices(rsolver%ilevel), &
            rsolver%rspaceTemp2, rsolver%rspaceTemp1, -1.0_DP, 1.0_DP)
        
        ! Apply the space solver
        call linsol_initData(rsolver%p_rspaceSolver,ierror)
        if (ierror .ne. 0) then
          call output_line ("Spatial system cannot be factorised in timestep "//&
              trim(sys_siL(istep,10)))
              
          call matio_writeBlockMatrixHR (rsolver%p_RspaceMatrices(rsolver%ilevel), &
              "matrix", .true., 0, "matrix.txt", "(E10.3)") 
              
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
        
        if (istep .lt. rd%NEQtime) then
          ! Load the previous timestep.
          call sptivec_getTimestepData(rsolver%rspaceTimeTemp1,istep+1,rsolver%rspaceTemp2)
          
          ! Subtract the primal.
          call stmv_getSubmatrix (rsolver%p_rmatrix, istep, istep+1, &
              rsolver%p_RspaceMatrices(rsolver%ilevel))

          call stmv_implementDefBCSubmatrix (rsolver%p_rmatrix%p_rboundaryCond, &
              rsolver%p_RspaceMatrices(rsolver%ilevel), istep, istep+1, &
              rsolver%p_RdiscreteBC(rsolver%ilevel))

          call lsysbl_blockMatVec (rsolver%p_RspaceMatrices(rsolver%ilevel), &
              rsolver%rspaceTemp2, rsolver%rspaceTemp1, -1.0_DP, 1.0_DP)
        end if
        
        ! Assemble diagonal submatrix of that timestep.
        call stmv_getSubmatrix (rsolver%p_rmatrix, istep, istep, &
            rsolver%p_RspaceMatrices(rsolver%ilevel))
            
        ! Apply the boundary conditions to the matrix
        call stmv_implementDefBCSubmatrix (rsolver%p_rmatrix%p_rboundaryCond, &
            rsolver%p_RspaceMatrices(rsolver%ilevel), istep, istep, &
            rsolver%p_RdiscreteBC(rsolver%ilevel))
        
        ! Subtract the diagonal.
        call sptivec_getTimestepData(rsolver%rspaceTimeTemp1,istep,rsolver%rspaceTemp2)
        
        call lsysbl_blockMatVec (rsolver%p_RspaceMatrices(rsolver%ilevel), &
            rsolver%rspaceTemp2, rsolver%rspaceTemp1, -1.0_DP, 1.0_DP)
        
        ! Apply the space solver
        call linsol_initData(rsolver%p_rspaceSolver,ierror)
        if (ierror .ne. 0) then
          call output_line ("Spatial system cannot be factorised in timestep "//&
              trim(sys_siL(istep,10)))
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
    call stls_precondMultigridInternal (rsolver, rsolver%ilevel)

    ! Put the (weighted) solution into rd as returm value
    call sptivec_vectorLinearComb (rsolver%p_Rvectors1(rsolver%ilevel),rd,&
        rsolver%domega,0.0_DP)

  end subroutine

  ! ***************************************************************************

  recursive subroutine stls_precondMultigridInternal (rsolver, ilevel)
  
  ! General preconditioning to a defect vector in rsolver%p_Rvectors2.
  ! Result written into rsolver%p_Rvectors1.
  
  ! Solver structure
  type(t_spacetimelinsol), intent(inout) :: rsolver
  
  ! Current level
  integer, intent(in) :: ilevel
  character (len=5) :: sstring
  real(DP) :: dfactor
  
    ! local variables
    integer :: ite,nite
    real(DP) :: dresInit, dresCurrent, drho, drhoAsymp
    real(DP), dimension(3) :: DlastResiduals
  
    ! Are we on the level of the coarse grid solver?
    if (ilevel .eq. rsolver%p_rcoarsegridsolver%ilevel) then
      ! Call the coarse grid solver.
      call stls_precondDefect (rsolver%p_rcoarsegridsolver,rsolver%p_Rvectors2(ilevel))
      
      ! Put the result to the solution.
      call sptivec_copyVector (rsolver%p_Rvectors2(ilevel),rsolver%p_Rvectors1(ilevel))
    else
      ! Do the iteration.
      
      ! p_Rvectors1 = solution
      ! p_Rvectors2 = rhs
      ! p_Rvectors3 = temp vector
      
      ! Clear the solution vector p_Rvectors1.
      call sptivec_clearVector (rsolver%p_Rvectors1(ilevel))

      ! Max. number of iterations.
      nite = rsolver%nmaxiterations
      
      ! If we are on a lower level, do as many iterations as prescribed by the cycle
      ! Otherwise, calculate the initial residuum for residuum checks.
      if (ilevel .lt. rsolver%ilevel) then
        nite = 1
      else
        ! Initial residuum.
        dresInit = sptivec_vectorNorm (rsolver%p_Rvectors2(ilevel),LINALG_NORML2)
        dresCurrent = dresInit
        DlastResiduals (:) = dresInit
      end if
      
      do ite = 1,nite
      
        if (ilevel .eq. rsolver%ilevel) then
          
          ! Stopping criterion check on the max. level.
        
          if (rsolver%ioutputlevel .ge. 2) then
            call output_line("Space-Time Multigrid: Step "//trim(sys_siL(ite-1,10))//", Residuum = "//&
                trim(sys_sdEP(dresCurrent,15,7)))
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
          end if
          
        end if

!        ! RHS      
!        call stpp_postproc (rsolver%p_rmatrix%p_rphysics,rsolver%p_Rvectors2(ilevel))
!        print *,"RHS"
!        read (*,*)

!        ! Temp. defect
!        call sptivec_copyVector (rsolver%p_Rvectors2(ilevel),rsolver%p_Rvectors3(ilevel))
!        call stmv_matvec (rsolver%p_rmatrix, &
!            rsolver%p_Rvectors1(ilevel),rsolver%p_Rvectors3(ilevel), -1.0_DP, 1.0_DP)
!        call spop_applyBC (rsolver%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, rsolver%p_Rvectors3(ilevel))
!        call stpp_postproc (rsolver%p_rmatrix%p_rphysics,rsolver%p_Rvectors3(ilevel))
!        call stpp_printDefectSubnormsDirect(rsolver%p_Rvectors3(ilevel))
!        print *,"current defect"
!        read (*,*)
                
        ! Presmoothing
        if (associated(rsolver%p_Rpresmoothers)) then
          if (rsolver%p_Rpresmoothers(ilevel)%csolverType .ne. STLS_TYPE_NONE) then
            call stls_solveAdaptively (rsolver%p_rpresmoothers(ilevel), &
                rsolver%p_Rvectors1(ilevel),rsolver%p_Rvectors2(ilevel),rsolver%p_Rvectors3(ilevel))
          end if
        end if

!        ! Temp. defect
!        call sptivec_copyVector (rsolver%p_Rvectors2(ilevel),rsolver%p_Rvectors3(ilevel))
!        call stmv_matvec (rsolver%p_rmatrix, &
!            rsolver%p_Rvectors1(ilevel),rsolver%p_Rvectors3(ilevel), -1.0_DP, 1.0_DP)
!        call spop_applyBC (rsolver%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, rsolver%p_Rvectors3(ilevel))
!        call stpp_postproc (rsolver%p_rmatrix%p_rphysics,rsolver%p_Rvectors3(ilevel))
!        call stpp_printDefectSubnormsDirect(rsolver%p_Rvectors3(ilevel))
!        print *,"current defect after sm."
!        read (*,*)
            
!        ! DEBUG: Real solution.
!        call sptivec_copyVector (rsolver%p_Rvectors2(ilevel),rsolver%p_Rvectors3(ilevel))
!        call stmv_matvec (rsolver%p_rmatrix, &
!            rsolver%p_Rvectors1(ilevel),rsolver%p_Rvectors3(ilevel), -1.0_DP, 1.0_DP)
!        call spop_applyBC (rsolver%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, rsolver%p_Rvectors3(ilevel))
!        call stls_precondDefect (rsolver%p_Rpreconditioners(ilevel),rsolver%p_Rvectors3(ilevel))
!        print *,"real sol."
!        call stpp_postproc (rsolver%p_rmatrix%p_rphysics,rsolver%p_Rvectors3(ilevel))
!        read (*,*)

        ! Build the defect into the temp vector.
        call sptivec_copyVector (rsolver%p_Rvectors2(ilevel),rsolver%p_Rvectors3(ilevel))
        call stmv_matvec (rsolver%p_rmatrix, &
            rsolver%p_Rvectors1(ilevel),rsolver%p_Rvectors3(ilevel), -1.0_DP, 1.0_DP)
        call spop_applyBC (rsolver%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, rsolver%p_Rvectors3(ilevel))
        
        !call stpp_postproc (rsolver%p_rmatrix%p_rphysics,rsolver%p_Rvectors3(ilevel))

        ! Restriction.
        call sptipr_performRestriction (rsolver%p_rspaceTimeProjection,ilevel,&
            rsolver%p_Rvectors2(ilevel-1),rsolver%p_Rvectors3(ilevel),&
            rsolver%p_RspaceVectors(ilevel-1),rsolver%p_RspaceVectors(ilevel))
        
        ! Boundary conditions.
        call spop_applyBC (rsolver%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, rsolver%p_Rvectors2(ilevel))

        ! Solve on the lower level
        call stls_precondMultigridInternal (rsolver, ilevel-1)
        
        ! Prolongation
        call sptipr_performProlongation (rsolver%p_rspaceTimeProjection,ilevel,&
            rsolver%p_Rvectors1(ilevel-1),rsolver%p_Rvectors3(ilevel),&
            rsolver%p_RspaceVectors(ilevel-1),rsolver%p_RspaceVectors(ilevel))
        
        ! Boundary conditions.
        call spop_applyBC (rsolver%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, rsolver%p_Rvectors3(ilevel))
        
!        print *,"cgr sol."
!        call stpp_postproc (rsolver%p_rmatrix%p_rphysics,rsolver%p_Rvectors3(ilevel))
!        read (*,*) 
        dfactor = 1.0_DP
!        if (sstring .eq. "y") then
!          print *,"ok"
!          dfactor = 2.0_DP
!        end if
        
        ! Coarse grid correction
        call sptivec_vectorLinearComb (rsolver%p_Rvectors3(ilevel),rsolver%p_Rvectors1(ilevel),&
          dfactor,1.0_DP)

        !call stpp_postproc (rsolver%p_rmatrix%p_rphysics,rsolver%p_Rvectors1(ilevel))

!          call sptivec_copyVector (rsolver%p_Rvectors2(ilevel),rsolver%p_Rvectors3(ilevel))
!          call stmv_matvec (rsolver%p_rmatrix, &
!              rsolver%p_Rvectors1(ilevel),rsolver%p_Rvectors3(ilevel), -1.0_DP, 1.0_DP)
!          call spop_applyBC (rsolver%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, rsolver%p_Rvectors3(ilevel))
!
!        call stpp_postproc (rsolver%p_rmatrix%p_rphysics,rsolver%p_Rvectors3(ilevel))
        
!        if (ilevel .eq. rsolver%ilevel) then
!          call output_line ("Defect before smoothing:")
!          call sptivec_copyVector (rsolver%p_Rvectors2(ilevel),rsolver%p_Rvectors3(ilevel))
!          call stmv_matvec (rsolver%p_rmatrix, &
!              rsolver%p_Rvectors1(ilevel),rsolver%p_Rvectors3(ilevel), -1.0_DP, 1.0_DP)
!          call spop_applyBC (rsolver%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, rsolver%p_Rvectors3(ilevel))
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
          end if
        end if
        
        !call stpp_postproc (rsolver%p_rmatrix%p_rphysics,rsolver%p_Rvectors1(ilevel))

        ! New defect on the max. level for residuum checks.
        if (ilevel .eq. rsolver%ilevel) then
          call sptivec_copyVector (rsolver%p_Rvectors2(ilevel),rsolver%p_Rvectors3(ilevel))
          call stmv_matvec (rsolver%p_rmatrix, &
              rsolver%p_Rvectors1(ilevel),rsolver%p_Rvectors3(ilevel), -1.0_DP, 1.0_DP)
          call spop_applyBC (rsolver%p_rmatrix%p_rboundaryCond, SPOP_DEFECT, rsolver%p_Rvectors3(ilevel))
          
!          call output_line ("Defect after smoothing:")
!          call stpp_printDefectSubnormsDirect (rsolver%p_Rvectors3(ilevel))
          
          dresCurrent = sptivec_vectorNorm (rsolver%p_Rvectors3(ilevel),LINALG_NORML2)
          
          ! Shift the residuals
          DlastResiduals (1:size(DlastResiduals)-1) = DlastResiduals(2:)
          DlastResiduals (size(DlastResiduals)) = dresCurrent
        end if
        
      end do
      
      if (ilevel .eq. rsolver%ilevel) then
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
          call output_lbrk()
          call output_line("Iterations/Conv. rate: "//&
              trim(sys_siL(ite-1,10))//"/"//trim(sys_sdEP(drho,15,7)))
        end if
      end if
      
    end if

  end subroutine

end module
