!##############################################################################
!# ****************************************************************************
!# <name> newtoniterationlinear </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module realises the linear subsolver of the Newton iteration in the 
!# control space of the optimal control problem.
!# </purpose>
!##############################################################################

module newtoniterationlinear

  use fsystem
  use genoutput
  use paramlist
  
  use spatialdiscretisation
  use timediscretisation
  use linearsystemscalar
  use linearsystemblock
  
  use scalarpde
  use linearformevaluation
  use bilinearformevaluation
  use feevaluation2
  use blockmatassemblybase
  use blockmatassembly
  use collection
  
  use spacetimevectors
  use analyticsolution
  
  use structuresdiscretisation
  use structuresoptcontrol
  use structuresgeneral
  use structuresoptflow
  use structuresoperatorasm
  use assemblytemplates
  
  use spacematvecassembly
  use spacelinearsolver
  use spacesolver
  
  use kktsystemspaces
  use kktsystem
  use kktsystemhierarchy
  
  use spacetimehierarchy
  use spacetimeinterlevelprojection
  
  use structuresnewton
  
  implicit none
  
  private

!<constants>

!<constantblock description = "Solver identifier">

  ! Richardson iteration
  integer, parameter, public :: NLIN_SOLVER_RICHARDSON = 0
  
  ! Multigrid iteration
  integer, parameter, public :: NLIN_SOLVER_MULTIGRID = 1

!</constantblock>

!</constants>

  public :: t_linsolMultigrid
  public :: t_linsolParameters

!<type>

  !<typeblock>
  
  ! Multigrid parameters and structures
  type t_linsolMultigrid

    ! Identifies the cycle.
    ! =0: F-cycle,
    ! =1: V-cycle,
    ! =2: W-cycle
    integer :: ccycle = 0

    ! Specifies the smoother.
    ! =0: Damped Richardson iteration
    integer :: csmoother = 0

    ! Specifies the coarse grid solver.
    ! =0: Damped Richardson iteration
    integer :: ccoarseGridSolver = 0

    ! Number of presmoothing steps
    integer :: nsmPre = 0

    ! Number of postsmoothing steps
    integer :: nsmPost = 4
    
    ! Coarse grid correction damping parameter
    real(DP) :: dcoarseGridCorrectionWeight = 1.0_DP
    
    ! Section in the parameter list with parameters for the smoother
    character(LEN=PARLST_MLSECTION) :: ssectionSmoother
    
    ! Section in the parameter list with parameters for the coarse grid solver
    character(LEN=PARLST_MLSECTION) :: ssectionCoarseGridSolver
    
    ! Link to the parameter list
    type(t_parlist), pointer :: p_rparList => null()
    
    ! For every level in the hierarchy, a linear solver structure
    ! that configures the smoother. The solver structure at level 1
    ! defines the coarse grid solver
    type(t_linsolParameters), dimension(:), pointer :: p_rsubSolvers => null()
  
    ! KKT system hierarchy for solutions on the different levels
    type(t_kktsystemDirDerivHierarchy) :: rsolutionHier

    ! KKT system hierarchy for defect vectors on the different levels
    type(t_kktsystemDirDerivHierarchy) :: rdefectHier

    ! KKT system hierarchy for RHS vectors on the different levels
    type(t_kktsystemDirDerivHierarchy) :: rrhsHier
  
    ! Cycle counter
    integer, dimension(:), pointer :: p_IcycleCount => null()

  end type
  
  !</typeblock>
  
  !<typeblock>
  
  ! Contains the parameters of the Newton iteration.
  type t_linsolParameters
  
    ! <!-- --------------------------------------- -->
    ! <!-- GENRERAL PARAMETERS AND SOLVER SETTINGS -->
    ! <!-- --------------------------------------- -->

    ! Specifies the type of the solver.
    ! =0: Damped Richardson iteration with damping parameter omega.
    ! =1: Space-time multigrid method applied to the underlying 
    !     space-time hierarchy
    integer :: csolverType = 0

    ! General newton parameters
    type(t_newtonPrecParameters) :: rprecParameters

    ! Parameters of the OptFlow solver
    type(t_settings_optflow), pointer :: p_rsettingsSolver => null()
  
    ! <!-- ----------------------------- -->
    ! <!-- SUBSOLVERS AND OTHER SETTINGS -->
    ! <!-- ----------------------------- -->

    ! Hierarchy of solvers in space for all levels.
    ! Linearised primal equation.
    type(t_spaceSolverHierarchy), pointer :: p_rsolverHierPrimalLin => null()

    ! Hierarchy of solvers in space for all levels.
    ! Linearised dual equation.
    type(t_spaceSolverHierarchy), pointer :: p_rsolverHierDualLin => null()

    ! Defines a policy how to generate the initial condition of a timestep.
    ! =0: Always take zero
    ! =1: Propagate the solution of the previous/next timestep to the
    !     current one. (Default)
    ! =2: Take the solution of the last space-time iteration
    ! Warning: Avoid the setting cspatialInitCondPolicy = 0/1 if the linear
    ! solver in space does not solve up to a high accuracy.
    ! Otherwise, the solution of the linearised primal equation in the first
    ! timestep is not better than the stopping criterion of the space solver
    ! and thus, the complete space-time solution will not be better!
    ! Remember, there is no nonlinear loop around each timestep!
    ! That means, the stopping criterion of the linear space-time solver is
    ! overwritten by the stopping criterion of the space-solver, which
    ! is usually an undesired behaviour: Although the linear space-time solver
    ! solves up to, e.g., 1E-15, the total solution of the space-time solver
    ! is not better than depsrel(linear space-solver)!!!
    integer :: cspatialInitCondPolicy = SPINITCOND_PREVITERATE
    
    ! Parameters for the multigrid solver
    type(t_linsolMultigrid), pointer :: p_rlinsolMultigrid => null()
    
    ! <!-- -------------- -->
    ! <!-- TEMPORARY DATA -->
    ! <!-- -------------- -->
    
    ! Underlying KKT system hierarchy that defines the shape of the solutions.
    type(t_kktsystemHierarchy), pointer :: p_rkktsystemHierarchy => null()

    ! Projection hierarchy for the interlevel projection in space/time, primal space.
    type(t_sptiProjHierarchyBlock), pointer :: p_rprjHierSpaceTimePrimal => null()

    ! Projection hierarchy for the interlevel projection in space/time, dual space.
    type(t_sptiProjHierarchyBlock), pointer :: p_rprjHierSpaceTimeDual => null()

    ! Projection hierarchy for the interlevel projection in space/time, control space.
    type(t_sptiProjHierarchyBlock), pointer :: p_rprjHierSpaceTimeControl => null()
    
    ! Temporary memory on all levels for calculations
    type(t_controlSpace), dimension(:), pointer :: p_RtempVectors => null()

  end type
  
  !</typeblock>
  
!</types>

  ! Apply preconditioning of a defect with the Newton preconditioner
  ! in the control space.
  public :: newtonlin_precond
  
  ! Basic initialisation of the Newton solver
  public :: newtonlin_init
  
  ! Structural initialisation
  public :: newtonlin_initStructure
  
  ! Final initialisation
  public :: newtonlin_initData
  
  ! Cleanup of data
  public :: newtonlin_doneData

  ! Cleanup of structures
  public :: newtonlin_doneStructure

  ! Final cleanup
  public :: newtonlin_done

contains

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_getResidual (&
      rlinsolParam,rkktsystemDirDeriv,rrhs,rresidual,dres,iresnorm)
  
!<description>
  ! Calculates the residual in the linearised control equation.
!</description>
  
!<inputoutput>
  ! Parameters for the iteration.
  type(t_linsolParameters), intent(inout) :: rlinsolParam
  
  ! Structure defining the linearised KKT system.
  type(t_kktsystemDirDeriv), intent(inout), target :: rkktsystemDirDeriv
  
  ! Right-hand side of the linearised control equation
  type(t_controlSpace), intent(inout), target :: rrhs
  
  ! On output, this structure receives a representation of the search
  ! direction / residual in the Newton iteration
  type(t_controlSpace), intent(inout) :: rresidual

  ! type of norm. A LINALG_NORMxxxx constant.
  integer, intent(in) :: iresnorm
!</inputoutput>

!<output>
  ! L2-Norm of the residual
  real(DP), intent(out) :: dres
!</output>

!</subroutine>

    ! -------------------------------------------------------------
    ! Step 1: Solve the primal and dual system.
    ! -------------------------------------------------------------

    if (rlinsolParam%rprecParameters%ioutputLevel .ge. 2) then
      call output_line ("Linear space-time Residual: Solving the primal equation")
    end if

    ! Solve the primal equation, update the primal solution.
    call kkt_solvePrimalDirDeriv (rkktsystemDirDeriv,&
        rlinsolParam%p_rsolverHierPrimalLin,rlinsolParam%cspatialInitCondPolicy)
    
    if (rlinsolParam%rprecParameters%ioutputLevel .ge. 2) then
      call output_line ("Linear space-time Residual: Solving the dual equation")
    end if

    ! Solve the dual equation, update the dual solution.
    call kkt_solveDualDirDeriv (rkktsystemDirDeriv,&
        rlinsolParam%p_rsolverHierDualLin,rlinsolParam%cspatialInitCondPolicy)

    ! -------------------------------------------------------------
    ! Step 2: Calculate the residual
    ! -------------------------------------------------------------

    ! Take the solution of the linearised primal/dual system and
    ! calculate the residual in the control space.
    call kkt_calcControlResDirDeriv (rkktsystemDirDeriv,rrhs,rresidual,dres,iresnorm)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_smoothCorrection (rsmootherParams,nsm,&
      rkktsystemDirDeriv,rrhs,rtempVector)
  
!<description>
  ! Applies the space-time smoother.
!</description>

!<input>
  ! Number of smoothing steps
  integer, intent(in) :: nsm
!</input>

!<inputoutput>
  ! Parameters of the smoother
  type(t_linsolParameters), intent(inout) :: rsmootherParams
  
  ! Structure defining the linearised KKT system.
  ! The linearised control in this structure is used as
  ! initial and target vector.
  ! The control on the maximum level receives the result.
  type(t_kktsystemDirDeriv), intent(inout), target :: rkktsystemDirDeriv
  
  ! Right-hand side of the linearised control equation
  type(t_controlSpace), intent(inout), target :: rrhs

  ! Temporary control vector
  type(t_controlSpace), intent(inout) :: rtempVector
!</inputoutput>

!</subroutine>

    ! Configure the solver to apply exactly nsm steps
    if (nsm .le. 0) return
    rsmootherParams%rprecParameters%nminiterations = nsm
    rsmootherParams%rprecParameters%nmaxiterations = nsm
    
    ! Call the smoother.
    select case (rsmootherParams%csolverType)
    
    ! --------------------------------------
    ! Richardson iteration
    ! --------------------------------------
    case (NLIN_SOLVER_RICHARDSON)
      ! Call the Richardson iteration to calculate an update
      ! for the control.
      call newtonlin_richardson (rsmootherParams,rkktsystemDirDeriv,&
          rrhs,rtempVector)
          
    end select
    

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_richardson (rlinsolParam,rkktsystemDirDeriv,rrhs,rtempVector)
  
!<description>
  ! Applies a Richardson iteration in the control space
  ! for the linearised KKT system.
!</description>

!<inputoutput>
  ! Parameters for the iteration.
  ! The output parameters are changed according to the iteration.
  type(t_linsolParameters), intent(inout) :: rlinsolParam

  ! Structure defining the linearised KKT system.
  ! The linearised control in this structure is used as
  ! initial and target vector.
  ! The control on the maximum level receives the result.
  type(t_kktsystemDirDeriv), intent(inout), target :: rkktsystemDirDeriv
  
  ! Right-hand side of the linearised control equation
  type(t_controlSpace), intent(inout), target :: rrhs

  ! Temporary control vector
  type(t_controlSpace), intent(inout) :: rtempVector
!</inputoutput>

!</subroutine>

    ! Richardson is a do loop...
    rlinsolParam%rprecParameters%niterations = 0

    do while (.true.)
    
      ! -------------------------------------------------------------
      ! Get the current residual / search direction
      ! -------------------------------------------------------------

      ! Compute the residual and its norm.
      call newtonlin_getResidual (rlinsolParam,rkktsystemDirDeriv,rrhs,&
          rtempVector,rlinsolParam%rprecParameters%dresFinal,&
          rlinsolParam%rprecParameters%iresnorm)

      if (rlinsolParam%rprecParameters%niterations .eq. 0) then
        ! Remember the initial residual
        rlinsolParam%rprecParameters%dresInit = rlinsolParam%rprecParameters%dresFinal
      end if

      if (rlinsolParam%rprecParameters%ioutputLevel .ge. 2) then
        call output_line ("Space-time Richardson: Iteration "// &
            trim(sys_siL(rlinsolParam%rprecParameters%niterations,10))// &
            ", ||res(u)|| = "// &
            trim(sys_sdEL(rlinsolParam%rprecParameters%dresFinal,10)))
      end if

      ! -------------------------------------------------------------
      ! Check for convergence
      ! -------------------------------------------------------------
      if (newtonlin_checkConvergence(rlinsolParam%rprecParameters)) exit
      
      ! -------------------------------------------------------------
      ! Check for divergence
      ! -------------------------------------------------------------
      if (newtonlin_checkDivergence(rlinsolParam%rprecParameters)) exit
      
      ! -------------------------------------------------------------
      ! Check other stopping criteria
      ! -------------------------------------------------------------
      if (newtonlin_checkIterationStop(rlinsolParam%rprecParameters)) exit
      
      ! -------------------------------------------------------------
      ! Update of the solution
      ! -------------------------------------------------------------

      ! Add the residual (damped) to the current control.
      ! This updates the current control.
      call kktsp_controlLinearComb (rtempVector,&
          rlinsolParam%rprecParameters%domega,&
          rkktsystemDirDeriv%p_rcontrolLin,1.0_DP)
    
      ! -------------------------------------------------------------
      ! Proceed with the iteration
      ! -------------------------------------------------------------
      ! Next iteration
      rlinsolParam%rprecParameters%niterations = &
          rlinsolParam%rprecParameters%niterations + 1
    
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_multigrid (rlinsolParam,rkktsysDirDerivHierarchy,rrhs,ilevel)
  
!<description>
  ! Applies a Multigrid iteration in the control space
  ! for the linearised KKT system.
!</description>

!<input>
  ! Level in the hierarchy.
  integer, intent(in) :: ilevel
!</input>

!<inputoutput>
  ! Parameters for the iteration.
  ! The output parameters are changed according to the iteration.
  type(t_linsolParameters), intent(inout) :: rlinsolParam

  ! Hierarchy of directional derivatives of the KKT system.
  ! The control on the maximum level receives the result.
  type(t_kktsystemDirDerivHierarchy), intent(inout) :: rkktsysDirDerivHierarchy

  ! Right-hand side of the linearised control equation
  type(t_controlSpace), intent(inout), target :: rrhs
!</inputoutput>

!</subroutine>

    ! Dall the MG driver on the maximum level.
    call newtonlin_multigridDriver (&
        rlinsolParam,ilevel,rlinsolParam%p_rkktsystemHierarchy%nlevels,&
        rkktsysDirDerivHierarchy,rrhs)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  recursive subroutine MG_initCycles (ccycle,ilevel,nlevels,IcycleCount)
  
!<description>
  ! Initialises the cycle counter on level 1..ilevel.
!</description>

!<input>
  ! Cycle type.
  ! =0: F-cycle, =1: V-cycle, =2: W-cycle
  integer, intent(in) :: ccycle

  ! Current level
  integer, intent(in) :: ilevel
  
  ! Maximum level
  integer, intent(in) :: nlevels
!</input>

!<inputoutput>
  ! Cycle counter.
  integer, dimension(:), intent(inout) :: IcycleCount
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i

    if (ccycle .eq. 0) then
      ! F-cycle. Cycle counter is =2 on all levels
      ! if ilevel=nlevels and =1 otherwise.
      if (ilevel .eq. nlevels) then
        do i=1,ilevel
          IcycleCount(i) = 2
        end do
      else
        do i=1,ilevel
          IcycleCount(i) = 1
        end do
      end if
      
    else
      ! V-cyclke, W-cycle,...
      do i=1,ilevel
        IcycleCount(i) = ccycle
      end do
    end if    

  end subroutine

  ! ***************************************************************************

!<subroutine>

  recursive subroutine newtonlin_multigridDriver (rlinsolParam,ilevel,nlevels,&
      rkktsysDirDerivHierarchy,rrhs)
  
!<description>
  ! This is the actual multigrid driver which applies the
  ! multilevel iteration.
!</description>

!<input>
  ! Current level
  integer, intent(in) :: ilevel
  
  ! Maximum level in the hierarchy
  integer, intent(in) :: nlevels
!</input>

!<inputoutput>
  ! Parameters for the iteration.
  ! The output parameters are changed according to the iteration.
  type(t_linsolParameters), intent(inout), target :: rlinsolParam

  ! Hierarchy of directional derivatives of the KKT system.
  ! The control on level ilevel receives the result.
  type(t_kktsystemDirDerivHierarchy), intent(inout) :: rkktsysDirDerivHierarchy
  
  ! Right-hand side of the linearised control equation
  type(t_controlSpace), intent(inout), target :: rrhs
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_linsolMultigrid), pointer :: p_rlinsolMultigrid
    integer :: nminIterations,nmaxIterations,ite
    real(DP) :: dresInit,dresFinal
    type(t_kktsystemDirDeriv), pointer :: p_rkktSysRhsCoarse, p_rkktSysSolCoarse, p_rkktSysSolFine
    type(t_kktsystemDirDeriv), pointer :: p_rkktSysSolution, p_rkktSysDefect
    
    ! Get multigrid parameters.
    p_rlinsolMultigrid => rlinsolParam%p_rlinsolMultigrid

    ! On level 1 apply the coarse grid solver
    if (ilevel .eq. 1) then
      
      if (nlevels .eq. 1) then
        if (rlinsolParam%rprecParameters%ioutputLevel .ge. 1) then
          call output_line (&
              "Space-time Multigrid: Only one level. Switching back to 1-level solver.")
        end if
      end if

      ! Call the coarse grid solver.
      output_iautoOutputIndent = output_iautoOutputIndent + 2
      call newtonlin_precond (p_rlinsolMultigrid%p_rsubSolvers(ilevel),&
          rkktsysDirDerivHierarchy,rrhs,ilevel)
      output_iautoOutputIndent = output_iautoOutputIndent - 2
    
      rlinsolParam%rprecParameters%dresInit = &
          p_rlinsolMultigrid%p_rsubSolvers(ilevel)%rprecParameters%dresInit
      rlinsolParam%rprecParameters%dresFinal = &
          p_rlinsolMultigrid%p_rsubSolvers(ilevel)%rprecParameters%dresFinal
      rlinsolParam%rprecParameters%niterations = &
          p_rlinsolMultigrid%p_rsubSolvers(ilevel)%rprecParameters%niterations
    
      ! Finish
      return
    end if
    
    ! On level 2..nlevels, apply the smoothing and coarse grid correction process.

    if (ilevel .eq. nlevels) then
      ! Initialise the cycle counter for all levels.
      call MG_initCycles (&
          p_rlinsolMultigrid%ccycle,nlevels,nlevels,p_rlinsolMultigrid%p_IcycleCount)

      ! On level nlevels, apply as many iterations as configured.
      nminIterations = rlinsolParam%rprecParameters%nminIterations
      nmaxIterations = rlinsolParam%rprecParameters%nmaxIterations
    else
      ! On all the other levels, apply as many iterations as configured in the cycle.
      nminIterations = p_rlinsolMultigrid%p_IcycleCount(ilevel)
      nmaxIterations = p_rlinsolMultigrid%p_IcycleCount(ilevel)
    end if
    
    ! Apply the multigrid iteration

    call kkth_getKKTsystemDirDeriv (rkktsysDirDerivHierarchy,ilevel,p_rkktSysSolution)
    call kkth_getKKTsystemDirDeriv (p_rlinsolMultigrid%rrhsHier,ilevel-1,p_rkktSysRhsCoarse)
    call kkth_getKKTsystemDirDeriv (p_rlinsolMultigrid%rsolutionHier,ilevel-1,p_rkktSysSolCoarse)
    call kkth_getKKTsystemDirDeriv (p_rlinsolMultigrid%rsolutionHier,ilevel,p_rkktSysSolFine)
    call kkth_getKKTsystemDirDeriv (p_rlinsolMultigrid%rdefectHier,ilevel,p_rkktSysDefect)
    
    ite = 0
    call kkt_clearDirDeriv (p_rkktSysSolution)

    do
    
      if (ilevel .lt. nlevels) then
        ! On level < nlevels, do not calculate the next residual but
        ! exit immediately if the maximum number of iterations
        ! (= cycle count) is reached.
        if (ite .ge. nmaxIterations) exit
      end if
    
      ! Get the residual
      call newtonlin_getResidual (rlinsolParam,&
          p_rkktSysSolution,rrhs,&
          p_rkktSysDefect%p_rcontrolLin,&
          dresFinal,rlinsolParam%rprecParameters%iresnorm)

      if (ite .eq. 0) dresInit = dresFinal

      ! Some checks and output on the maximum level
      if (ilevel .eq. nlevels) then
      
        if (rlinsolParam%rprecParameters%ioutputLevel .ge. 2) then
          call output_line ("Space-time MG: Iteration "// &
              trim(sys_siL(ite,10))// &
              ", ||res(u)|| = "// &
              trim(sys_sdEL(dresFinal,10)))
        end if
        
        ! Save statistics on the maximum level
        rlinsolParam%rprecParameters%dresInit = dresInit
        rlinsolParam%rprecParameters%dresFinal = dresFinal
        rlinsolParam%rprecParameters%niterations = ite

        ! -------------------------------------------------------------
        ! Check for convergence
        ! -------------------------------------------------------------
        if (newtonlin_checkConvergence(rlinsolParam%rprecParameters)) exit
        
        ! -------------------------------------------------------------
        ! Check for divergence
        ! -------------------------------------------------------------
        if (newtonlin_checkDivergence(rlinsolParam%rprecParameters)) exit
        
        ! -------------------------------------------------------------
        ! Check other stopping criteria
        ! -------------------------------------------------------------
        if (newtonlin_checkIterationStop(rlinsolParam%rprecParameters)) exit
      
      end if
      
      ! Apply presmoothing
      if (p_rlinsolMultigrid%nsmpre .ne. 0) then
        if (rlinsolParam%rprecParameters%ioutputLevel .ge. 3) then
          call output_line (&
              "Space-time Multigrid: Invoking pre-smoother.")
        end if

        output_iautoOutputIndent = output_iautoOutputIndent + 2
        call newtonlin_smoothCorrection (p_rlinsolMultigrid%p_rsubSolvers(ilevel),&
              p_rlinsolMultigrid%nsmpre,&
              p_rkktSysSolution,rrhs,&
              rlinsolParam%p_rtempVectors(ilevel))
        output_iautoOutputIndent = output_iautoOutputIndent - 2

        ! Get the new residual -- it has changed due to the smoothing.
        call newtonlin_getResidual (rlinsolParam,&
            p_rkktSysSolution,rrhs,&
            p_rkktSysDefect%p_rcontrolLin,&
            dresFinal,rlinsolParam%rprecParameters%iresnorm)
      end if
      
      if (rlinsolParam%rprecParameters%ioutputLevel .ge. 3) then
        call output_line (&
            "Space-time Multigrid: Switching to level "//trim(sys_siL(ilevel-1,10))//".")
      end if

      ! Restriction of the defect
      call sptipr_performRestriction (&
          rlinsolParam%p_rprjHierSpaceTimeControl,ilevel,CCSPACE_CONTROL,&
          p_rkktSysRhsCoarse%p_rcontrolLin%p_rvectorAccess, &
          p_rkktSysDefect%p_rcontrolLin%p_rvectorAccess)
      
      if (ilevel-1 .eq. 1) then
        if (rlinsolParam%rprecParameters%ioutputLevel .ge. 3) then
          call output_line (&
              "Space-time Multigrid: Invoking coarse-grid solver.")
        end if
      end if

      ! Coarse grid correction
      output_iautoOutputIndent = output_iautoOutputIndent + 2
      call newtonlin_multigridDriver (rlinsolParam,ilevel-1,nlevels,&
          p_rlinsolMultigrid%rsolutionHier,&
          p_rkktSysRhsCoarse%p_rcontrolLin)
      output_iautoOutputIndent = output_iautoOutputIndent - 2
      
      if (rlinsolParam%rprecParameters%ioutputLevel .ge. 3) then
        call output_line (&
            "Space-time Multigrid: Switching to level "//trim(sys_siL(ilevel,10))//".")
      end if

      ! Prolongation of the correction.
      ! Prolongate to the "defect" vector on the current level and use
      ! this for the coarse grid correction. 
      ! Do not use the "solution" vector as target of the prolongation since
      ! that may be in use due to the multigrid iteration.
      call sptipr_performProlongation (rlinsolParam%p_rprjHierSpaceTimeControl,ilevel,&
          p_rkktSysSolCoarse%p_rcontrolLin%p_rvectorAccess,&
          p_rkktSysDefect%p_rcontrolLin%p_rvectorAccess)

      ! And the actual correction...
      call kktsp_controlLinearComb (&
          p_rkktSysDefect%p_rcontrolLin,&
          p_rlinsolMultigrid%dcoarseGridCorrectionWeight,&
          p_rkktSysSolution%p_rcontrolLin,&
          1.0_DP)

      ! Apply postsmoothing
      if (p_rlinsolMultigrid%nsmpost .ne. 0) then
        if (rlinsolParam%rprecParameters%ioutputLevel .ge. 3) then
          call output_line (&
              "Space-time Multigrid: Invoking post-smoother.")
        end if

        output_iautoOutputIndent = output_iautoOutputIndent + 2
        call newtonlin_smoothCorrection (p_rlinsolMultigrid%p_rsubSolvers(ilevel),&
              p_rlinsolMultigrid%nsmpost,&
              p_rkktSysSolution,rrhs,&
              rlinsolParam%p_rtempVectors(ilevel))
        output_iautoOutputIndent = output_iautoOutputIndent - 2
              
        ! Calculation of the new residual is done in the beginning of the loop.
      end if

      ! Next iteration
      ite = ite + 1

    end do
    
    ! Initialise the cycle counter for the next run.
    call MG_initCycles (&
        p_rlinsolMultigrid%ccycle,ilevel,nlevels,p_rlinsolMultigrid%p_IcycleCount)
        
  end subroutine

  ! ***************************************************************************

!<subroutine>

  recursive subroutine newtonlin_precond (rlinsolParam,rkktsysDirDerivHierarchy,rnewtonDir,ilevel)
  
!<description>
  ! Calculates the Newton search direction by applying the Newton
  ! preconditioner to the given residual rnewtonDir.
  ! This calculates a correction g with
  !      J''(u) g = d
  ! where d=rnewtonDir. On return of this routine, there is rnewtonDir=g.
!</description>

!<input>
  ! Level in the hierarchy. If not present, the maximum level is assumed.
  integer, intent(in), optional :: ilevel
!</input>

!<inputoutput>
  ! Parameters for the iteration.
  ! The output parameters are changed according to the iteration.
  type(t_linsolParameters), intent(inout) :: rlinsolParam

  ! Hierarchy of directional derivatives of the KKT system.
  ! Used for temporary calculations.
  type(t_kktsystemDirDerivHierarchy), intent(inout) :: rkktsysDirDerivHierarchy

  ! This structure contains the search direction to be
  ! preconditioned. Will be replaced by the precoditioned search direction.
  type(t_controlSpace), intent(inout) :: rnewtonDir
!</inputoutput>

!</subroutine>

    integer :: ilv
    type(t_kktsystemDirDeriv), pointer :: p_rkktsysDirDeriv
    
    ! The calculation is done on the topmost level
    ilv = rkktsysDirDerivHierarchy%nlevels
    if (present(ilevel)) ilv = ilevel
    
    p_rkktsysDirDeriv => rkktsysDirDerivHierarchy%p_RkktSysDirDeriv(ilv)

    ! For the calculation of the Newon search direction, we have to solve
    ! the linear system
    !
    !    J''(u) g = d
    !
    ! On the entry of this routine, there is d=rnewtonDir. When the routine
    ! is left, rnewtonDir receives g.
    !
    ! For this purpose, we need a linear solver based on J''(u) in the
    ! control space. For distributed control, this is multigrid-based,
    ! for boundary control (or general control with a limited number of
    ! control variables) a single-grid solver.
    !
    ! The system to solve is rather similar to the nonlinear system,
    ! but involves a couple of additional terms. It is still an iteration
    ! in u and thus, involves a forward and a backward solve for every 
    ! update -- probably on different refinement levels.
    !
    ! In the following, rnewtonDir is used as right-hand side of the
    ! equation. The solution is calculated in 
    !     rkktsystemDirDeriv%p_rcontrolLin
    ! which starts with zero.
    
    call kkt_clearDirDeriv (p_rkktsysDirDeriv)

    select case (rlinsolParam%csolverType)
    
    ! --------------------------------------
    ! Richardson iteration
    ! --------------------------------------
    case (NLIN_SOLVER_RICHARDSON)
      ! Call the Richardson iteration to calculate an update
      ! for the control.
      call newtonlin_richardson (rlinsolParam,p_rkktsysDirDeriv,&
          rnewtonDir,rlinsolParam%p_rtempVectors(ilv))
          
    ! --------------------------------------
    ! Multigrid iteration
    ! --------------------------------------
    case (NLIN_SOLVER_MULTIGRID)
      ! Invoke multigrid.
      call newtonlin_multigrid (rlinsolParam,rkktsysDirDerivHierarchy,rnewtonDir,ilv)
      
    end select
    
    ! Overwrite the rnewtonDir with the update.
    call kktsp_controlLinearComb (&
        p_rkktsysDirDeriv%p_rcontrolLin,1.0_DP,rnewtonDir,0.0_DP)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_init (rlinsolParam,rsettingsSolver,&
      rsolverHierPrimalLin,rsolverHierDualLin,rparlist,ssection)
  
!<description>
  ! Basic initialisation of the preconditioner.
!</description>

!<input>
  ! Parameters of the OptFlow solver
  type(t_settings_optflow), intent(in), target :: rsettingsSolver

  ! Hierarchy of solvers in space for all levels.
  ! Linearised primal equation.
  type(t_spaceSolverHierarchy), target :: rsolverHierPrimalLin
  
  ! Hierarchy of solvers in space for all levels.
  ! Linearised dual equation.
  type(t_spaceSolverHierarchy), target :: rsolverHierDualLin

  ! Parameter list that contains the data for the preconditioning.
  type(t_parlist), intent(in) :: rparlist
  
  ! Entry Section in the parameter list containing the data of the 
  ! preconditioner.
  character(len=*), intent(in) :: ssection
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_linsolParameters), intent(out) :: rlinsolParam
!</inputoutput>

!</subroutine>

    character(len=SYS_STRLEN) :: ssectionSpaceTimeMG
   
    ! Remember the settings for later use
    rlinsolParam%p_rsettingsSolver => rsettingsSolver
    rlinsolParam%p_rsolverHierPrimalLin => rsolverHierPrimalLin
    rlinsolParam%p_rsolverHierDualLin => rsolverHierDualLin
   
    ! Get the corresponding parameters from the parameter list.
    call newtonlin_initBasicParams (rlinsolParam%rprecParameters,ssection,rparlist)

    ! Solver type
    call parlst_getvalue_int (rparlist, ssection, "csolverType", &
        rlinsolParam%csolverType, rlinsolParam%csolverType)

    ! Generation of the initial condition
    call parlst_getvalue_int (rparlist, ssection, "cspatialInitCondPolicy", &
        rlinsolParam%cspatialInitCondPolicy, rlinsolParam%cspatialInitCondPolicy)
        
    select case (rlinsolParam%csolverType)
    
    ! -----------------------------------------------------
    ! Multigrid inítialisation
    ! -----------------------------------------------------
    case (NLIN_SOLVER_MULTIGRID)
    
      ! Get the section configuring the space-time MG
      call parlst_getvalue_string (rparlist, ssection, &
          "ssectionSpaceTimeMG", ssectionSpaceTimeMG, "SPACETIME-LINSOLVER",bdequote=.true.)
          
      ! Initialise the MG parameters
      allocate(rlinsolParam%p_rlinsolMultigrid)
      call newtonlin_initMG (rlinsolParam%p_rlinsolMultigrid,rparlist,ssectionSpaceTimeMG)
    
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_initMG (rlinsolMultigrid,rparlist,ssection)
  
!<description>
  ! Basic initialisation of the multigrid preconditioner.
!</description>

!<input>
  ! Parameter list that contains the data for the preconditioning.
  type(t_parlist), intent(in), target :: rparlist
  
  ! Entry Section in the parameter list containing the data of the 
  ! preconditioner.
  character(len=*), intent(in) :: ssection
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_linsolMultigrid), intent(out) :: rlinsolMultigrid
!</inputoutput>

!</subroutine>

    ! Read parameters from the DAT file.
    call parlst_getvalue_int (rparlist, ssection, "ccycle", &
        rlinsolMultigrid%ccycle, rlinsolMultigrid%ccycle)
    
    call parlst_getvalue_int (rparlist, ssection, "csmoother", &
        rlinsolMultigrid%csmoother, rlinsolMultigrid%csmoother)

    call parlst_getvalue_int (rparlist, ssection, "ccoarseGridSolver", &
        rlinsolMultigrid%ccoarseGridSolver, rlinsolMultigrid%ccoarseGridSolver)

    call parlst_getvalue_int (rparlist, ssection, "nsmPre", &
        rlinsolMultigrid%nsmPre, rlinsolMultigrid%nsmPre)
        
    call parlst_getvalue_int (rparlist, ssection, "nsmPost", &
        rlinsolMultigrid%nsmPost, rlinsolMultigrid%nsmPost)

    call parlst_getvalue_double (rparlist, ssection, "dcoarseGridCorrectionWeight", &
        rlinsolMultigrid%dcoarseGridCorrectionWeight, &
        rlinsolMultigrid%dcoarseGridCorrectionWeight)

    call parlst_getvalue_string (rparlist, ssection, "ssectionSmoother", &
        rlinsolMultigrid%ssectionSmoother,"SPACETIME-SMOOTHER",&
        bdequote=.true.)

    call parlst_getvalue_string (rparlist, ssection, "ssectionCoarseGridSolver", &
        rlinsolMultigrid%ssectionCoarseGridSolver,"SPACETIME-COARSEGRIDSOLVER",&
        bdequote=.true.)
        
    ! Remember the parameter list for later initialisations
    ! of the subsolvers.
    rlinsolMultigrid%p_rparList => rparList

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_initStructure (rlinsolParam,rkktsystemHierarchy,&
      rprjHierSpaceTimePrimal,rprjHierSpaceTimeDual,rprjHierSpaceTimeControl)
  
!<description>
  ! Structural initialisation of the preconditioner
!</description>

!<input>
  ! Defines the basic hierarchy of the solutions of the KKT system.
  ! This can be a 'template' structure, i.e., memory for the solutions
  ! in rkktsystemHierarchy does not have to be allocated.
  type(t_kktsystemHierarchy), intent(in), target :: rkktsystemHierarchy

  ! Projection hierarchy for the interlevel projection in space/time, primal space.
  type(t_sptiProjHierarchyBlock), intent(in), target :: rprjHierSpaceTimePrimal

  ! Projection hierarchy for the interlevel projection in space/time, dual space.
  type(t_sptiProjHierarchyBlock), intent(in), target :: rprjHierSpaceTimeDual

  ! Projection hierarchy for the interlevel projection in space/time, control space.
  type(t_sptiProjHierarchyBlock), intent(in), target :: rprjHierSpaceTimeControl
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_linsolParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>

    integer :: i,nlevels

    ! Remember the structure of the solutions.
    rlinsolParam%p_rkktsystemHierarchy => rkktsystemHierarchy
    
    ! Save the projection hierarchies which allows us to apply
    ! prolongation/restriction.
    rlinsolParam%p_rprjHierSpaceTimePrimal => rprjHierSpaceTimePrimal
    rlinsolParam%p_rprjHierSpaceTimeDual => rprjHierSpaceTimeDual
    rlinsolParam%p_rprjHierSpaceTimeControl => rprjHierSpaceTimeControl

    ! Prepare temporary vectors.
    nlevels = rlinsolParam%p_rkktsystemHierarchy%nlevels
    allocate(rlinsolParam%p_RtempVectors(nlevels))
    
    do i=1,nlevels
      call kkth_initControl (&
          rlinsolParam%p_RtempVectors(i),rkktsystemHierarchy,i)
    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_initDataMultigrid (rsolver,rlinsolMultigrid,rsolutionHierarchy)
  
!<description>
  ! Initialises a multigrid substructure.
!</description>

!<input>
  ! Solver structure of the main solver.
  type(t_linsolParameters), intent(in) :: rsolver

  ! Defines a hierarchy of the solutions of the KKT system.
  type(t_kktsystemHierarchy), intent(in), target :: rsolutionHierarchy
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_linsolMultigrid), intent(inout) :: rlinsolMultigrid
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ilevel

    ! Allocate the smoothers / coarse grid solvers
    allocate (rlinsolMultigrid%p_rsubSolvers(rsolutionHierarchy%nlevels))
    
    ! Propagate structures
    do ilevel = 1,rsolutionHierarchy%nlevels
      rlinsolMultigrid%p_rsubSolvers(ilevel)%p_rsolverHierPrimalLin => rsolver%p_rsolverHierPrimalLin
      rlinsolMultigrid%p_rsubSolvers(ilevel)%p_rsolverHierDualLin => rsolver%p_rsolverHierDualLin
      rlinsolMultigrid%p_rsubSolvers(ilevel)%p_rkktsystemHierarchy => rsolver%p_rkktsystemHierarchy
      rlinsolMultigrid%p_rsubSolvers(ilevel)%p_rprjHierSpaceTimePrimal => rsolver%p_rprjHierSpaceTimePrimal
      rlinsolMultigrid%p_rsubSolvers(ilevel)%p_rprjHierSpaceTimeDual => rsolver%p_rprjHierSpaceTimeDual
      rlinsolMultigrid%p_rsubSolvers(ilevel)%p_rprjHierSpaceTimeControl => rsolver%p_rprjHierSpaceTimeControl
      rlinsolMultigrid%p_rsubSolvers(ilevel)%p_RtempVectors => rsolver%p_RtempVectors
    end do
    
    ! Level 1 is the coarse grid solver
    select case (rlinsolMultigrid%ccoarseGridSolver)
    
    ! -------------------------------
    ! Richardson solver
    ! -------------------------------
    case (NLIN_SOLVER_RICHARDSON)
    
      ! Partial initialisation of the corresponding linear solver structure
      rlinsolMultigrid%p_rsubSolvers(1)%csolverType = NLIN_SOLVER_RICHARDSON
      call newtonlin_initBasicParams (&
          rlinsolMultigrid%p_rsubSolvers(1)%rprecParameters,&
          rlinsolMultigrid%ssectionCoarseGridSolver,rlinsolMultigrid%p_rparList)

    end select
    
    ! Level 2..nlevels are smoothers
    do ilevel = 2,rsolutionHierarchy%nlevels

      ! Level 1 is the coarse grid solver
      select case (rlinsolMultigrid%ccoarseGridSolver)
      
      ! -------------------------------
      ! Richardson smoother
      ! -------------------------------
      case (NLIN_SOLVER_RICHARDSON)

      ! Partial initialisation of the corresponding linear solver structure
      rlinsolMultigrid%p_rsubSolvers(ilevel)%csolverType = NLIN_SOLVER_RICHARDSON
        call newtonlin_initBasicParams (&
            rlinsolMultigrid%p_rsubSolvers(ilevel)%rprecParameters,&
            rlinsolMultigrid%ssectionSmoother,rlinsolMultigrid%p_rparList)
            
      end select
    
    end do
    
    ! Allocate the cycle counter.
    allocate(rlinsolMultigrid%p_IcycleCount(rsolutionHierarchy%nlevels))
    
    ! Create temp vectors for the iteration based on the given KKT system hierarchy
    call kkth_initHierarchyDirDeriv (rlinsolMultigrid%rsolutionHier,rsolutionHierarchy)
    call kkth_initHierarchyDirDeriv (rlinsolMultigrid%rdefectHier,rsolutionHierarchy)
    call kkth_initHierarchyDirDeriv (rlinsolMultigrid%rrhsHier,rsolutionHierarchy)
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_initData (rlinsolParam,rsolutionHierarchy)
  
!<description>
  ! Final preparation of the preconditioner for preconditioning.
!</description>

!<inputoutput>
  ! Defines a hierarchy of the solutions of the KKT system.
  ! On the highest level, the solution is expected.
  ! The solution is propagated to all lower levels.
  type(t_kktsystemHierarchy), intent(inout), target :: rsolutionHierarchy

  ! Structure to be initialised.
  type(t_linsolParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>

    integer :: ilevel
    type(t_kktsystem), pointer :: p_rkktSystem,p_rkktSystemCoarse

    ! Initialise subsolvers
    select case (rlinsolParam%csolverType)
    
    ! -----------------------------------------------------
    ! Multigrid inítialisation
    ! -----------------------------------------------------
    case (NLIN_SOLVER_MULTIGRID)

      call newtonlin_initDataMultigrid (&
          rlinsolParam,rlinsolParam%p_rlinsolMultigrid,rsolutionHierarchy)
      
    end select
    
    ! Propagate the solution to the lower levels.
    do ilevel = rsolutionHierarchy%nlevels,2,-1

      call kkth_getKKTsystem (rsolutionHierarchy,ilevel-1,p_rkktSystemCoarse)
      call kkth_getKKTsystem (rsolutionHierarchy,ilevel,p_rkktSystem)
      
      call sptipr_performInterpolation (rlinsolParam%p_RprjHierSpaceTimePrimal,ilevel,&
          p_rkktSystemCoarse%p_rprimalSol%p_rvectorAccess, &
          p_rkktSystem%p_rprimalSol%p_rvectorAccess)

      call sptipr_performInterpolation (rlinsolParam%p_RprjHierSpaceTimeDual,ilevel,&
          p_rkktSystemCoarse%p_rdualSol%p_rvectorAccess, &
          p_rkktSystem%p_rdualSol%p_rvectorAccess)

      call sptipr_performInterpolation (rlinsolParam%p_RprjHierSpaceTimeControl,ilevel,&
          p_rkktSystemCoarse%p_rcontrol%p_rvectorAccess, &
          p_rkktSystem%p_rcontrol%p_rvectorAccess)
    
    end do
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_doneData (rlinsolParam)
  
!<description>
  ! Cleanup of the data initalised in newtonlin_initData.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_linsolParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>
   
    ! Clean up subsolvers
    select case (rlinsolParam%csolverType)
    
    ! -----------------------------------------------------
    ! Multigrid
    ! -----------------------------------------------------
    case (NLIN_SOLVER_MULTIGRID)

      call newtonlin_doneDataMultigrid (rlinsolParam%p_rlinsolMultigrid)
      
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_doneDataMultigrid (rlinsolMultigrid)
  
!<description>
  ! Cleanup of the data initalised in newtonlin_initStructure.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_linsolMultigrid), intent(inout) :: rlinsolMultigrid
!</inputoutput>

!</subroutine>

    ! Release temp data
    call kkth_doneHierarchyDirDeriv (rlinsolMultigrid%rsolutionHier)
    call kkth_doneHierarchyDirDeriv (rlinsolMultigrid%rdefectHier)
    call kkth_doneHierarchyDirDeriv (rlinsolMultigrid%rrhsHier)

    ! Release memory
    deallocate(rlinsolMultigrid%p_IcycleCount)
    deallocate(rlinsolMultigrid%p_rsubSolvers)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_doneStructure (rlinsolParam)
  
!<description>
  ! Cleanup of the data initalised in newtonlin_initStructure.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_linsolParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>

    integer :: i

    ! Release temp vectors
    do i=rlinsolParam%p_rkktsystemHierarchy%nlevels,1,-1
      call kktsp_doneControlVector (rlinsolParam%p_RtempVectors(i))
    end do
    deallocate(rlinsolParam%p_RtempVectors)
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_done (rlinsolParam)
  
!<description>
  ! Clean up a preconditioner.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_linsolParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>

    ! Release memory
    if (associated(rlinsolParam%p_rlinsolMultigrid)) then
      deallocate(rlinsolParam%p_rlinsolMultigrid)
    end if
   
  end subroutine

end module
