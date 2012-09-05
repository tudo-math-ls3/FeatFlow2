
!##############################################################################
!# ****************************************************************************
!# <name> spacesolver </name>
!# ****************************************************************************
!#
!# <purpose>
!# Contains routines for solving nonlinear/linear systems in space.
!#
!# For a forward equation, the space solver applies a nonlinear loop.
!# For a backward equation or a linearised forward/backward equation,
!# the space solver directly delegates the solution task to a linear
!# subsolver.
!# </purpose>
!##############################################################################

module spacesolver

  use fsystem
  use genoutput
  
  use spatialdiscretisation
  use timediscretisation
  use linearalgebra
  use linearsystemscalar
  use linearsystemblock
  use vectorfilters
  use iterationcontrol
  
  use scalarpde
  use linearformevaluation
  use bilinearformevaluation
  use feevaluation2
  use blockmatassemblybase
  use blockmatassembly
  use collection
  
  use spacetimevectors
  use analyticsolution
  use spacetimehierarchy
  
  use structuresdiscretisation
  use structuresoptcontrol
  use structuresgeneral
  use structuresoptflow
  use structuresboundaryconditions
  use assemblytemplates
  
  use spatialbc
  
  use structuresoperatorasm
  use spacematvecassembly
  
  use kktsystemspaces
  
  use spacelinearsolver
  use structuresnewton

  implicit none
  
  private

!<types>

!<typeblock>

  ! Linear solver statistics
  type t_spaceslSolverStat
  
    ! Number of iterations necessary for the solver
    integer :: niterations = 0
    
    ! Total time necessary for the linear solver
    type(t_timer) :: rtotalTime

    ! Time necessary for the creation of nonlinear defects
    type(t_timer) :: rtimeDefect

    ! Time necessary for the creation of RHS vectors in space
    type(t_timer) :: rtimeRhs

    ! Time necessary for the creation of matrices
    type(t_timer) :: rtimeMatrixAssembly

    ! Statistics of the linear subsolver
    type(t_lssSolverStat) :: rlssSolverStat
    
  end type

!</typeblock>

  public :: t_spaceslSolverStat

!<typeblock>

  ! The space solver hierarchy structure is a solver which is able
  ! to solve in every timestep of a space-time hierarchy the corresponding
  ! nonlinear forward system or the corresponding linear forward/backward
  ! systems.
  type t_spaceSolverHierarchy
  
    ! <!-- --------------------------- -->
    ! <!-- NONLINEAR SOLVER PARAMETERS -->
    ! <!-- --------------------------- -->

    ! General Newton parameters for nonlinear iterations.
    ! This is level independent.
    type(t_newtonParameters) :: rnewtonParams
    
    ! Iteration control parameters
    type(t_iterationControl) :: riter
    
    ! <!-- ------------------------ -->
    ! <!-- LINEAR SOLVER STRUCTURES -->
    ! <!-- ------------------------ -->

    ! Hierarchy of linear solvers in space used for solving
    ! auxiliary linear subproblems. The solver on level ilevel
    ! will be used to solve linear subproblems.
    type(t_linsolHierarchySpace), pointer :: p_rlssHierarchy => null()

    ! Hierarchy of linear solvers in space used for solving
    ! auxiliary linear subproblems. The solver on level ilevel
    ! will be used to solve linear subproblems.
    ! These sovlers are used as fallback if the standard solver above fails.
    type(t_linsolHierarchySpace), pointer :: p_rlssHierarchy2 => null()

    ! <!-- ------------------------------ -->
    ! <!-- BOUNDARY CONDITIONS STRUCTURES -->
    ! <!-- ------------------------------ -->

    ! Boundary condition hierarchy for all space levels, primal and dual space
    type(t_optcBDCSpaceHierarchy), pointer :: p_roptcBDCSpaceHierarchy => null()

    ! <!-- ---------------- -->
    ! <!-- HIERARCHIES      -->
    ! <!-- ---------------- -->
    
    ! A hierarchy of operator assembly structures for all levels.
    ! This is a central discretisation structure passed to all assembly routines.
    type(t_spacetimeOpAsmHierarchy), pointer :: p_roperatorAsmHier => null()
    
    ! <!-- ---------------- -->
    ! <!-- GENERAL SETTINGS -->
    ! <!-- ---------------- -->

    ! Type of equation to be solved here. This can be 
    ! OPTP_PRIMAL, OPTP_DUAL, OPTP_PRIMALLIN or OPTP_DUALLIN,
    ! depending on which equation to solve.
    integer :: copType = 0

    ! <!-- -------------- -->
    ! <!-- TEMPORARY DATA -->
    ! <!-- -------------- -->
    
    ! Current space level, this solver structure is configured for.
    integer :: ispacelevel = 0
    
    ! Current time level, this solver structure is configured for.
    integer :: itimelevel = 0
    
    ! Temporary vectors for nonlinear iterations
    type(t_vectorBlock), pointer :: p_rd => null()
    type(t_vectorBlock), pointer :: p_rd2 => null()

    ! Hierarchy of system matrices
    type(t_matrixBlock), dimension(:), pointer :: p_Rmatrices => null()
    
    ! Hierarchy of assembly flags.
    type(t_assemblyFlags), dimension(:), pointer :: p_RassemblyFlags => null()
    
    ! Temporary assembly data
    type(t_assemblyTempDataSpace) :: rtempData
    
  end type

!</typeblock>

  public :: t_spaceSolverHierarchy

!</types>

!<constants>

!<constantblock description = "Specification constants for vectors in a system in space.">

  ! Solution vector
  integer, parameter :: SPACESLH_VEC_SOLUTION = 0
  
  ! RHS vector
  integer, parameter :: SPACESLH_VEC_RHS      = 1
  
  ! Defect vector
  integer, parameter :: SPACESLH_VEC_DEFECT   = 2

!</constantblock>

!<constantblock description = "Specifies some flags that modify the underlying equation to solve.">

  ! Assemble / solve the full system
  integer, parameter, public :: SPACESLH_EQNF_DEFAULT = 0
  
  ! Simplify the system: No Newton term in the linearised dual equation
  ! This gives the same result as solving the system with the "simple" solver
  ! OPTP_PRIMALLIN_SIMPLE.
  integer, parameter, public :: SPACESLH_EQNF_NONEWTONPRIMAL = 1

  ! Simplify the system: No Newton term in the linearised dual equation
  ! This gives the same result as solving the system with the "simple" solver
  ! OPTP_DUALLIN_SIMPLE.
  integer, parameter, public :: SPACESLH_EQNF_NONEWTONDUAL = 2
  
  ! Simplify the system: No Newton term in the linearised primal and dual equation.
  ! This gives the same result as solving the system with the "simple" solver
  ! OPTP_PRIMALLIN_SIMPLE / OPTP_DUALLIN_SIMPLE.
  integer, parameter, public :: SPACESLH_EQNF_NONEWTON = 3

!</constantblock>

!</constants>

  ! Basic initialisation of a space solver
  public :: spaceslh_init

  ! Structural initialisation
  public :: spaceslh_initStructure
  
  ! Solves the spatial system
  public :: spaceslh_solve
  
  ! Cleanup of structures
  public :: spaceslh_doneStructure

  ! Final cleanup
  public :: spaceslh_done
  
  ! Set the stopping criterion for the adaptive Newton algorithm
  public :: spaceslh_setEps

  ! Clears a statistics block
  public :: spacesl_clearStatistics
  
  ! Sums up two statistic blocks
  public :: spacesl_sumStatistics
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_init (rsolver,copType,rlssHierarchy,rlssHierarchy2,ssection,rparamList)
  
!<description>
  ! Initialises the solver parameters according to a parameter list.
  ! For each level in rlssHierarchy, a space solver will be created.
!</description>
  
!<input>
  ! Type of equation to be solved here. This can be 
  ! OPTP_PRIMAL, OPTP_DUAL, OPTP_PRIMALLIN or OPTP_DUALLIN,
  ! depending on which equation to solve.
  integer, intent(in) :: copType

  ! Hierarchy of linear solvers in space used for solving
  ! auxiliary linear subproblems. The solver on level ispacelevel
  ! will be used to solve linear subproblems.
  type(t_linsolHierarchySpace), intent(in), target :: rlssHierarchy

  ! Hierarchy of linear solvers in space used for solving
  ! auxiliary linear subproblems. The solver on level ispacelevel
  ! will be used to solve linear subproblems.
  ! THese sovlers are used as fallback if the standard solver fails.
  type(t_linsolHierarchySpace), intent(in), target :: rlssHierarchy2

  ! OPTIONAL: Name of the section in the parameter list containing the parameters
  ! of the nonlinear solver.
  ! Can be omitted for linear solvers in space.
  character(LEN=*), intent(in), optional :: ssection

  ! OPTIONAL: Parameter list with the parameters configuring the nonlinear solver
  ! Can be omitted for linear solvers in space.
  type(t_parlist), intent(in), optional :: rparamList
!</input>

!<output>
  ! Solver structure receiving the parameters
  type(t_spaceSolverHierarchy), intent(inout) :: rsolver
!</output>

!</subroutine>

    ! Remember the solver settings for later use
    rsolver%p_rlssHierarchy => rlssHierarchy
    rsolver%p_rlssHierarchy2 => rlssHierarchy2
    rsolver%copType = copType

    if (present(ssection) .and. present(rparamList)) then
    
      ! Initialise basic parameters of the nonlinear solver
      call newtonit_initBasicParams (rsolver%rnewtonParams,rsolver%riter,ssection,rparamList)
    
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_initStructure (rsolver, ispacelevel, itimelevel,&
      roperatorAsmHier,rstatistics,ierror)
  
!<description>
  ! Structural initialisation of the Newton solver.
!</description>

!<input>
  ! Space-level in the global space-time hierarchy, the solver should be applied to.
  integer, intent(in) :: ispacelevel

  ! Time-level in the global space-time hierarchy, the solver should be applied to.
  integer, intent(in) :: itimelevel

  ! A hierarchy of operator assembly structures for all levels.
  ! This is a central discretisation structure passed to all assembly routines.
  type(t_spacetimeOpAsmHierarchy), intent(in), target :: roperatorAsmHier
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_spaceSolverHierarchy), intent(inout) :: rsolver
!</inputoutput>

!<output>
  ! Statistics structure
  type(t_spaceslSolverStat), intent(out) :: rstatistics

  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!</subroutine>

    integer :: ilev
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    type(t_spacetimeOperatorAsm) :: roperatorAsm
    type(t_lssSolverStat) :: rlocalStatLss
    
    if (rsolver%p_rlssHierarchy%nlmin .gt. ispaceLevel) then
      call output_line ("Invalid space level.", &
          OU_CLASS_ERROR,OU_MODE_STD,"spaceslh_initStructure")
      call sys_halt()
    end if
    
    ! Remember the level and the hierarchy structures. 
    rsolver%ispacelevel = ispacelevel
    rsolver%itimelevel = itimelevel
    rsolver%p_roperatorAsmHier => roperatorAsmHier

    ! Prepare boundary condition structures
    allocate(rsolver%p_roptcBDCSpaceHierarchy)
    
    call sbch_initBDCHierarchy (rsolver%p_roptcBDCSpaceHierarchy,&
        rsolver%p_rlssHierarchy%nlmin,rsolver%ispaceLevel)
    
    ! Allocate temporary vectors / matrices for all space levels.
    allocate(rsolver%p_Rmatrices(rsolver%ispaceLevel))
    allocate(rsolver%p_RassemblyFlags(rsolver%ispaceLevel))
    
    do ilev = 1,rsolver%ispaceLevel
    
      ! Initialise a temporary vector for the nonlinearity
      select case (rsolver%copType)
      case (OPTP_PRIMAL,OPTP_PRIMALLIN,OPTP_PRIMALLIN_SIMPLE)
        p_rdiscretisation => &
            roperatorAsmHier%p_rfeHierarchyPrimal%p_rfeSpaces(ilev)%p_rdiscretisation

      case (OPTP_DUAL,OPTP_DUALLIN,OPTP_DUALLIN_SIMPLE)
        p_rdiscretisation => &
            roperatorAsmHier%p_rfeHierarchyDual%p_rfeSpaces(ilev)%p_rdiscretisation
      
      case default
        call output_line ("Unknown space.", &
            OU_CLASS_ERROR,OU_MODE_STD,"spaceslh_initStructure")
        call sys_halt()
      end select
      
      ! Take a look at the solver type to configure the assembly flags.
      select case (rsolver%p_rlssHierarchy%p_RlinearSolvers(ilev)%isolverType)
      
      ! -------------------------------
      ! UMFPACK
      ! -------------------------------
      case (LSS_LINSOL_UMFPACK)
      
        rsolver%p_RassemblyFlags(ilev)%bumfpackSolver = .true.

      ! -------------------------------
      ! MULTIGRID
      ! -------------------------------
      case (LSS_LINSOL_MG)
      
        ! Check if UMFPACK is the coarse grid solver.
        if ((ilev .eq. 1) .and. &
            (rsolver%p_rlssHierarchy%p_RlinearSolvers(ilev)%icoarseGridSolverType .eq. 0)) then
          rsolver%p_RassemblyFlags(ilev)%bumfpackSolver = .true.
        end if
      
      end select

      ! Take a look at the solver type to configure the assembly flags.
      select case (rsolver%p_rlssHierarchy2%p_RlinearSolvers(ilev)%isolverType)
      
      ! -------------------------------
      ! UMFPACK
      ! -------------------------------
      case (LSS_LINSOL_UMFPACK)
      
        rsolver%p_RassemblyFlags(ilev)%bumfpackSolver = .true.

      ! -------------------------------
      ! MULTIGRID
      ! -------------------------------
      case (LSS_LINSOL_MG)
      
        ! Check if UMFPACK is the coarse grid solver.
        if ((ilev .eq. 1) .and. &
            (rsolver%p_rlssHierarchy2%p_RlinearSolvers(ilev)%icoarseGridSolverType .eq. 0)) then
          rsolver%p_RassemblyFlags(ilev)%bumfpackSolver = .true.
        end if
      
      end select
      
      ! Get the corresponding operator assembly structure and
      ! initialise a system matrix
      call stoh_getOpAsm_slvtlv (&
          roperatorAsm,roperatorAsmHier,ilev,rsolver%itimelevel)
          
      call smva_allocSystemMatrix (rsolver%p_Rmatrices(ilev),&
          roperatorAsmHier%ranalyticData%p_rphysics,&
          roperatorAsmHier%ranalyticData%p_rsettingsOptControl,&
          roperatorAsmHier%ranalyticData%p_rsettingsSpaceDiscr,&
          p_rdiscretisation,roperatorAsm%p_rasmTemplates,&
          rsolver%p_RassemblyFlags(ilev))
          
      ! Provide the matrix to the linear solver
      call lssh_setMatrix(rsolver%p_rlsshierarchy,ilev,rsolver%p_Rmatrices(ilev))
      call lssh_setMatrix(rsolver%p_rlsshierarchy2,ilev,rsolver%p_Rmatrices(ilev))
      
      if (ilev .eq. rsolver%ispaceLevel) then
      
        ! On the highest space level, allocate some temp vectors on the topmost level

        allocate (rsolver%p_rd)
        call lsysbl_createVectorBlock (p_rdiscretisation,rsolver%p_rd)

        allocate (rsolver%p_rd2)
        call lsysbl_createVectorBlock (p_rdiscretisation,rsolver%p_rd2)
        
        select case (rsolver%copType)
        case (OPTP_PRIMAL,OPTP_PRIMALLIN,OPTP_PRIMALLIN_SIMPLE)
          call smva_allocTempData (rsolver%rtempData,&
              roperatorAsmHier%ranalyticData%p_rphysics,&
              roperatorAsmHier%ranalyticData%p_rsettingsOptControl,&
              roperatorAsmHier%ranalyticData%p_rsettingsSpaceDiscr,&
              ilev,roperatorAsmHier,rsolver%p_RassemblyFlags(ilev))

        case (OPTP_DUAL,OPTP_DUALLIN,OPTP_DUALLIN_SIMPLE)
          call smva_allocTempData (rsolver%rtempData,&
              roperatorAsmHier%ranalyticData%p_rphysics,&
              roperatorAsmHier%ranalyticData%p_rsettingsOptControl,&
              roperatorAsmHier%ranalyticData%p_rsettingsSpaceDiscr,&
              ilev,roperatorAsmHier,rsolver%p_RassemblyFlags(ilev))
        case default
          call output_line ("Unknown space.", &
              OU_CLASS_ERROR,OU_MODE_STD,"spaceslh_initStructure")
          call sys_halt()
        end select

        ! Initialise the structures of the associated linear subsolver
        call lssh_initStructure (rsolver%p_rlsshierarchy,ilev,rlocalStatLss,ierror)
        call lssh_initStructure (rsolver%p_rlsshierarchy2,ilev,rlocalStatLss,ierror)
        
        call lss_sumStatistics(rlocalStatLss,rstatistics%rlssSolverStat)

      end if
    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_initData_bdc (rsolver, idofTime, rcontrol)
  
!<description>
  ! Initialises the boundary condition structures for time DOF idofTime
!</description>

!<input>
  ! Number of the DOF in time.
  integer, intent(in) :: idofTime

  ! OPTIONAL: Structure that defines the current control.
  type(t_controlSpace), intent(inout), optional :: rcontrol
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_spaceSolverHierarchy), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ilev
    type(t_vectorBlock), pointer :: p_rcontrolVec
    
    ! Clean up the boundary conditions.
    call sbch_resetBCstructure (rsolver%p_roptcBDCSpaceHierarchy)
    
    ! Loop over all levels, from the highest to the lowest,
    ! and set up the boundary conditions.
    do ilev = rsolver%ispacelevel,rsolver%p_rlssHierarchy%nlmin,-1
    
      if (.not. present(rcontrol)) then
      
        ! Assemble Dirichlet/Neumann boundary conditions
        call smva_initDirichletNeumannBC (&
            rsolver%p_roptcBDCSpaceHierarchy%p_RoptcBDCspace(ilev),&
            rsolver%p_roperatorAsmHier%ranalyticData%p_roptcBDC,&
            rsolver%copType,&
            rsolver%p_roperatorAsmHier,ilev,rsolver%itimelevel,idoftime,&
            rsolver%p_roperatorAsmHier%ranalyticData%p_rglobalData)
      
      else
      
        call sptivec_getVectorFromPool (rcontrol%p_rvectorAccess,idofTime,p_rcontrolVec)

        ! Assemble Dirichlet/Neumann boundary conditions
        call smva_initDirichletNeumannBC (&
            rsolver%p_roptcBDCSpaceHierarchy%p_RoptcBDCspace(ilev),&
            rsolver%p_roperatorAsmHier%ranalyticData%p_roptcBDC,&
            rsolver%copType,&
            rsolver%p_roperatorAsmHier,ilev,rsolver%itimelevel,idoftime,&
            rsolver%p_roperatorAsmHier%ranalyticData%p_rglobalData,p_rcontrolVec)
      end if
    
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_implementBDC (rvector, cvectype, ilevel, rsolver)
  
!<description>
  ! Implementa boundary conditions into a vector rvector
!</description>

!<input>
  ! Type of the vector. One of the SPACESLH_VEC_xxxx constants
  integer, intent(in) :: cvectype
  
  ! Level corresponding to rvector
  integer, intent(in) :: ilevel

  ! Solver structure
  type(t_spaceSolverHierarchy), intent(in) :: rsolver
!</input>

!<inputoutput>
  ! Vector that receives the boundary conditions.
  type(t_vectorBlock), intent(inout) :: rvector
!</inputoutput>

!</subroutine>
    
    ! local variables
    type(t_optcBDCSpace), pointer :: p_roptcBDCspace

    ! Get the boundary conditions structure
    p_roptcBDCspace => rsolver%p_roptcBDCSpaceHierarchy%p_RoptcBDCspace(ilevel)

    ! Type of vector?
    select case (cvectype)
    
    ! ---------------------------------
    ! Solution vector
    ! ---------------------------------
    case (SPACESLH_VEC_SOLUTION)
    
      ! Dirichlet/Neumann BC
      call vecfil_discreteBCsol (rvector,p_roptcBDCspace%rdiscreteBC)
      
      ! Which equation do we have?    
      select case (rsolver%p_roperatorAsmHier%ranalyticData%p_rphysics%cequation)

      ! *************************************************************
      ! Stokes/Navier Stokes.
      ! *************************************************************
      case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)

        ! Pure Dirichlet problem?

        ! -----------------------------------------------
        ! Integral-mean-value-zero filter for pressure
        ! -----------------------------------------------
        if (p_roptcBDCSpace%rneumannBoundary%nregions .eq. 0) then
          call vecfil_subvectorToL20 (rvector,3)
        end if

      ! *************************************************************
      ! Heat equation
      ! *************************************************************
      case (CCEQ_HEAT2D)

        ! Pure Neumann problem?

        ! -----------------------------------------------
        ! Integral-mean-value-zero filter for the vector
        ! -----------------------------------------------
        if (p_roptcBDCSpace%rdirichletBoundary%nregions .eq. 0) then
          call vecfil_subvectorToL20 (rvector,3)
        end if

      end select      

    ! ---------------------------------
    ! RHS vector
    ! ---------------------------------
    case (SPACESLH_VEC_RHS)
      ! Dirichlet/Neumann BC
      call vecfil_discreteBCrhs (rvector,p_roptcBDCspace%rdiscreteBC)

      ! Which equation do we have?    
      select case (rsolver%p_roperatorAsmHier%ranalyticData%p_rphysics%cequation)

      ! *************************************************************
      ! Stokes/Navier Stokes.
      ! *************************************************************
      case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)

        ! Pure Dirichlet problem?

        ! -----------------------------------------------
        ! Integral-mean-value-zero filter for pressure
        ! -----------------------------------------------
        if (p_roptcBDCSpace%rneumannBoundary%nregions .eq. 0) then
          call vecfil_subvectorToL20 (rvector,3)
        end if

      ! *************************************************************
      ! Heat equation
      ! *************************************************************
      case (CCEQ_HEAT2D)

        ! Pure Neumann problem?

        ! -----------------------------------------------
        ! Integral-mean-value-zero filter for the vector
        ! -----------------------------------------------
        if (p_roptcBDCSpace%rdirichletBoundary%nregions .eq. 0) then
          call vecfil_subvectorToL20 (rvector,3)
        end if

      end select      

    ! ---------------------------------
    ! Defect vector
    ! ---------------------------------
    case (SPACESLH_VEC_DEFECT)
      ! Dirichlet/Neumann BC
      call vecfil_discreteBCdef (rvector,p_roptcBDCspace%rdiscreteBC)
    
      ! Which equation do we have?    
      select case (rsolver%p_roperatorAsmHier%ranalyticData%p_rphysics%cequation)

      ! *************************************************************
      ! Stokes/Navier Stokes.
      ! *************************************************************
      case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)

        ! Pure Dirichlet problem?

        ! -----------------------------------------------
        ! Integral-mean-value-zero filter for pressure
        ! -----------------------------------------------
        if (p_roptcBDCSpace%rneumannBoundary%nregions .eq. 0) then
          call vecfil_subvectorToL20 (rvector,3)
        end if

      ! *************************************************************
      ! Heat equation
      ! *************************************************************
      case (CCEQ_HEAT2D)

        ! Pure Neumann problem?

        ! -----------------------------------------------
        ! Integral-mean-value-zero filter for the vector
        ! -----------------------------------------------
        if (p_roptcBDCSpace%rdirichletBoundary%nregions .eq. 0) then
          call vecfil_subvectorToL20 (rvector,3)
        end if

      end select      

    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_initData_primal (rsolver, ierror, idofTime, &
      rprimalSol, isollevelSpace, rstatistics, coptype, bfallbackSolver)
  
!<description>
  ! Final preparation of the Newton solver. Primal space.
!</description>

!<input>
  ! Number of the DOF in time which should be calculated into rdest.
  integer, intent(in) :: idofTime

  ! Current solution of the primal equation. Defines the nonlinearity.
  ! The space level of the solution is specified by isollevelSpace.
  ! The time level must match rsolver%itimelevel.
  type(t_primalSpace), intent(inout), target :: rprimalSol
  
  ! Space level corresponding to rprimalSol.
  integer, intent(in) :: isollevelSpace
  
  ! Operator type. Either OPTP_PRIMALLIN or OPTP_PRIMALLIN_SIMPLE.
  integer, intent(in) :: coptype

  ! Whether or not to initialise the fallback solver.
  ! =FALSE: Create matrices and initialise the default solver.
  ! =TRUE: Use the matrices calculated for FALSE and initialise the 
  ! fallback solver. THe matrices must have been initialised using FALSE.
  logical, intent(in) :: bfallbackSolver
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_spaceSolverHierarchy), intent(inout) :: rsolver
!</inputoutput>

!<output>
  ! Statistics structure
  type(t_spaceslSolverStat), intent(out) :: rstatistics

  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!</subroutine>

    ! local variables
    integer :: ilev,iprevlv
    type(t_lssSolverStat) :: rlocalStatLss
    
    ! Calculate matrices and initialise the default solver?
    if (.not. bfallbackSolver) then
    
      ! iprevlv contains the last assembled level.
      iprevlv = 0

      ! Loop over all levels, from the highest to the lowest,
      ! and set up the system matrices.
      call stat_startTimer (rstatistics%rtimeMatrixAssembly)
      
      do ilev = rsolver%ispacelevel,rsolver%p_rlssHierarchy%nlmin,-1
      
        ! Assemble the matrix
        call smva_assembleMatrix_primal (rsolver%p_Rmatrices(ilev),ilev,idofTime,&
            rsolver%p_roperatorAsmHier,&
            rprimalSol,isollevelSpace,rsolver%itimelevel,coptype,&
            rsolver%p_roptcBDCSpaceHierarchy%p_RoptcBDCspace(ilev),&
            rsolver%p_RassemblyFlags(ilev),rsolver%rtempData,iprevlv)
            
        iprevlv = ilev
      
      end do
      
      call stat_stopTimer (rstatistics%rtimeMatrixAssembly)

      ! Initialise the structures of the associated linear subsolver
      rsolver%p_rlsshierarchy%p_rdebugFlags%sstringTag = &
          "primal_"//trim(sys_siL(idofTime,10))//"_"//trim(sys_siL(isollevelSpace,10))

      call lssh_initData (&
          rsolver%p_rlsshierarchy,rsolver%p_roptcBDCSpaceHierarchy,&
          rsolver%ispacelevel,rlocalStatLss,ierror)
      
      call lss_sumStatistics(rlocalStatLss,rstatistics%rlssSolverStat)

    else
    
      ! Matrices have already been initialised.
      ! Initialise the fallback solver.
      rsolver%p_rlsshierarchy%p_rdebugFlags%sstringTag = &
          "primal2_"//trim(sys_siL(idofTime,10))//"_"//trim(sys_siL(isollevelSpace,10))

      call lssh_initData (&
          rsolver%p_rlsshierarchy2,rsolver%p_roptcBDCSpaceHierarchy,&
          rsolver%ispacelevel,rlocalStatLss,ierror)
      
      call lss_sumStatistics(rlocalStatLss,rstatistics%rlssSolverStat)
    
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_initData_primallin (rsolver, ierror, idofTime, &
      rprimalSol, isollevelSpace, coptype, rstatistics, bfallbackSolver)
  
!<description>
  ! Final preparation of the Newton solver. Linearised primal space.
!</description>

!<input>
  ! Number of the DOF in time which should be calculated into rdest.
  integer, intent(in) :: idofTime

  ! Current solution of the primal equation. Defines the nonlinearity.
  ! The space level of the solution is specified by isollevelSpace.
  ! The time level must match rsolver%itimelevel.
  type(t_primalSpace), intent(inout), target :: rprimalSol
  
  ! Space level corresponding to rprimalSol.
  integer, intent(in) :: isollevelSpace
  
  ! Operator type. Either OPTP_PRIMALLIN or OPTP_PRIMALLIN_SIMPLE.
  integer, intent(in) :: coptype

  ! Whether or not to initialise the fallback solver.
  ! =FALSE: Create matrices and initialise the default solver.
  ! =TRUE: Use the matrices calculated for FALSE and initialise the 
  ! fallback solver. THe matrices must have been initialised using FALSE.
  logical, intent(in) :: bfallbackSolver
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_spaceSolverHierarchy), intent(inout) :: rsolver
!</inputoutput>

!<output>
  ! Statistics structure
  type(t_spaceslSolverStat), intent(out) :: rstatistics

  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!</subroutine>

    ! local variables
    integer :: ilev,iprevlv
    type(t_lssSolverStat) :: rlocalStatLss
    
    ! Calculate matrices and initialise the default solver?
    if (.not. bfallbackSolver) then
    
      ! iprevlv contains the last assembled level.
      iprevlv = 0

      ! Loop over all levels, from the highest to the lowest,
      ! and set up the system matrices.
      call stat_startTimer (rstatistics%rtimeMatrixAssembly)
      
      do ilev = rsolver%ispacelevel,rsolver%p_rlssHierarchy%nlmin,-1
      
        ! Assemble the matrix
        call smva_assembleMatrix_primal (rsolver%p_Rmatrices(ilev),ilev,idofTime,&
            rsolver%p_roperatorAsmHier,&
            rprimalSol,isollevelSpace,rsolver%itimelevel,coptype,&
            rsolver%p_roptcBDCSpaceHierarchy%p_RoptcBDCspace(ilev),&
            rsolver%p_RassemblyFlags(ilev),rsolver%rtempData,iprevlv)
            
        iprevlv = ilev
      
      end do
      
      call stat_stopTimer (rstatistics%rtimeMatrixAssembly)

      ! Initialise the structures of the associated linear subsolver
      rsolver%p_rlsshierarchy%p_rdebugFlags%sstringTag = &
          "primallin_"//trim(sys_siL(idofTime,10))//"_"//trim(sys_siL(isollevelSpace,10))
      
      call lssh_initData (&
          rsolver%p_rlsshierarchy,rsolver%p_roptcBDCSpaceHierarchy,&
          rsolver%ispacelevel,rlocalStatLss,ierror)
      
      call lss_sumStatistics(rlocalStatLss,rstatistics%rlssSolverStat)

    else
    
      ! Matrices have already been initialised.
      ! Initialise the fallback solver.

      rsolver%p_rlsshierarchy%p_rdebugFlags%sstringTag = &
          "primallin2_"//trim(sys_siL(idofTime,10))//"_"//trim(sys_siL(isollevelSpace,10))
      
      call lssh_initData (&
          rsolver%p_rlsshierarchy2,rsolver%p_roptcBDCSpaceHierarchy,&
          rsolver%ispacelevel,rlocalStatLss,ierror)
      
      call lss_sumStatistics(rlocalStatLss,rstatistics%rlssSolverStat)

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_initData_dual (rsolver, ierror, idofTime, &
      rprimalSol, isollevelSpace, rstatistics, bfallbackSolver)
  
!<description>
  ! Final preparation of the Newton solver. Dual space.
!</description>

!<input>
  ! Number of the DOF in time which should be calculated into rdest.
  integer, intent(in) :: idofTime

  ! Current solution of the primal equation. Defines the nonlinearity.
  ! The space level of the solution is specified by isollevelSpace.
  ! The time level must match rsolver%itimelevel.
  type(t_primalSpace), intent(inout), target :: rprimalSol
  
  ! Space level corresponding to rprimalSol.
  integer, intent(in) :: isollevelSpace

  ! Whether or not to initialise the fallback solver.
  ! =FALSE: Create matrices and initialise the default solver.
  ! =TRUE: Use the matrices calculated for FALSE and initialise the 
  ! fallback solver. THe matrices must have been initialised using FALSE.
  logical, intent(in) :: bfallbackSolver
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_spaceSolverHierarchy), intent(inout) :: rsolver
!</inputoutput>

!<output>
  ! Statistics structure
  type(t_spaceslSolverStat), intent(inout) :: rstatistics

  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!</subroutine>

    ! local variables
    integer :: ilev,iprevlv
    type(t_lssSolverStat) :: rlocalStatLss
    
    ! Calculate matrices and initialise the default solver?
    if (.not. bfallbackSolver) then

      ! iprevlv contains the last assembled level.
      iprevlv = 0

      ! Loop over all levels, from the highest to the lowest,
      ! and set up the system matrices.
      call stat_startTimer (rstatistics%rtimeMatrixAssembly)
      
      do ilev = rsolver%ispacelevel,rsolver%p_rlssHierarchy%nlmin,-1
      
        ! Assemble the matrix
        call smva_assembleMatrix_dual (rsolver%p_Rmatrices(ilev),ilev,idofTime,&
            rsolver%p_roperatorAsmHier,rprimalSol,isollevelSpace,rsolver%itimelevel,&
            rsolver%p_roptcBDCSpaceHierarchy%p_RoptcBDCspace(ilev),&
            rsolver%p_RassemblyFlags(ilev),rsolver%rtempData,iprevlv)
            
        iprevlv = ilev
      
      end do
      
      call stat_stopTimer (rstatistics%rtimeMatrixAssembly)

      ! Initialise the structures of the associated linear subsolver
      rsolver%p_rlsshierarchy%p_rdebugFlags%sstringTag = &
          "dual_"//trim(sys_siL(idofTime,10))//"_"//trim(sys_siL(isollevelSpace,10))
      
      call lssh_initData (&
          rsolver%p_rlsshierarchy,rsolver%p_roptcBDCSpaceHierarchy,&
          rsolver%ispacelevel,rlocalStatLss,ierror)
      
      call lss_sumStatistics(rlocalStatLss,rstatistics%rlssSolverStat)

    else
    
      ! Matrices have already been initialised.
      ! Initialise the fallback solver.
      rsolver%p_rlsshierarchy%p_rdebugFlags%sstringTag = &
          "dual2_"//trim(sys_siL(idofTime,10))//"_"//trim(sys_siL(isollevelSpace,10))
      
      call lssh_initData (&
          rsolver%p_rlsshierarchy2,rsolver%p_roptcBDCSpaceHierarchy,&
          rsolver%ispacelevel,rlocalStatLss,ierror)
      
      call lss_sumStatistics(rlocalStatLss,rstatistics%rlssSolverStat)

    end if
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_initData_dualLin (rsolver, ierror, idofTime, &
      rprimalSol, isollevelSpace, rstatistics, bfallbackSolver)
  
!<description>
  ! Final preparation of the Newton solver. Linearised dual space.
!</description>

!<input>
  ! Number of the DOF in time which should be calculated into rdest.
  integer, intent(in) :: idofTime

  ! Current solution of the primal equation. Defines the nonlinearity.
  ! The space level of the solution is specified by isollevelSpace.
  ! The time level must match rsolver%itimelevel.
  type(t_primalSpace), intent(inout), target :: rprimalSol
  
  ! Space level corresponding to rprimalSol.
  integer, intent(in) :: isollevelSpace

  ! Whether or not to initialise the fallback solver.
  ! =FALSE: Create matrices and initialise the default solver.
  ! =TRUE: Use the matrices calculated for FALSE and initialise the 
  ! fallback solver. THe matrices must have been initialised using FALSE.
  logical, intent(in) :: bfallbackSolver
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_spaceSolverHierarchy), intent(inout) :: rsolver
!</inputoutput>

!<output>
  ! Statistics structure
  type(t_spaceslSolverStat), intent(out) :: rstatistics

  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!</subroutine>

    ! Implementation identical to spaceslh_initData_dual.

    ! local variables
    integer :: ilev,iprevlv
    type(t_lssSolverStat) :: rlocalStatLss
    
    ! Calculate matrices and initialise the default solver?
    if (.not. bfallbackSolver) then

      ! iprevlv contains the last assembled level.
      iprevlv = 0

      ! Loop over all levels, from the highest to the lowest,
      ! and set up the system matrices.
      call stat_startTimer (rstatistics%rtimeMatrixAssembly)
      
      do ilev = rsolver%ispacelevel,rsolver%p_rlssHierarchy%nlmin,-1
      
        ! Assemble the matrix
        call smva_assembleMatrix_dual (rsolver%p_Rmatrices(ilev),ilev,idofTime,&
            rsolver%p_roperatorAsmHier,rprimalSol,isollevelSpace,rsolver%itimelevel,&
            rsolver%p_roptcBDCSpaceHierarchy%p_RoptcBDCspace(ilev),&
            rsolver%p_RassemblyFlags(ilev),rsolver%rtempData,iprevlv)
            
        iprevlv = ilev
      
      end do
      
      call stat_stopTimer (rstatistics%rtimeMatrixAssembly)

      ! Initialise the structures of the associated linear subsolver
      rsolver%p_rlsshierarchy%p_rdebugFlags%sstringTag = &
          "duallin_"//trim(sys_siL(idofTime,10))//"_"//trim(sys_siL(isollevelSpace,10))

      call lssh_initData (&
          rsolver%p_rlsshierarchy,rsolver%p_roptcBDCSpaceHierarchy,&
          rsolver%ispacelevel,rlocalStatLss,ierror)
      
      call lss_sumStatistics(rlocalStatLss,rstatistics%rlssSolverStat)

    else
    
      ! Matrices have already been initialised.
      ! Initialise the fallback solver.
      rsolver%p_rlsshierarchy%p_rdebugFlags%sstringTag = &
          "duallin2_"//trim(sys_siL(idofTime,10))//"_"//trim(sys_siL(isollevelSpace,10))

      call lssh_initData (&
          rsolver%p_rlsshierarchy2,rsolver%p_roptcBDCSpaceHierarchy,&
          rsolver%ispacelevel,rlocalStatLss,ierror)
      
      call lss_sumStatistics(rlocalStatLss,rstatistics%rlssSolverStat)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_doneData (rsolver,bfallbackSolver)
  
!<description>
  ! Cleanup of the data initalised in spaceslh_initData.
!</description>

!<input>
  ! Whether or not to clean up the fallback solver.
  logical, intent(in) :: bfallbackSolver
!</input>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_spaceSolverHierarchy), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>

    ! Clean up the structures of the associated linear subsolver
    if (.not. bfallbackSolver) then
      call lssh_doneData (rsolver%p_rlsshierarchy,rsolver%ispacelevel)
    else
      call lssh_doneData (rsolver%p_rlsshierarchy2,rsolver%ispacelevel)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_doneStructure (rsolver)
  
!<description>
  ! Cleanup of the data initalised in spaceslh_initStructure.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_spaceSolverHierarchy), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>

    integer :: ilev

    ! Clean up boundary conditions   
    call sbch_doneBDCHierarchy (rsolver%p_roptcBDCSpaceHierarchy)
    deallocate(rsolver%p_roptcBDCSpaceHierarchy)

    ! Clean up the structures of the associated linear subsolver
    call lssh_doneStructure (rsolver%p_rlsshierarchy2,rsolver%ispacelevel)
    call lssh_doneStructure (rsolver%p_rlsshierarchy,rsolver%ispacelevel)
   
   ! Release assembly data
    call smva_releaseTempData (rsolver%rtempData)
   
    ! deallocate temporary vectors
    do ilev=rsolver%ispacelevel,rsolver%p_rlssHierarchy%nlmin,-1
      call lsysbl_releaseMatrix (rsolver%p_Rmatrices(ilev))
    end do

    deallocate(rsolver%p_Rmatrices)
    deallocate(rsolver%p_RassemblyFlags)
    
    ! Deallocate temporary data
    call lsysbl_releaseVector (rsolver%p_rd2)
    deallocate (rsolver%p_rd2)

    call lsysbl_releaseVector (rsolver%p_rd)
    deallocate (rsolver%p_rd)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_done (rsolver)
  
!<description>
  ! Clean up the Newton iteration.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_spaceSolverHierarchy), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>

    type (t_spaceSolverHierarchy) :: rtemplate

    ! Initialise with standard parameters    
    rsolver = rtemplate

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_solve (rsolver,idofTime,csolgeneration,ceqnflags,rstatistics,&
      isollevelSpace,rprimalSol,rdualSol,rcontrol,&
      rprimalSolLin,rdualSolLin,rcontrolLin)
  
!<description>
  ! Solves the spatial linear/nonlinear system.
!</description>

!<inputoutput>
  ! Parameters for the iteration.
  ! The output parameters are changed according to the iteration.
  type(t_spaceSolverHierarchy), intent(inout) :: rsolver

  ! Number of the DOF in time which should be calculated into rdest.
  integer, intent(in) :: idofTime

  ! Defines a policy how to generate the solution to take the defect from.
  ! =0: Always take zero
  ! =1: Propagate the solution of the previous/next timestep to the
  !     current one.
  ! =2: Take the solution of the last space-time iteration
  integer, intent(in) :: csolgeneration
  
  ! Equation flags that specify modifications to the equation to solve.
  ! One of the SPACESLH_EQNF_xxxx constants.
  integer, intent(in) :: ceqnflags

  ! Space level corresponding to the solution structures.
  integer, intent(in) :: isollevelSpace

  ! Current solution of the primal equation.
  type(t_primalSpace), intent(inout), target :: rprimalSol

  ! Current solution of the dual equation.
  type(t_dualSpace), intent(inout), optional, target :: rdualSol

  ! Current solution of the control equation
  type(t_controlSpace), intent(inout), optional, target :: rcontrol

  ! Current solution of the primal equation. Defines the nonlinearity.
  ! The space level of the solution is specified by isollevelSpace.
  ! The time level must match rsolver%itimelevel.
  type(t_primalSpace), intent(inout), optional, target :: rprimalSolLin

  ! Current solution of the linearised dual equation.
  ! The space level of the solution is specified by isollevelSpace.
  ! The time level must match rsolver%itimelevel.
  type(t_dualSpace), intent(inout), optional, target :: rdualSolLin

  ! Current solution of the linearised control equation
  ! The space level of the solution is specified by isollevelSpace.
  ! The time level must match rsolver%itimelevel.
  type(t_controlSpace), intent(inout), optional, target :: rcontrolLin
!</inputoutput>

!<output>
  ! Statistics structure
  type(t_spaceslSolverStat), intent(out) :: rstatistics
!</output>

!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rd, p_rd2, p_rx
    integer :: ierror, coptype
    real(DP) :: dres
    type(t_linsolNode), pointer :: p_rsolverNode
    real(DP), dimension(:), pointer :: p_Dd, p_Dx
    type(t_lssSolverStat) :: rlocalStatLss
    type(t_spaceslSolverStat) :: rlocalStat
    type(t_timer) :: rtimer
   
    ! Clear statistics
    call stat_startTimer (rstatistics%rtotalTime)
   
    ! At first, get a temp vector we can use for creating defects,
    p_rd => rsolver%p_rd
    p_rd2 => rsolver%p_rd2
    
    ! DEBUG!!!
    call lsysbl_getbase_double (p_rd,p_Dd)

    ! Which type of operator is to be solved?
    select case (rsolver%coptype)

    ! ---------------------------------------------------------------
    ! Primal equation. Nonlinear loop.
    ! ---------------------------------------------------------------
    case (OPTP_PRIMAL)

      ! This is the most complicated part. On level ispacelevel, we have
      ! to apply a nonlinear loop.
      
      ! Apply the Newton iteration
      call itc_initIteration(rsolver%riter)

      do while (.true.)
      
        ! -------------------------------------------------------------
        ! Initialise boundary conditions
        ! -------------------------------------------------------------
        ! No filter for the initial condition.
        if (idoftime .gt. 1) then
          call stat_startTimer (rtimer)
          call spaceslh_initData_bdc (rsolver, idofTime, rcontrol)      
          call stat_stopTimer (rtimer)
        end if
      
        ! -------------------------------------------------------------
        ! Get the nonlinear defect
        ! -------------------------------------------------------------
        call stat_startTimer (rstatistics%rtimeDefect)
        
        ! Compute the basic (unpreconditioned) search direction in rd.
        if (rsolver%riter%niterations .eq. 0) then
          call smva_getDef_primal (p_rd,&
              rsolver%ispacelevel,rsolver%itimelevel,idofTime,&
              rsolver%p_roperatorAsmHier,rprimalSol,rcontrol,&
              rsolver%p_roptcBDCSpaceHierarchy%p_RoptcBDCspace(rsolver%ispacelevel),&
              rsolver%rtempData,csolgeneration,rtimer)
        else
          ! Start from the currently available iterate and continue the nonlinear
          ! iteration
          call smva_getDef_primal (p_rd,&
              rsolver%ispacelevel,rsolver%itimelevel,idofTime,&
              rsolver%p_roperatorAsmHier,rprimalSol,rcontrol,&
              rsolver%p_roptcBDCSpaceHierarchy%p_RoptcBDCspace(rsolver%ispacelevel),&
              rsolver%rtempData,SPINITCOND_PREVITERATE,rtimer)
        end if    
        
        call spaceslh_implementBDC (p_rd, SPACESLH_VEC_DEFECT, &
            rsolver%ispacelevel, rsolver)

        call stat_stopTimer (rstatistics%rtimeDefect)    
        call stat_addTimers (rtimer,rstatistics%rtimeRHS)
        
        dres = lsysbl_vectorNorm(p_rd,rsolver%rnewtonParams%iresNorm)
            
        if (rsolver%riter%cstatus .eq. ITC_STATUS_UNDEFINED) then
          ! Remember the initial residual
          call itc_initResidual(rsolver%riter,dres)
        else
          ! Push the residual, increase the iteration counter
          call itc_pushResidual(rsolver%riter,dres)
        end if
        
        if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
          call output_line (&
              trim(sys_si(idofTime-1,8)) // &
              " " // trim(sys_si(rsolver%riter%niterations,4)) // &
              " " // trim(sys_sdEL(dres,1)) )
        end if

        ! -------------------------------------------------------------
        ! Check for convergence / divergence / ...
        ! -------------------------------------------------------------
        if (rsolver%riter%cstatus .ne. ITC_STATUS_CONTINUE) exit

        ! -------------------------------------------------------------
        ! Adaptive Newton for the next iteration?
        ! -------------------------------------------------------------
        if (rsolver%rnewtonParams%ctypeIteration .eq. 3) then
          call spaceslh_adNewton_setEps (&
              rsolver,rsolver%rnewtonParams%radaptiveNewton,&
              rsolver%riter%dresInitial,dres)
        end if
      
        ! -------------------------------------------------------------
        ! Partial Newton for the next iteration?
        ! -------------------------------------------------------------

        ! By default use partial Newton.
        select case (rsolver%rnewtonParams%radaptiveNewton%cpartialNewton)
        case (NEWTN_PN_FULLNEWTON)
          coptype = OPTP_PRIMALLIN

        case (NEWTN_PN_PARTIALNEWTON)
          coptype = OPTP_PRIMALLIN_SIMPLE
        end select
        
        ! Should we switch to full Newton?
        if (rsolver%rnewtonParams%radaptiveNewton%cpartialNewton .ne. NEWTN_PN_FULLNEWTON) then
        
          if (rsolver%riter%niterations .gt. &
              rsolver%rnewtonParams%radaptiveNewton%nmaxPartialNewtonIterations) then
            ! Full Newton
            coptype = OPTP_PRIMALLIN
            
          else if (rsolver%riter%niterations .ge. &
              rsolver%rnewtonParams%radaptiveNewton%nminPartialNewtonIterations) then
            ! We are allowed to switch to the full Newton if some conditions
            ! are met.
            if (rsolver%rnewtonParams%radaptiveNewton%dtolAbsPartialNewton .gt. 0.0_DP) then
              ! Check the absolute residual
              if (dres .le. rsolver%rnewtonParams%radaptiveNewton%dtolAbsPartialNewton) then
                coptype = OPTP_PRIMALLIN
              end if
            end if
            if (rsolver%rnewtonParams%radaptiveNewton%dtolRelPartialNewton .gt. 0.0_DP) then
              ! Check the relative residual
              if (dres .le. &
                  rsolver%rnewtonParams%radaptiveNewton%dtolRelPartialNewton*rsolver%riter%dresInitial) then
                coptype = OPTP_PRIMALLIN
              end if
            end if
              
          end if
        
          ! Print out the choice
          if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
            select case (coptype)
            case (OPTP_PRIMALLIN)
              call output_line ("Space Newton: Full Newton selected.")
            case (OPTP_PRIMALLIN_SIMPLE)
              call output_line ("Space Newton: Partial Newton selected.")
            end select
          end if

        end if      

        ! -------------------------------------------------------------
        ! Preconditioning with the Newton matrix
        ! -------------------------------------------------------------

        ! -------------------------------------------------------------
        ! Assemble the matrices/boundary conditions on all levels.
        ! -------------------------------------------------------------
        ! Copy the defect vector in case we have to use the backup solver.
        call lsysbl_copyVector (p_rd,p_rd2)

        ! Assemble the matrices/boundary conditions on all levels.
        call spaceslh_initData_primal (rsolver, ierror, idofTime, &
            rprimalSol,rsolver%ispacelevel,rlocalStat,coptype,.false.)
            
        call spacesl_sumStatistics (rlocalStat,rstatistics,.false.)
        
        ! -------------------------------------------------------------
        ! Call the linear solver
        ! -------------------------------------------------------------
        output_iautoOutputIndent = output_iautoOutputIndent + 2
        call lssh_precondDefect (&
            rsolver%p_rlsshierarchy,rsolver%ispacelevel,p_rd,rlocalStatLss,p_rsolverNode)
        output_iautoOutputIndent = output_iautoOutputIndent - 2
        
        call lss_sumStatistics (rlocalStatLss,rstatistics%rlssSolverStat)

        ! Cleanup
        call spaceslh_doneData (rsolver,.false.)
        
        ! -------------------------------------------------------------
        ! Call the fallback solver if the previous solver failed
        ! -------------------------------------------------------------
        if (p_rsolverNode%iresult .ne. 0) then

          if (rsolver%rnewtonParams%ioutputLevel .ge. 1) then
            call output_line ("Linear solver in space failed in timestep "//&
                trim(sys_siL(idofTime,10))//". Invoking fallback solver.", &
                OU_CLASS_WARNING,OU_MODE_STD,"spaceslh_solve")
          end if

          ! Reconstruct the defect vector
          call lsysbl_copyVector (p_rd2,p_rd)
          
          call spaceslh_initData_primal (rsolver, ierror, idofTime, &
              rprimalSol,rsolver%ispacelevel,rlocalStat,coptype,.true.)
              
          call spacesl_sumStatistics (rlocalStat,rstatistics,.false.)
          
          ! Solve the system
          output_iautoOutputIndent = output_iautoOutputIndent + 2
          call lssh_precondDefect (&
              rsolver%p_rlsshierarchy2,rsolver%ispacelevel,p_rd,rlocalStatLss,p_rsolverNode)
          output_iautoOutputIndent = output_iautoOutputIndent - 2
          
          call lss_sumStatistics (rlocalStatLss,rstatistics%rlssSolverStat)

          ! Cleanup
          call spaceslh_doneData (rsolver,.true.)
          
          ! Ignore the solution if everything fails.
          if (p_rsolverNode%iresult .ne. 0) then
            call output_line ("Linear solver in space failed. Solution ignored.", &
                OU_CLASS_WARNING,OU_MODE_STD,"spaceslh_solve")
            call lsysbl_clearVector (p_rd)
          end if
        end if
        
        ! -------------------------------------------------------------
        ! Update of the solution
        ! -------------------------------------------------------------

        ! Update the control according to the search direction:
        !
        !    u_n+1  =  u_n  +  g_n
        !
        ! or to any configured step-length control rule.
        call stat_startTimer (rtimer)
        call sptivec_getVectorFromPool(rprimalSol%p_rvectorAccess,idofTime,p_rx)
        
        ! DEBUG!!!
        call lsysbl_getbase_double (p_rx,p_Dx)
        
        call lsysbl_vectorLinearComb (p_rd,p_rx,rsolver%rnewtonParams%domega,1.0_DP)

        ! Implement boundary conditions
        call spaceslh_implementBDC (p_rx, SPACESLH_VEC_SOLUTION, &
            rsolver%ispacelevel, rsolver)
        
        call sptivec_commitVecInPool(rprimalSol%p_rvectorAccess,idofTime)
        call stat_stopTimer (rtimer)        
        ! -------------------------------------------------------------
        ! Clean up boundary conditions
        ! -------------------------------------------------------------
        call sbch_resetBCstructure (rsolver%p_roptcBDCSpaceHierarchy)

        ! -------------------------------------------------------------
        ! Proceed with the next iteration
        ! -------------------------------------------------------------
      
      end do
      
      rstatistics%niterations = rsolver%riter%niterations
      
    ! ---------------------------------------------------------------
    ! Dual equation. Linear loop.
    ! ---------------------------------------------------------------
    case (OPTP_DUAL)
      ! We apply one defect correction, so to pass a defect to the 
      ! linear solver.

      ! -------------------------------------------------------------
      ! Initialise boundary conditions
      ! -------------------------------------------------------------
      call spaceslh_initData_bdc (rsolver, idofTime)      

      ! -------------------------------------------------------------
      ! Create a defect
      ! -------------------------------------------------------------
      call stat_startTimer (rstatistics%rtimeDefect)
      
      call smva_getDef_dual (p_rd,&
          rsolver%ispacelevel,rsolver%itimelevel,idofTime,&
          rsolver%p_roperatorAsmHier,rprimalSol,rdualSol,&
          rsolver%p_roptcBDCSpaceHierarchy%p_RoptcBDCspace(rsolver%ispacelevel),&
          rsolver%rtempData,csolgeneration,rtimer)

      call spaceslh_implementBDC (p_rd, SPACESLH_VEC_DEFECT, &
          rsolver%ispacelevel, rsolver)
      
      call stat_stopTimer (rstatistics%rtimeDefect)
      call stat_addTimers (rtimer,rstatistics%rtimeRHS)
      
      ! Get the residual
      dres = lsysbl_vectorNorm(p_rd,rsolver%rnewtonParams%iresNorm)
      call itc_initResidual (rsolver%riter,dres)
      
      ! Print the initial residual
      if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
        call output_line (&
            trim(sys_si(idofTime-1,8)) // &
            " " // trim(sys_si(0,4)) // &
            " " // trim(sys_sdEL(dres,1)) )
      end if

      ! Cancel if the defect is zero.
      if (dres .gt. 10.0_DP*SYS_EPSREAL_DP) then

        ! -------------------------------------------------------------
        ! Assemble the matrices/boundary conditions on all levels.
        ! -------------------------------------------------------------
        ! Copy the defect vector in case we have to use the backup solver.
        call lsysbl_copyVector (p_rd,p_rd2)

        call spaceslh_initData_dual (rsolver, ierror, idofTime, &
            rprimalSol,rsolver%ispacelevel,rlocalStat,.false.)
            
        call spacesl_sumStatistics (rlocalStat,rstatistics,.false.)

        ! -------------------------------------------------------------
        ! Call the linear solver, it does the job for us.
        ! -------------------------------------------------------------
        call lssh_precondDefect (&
            rsolver%p_rlsshierarchy,rsolver%ispacelevel,p_rd,rlocalStatLss,p_rsolverNode)
        
        call lss_sumStatistics (rlocalStatLss,rstatistics%rlssSolverStat)

        ! -------------------------------------------------------------
        ! Solver-Cleanup
        ! -------------------------------------------------------------
        call spaceslh_doneData (rsolver,.false.)
        
        ! -------------------------------------------------------------
        ! Call the fallback solver if the previous solver failed
        ! -------------------------------------------------------------
        if (p_rsolverNode%iresult .ne. 0) then

          if (rsolver%rnewtonParams%ioutputLevel .ge. 1) then
            call output_line ("Linear solver in space failed in timestep "//&
                trim(sys_siL(idofTime,10))//". Invoking fallback solver.", &
                OU_CLASS_WARNING,OU_MODE_STD,"spaceslh_solve")
          end if

          ! Reconstruct the defect vector
          call lsysbl_copyVector (p_rd2,p_rd)

          call spaceslh_initData_dual (rsolver, ierror, idofTime, &
              rprimalSol,rsolver%ispacelevel,rlocalStat,.true.)
              
          call spacesl_sumStatistics (rlocalStat,rstatistics,.false.)

          call lssh_precondDefect (&
              rsolver%p_rlsshierarchy2,rsolver%ispacelevel,p_rd,rlocalStatLss,p_rsolverNode)
          
          call lss_sumStatistics (rlocalStatLss,rstatistics%rlssSolverStat)

          call spaceslh_doneData (rsolver,.true.)

          ! Ignore the solution if everything fails.
          if (p_rsolverNode%iresult .ne. 0) then
            call output_line ("Linear solver in space failed. Solution ignored.", &
                OU_CLASS_WARNING,OU_MODE_STD,"spaceslh_solve")
            call lsysbl_clearVector (p_rd)
          end if
        end if      

        ! -------------------------------------------------------------
        ! Statistics
        ! -------------------------------------------------------------

        ! Get the final residual
        dres = p_rsolverNode%dfinalDefect
        call itc_pushResidual (rsolver%riter,dres)
      
        ! Print the final residual
        if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
          call output_line (&
              trim(sys_si(idofTime-1,8)) // &
              " " // trim(sys_si(1,4)) // &
              " " // trim(sys_sdEL(dres,1)) )
        end if

        ! -------------------------------------------------------------
        ! Update the control according to the defect:
        !
        !    u_new  =  u_old  +  g_n
        ! -------------------------------------------------------------

        call sptivec_getVectorFromPool(rdualSol%p_rvectorAccess,idofTime,p_rx)
        
        ! DEBUG!!!
        call lsysbl_getbase_double (p_rx,p_Dx)
        
        call lsysbl_vectorLinearComb (p_rd,p_rx,1.0_DP,1.0_DP)

        ! Implement boundary conditions
        call spaceslh_implementBDC (p_rx, SPACESLH_VEC_SOLUTION, &
              rsolver%ispacelevel, rsolver)

        call sptivec_commitVecInPool(rdualSol%p_rvectorAccess,idofTime)

      end if
    
      ! -------------------------------------------------------------
      ! Clean up boundary conditions
      ! -------------------------------------------------------------
      call sbch_resetBCstructure (rsolver%p_roptcBDCSpaceHierarchy)
      
      ! Only one iteration
      rsolver%riter%niterations = 1
      rstatistics%niterations = 1

    ! ---------------------------------------------------------------
    ! Linearised primal equation. Linear loop.
    ! ---------------------------------------------------------------
    case (OPTP_PRIMALLIN,OPTP_PRIMALLIN_SIMPLE)
      ! We apply one defect correction, so to pass a defect to the 
      ! linear solver.

      ! -------------------------------------------------------------
      ! Initialise boundary conditions
      ! -------------------------------------------------------------
      call stat_startTimer (rtimer)
      call spaceslh_initData_bdc (rsolver, idofTime, rcontrolLin)      
      call stat_stopTimer (rtimer)

      ! -------------------------------------------------------------
      ! Create a defect
      ! -------------------------------------------------------------
      call stat_startTimer (rstatistics%rtimeDefect)
      
      ! Type of the operator.
      ! If the flags disable the primal Newton, force PRIMALLIN_SIMPLE.
      coptype = rsolver%coptype
      if (iand(ceqnflags,SPACESLH_EQNF_NONEWTONPRIMAL) .ne. 0) then
        coptype = OPTP_PRIMALLIN_SIMPLE
      end if
      
      call smva_getDef_primalLin (p_rd,&
          rsolver%ispacelevel,rsolver%itimelevel,idofTime,&
          rsolver%p_roperatorAsmHier,rprimalSol,rprimalSolLin,rcontrolLin,&
          coptype,&
          rsolver%p_roptcBDCSpaceHierarchy%p_RoptcBDCspace(rsolver%ispacelevel),&
          rsolver%rtempData,csolgeneration,rtimer)
          
      call spaceslh_implementBDC (p_rd, SPACESLH_VEC_DEFECT, &
          rsolver%ispacelevel, rsolver)

      call stat_stopTimer (rstatistics%rtimeDefect)
      call stat_addTimers (rtimer,rstatistics%rtimeRHS)
          
      ! Get the residual
      dres = lsysbl_vectorNorm(p_rd,rsolver%rnewtonParams%iresNorm)
      call itc_initResidual (rsolver%riter,dres)

      ! Print the initial residual
      if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
        call output_line (&
            trim(sys_si(idofTime-1,8)) // &
            " " // trim(sys_si(0,4)) // &
            " " // trim(sys_sdEL(dres,1)) )
      end if

      ! Cancel if the defect is zero.
      if (dres .gt. 10.0_DP*SYS_EPSREAL_DP) then

        ! -------------------------------------------------------------
        ! Assemble the matrices/boundary conditions on all levels.
        ! -------------------------------------------------------------
        ! Copy the defect vector in case we have to use the backup solver.
        call lsysbl_copyVector (p_rd,p_rd2)

        call spaceslh_initData_primalLin (rsolver, ierror, idofTime, &
            rprimalSol,rsolver%ispacelevel,coptype,rlocalStat,.false.)
            
        call spacesl_sumStatistics (rlocalStat,rstatistics,.false.)

        ! -------------------------------------------------------------
        ! Call the linear solver, it does the job for us.
        ! -------------------------------------------------------------
        call lssh_precondDefect (&
            rsolver%p_rlsshierarchy,rsolver%ispacelevel,p_rd,rlocalStatLss,p_rsolverNode)
        
        call lss_sumStatistics (rlocalStatLss,rstatistics%rlssSolverStat)

        ! -------------------------------------------------------------
        ! Solver-Cleanup
        ! -------------------------------------------------------------
        call spaceslh_doneData (rsolver,.false.)

        ! -------------------------------------------------------------
        ! Call the fallback solver if the previous solver failed
        ! -------------------------------------------------------------
        if (p_rsolverNode%iresult .ne. 0) then

          if (rsolver%rnewtonParams%ioutputLevel .ge. 1) then
            call output_line ("Linear solver in space failed in timestep "//&
                trim(sys_siL(idofTime,10))//". Invoking fallback solver.", &
                OU_CLASS_WARNING,OU_MODE_STD,"spaceslh_solve")
          end if

          ! Reconstruct the defect vector
          call lsysbl_copyVector (p_rd2,p_rd)

          call spaceslh_initData_primalLin (rsolver, ierror, idofTime, &
              rprimalSol,rsolver%ispacelevel,coptype,rlocalStat,.true.)
              
          call spacesl_sumStatistics (rlocalStat,rstatistics,.false.)

          call lssh_precondDefect (&
              rsolver%p_rlsshierarchy2,rsolver%ispacelevel,p_rd,rlocalStatLss,p_rsolverNode)
          
          call lss_sumStatistics (rlocalStatLss,rstatistics%rlssSolverStat)
          
          call spaceslh_doneData (rsolver,.true.)

          ! Ignore the solution if everything fails.
          if (p_rsolverNode%iresult .ne. 0) then
            call output_line ("Linear solver in space failed. Solution ignored.", &
                OU_CLASS_WARNING,OU_MODE_STD,"spaceslh_solve")
            call lsysbl_clearVector (p_rd)
          end if
        end if
        
        ! -------------------------------------------------------------
        ! Statistics
        ! -------------------------------------------------------------

        ! Get the final residual
        dres = p_rsolverNode%dfinalDefect
        call itc_pushResidual (rsolver%riter,dres)
      
        ! Print the final residual
        if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
          call output_line (&
              trim(sys_si(idofTime-1,8)) // &
              " " // trim(sys_si(1,4)) // &
              " " // trim(sys_sdEL(dres,1)) )
        end if

        ! -------------------------------------------------------------
        ! Update the control according to the defect
        !
        !    u_new  =  u_old  +  g_n
        ! -------------------------------------------------------------

        call stat_startTimer (rtimer)
        call sptivec_getVectorFromPool(rprimalSolLin%p_rvectorAccess,idofTime,p_rx)
        
        ! DEBUG!!!
        call lsysbl_getbase_double (p_rx,p_Dx)
        
        call lsysbl_vectorLinearComb (p_rd,p_rx,1.0_DP,1.0_DP)

        ! Implement boundary conditions.
        call spaceslh_implementBDC (p_rx, SPACESLH_VEC_SOLUTION, &
              rsolver%ispacelevel, rsolver)

        call sptivec_commitVecInPool(rprimalSolLin%p_rvectorAccess,idofTime)
        call stat_startTimer (rtimer)

      end if

      ! -------------------------------------------------------------
      ! Clean up boundary conditions
      ! -------------------------------------------------------------
      call sbch_resetBCstructure (rsolver%p_roptcBDCSpaceHierarchy)

      ! Only one iteration
      rsolver%riter%niterations = 1
      rstatistics%niterations = 1

    ! ---------------------------------------------------------------
    ! Linearised dual equation. Linear loop.
    ! ---------------------------------------------------------------
    case (OPTP_DUALLIN,OPTP_DUALLIN_SIMPLE)
      ! We apply one defect correction, so to pass a defect to the 
      ! linear solver.

      ! -------------------------------------------------------------
      ! Initialise boundary conditions
      ! -------------------------------------------------------------
      call spaceslh_initData_bdc (rsolver, idofTime)      

      ! -------------------------------------------------------------
      ! Create a defect
      ! -------------------------------------------------------------
      call stat_startTimer (rstatistics%rtimeDefect)
      
      ! Type of the operator.
      ! If the flags disable the primal Newton, force DUALLIN_SIMPLE.
      coptype = rsolver%coptype
      if ((iand(ceqnflags,SPACESLH_EQNF_NONEWTONPRIMAL) .ne. 0) .or. &
          (iand(ceqnflags,SPACESLH_EQNF_NONEWTONDUAL) .ne. 0)) then
        coptype = OPTP_DUALLIN_SIMPLE
      end if

      call smva_getDef_dualLin (p_rd,&
          rsolver%ispacelevel,rsolver%itimelevel,idofTime,&
          rsolver%p_roperatorAsmHier,rprimalSol,rdualSol,rprimalSolLin,rdualSolLin,&
          coptype,&
          rsolver%p_roptcBDCSpaceHierarchy%p_RoptcBDCspace(rsolver%ispacelevel),&
          rsolver%rtempData,csolgeneration,rtimer)
          
      call spaceslh_implementBDC (p_rd, SPACESLH_VEC_DEFECT, &
          rsolver%ispacelevel, rsolver)

      call stat_stopTimer (rstatistics%rtimeDefect)
      call stat_addTimers (rtimer,rstatistics%rtimeRHS)
          
      ! Get the residual
      dres = lsysbl_vectorNorm(p_rd,rsolver%rnewtonParams%iresNorm)
      call itc_initResidual (rsolver%riter,dres)

      ! Print the initial residual
      if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
        call output_line (&
            trim(sys_si(idofTime-1,8)) // &
            " " // trim(sys_si(0,4)) // &
            " " // trim(sys_sdEL(dres,1)) )
      end if

      ! Cancel if the defect is zero.
      if (dres .gt. 10.0_DP*SYS_EPSREAL_DP) then

        ! -------------------------------------------------------------
        ! Assemble the matrices/boundary conditions on all levels.
        ! -------------------------------------------------------------
        ! Copy the defect vector in case we have to use the backup solver.
        call lsysbl_copyVector (p_rd,p_rd2)

        call spaceslh_initData_dualLin (rsolver, ierror, idofTime, &
            rprimalSol,rsolver%ispacelevel,rlocalStat,.false.)
        
        call spacesl_sumStatistics (rlocalStat,rstatistics,.false.)

        ! -------------------------------------------------------------
        ! Call the linear solver, it does the job for us.
        ! -------------------------------------------------------------
        call lssh_precondDefect (&
            rsolver%p_rlsshierarchy,rsolver%ispacelevel,p_rd,rlocalStatLss,p_rsolverNode)
        
        call lss_sumStatistics (rlocalStatLss,rstatistics%rlssSolverStat)

        ! -------------------------------------------------------------
        ! Solver-Cleanup
        ! -------------------------------------------------------------
        call spaceslh_doneData (rsolver,.false.)

        ! -------------------------------------------------------------
        ! Call the fallback solver if the previous solver failed
        ! -------------------------------------------------------------
        if (p_rsolverNode%iresult .ne. 0) then

          if (rsolver%rnewtonParams%ioutputLevel .ge. 1) then
            call output_line ("Linear solver in space failed in timestep "//&
                trim(sys_siL(idofTime,10))//". Invoking fallback solver.", &
                OU_CLASS_ERROR,OU_MODE_STD,"spaceslh_solve")
          end if

          ! Reconstruct the defect vector
          call lsysbl_copyVector (p_rd2,p_rd)

          call spaceslh_initData_dualLin (rsolver, ierror, idofTime, &
              rprimalSol,rsolver%ispacelevel,rlocalStat,.true.)
          
          call spacesl_sumStatistics (rlocalStat,rstatistics,.false.)

          call lssh_precondDefect (&
              rsolver%p_rlsshierarchy2,rsolver%ispacelevel,p_rd,rlocalStatLss,p_rsolverNode)
          
          call lss_sumStatistics (rlocalStatLss,rstatistics%rlssSolverStat)

          call spaceslh_doneData (rsolver,.true.)

          ! Ignore the solution if everything fails.
          if (p_rsolverNode%iresult .ne. 0) then
            call output_line ("Linear solver in space failed. Solution ignored.", &
                OU_CLASS_ERROR,OU_MODE_STD,"spaceslh_solve")
            call lsysbl_clearVector (p_rd)
          end if

        end if      

        ! -------------------------------------------------------------
        ! Statistics
        ! -------------------------------------------------------------

        ! Get the final residual
        dres = p_rsolverNode%dfinalDefect
        call itc_pushResidual (rsolver%riter,dres)
      
        ! Print the final residual
        if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
          call output_line (&
              trim(sys_si(idofTime-1,8)) // &
              " " // trim(sys_si(1,4)) // &
              " " // trim(sys_sdEL(dres,1)) )
        end if

        ! -------------------------------------------------------------
        ! Update the control according to the defect
        !
        !    u_new  =  u_old  +  g_n
        ! -------------------------------------------------------------

        call sptivec_getVectorFromPool(rdualSolLin%p_rvectorAccess,idofTime,p_rx)

        ! DEBUG!!!
        call lsysbl_getbase_double (p_rx,p_Dx)
        
        call lsysbl_vectorLinearComb (p_rd,p_rx,1.0_DP,1.0_DP)

        ! Implement boundary conditions.
        ! Use the defect filter since we are in the defect space.
        call spaceslh_implementBDC (p_rx, SPACESLH_VEC_SOLUTION, &
              rsolver%ispacelevel, rsolver)
        
        call sptivec_commitVecInPool(rdualSolLin%p_rvectorAccess,idofTime)

      end if

      ! -------------------------------------------------------------
      ! Clean up boundary conditions
      ! -------------------------------------------------------------
      call sbch_resetBCstructure (rsolver%p_roptcBDCSpaceHierarchy)

      ! Only one iteration
      rsolver%riter%niterations = 1
      rstatistics%niterations = 1

    end select

    ! Measure the total time
    call stat_stopTimer (rstatistics%rtotalTime)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_setEps (rsolver,depsAbs,depsRel)
  
!<description>
  ! Sets the stopping criterion for the use application of the adaptive
  ! Newton algorithm.
  !
  ! WARNING!!! THIS ROUTINE HAS A SIDE EFFECT!
  ! IT SETS THE STOPPING CRITERION OF ALL SOLVERS IN SPACE TO AN APPROPRIATE
  ! VALUE! IF THE SPACE SOVLERS AE USED SOMEWHERE ELSE, THE STOPPING CRITERION
  ! IS LOST AND THUS, THEY MAY BEHAVE NOT AS EXPECTED!!!
!</description>

!<inputoutput>
  ! Parameters for the iteration.
  ! The output parameters are changed according to the iteration.
  type(t_spaceSolverHierarchy), intent(inout) :: rsolver
  
  ! New absolute stopping criterion. Absolute residual.
  ! =0.0: Switch off the check.
  real(DP), intent(in) :: depsAbs

  ! New absolute stopping criterion. Relative residual.
  ! =0.0: Switch off the check.
  real(DP), intent(in) :: depsRel
!</inputoutput>

!</subroutine>

    if (rsolver%copType .eq. OPTP_PRIMAL) then
      
      ! Forward solver is nonlinear. 
      ! Modify the stopping criteria of the nonlinear solver.
      rsolver%riter%dtolRel = depsRel
      rsolver%riter%dtolAbs = depsAbs
      if (depsAbs .eq. 0.0_DP) then
        rsolver%riter%ctolMode = ITC_STOP_MODE_REL
      else
        rsolver%riter%ctolMode = ITC_STOP_MODE_ABS
      end if
    
    else

      ! Set the eps of the preconditioner.
      call spaceslh_setEps_prec (rsolver,depsAbs,depsRel)

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_setEps_prec (rsolver,depsAbs,depsRel)
  
!<description>
  ! Sets the stopping criterion of the preconditioner for the use 
  ! application of the adaptive Newton algorithm.
  !
  ! WARNING!!! THIS ROUTINE HAS A SIDE EFFECT!
  ! IT SETS THE STOPPING CRITERION OF ALL SOLVERS IN SPACE TO AN APPROPRIATE
  ! VALUE! IF THE SPACE SOVLERS AE USED SOMEWHERE ELSE, THE STOPPING CRITERION
  ! IS LOST AND THUS, THEY MAY BEHAVE NOT AS EXPECTED!!!
!</description>

!<inputoutput>
  ! Parameters for the iteration.
  ! The output parameters are changed according to the iteration.
  type(t_spaceSolverHierarchy), intent(inout) :: rsolver
  
  ! New absolute stopping criterion. Absolute residual.
  ! =0.0: Switch off the check.
  real(DP), intent(in) :: depsAbs

  ! New absolute stopping criterion. Relative residual.
  ! =0.0: Switch off the check.
  real(DP), intent(in) :: depsRel
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ilev
    type(t_linsolSpace), pointer :: p_rlinsol

    ! Loop through all levels
    do ilev = rsolver%p_rlssHierarchy%nlmin,rsolver%p_rlssHierarchy%nlmax
    
      ! Get the space solver
      p_rlinsol => rsolver%p_rlssHierarchy%p_RlinearSolvers(ilev)
      
      ! Solver type?
      select case (p_rlinsol%isolverType)
      
      ! ------------------------------
      ! UMFPACK. Nothing to do.
      ! ------------------------------
      case (LSS_LINSOL_UMFPACK)

      ! -------------------------------
      ! Multigrid
      ! -------------------------------
      case (LSS_LINSOL_MG)
      
        ! Modify the stopping criterion of the solver.
        p_rlinsol%p_rsolverNode%depsRel = depsRel
        p_rlinsol%p_rsolverNode%depsAbs = depsAbs
        
      case default
        call output_line ("Unknown solver in space.", &
            OU_CLASS_ERROR,OU_MODE_STD,"spaceslh_adNewton_setEps")
        call sys_halt()
        
      end select
    
    end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_adNewton_setEps (rsolver,radNewtonParams,dresInitial,dresLastIte)
  
!<description>
  ! Realises the adaptive Newton algorithm. Sets the stopping criterions of
  ! all solvers to appropriate values.
  !
  ! WARNING!!! THIS ROUTINE HAS A SIDE EFFECT!
  ! IT SETS THE STOPPING CRITERION OF ALL SOLVERS IN SPACE TO AN APPROPRIATE
  ! VALUE! IF THE SPACE SOVLERS AE USED SOMEWHERE ELSE, THE STOPPING CRITERION
  ! IS LOST AND THUS, THEY MAY BEHAVE NOT AS EXPECTED!!!
!</description>

!<input>
  ! Parameters of the adaptive Newton algotithm.
  type(t_ccDynamicNewtonControl), intent(in) :: radNewtonParams
  
  ! Initial residual.
  ! May be set to 0.0 if there is no initial residual.
  real(DP), intent(in) :: dresInitial
  
  ! Residual obtained in the last nonlinear iteration.
  ! In the first call, this should be set to dresInitial.
  real(DP), intent(in) :: dresLastIte
!</input>

!<inputoutput>
  ! Parameters for the Newton iteration.
  ! The output parameters are changed according to the iteration.
  type(t_spaceSolverHierarchy), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>

    real(DP) :: depsAbs,depsRel,ddigitsGained, ddigitsToGain
    
    if ((dresInitial .eq. 0.0_DP) .and. (dresLastIte .eq. 0.0_DP)) then
      
      ! Gain two digits for the initial residual, that is enough
      depsAbs = 0.0_DP
      depsRel = 1.0E-2_DP
      
      if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
        call output_line ("Adaptive Newton: New stopping criterion. ||res_rel|| < "//&
            trim(sys_sdEL(depsRel,10)))
      end if

    else
    
      ! We have to determine a new dresAbs for all solver components:
      ! - At least, gain as many digits as configured in the adaptive-Newton 
      !   structure has to be gained
      ! - The number of digits to gain in the next iteration has to be an
      !   appropriate multiple of the digits already gained.
      ! So...
      
      ddigitsGained = dresLastIte/dresInitial
      
      ddigitsToGain = min(radNewtonParams%dinexactNewtonTolRel*ddigitsGained,&
          ddigitsGained ** radNewtonParams%dinexactNewtonExponent)
          
      depsRel = 0.0_DP
      depsAbs = max(dresInitial * ddigitsToGain,radNewtonParams%dinexactNewtonTolAbs)
      
      ! Do not gain too much.
      depsAbs = max(depsAbs,&
          max(dresInitial * rsolver%riter%dtolrel * radNewtonParams%dinexactNewtonTolRel,&
              rsolver%riter%dtolabs * radNewtonParams%dinexactNewtonTolRel))

      if (rsolver%rnewtonParams%ioutputLevel .ge. 3) then
        call output_line ("Adaptive Newton: New stopping criterion. ||res|| < "//&
            trim(sys_sdEL(depsAbs,10)))
      end if
    end if
    
    ! Initialise the nonlinear and linear solver in space which are used
    ! for the calvulation of the current residual.
    call spaceslh_setEps_prec (rsolver,depsAbs,depsRel)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spacesl_clearStatistics(rstatistics)
  
!<description>
  ! Resets a statistic structure.
!</description>

!<inputoutput>
  ! Structure to be reset.
  type(t_spaceslSolverStat), intent(inout) :: rstatistics
!</inputoutput>

!</subroutine>

    rstatistics%niterations = 0
    
    call stat_clearTimer(rstatistics%rtotalTime)
    call stat_clearTimer(rstatistics%rtimeDefect)
    call stat_clearTimer(rstatistics%rtimeMatrixAssembly)
    
    call lss_clearStatistics (rstatistics%rlssSolverStat)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spacesl_sumStatistics(rstatistics1,rstatistics2,btotalTime)
  
!<description>
  ! Sums up the data of rstatistics1 to the data in rstatistics2.
!</description>

!<input>
  ! Source structure
  type(t_spaceslSolverStat), intent(in) :: rstatistics1
  
  ! OPTIONAL: Whether or not to sum up the total time.
  ! If not present, TRUE is assumed.
  logical, intent(in), optional :: btotalTime
!</input>

!<inputoutput>
  ! Destination structure.
  type(t_spaceslSolverStat), intent(inout) :: rstatistics2
!</inputoutput>

!</subroutine>

    rstatistics2%niterations = rstatistics2%niterations + rstatistics1%niterations
    
    if (.not. present(btotalTime)) then
      call stat_addTimers(rstatistics1%rtotalTime,rstatistics2%rtotalTime)
    else if (btotalTime) then
      call stat_addTimers(rstatistics1%rtotalTime,rstatistics2%rtotalTime)
    end if
    
    call stat_addTimers(rstatistics1%rtimeDefect,rstatistics2%rtimeDefect)
    call stat_addTimers(rstatistics1%rtimeRHS,rstatistics2%rtimeRHS)
    call stat_addTimers(rstatistics1%rtimeMatrixAssembly,rstatistics2%rtimeMatrixAssembly)
    
    call lss_sumStatistics (rstatistics1%rlssSolverStat,rstatistics2%rlssSolverStat)

  end subroutine

end module
