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
    
    ! <!-- ------------------------ -->
    ! <!-- LINEAR SOLVER STRUCTURES -->
    ! <!-- ------------------------ -->

    ! Hierarchy of linear solvers in space used for solving
    ! auxiliary linear subproblems. The solver on level ilevel
    ! will be used to solve linear subproblems.
    type(t_linsolHierarchySpace), pointer :: p_rlssHierarchy => null()

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

    ! Hierarchy of system matrices
    type(t_matrixBlock), dimension(:), pointer :: p_Rmatrices => null()
    
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

contains

  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_init (rsolver,copType,rlssHierarchy,ssection,rparamList)
  
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
  type(t_linsolHierarchySpace), target :: rlssHierarchy

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
    rsolver%copType = copType

    if (present(ssection) .and. present(rparamList)) then
    
      ! Initialise basic parameters of the nonlinear solver
      call newtonit_initBasicParams (rsolver%rnewtonParams,ssection,rparamList)
    
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_initStructure (rsolver, ispacelevel, itimelevel,&
      roperatorAsmHier,ierror)
  
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
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!</subroutine>

    integer :: ilev
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    type(t_spacetimeOperatorAsm) :: roperatorAsm
    type(t_discreteBC), pointer :: p_rdiscreteBC

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
    
    do ilev = rsolver%p_rlssHierarchy%nlmin,rsolver%ispaceLevel
    
      ! Initialise a temporary vector for the nonlinearity
      select case (rsolver%copType)
      case (OPTP_PRIMAL,OPTP_PRIMALLIN)
        p_rdiscretisation => &
            roperatorAsmHier%p_rfeHierarchyPrimal%p_rfeSpaces(ilev)%p_rdiscretisation

      case (OPTP_DUAL,OPTP_DUALLIN)
        p_rdiscretisation => &
            roperatorAsmHier%p_rfeHierarchyDual%p_rfeSpaces(ilev)%p_rdiscretisation
      
      case default
        call output_line ("Unknown space.", &
            OU_CLASS_ERROR,OU_MODE_STD,"spaceslh_initStructure")
        call sys_halt()
      end select
      
      ! Get the corresponding operator assembly structure and
      ! initialise a system matrix
      call stoh_getOpAsm_slvtlv (&
          roperatorAsm,roperatorAsmHier,ilev,rsolver%itimelevel)
          
      call smva_allocSystemMatrix (rsolver%p_Rmatrices(ilev),&
          roperatorAsmHier%ranalyticData%p_rphysics,&
          roperatorAsmHier%ranalyticData%p_rsettingsOptControl,&
          roperatorAsmHier%ranalyticData%p_rsettingsSpaceDiscr,&
          p_rdiscretisation,roperatorAsm%p_rasmTemplates)
          
      ! Bind discrete boundary conditions to the matrix
      call sbch_getDiscreteBC (rsolver%p_roptcBDCSpaceHierarchy,ilev,p_rdiscreteBC)
      call lsysbl_assignDiscreteBC (rsolver%p_Rmatrices(ilev),p_rdiscreteBC)
          
      ! Provide the matrix to the linear solver
      call lssh_setMatrix(rsolver%p_rlsshierarchy,ilev,rsolver%p_Rmatrices(ilev))
      
      if (ilev .eq. rsolver%ispaceLevel) then
      
        ! On the highest space level, allocate some temp vectors on the topmost level

        allocate (rsolver%p_rd)
        call lsysbl_createVectorBlock (p_rdiscretisation,rsolver%p_rd)
        
        ! Bind the boundary conditions to that vector
        call lsysbl_assignDiscreteBC (rsolver%p_Rmatrices(ilev),p_rdiscreteBC)
            
        select case (rsolver%copType)
        case (OPTP_PRIMAL,OPTP_PRIMALLIN)
          call smva_allocTempData (rsolver%rtempData,&
              roperatorAsmHier%ranalyticData%p_rphysics,&
              roperatorAsmHier%ranalyticData%p_rsettingsOptControl,&
              roperatorAsmHier%ranalyticData%p_rsettingsSpaceDiscr,&
              ilev,roperatorAsmHier)

        case (OPTP_DUAL,OPTP_DUALLIN)
          call smva_allocTempData (rsolver%rtempData,&
              roperatorAsmHier%ranalyticData%p_rphysics,&
              roperatorAsmHier%ranalyticData%p_rsettingsOptControl,&
              roperatorAsmHier%ranalyticData%p_rsettingsSpaceDiscr,&
              ilev,roperatorAsmHier)
        end select

        ! Initialise the structures of the associated linear subsolver
        call lssh_initStructure (rsolver%p_rlsshierarchy,ilev,ierror)

      end if
    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_initData_bdc (rsolver, idofTime)
  
!<description>
  ! Initialises the boundary condition structures for time DOF idofTime
!</description>

!<input>
  ! Number of the DOF in time.
  integer, intent(in) :: idofTime
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_spaceSolverHierarchy), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ilev
    
    ! Clean up the boundary conditions.
    call sbch_resetBCstructure (rsolver%p_roptcBDCSpaceHierarchy)
    
    ! Loop over all levels, from the highest to the lowest,
    ! and set up the boundary conditions.
    do ilev = rsolver%ispacelevel,rsolver%p_rlssHierarchy%nlmin,-1
    
      ! Assemble Dirichlet/Neumann boundary conditions
      call smva_initDirichletNeumannBC (&
          rsolver%p_roptcBDCSpaceHierarchy%p_RoptcBDCspace(ilev),&
          rsolver%p_roperatorAsmHier%ranalyticData%p_roptcBDC,&
          rsolver%copType,&
          rsolver%p_roperatorAsmHier,ilev,rsolver%itimelevel,idoftime,&
          rsolver%p_roperatorAsmHier%ranalyticData%p_rglobalData)
    
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
      rprimalSol, isollevelSpace)
  
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
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_spaceSolverHierarchy), intent(inout) :: rsolver
!</inputoutput>

!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!</subroutine>

    ! local variables
    integer :: ilev,iprevlv
    logical :: bfull
    
    ! Full Newton?
    bfull = rsolver%rnewtonParams%ctypeIteration .eq. 2
    
    ! iprevlv contains the last assembled level.
    iprevlv = 0

    ! Loop over all levels, from the highest to the lowest,
    ! and set up the system matrices.
    do ilev = rsolver%ispacelevel,rsolver%p_rlssHierarchy%nlmin,-1
    
      ! Assemble the matrix
      call smva_assembleMatrix_primal (rsolver%p_Rmatrices(ilev),ilev,idofTime,&
          rsolver%p_roperatorAsmHier,&
          rprimalSol,isollevelSpace,rsolver%itimelevel,bfull,&
          rsolver%p_roptcBDCSpaceHierarchy%p_RoptcBDCspace(ilev),&
          rsolver%rtempData,iprevlv)
          
      iprevlv = ilev
    
    end do

    ! Initialise the structures of the associated linear subsolver
    rsolver%p_rlsshierarchy%p_rdebugFlags%sstringTag = &
        "primal_"//trim(sys_siL(idofTime,10))//"_"//trim(sys_siL(isollevelSpace,10))
    call lssh_initData (&
        rsolver%p_rlsshierarchy,rsolver%p_roptcBDCSpaceHierarchy,&
        rsolver%ispacelevel,ierror)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_initData_primallin (rsolver, ierror, idofTime, &
      rprimalSol, isollevelSpace, bfull)
  
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
  
  ! Whether or not to apply the full Newton operator.
  logical, intent(in) :: bfull
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_spaceSolverHierarchy), intent(inout) :: rsolver
!</inputoutput>

!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!</subroutine>

    ! local variables
    integer :: ilev,iprevlv
    
    ! iprevlv contains the last assembled level.
    iprevlv = 0

    ! Loop over all levels, from the highest to the lowest,
    ! and set up the system matrices.
    do ilev = rsolver%ispacelevel,rsolver%p_rlssHierarchy%nlmin,-1
    
      ! Assemble the matrix
      call smva_assembleMatrix_primal (rsolver%p_Rmatrices(ilev),ilev,idofTime,&
          rsolver%p_roperatorAsmHier,&
          rprimalSol,isollevelSpace,rsolver%itimelevel,bfull,&
          rsolver%p_roptcBDCSpaceHierarchy%p_RoptcBDCspace(ilev),&
          rsolver%rtempData,iprevlv)
          
      iprevlv = ilev
    
    end do

    ! Initialise the structures of the associated linear subsolver
    rsolver%p_rlsshierarchy%p_rdebugFlags%sstringTag = &
        "primallin_"//trim(sys_siL(idofTime,10))//"_"//trim(sys_siL(isollevelSpace,10))
    call lssh_initData (&
        rsolver%p_rlsshierarchy,rsolver%p_roptcBDCSpaceHierarchy,&
        rsolver%ispacelevel,ierror)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_initData_dual (rsolver, ierror, idofTime, &
      rprimalSol, isollevelSpace)
  
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
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_spaceSolverHierarchy), intent(inout) :: rsolver
!</inputoutput>

!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!</subroutine>

    ! local variables
    integer :: ilev,iprevlv
    
    ! iprevlv contains the last assembled level.
    iprevlv = 0

    ! Loop over all levels, from the highest to the lowest,
    ! and set up the system matrices.
    do ilev = rsolver%ispacelevel,rsolver%p_rlssHierarchy%nlmin,-1
    
      ! Assemble the matrix
      call smva_assembleMatrix_dual (rsolver%p_Rmatrices(ilev),ilev,idofTime,&
          rsolver%p_roperatorAsmHier,rprimalSol,isollevelSpace,rsolver%itimelevel,&
          rsolver%p_roptcBDCSpaceHierarchy%p_RoptcBDCspace(ilev),&
          rsolver%rtempData,iprevlv)
          
      iprevlv = ilev
    
    end do

    ! Initialise the structures of the associated linear subsolver
    rsolver%p_rlsshierarchy%p_rdebugFlags%sstringTag = &
        "dual_"//trim(sys_siL(idofTime,10))//"_"//trim(sys_siL(isollevelSpace,10))
    call lssh_initData (&
        rsolver%p_rlsshierarchy,rsolver%p_roptcBDCSpaceHierarchy,&
        rsolver%ispacelevel,ierror)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_initData_dualLin (rsolver, ierror, idofTime, &
      rprimalSol, isollevelSpace, bfull)
  
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

  ! Whether or not to apply the full Newton operator.
  logical, intent(in) :: bfull
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_spaceSolverHierarchy), intent(inout) :: rsolver
!</inputoutput>

!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!</subroutine>

    ! Implementation identical to spaceslh_initData_dual.
    ! bfull does not have to be respected.

    ! local variables
    integer :: ilev,iprevlv
    
    ! iprevlv contains the last assembled level.
    iprevlv = 0

    ! Loop over all levels, from the highest to the lowest,
    ! and set up the system matrices.
    do ilev = rsolver%ispacelevel,rsolver%p_rlssHierarchy%nlmin,-1
    
      ! Assemble the matrix
      call smva_assembleMatrix_dual (rsolver%p_Rmatrices(ilev),ilev,idofTime,&
          rsolver%p_roperatorAsmHier,rprimalSol,isollevelSpace,rsolver%itimelevel,&
          rsolver%p_roptcBDCSpaceHierarchy%p_RoptcBDCspace(ilev),&
          rsolver%rtempData,iprevlv)
          
      iprevlv = ilev
    
    end do

    ! Initialise the structures of the associated linear subsolver
    rsolver%p_rlsshierarchy%p_rdebugFlags%sstringTag = &
        "duallin_"//trim(sys_siL(idofTime,10))//"_"//trim(sys_siL(isollevelSpace,10))
    call lssh_initData (&
        rsolver%p_rlsshierarchy,rsolver%p_roptcBDCSpaceHierarchy,&
        rsolver%ispacelevel,ierror)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_doneData (rsolver)
  
!<description>
  ! Cleanup of the data initalised in spaceslh_initData.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_spaceSolverHierarchy), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>

    ! Clean up the structures of the associated linear subsolver
    call lssh_doneData (rsolver%p_rlsshierarchy,rsolver%ispacelevel)

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
    call lssh_doneStructure (rsolver%p_rlsshierarchy,rsolver%ispacelevel)
   
   ! Release assembly data
    call smva_releaseTempData (rsolver%rtempData)
   
    ! deallocate temporary vectors
    do ilev=rsolver%ispacelevel,rsolver%p_rlssHierarchy%nlmin,-1
      call lsysbl_releaseMatrix (rsolver%p_Rmatrices(ilev))
    end do

    deallocate(rsolver%p_Rmatrices)
    
    ! Deallocate temporary data
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

  subroutine spaceslh_solve (rsolver,idofTime,csolgeneration,isollevelSpace,&
      rprimalSol,rdualSol,rcontrol,rprimalSolLin,rdualSolLin,rcontrolLin)
  
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

!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rd, p_rx
    integer :: ierror
    type(t_linsolNode), pointer :: p_rsolverNode
    real(DP), dimension(:), pointer :: p_Dd, p_Dx
    type(t_discreteBC), pointer :: p_rdiscreteBC
   
    ! At first, get a temp vector we can use for creating defects,
    p_rd => rsolver%p_rd
    
    ! Get the boundary conditions for the current timestep.
    call sbch_getDiscreteBC (&
        rsolver%p_roptcBDCSpaceHierarchy,rsolver%ispacelevel,p_rdiscreteBC)
        
    ! Assign the boundary conditions to the temp vector.
    ! Necessary for the linear solver.
    call lsysbl_assignDiscreteBC (p_rd,p_rdiscreteBC)
    
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
      rsolver%rnewtonParams%nnonlinearIterations = 0
      do while (.true.)
      
        ! -------------------------------------------------------------
        ! Initialise boundary conditions
        ! -------------------------------------------------------------
        ! No filter for the initial condition.
        if (idoftime .gt. 1) then
          call spaceslh_initData_bdc (rsolver, idofTime)      
        end if
      
        ! -------------------------------------------------------------
        ! Get the nonlinear defect
        ! -------------------------------------------------------------
        ! Compute the basic (unpreconditioned) search direction in rd.
        if (rsolver%rnewtonParams%nnonlinearIterations .eq. 0) then
          call smva_getDef_primal (p_rd,&
              rsolver%ispacelevel,rsolver%itimelevel,idofTime,&
              rsolver%p_roperatorAsmHier,rprimalSol,rcontrol,&
              rsolver%p_roptcBDCSpaceHierarchy%p_RoptcBDCspace(rsolver%ispacelevel),&
              rsolver%rtempData,csolgeneration)
        else
          ! Start from the currently available iterate and continue the nonlinear
          ! iteration
          call smva_getDef_primal (p_rd,&
              rsolver%ispacelevel,rsolver%itimelevel,idofTime,&
              rsolver%p_roperatorAsmHier,rprimalSol,rcontrol,&
              rsolver%p_roptcBDCSpaceHierarchy%p_RoptcBDCspace(rsolver%ispacelevel),&
              rsolver%rtempData,SPINITCOND_PREVITERATE)
        end if        
        
        rsolver%rnewtonParams%dresFinal = &
            lsysbl_vectorNorm(p_rd,rsolver%rnewtonParams%iresNorm)
            
        if (rsolver%rnewtonParams%nnonlinearIterations .eq. 0) then
          ! Remember the initial residual
          rsolver%rnewtonParams%dresInit = rsolver%rnewtonParams%dresFinal
        end if
        
        if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
          call output_line (&
              trim(sys_si(idofTime-1,8)) // &
              " " // trim(sys_si(rsolver%rnewtonParams%nnonlinearIterations,4)) // &
              " " // trim(sys_sdEL(rsolver%rnewtonParams%dresFinal,1)) )
        end if

        ! -------------------------------------------------------------
        ! Check for convergence
        ! -------------------------------------------------------------
        if (newtonit_checkConvergence (rsolver%rnewtonParams)) exit
        
        ! -------------------------------------------------------------
        ! Check for divergence
        ! -------------------------------------------------------------
        if (newtonit_checkDivergence (rsolver%rnewtonParams)) exit
        
        ! -------------------------------------------------------------
        ! Check other stopping criteria
        ! -------------------------------------------------------------
        if (newtonit_checkIterationStop (rsolver%rnewtonParams)) exit
        
        ! -------------------------------------------------------------
        ! Preconditioning with the Newton matrix
        ! -------------------------------------------------------------

        ! Assemble the matrices/boundary conditions on all levels.
        call spaceslh_initData_primal (rsolver, ierror, idofTime, &
            rprimalSol,rsolver%ispacelevel)
        
        ! Solve the system
        output_iautoOutputIndent = output_iautoOutputIndent + 2
        call lssh_precondDefect (rsolver%p_rlsshierarchy,rsolver%ispacelevel,p_rd)
        output_iautoOutputIndent = output_iautoOutputIndent - 2
        
        ! -------------------------------------------------------------
        ! Solver-Cleanup
        ! -------------------------------------------------------------
        call spaceslh_doneData (rsolver)
        
        ! -------------------------------------------------------------
        ! Update of the solution
        ! -------------------------------------------------------------

        ! Update the control according to the search direction:
        !
        !    u_n+1  =  u_n  +  g_n
        !
        ! or to any configured step-length control rule.
        call sptivec_getVectorFromPool(rprimalSol%p_rvectorAccess,idofTime,p_rx)
        
        ! DEBUG!!!
        call lsysbl_getbase_double (p_rx,p_Dx)
        
        call lsysbl_vectorLinearComb (p_rd,p_rx,rsolver%rnewtonParams%domega,1.0_DP)

        ! Implement boundary conditions
        call spaceslh_implementBDC (p_rx, SPACESLH_VEC_SOLUTION, &
            rsolver%ispacelevel, rsolver)
        
        call sptivec_commitVecInPool(rprimalSol%p_rvectorAccess,idofTime)
        
        ! -------------------------------------------------------------
        ! Clean up boundary conditions
        ! -------------------------------------------------------------
        call sbch_resetBCstructure (rsolver%p_roptcBDCSpaceHierarchy)

        ! -------------------------------------------------------------
        ! Proceed with the next iteration
        ! -------------------------------------------------------------
        ! Next iteration
        rsolver%rnewtonParams%nnonlinearIterations = &
            rsolver%rnewtonParams%nnonlinearIterations + 1
      
      end do
      
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
      call smva_getDef_dual (p_rd,&
          rsolver%ispacelevel,rsolver%itimelevel,idofTime,&
          rsolver%p_roperatorAsmHier,rprimalSol,rdualSol,&
          rsolver%p_roptcBDCSpaceHierarchy%p_RoptcBDCspace(rsolver%ispacelevel),&
          rsolver%rtempData,csolgeneration)
      
      rsolver%rnewtonParams%dresFinal = &
          lsysbl_vectorNorm(p_rd,rsolver%rnewtonParams%iresNorm)
      
      ! Remember the initial residual
      rsolver%rnewtonParams%dresInit = rsolver%rnewtonParams%dresFinal

      ! Print the initial residual
      if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
        call output_line (&
            trim(sys_si(idofTime-1,8)) // &
            " " // trim(sys_si(0,4)) // &
            " " // trim(sys_sdEL(rsolver%rnewtonParams%dresFinal,1)) )
      end if

      ! -------------------------------------------------------------
      ! Assemble the matrices/boundary conditions on all levels.
      ! -------------------------------------------------------------
      call spaceslh_initData_dual (rsolver, ierror, idofTime, &
          rprimalSol,rsolver%ispacelevel)

      ! -------------------------------------------------------------
      ! Call the linear solver, it does the job for us.
      ! -------------------------------------------------------------
      call lssh_precondDefect (rsolver%p_rlsshierarchy,rsolver%ispacelevel,p_rd,p_rsolverNode)
      
      ! -------------------------------------------------------------
      ! Solver-Cleanup
      ! -------------------------------------------------------------
      call spaceslh_doneData (rsolver)

      ! Get the final residual
      rsolver%rnewtonParams%dresFinal = p_rsolverNode%dfinalDefect
    
      ! Print the final residual
      if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
        call output_line (&
            trim(sys_si(idofTime-1,8)) // &
            " " // trim(sys_si(1,4)) // &
            " " // trim(sys_sdEL(rsolver%rnewtonParams%dresFinal,1)) )
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
    
      ! -------------------------------------------------------------
      ! Clean up boundary conditions
      ! -------------------------------------------------------------
      call sbch_resetBCstructure (rsolver%p_roptcBDCSpaceHierarchy)

    ! ---------------------------------------------------------------
    ! Linearised primal equation. Linear loop.
    ! ---------------------------------------------------------------
    case (OPTP_PRIMALLIN,OPTP_PRIMALLIN_SIMPLE)
      ! We apply one defect correction, so to pass a defect to the 
      ! linear solver.

      ! -------------------------------------------------------------
      ! Initialise boundary conditions
      ! -------------------------------------------------------------
      call spaceslh_initData_bdc (rsolver, idofTime)      

      ! -------------------------------------------------------------
      ! Create a defect
      ! -------------------------------------------------------------
      call smva_getDef_primalLin (p_rd,&
          rsolver%ispacelevel,rsolver%itimelevel,idofTime,&
          rsolver%p_roperatorAsmHier,rprimalSol,rcontrol,rprimalSolLin,&
          rcontrolLin,rsolver%coptype .eq. OPTP_PRIMALLIN,&
          rsolver%p_roptcBDCSpaceHierarchy%p_RoptcBDCspace(rsolver%ispacelevel),&
          rsolver%rtempData,csolgeneration)
          
      rsolver%rnewtonParams%dresFinal = &
          lsysbl_vectorNorm(p_rd,rsolver%rnewtonParams%iresNorm)
      
      ! Remember the initial residual
      rsolver%rnewtonParams%dresInit = rsolver%rnewtonParams%dresFinal

      ! Print the initial residual
      if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
        call output_line (&
            trim(sys_si(idofTime-1,8)) // &
            " " // trim(sys_si(0,4)) // &
            " " // trim(sys_sdEL(rsolver%rnewtonParams%dresFinal,1)) )
      end if

      ! -------------------------------------------------------------
      ! Assemble the matrices/boundary conditions on all levels.
      ! -------------------------------------------------------------
      call spaceslh_initData_primalLin (rsolver, ierror, idofTime, &
          rprimalSol,rsolver%ispacelevel,rsolver%coptype .eq. OPTP_PRIMALLIN)

      ! -------------------------------------------------------------
      ! Call the linear solver, it does the job for us.
      ! -------------------------------------------------------------
      call lssh_precondDefect (rsolver%p_rlsshierarchy,rsolver%ispacelevel,p_rd,p_rsolverNode)
      
      ! -------------------------------------------------------------
      ! Solver-Cleanup
      ! -------------------------------------------------------------
      call spaceslh_doneData (rsolver)

      ! Get the final residual
      rsolver%rnewtonParams%dresFinal = p_rsolverNode%dfinalDefect
    
      ! Print the final residual
      if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
        call output_line (&
            trim(sys_si(idofTime-1,8)) // &
            " " // trim(sys_si(1,4)) // &
            " " // trim(sys_sdEL(rsolver%rnewtonParams%dresFinal,1)) )
      end if

      ! -------------------------------------------------------------
      ! Update the control according to the defect
      !
      !    u_new  =  u_old  +  g_n
      ! -------------------------------------------------------------

      call sptivec_getVectorFromPool(rprimalSolLin%p_rvectorAccess,idofTime,p_rx)
      
      ! DEBUG!!!
      call lsysbl_getbase_double (p_rx,p_Dx)
      
      call lsysbl_vectorLinearComb (p_rd,p_rx,1.0_DP,1.0_DP)

      ! Implement boundary conditions.
      ! Use the defect filter since we are in the defect space.
      call spaceslh_implementBDC (p_rx, SPACESLH_VEC_DEFECT, &
            rsolver%ispacelevel, rsolver)

      call sptivec_commitVecInPool(rprimalSolLin%p_rvectorAccess,idofTime)

      ! -------------------------------------------------------------
      ! Clean up boundary conditions
      ! -------------------------------------------------------------
      call sbch_resetBCstructure (rsolver%p_roptcBDCSpaceHierarchy)

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
      call smva_getDef_dualLin (p_rd,&
          rsolver%ispacelevel,rsolver%itimelevel,idofTime,&
          rsolver%p_roperatorAsmHier,rprimalSol,rdualSol,rprimalSolLin,&
          rdualSolLin,rsolver%coptype .eq. OPTP_DUALLIN,&
          rsolver%p_roptcBDCSpaceHierarchy%p_RoptcBDCspace(rsolver%ispacelevel),&
          rsolver%rtempData,csolgeneration)
          
      ! Cleanup
      call spaceslh_doneData (rsolver)

      rsolver%rnewtonParams%dresFinal = &
          lsysbl_vectorNorm(p_rd,rsolver%rnewtonParams%iresNorm)
      
      ! Remember the initial residual
      rsolver%rnewtonParams%dresInit = rsolver%rnewtonParams%dresFinal

      ! Print the initial residual
      if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
        call output_line (&
            trim(sys_si(idofTime-1,8)) // &
            " " // trim(sys_si(0,4)) // &
            " " // trim(sys_sdEL(rsolver%rnewtonParams%dresFinal,1)) )
      end if

      ! -------------------------------------------------------------
      ! Assemble the matrices/boundary conditions on all levels.
      ! -------------------------------------------------------------
      call spaceslh_initData_dualLin (rsolver, ierror, idofTime, &
          rprimalSol,rsolver%ispacelevel,rsolver%coptype .eq. OPTP_DUALLIN)

      ! -------------------------------------------------------------
      ! Call the linear solver, it does the job for us.
      ! -------------------------------------------------------------
      call lssh_precondDefect (rsolver%p_rlsshierarchy,rsolver%ispacelevel,p_rd,p_rsolverNode)
      
      ! -------------------------------------------------------------
      ! Solver-Cleanup
      ! -------------------------------------------------------------
      call spaceslh_doneData (rsolver)

      ! Get the final residual
      rsolver%rnewtonParams%dresFinal = p_rsolverNode%dfinalDefect
    
      ! Print the final residual
      if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
        call output_line (&
            trim(sys_si(idofTime-1,8)) // &
            " " // trim(sys_si(1,4)) // &
            " " // trim(sys_sdEL(rsolver%rnewtonParams%dresFinal,1)) )
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
      call spaceslh_implementBDC (p_rx, SPACESLH_VEC_DEFECT, &
            rsolver%ispacelevel, rsolver)
      
      call sptivec_commitVecInPool(rdualSolLin%p_rvectorAccess,idofTime)

      ! -------------------------------------------------------------
      ! Clean up boundary conditions
      ! -------------------------------------------------------------
      call sbch_resetBCstructure (rsolver%p_roptcBDCSpaceHierarchy)

    end select

  end subroutine


end module
