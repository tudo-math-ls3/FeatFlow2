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

    ! Parameters of the OptFlow solver.
    ! The global space-time hierarchy in this structure defines
    ! the hierarchy, this solver is build upon.
    type(t_settings_optflow), pointer :: p_rsettingsSolver => null()

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

  subroutine spaceslh_init (rsolver,rsettingsSolver,copType,&
      rlssHierarchy,roptcBDCSpaceHierarchy,ssection,rparamList)
  
!<description>
  ! Initialises the solver parameters according to a parameter list.
  ! For each level in rlssHierarchy, a space solver will be created.
!</description>
  
!<input>
  ! Parameters of the OptFlow solver
  type(t_settings_optflow), intent(in), target :: rsettingsSolver
  
  ! Type of equation to be solved here. This can be 
  ! OPTP_PRIMAL, OPTP_DUAL, OPTP_PRIMALLIN or OPTP_DUALLIN,
  ! depending on which equation to solve.
  integer, intent(in) :: copType

  ! Hierarchy of linear solvers in space used for solving
  ! auxiliary linear subproblems. The solver on level ispacelevel
  ! will be used to solve linear subproblems.
  type(t_linsolHierarchySpace), target :: rlssHierarchy

  ! Boundary condition hierarchy for all space levels, primal and dual space
  type(t_optcBDCSpaceHierarchy), target :: roptcBDCSpaceHierarchy

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
    rsolver%p_rsettingsSolver => rsettingsSolver
    rsolver%p_rlssHierarchy => rlssHierarchy
    rsolver%p_roptcBDCSpaceHierarchy => roptcBDCSpaceHierarchy
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

    if (rsolver%p_rlssHierarchy%nlmin .gt. ispaceLevel) then
      call output_line ("Invalid space level.", &
          OU_CLASS_ERROR,OU_MODE_STD,"spaceslh_initStructure")
      call sys_halt()
    end if
    
    ! Remember the level and the hierarchy structures. 
    rsolver%ispacelevel = ispacelevel
    rsolver%itimelevel = itimelevel
    rsolver%p_roperatorAsmHier => roperatorAsmHier

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
          
      ! Provide the matrix to the linear solver
      call lssh_setMatrix(rsolver%p_rlsshierarchy,ilev,rsolver%p_Rmatrices(ilev))
      
      if (ilev .eq. rsolver%ispaceLevel) then
      
        ! On the highest space level, allocate some temp vectors on the topmost level

        allocate (rsolver%p_rd)
        call lsysbl_createVectorBlock (p_rdiscretisation,rsolver%p_rd)
            
        select case (rsolver%copType)
        case (OPTP_PRIMAL,OPTP_PRIMALLIN)
          call smva_allocTempData (rsolver%rtempData,&
              roperatorAsmHier%ranalyticData%p_rphysics,&
              roperatorAsmHier%ranalyticData%p_rsettingsOptControl,&
              roperatorAsmHier%ranalyticData%p_rsettingsSpaceDiscr,&
              ilev,roperatorAsmHier%p_rfeHierarchyPrimal,roperatorAsmHier)

        case (OPTP_DUAL,OPTP_DUALLIN)
          call smva_allocTempData (rsolver%rtempData,&
              roperatorAsmHier%ranalyticData%p_rphysics,&
              roperatorAsmHier%ranalyticData%p_rsettingsOptControl,&
              roperatorAsmHier%ranalyticData%p_rsettingsSpaceDiscr,&
              ilev,roperatorAsmHier%p_rfeHierarchyDual,roperatorAsmHier)
        end select

        ! Initialise the structures of the associated linear subsolver
        call lssh_initStructure (rsolver%p_rlsshierarchy,ilev,ierror)

      end if
    end do
    
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
    integer :: ilev
    
    ! A combination of LSS_SLFLAGS_xxxx constants defining the way,
    ! the solver performs filtering.
    integer(I32) :: cflags
    
    ! Loop over all levels, from the highest to the lowest,
    ! and set up the system matrices.
    do ilev = rsolver%ispacelevel,rsolver%p_rlssHierarchy%nlmin,-1
    
      call smva_assembleMatrix_primal (rsolver%p_Rmatrices(ilev),ilev,idofTime,&
          rsolver%p_roperatorAsmHier,&
          rprimalSol,isollevelSpace,rsolver%itimelevel,&
          rsolver%rnewtonParams%ctypeIteration .eq. 2,rsolver%rtempData)
    
    end do

    ! Initialise the structures of the associated linear subsolver
    call lssh_initData (&
        rsolver%p_rlsshierarchy,rsolver%ispacelevel,cflags,ierror)

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

  subroutine spaceslh_solve (rsolver,idofTime,isollevelSpace,&
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
   
    ! Which type of operator is to be solved?
    select case (rsolver%coptype)

    ! ---------------------------------------------------------------
    ! Primal equation. Nonlinear loop.
    ! ---------------------------------------------------------------
    case (OPTP_PRIMAL)

      ! This is the most complicated part. On level ispacelevel, we have
      ! to apply a nonlinear loop.
      !
      ! At first, get a temp vector we can use for creating the
      ! nonlinear defect.
      p_rd => rsolver%p_rd
      
      ! Apply the Newton iteration
      rsolver%rnewtonParams%nnonlinearIterations = 0
      do while (.true.)
      
        ! Compute the basic (unpreconditioned) search direction in rd.
        call smva_getDef_primal (p_rd,&
            rsolver%ispacelevel,rsolver%itimelevel,idofTime,&
            rsolver%p_roperatorAsmHier,&
            rprimalSol,rcontrol,rsolver%rtempData)

        if (rsolver%rnewtonParams%nnonlinearIterations .eq. 1) then
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

        ! Assemble the matrices on all levels.
        call spaceslh_initData_primal (rsolver, ierror, idofTime, &
            rprimalSol,rsolver%ispacelevel)
        
        ! Solve the system
        call lssh_precondDefect (rsolver%p_rlsshierarchy,rsolver%ispacelevel,p_rd)
        
        ! Cleanup
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
        call lsysbl_vectorLinearComb (p_rd,p_rx,rsolver%rnewtonParams%domega,1.0_DP)
        call sptivec_commitVecInPool(rprimalSol%p_rvectorAccess,idofTime)
        
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
      ! Call the linear solver, it does the job for us.
      !call lssh_precondDefect (rsolver%p_rlsshierarchy,ilevel,rd)
    
    ! ---------------------------------------------------------------
    ! Linearised primal equation. Linear loop.
    ! ---------------------------------------------------------------
    case (OPTP_PRIMALLIN,OPTP_PRIMALLIN_SIMPLE)
      ! Call the linear solver, it does the job for us.
      !call lssh_precondDefect (rsolver%p_rlsshierarchy,ilevel,rd)

    ! ---------------------------------------------------------------
    ! Linearised dual equation. Linear loop.
    ! ---------------------------------------------------------------
    case (OPTP_DUALLIN)
      ! Call the linear solver, it does the job for us.
      !call lssh_precondDefect (rsolver%p_rlsshierarchy,ilevel,rd)

    end select

  end subroutine


end module
