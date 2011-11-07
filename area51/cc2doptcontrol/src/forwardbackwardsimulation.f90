!##############################################################################
!# ****************************************************************************
!# <name> forwardbackwardsimulation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module implements a pure forward or backward simulation
!# acting on a primal/dual solution vector.
!#
!# The following routines can be found here:
!#
!# 1.) fbsim_init
!#     -> Initialise a forward-backward solver.
!#
!# 2.) fbsim_done
!#     -> Release a forward-backward solver
!#
!# 3.) fbsim_setMatrix
!#     -> Assign a space-time matrix to a forward-backward solver
!#
!# 4.) fbsim_simulate
!#      -> Performs a forward or backward simulation
!# </purpose>
!##############################################################################

module forwardbackwardsimulation

  use fsystem
  use linearsolver
  
  use basicstructures
  use spacepreconditioner
  use spacetimelinearsystem
  use spacepreconditioner
  
  use paramlist

  implicit none

!<types>

!<typeblock>

  ! This structure controls the Newton iteration -- i.e. the preconditioning
  ! with the Frechet derivative of the Navier--Stokes equation, which
  ! can lead to quadratic covergence of the nonlinear solver.
  ! As Newton works only in the basin of attraction of the solution,
  ! the parameters in this structure allow to define a switching criterion
  ! when to use Newton. In the first couple of iterations, defect correction
  ! is used, while the iteration switches to Newton if the residuum is small
  ! enough.
  ! This block is used if CCPREC_NEWTONDYNAMIC is used as preconditioner.
  type t_ccDynamicNewtonControl
  
    ! Minimum number of fix point iterations before switching to
    ! preconditioning with the Newton matrix. (IFIXMIN)

    integer :: nminFixPointIterations = 0

    ! Maximum number of fix point iterations before switching to
    ! preconditioning with the Newton matrix. (IFIXMAX)

    integer :: nmaxFixPointIterations = 999

    ! Norm of absolute residuum before applying Newton. 
    ! Newton is only applied
    ! if   ||absolute residuum|| < depsAbsNewton
    ! and  ||relative residuum|| < depsRelNewton.
    ! Otherwise, the usual fix point iteration is used.
    ! Stamndard value = 1E-5.

    real(DP) :: depsAbsNewton = 1.0E-5_DP

    ! Norm of relative residuum before applying Newton. 
    ! Newton is only applied
    ! if   ||absolute residuum|| < depsAbsNewton
    ! and  ||relative residuum|| < depsRelNewton.
    ! Otherwise, the usual fix point iteration is used.
    ! Standard value = 1E99 -> The absolute residuum counts.

    real(DP) :: depsRelNewton = 1.0E99_DP
  
    ! Whether to use the inexact Newton iteration or not.
    ! The inexact Newton controls the stopping criterion of the linear
    ! solver according to the nonlinear residual.
    
    integer :: cinexactNewton = 1

    ! Stopping criterion for the linear solver in the inexact Newton iteration.
    ! Controls the minimum number of digits to gain in the linear solver
    ! per Newton iteration. Only used if cinexactNewton = 1.
    
    real(dp) :: dinexactNewtonEpsRel = 1.0E-2_DP

    ! Exponent to control the stopping criterion for the linear solver in
    ! an inexact Newton iteration. =2 result in quadratic convergence,
    ! =1.5 in superlinear convergence. Only used if cinexactNewton = 1.
    
    real(dp) :: dinexactNewtonExponent = 2.0_DP
  
  end type
  
!</typeblock>

!<typeblock>

  ! This configuration block configures all parameters that are needed
  ! by the callback routines to perform the nonlinear iteration.
  ! On start of the program, a structure of this type is initialised.
  ! The entries of this structure are saved to the collection structure.
  ! When a callback routine is called, the structure is rebuild from
  ! the collection. When he nonlinear iteration is finished, the
  ! parameters are removed from the colletion again.
  type t_fbsimNonlinearIteration
  
    ! Norm of initial residuum for each subvector of the given block vector.
    ! Only valid, if the fcb_resNormCheck callback procedure is not specified
    ! in the call to the solver, otherwise undefined!
    real(DP), dimension(2)  :: DinitialDefect

    ! Norm of final residuum for each subvector of the given block vector.
    ! Only valid, if the fcb_resNormCheck callback procedure is not specified
    ! in the call to the solver, otherwise undefined!
    real(DP), dimension(2)  :: DfinalDefect

    ! Norm of initial residuum for the complete solution vector.
    ! Only valid, if the fcb_resNormCheck callback procedure is not specified
    ! in the call to the solver, otherwise undefined!
    real(DP) :: dinitialDefectTotal = 0.0_DP

    ! Norm of final residuum for the complete solution vector.
    ! Only valid, if the fcb_resNormCheck callback procedure is not specified
    ! in the call to the solver, otherwise undefined!
    real(DP) :: dfinalDefectTotal = 0.0_DP

    ! Type of stopping criterion to use for standard convergence test. One of the
    ! NLSOL_STOP_xxxx constants.
    ! Note: This parameter is only evaluated in the stanard convergence test.
    ! If the caller of the nonlinear solver specifies a callback routine fcb_resNormCheck
    ! for checking the convergence, that callback routine must implement its own
    ! logic to handle relative and absolute convrgence criteria!
    integer :: istoppingCriterion = NLSOL_STOP_STANDARD

    ! Type of norm to use in the residual checking of the total vector
    ! (cf. linearalgebra.f90).
    ! =0: euclidian norm, =1: l1-norm, =2: l2-norm, =3: MAX-norm
    integer :: iresNormTotal = 2
    
    ! Minimum number of iterations top perform
    integer :: nminIterations = 1

    ! Maximum number of iterations top perform
    integer :: nmaxIterations = 50
    
    ! Output level of the nonlinear iteration
    integer :: ioutputLevel = 2
    
    ! A filter chain that is used for implementing boundary conditions into
    ! (linear and nonlinear) defect vectors.
    type(t_filterChain), dimension(:), pointer :: p_RfilterChain
    
    ! Auxiliary variable: Saves the initial defect in the nonlinear iteration
    real(DP), dimension(2) :: DresidualInit = 0.0_DP
    
    ! Auxiliary variable: Saves the last defect in the nonlinear iteration
    real(DP), dimension(2) :: DresidualOld = 0.0_DP

    ! Auxiliary variable: Norm of the relative change = norm of the 
    ! preconditioned residual in the nonlinear iteration
    real(DP), dimension(3) :: DresidualCorr = 0.0_DP
    
    ! Auxiliary variable: Convergence criteria of the nonlinear solver
    real(DP), dimension(5) :: DepsNL = 0.0_DP
    
    ! Auxiliary variable: Last calculated damping parameter
    real(DP) :: domegaNL = 0.0_DP
    
    ! Auxiliary variable: Convergence rate of linear solver (if this is
    ! applied as preconditioner).
    real(DP) :: drhoLinearSolver = 0.0_DP
    
    ! Type of the iteration.
    ! =0: undefined
    ! =1: linear solver.
    ! =2: Newton solver
    ! =3: adaptive Newton with parameters in radaptiveNewton
    integer :: ctypeIteration = 0
    
    ! Output mode. Used for printing messages.
    ! =OU_MODE_STD: Print messages to the terminal and probably to a log 
    ! file (if a log file is opened).
    integer(I32)               :: coutputmode = OU_MODE_STD

    ! Parameters for the adaptive Newton iteration.
    type(t_ccDynamicNewtonControl) :: radaptiveNewton
    
  end type

!</typeblock>

!<typeblock>
  ! This type is used to save the preconditioner configuration and parameters
  ! Here it's noted, if and whose matrices exist and/or must be assmebled 
  ! transposed to be compatible with the preconditioner and more. It's more or less
  ! a collection if different flags plus the matrices/structures of the
  ! preconditioner which is used in space.
  type t_fbsimPreconditioner
  
    ! Type of the preconditioner.
    ! =0: undefined
    ! =1: linear solver.
    integer :: ctypePreconditioner = 0
  
    ! Solution space, the preconditioner is applied to.
    ! CCSPACE_PRIMAL: BC for the primal space.
    ! CCSPACE_DUAL: BC for the dual space.
    ! CCSPACE_PRIMALDUAL: BC for primal and dual space.
    integer :: cspace = CCSPACE_PRIMAL

    ! Whether to use 'adaptive matrices', i.e. set up coarse grid matrices
    ! with the help of fine grid matrices. This is used for very special
    ! discretisations only (e.g. Q1~/Q0). =0: deactivate
    integer :: iadaptiveMatrices    = 0
    
    ! A configuration parameter for adaptive matrices.
    real(DP) :: dadMatThreshold     = 0.0_DP

    ! If the preconditioner is a linear solver:
    ! Type of solver.
    ! =0: Gauss elimination (UMFPACK)
    ! =1: Multigrid solver
    integer :: isolverType = 0
    
    ! If the preconditioner is the linear multigrid solver:
    ! Type of smoother.
    ! =0: general VANKA (slow, but independent of the discretisation and of the problem)
    ! =1: general VANKA; 'direct' method, bypassing the defect correction approach.
    !     (-> specialised variant of 0, but slightly faster)
    ! =2: Simple Jacobi-like VANKA, 2D Navier Stokes problem, general discretisation
    !     (i.e. automatically chooses the best suitable VANKA variant).
    ! =3: Simple Jacobi-like VANKA, 2D Navier Stokes problem, general discretisation
    !     (i.e. automatically chooses the best suitable VANKA variant).
    !     'direct' method, bypassing the defect correction approach.
    !     (-> specialised variant of 8, but faster)
    ! =4: Full VANKA, 2D Navier Stokes problem, general discretisation
    !     (i.e. automatically chooses the best suitable VANKA variant).
    ! =5: Full VANKA, 2D Navier Stokes problem, general discretisation
    !     (i.e. automatically chooses the best suitable VANKA variant).
    !     'direct' method, bypassing the defect correction approach.
    !     (-> specialised variant of 10, but faster)
    integer :: ismootherType = 3
    
    ! If the preconditioner is the linear multigrid solver:
    ! Type of coarse grid solver.    
    ! =0: Gauss elimination (UMFPACK)
    ! =1: Defect correction with diagonal VANKA preconditioning.
    ! =2: BiCGStab with diagonal VANKA preconditioning
    integer :: icoarseGridSolverType = 1
        
    ! This flag is set to .TRUE. if there are no Neumann boundary
    ! components. In that case, the pressure matrices of direct
    ! solvers must be changed.
    logical :: bneedPressureDiagonalBlock = .false.
    
    ! Set to TRUE if the preconditioner needs virtually transposed B matrices
    ! as D matrices on all levels except for the coarse mesh.
    logical :: bneedVirtTransposedD = .false.
    
    ! Set to TRUE if the preconditioner needs virtually transposed B matrices
    ! as D matrices on the coarse mesh.
    logical :: bneedVirtTransposedDonCoarse = .false.
    
    !<!-- Parameters / structures for Preconditioner: Linear solver, e.g. MG -->
    
    ! Minimum refinement level
    integer :: NLMIN = 0
    
    ! Maximum refinement level
    integer :: NLMAX = 0
    
    ! An array of t_problem_lvl structures, each corresponding
    ! to one level of the discretisation. 
    type(t_problem_lvl), dimension(:), pointer :: p_RlevelInfo => null()
    
    ! A solver node that accepts parameters for the linear solver    
    type(t_linsolNode), pointer :: p_rsolverNode => null()

    ! If the linear solver contains a multigrid subsolver, this is a reference
    ! to the coarse grid solver.
    type(t_linsolNode), pointer :: p_rcoarseGridSolverNode => null()

    ! A filter chain that is used for implementing boundary conditions or other
    ! things when invoking the linear solver.
    type(t_filterChain), dimension(3) :: RfilterChain

    ! Interlevel projection hierarchy.
    type(t_interlevelProjectionHier) :: rprjHierarchy
    
    ! An array of matrices used for preconditioning in space.
    ! These matrices are 6x6 matrices for both, primal and dual space.
    ! The actual preconditioner matrices for either the primal or the
    ! dual space can be found in p_RmatrixPrecond and are shares submatrices
    ! of these matrices.
    type(t_matrixBlock), dimension(:), pointer :: p_RmatrixPrecondFullSpace

    ! An array of matrices used for preconditioning in space,
    ! only for primal or dual space.
    type(t_matrixBlock), dimension(:), pointer :: p_RmatrixPrecond
    
    ! Discrete boundary conditions on all levels corresponding to the
    ! discretisation of p_RmatrixPrecond.
    type(t_discreteBC), dimension(:), pointer :: p_RdiscreteBC
    
    ! Discrete fictitious BC`s on all levels corresponding to the
    ! discretisation of p_RmatrixPrecond.
    type(t_discreteFBC), dimension(:), pointer :: p_RdiscreteFBC
    
    ! An array of 3x#levels temp vectors.
    type(t_vectorBlock), dimension(:,:), pointer :: p_RtempVec
    
    ! Whether there are Neumann boundary components in the BC`s.
    logical :: bhasNeumann
    
  end type

!</typeblock>

!<typeblock>

  ! Configuration block for the simulation solver.
  type t_simSolver
  
    ! Type of simulation solver
    ! =0: Navier-Stokes forward solver on the primal solution vector.
    !     This iteration does not use roseenSolution.
    !     The dual solution influences the forward simulation; to prevent
    !     influence of the dual solution, it should be set to zero.
    ! =1: Oseen forward solver on the primal solution vector.
    !     The vector roseenSolution specifies the evaluation point of the 
    !     nonlinearity.
    ! =2: Oseen backward solver on the dual solution vector.
    !     The vector roseenSolution specifies the evaluation point of the 
    !     nonlinearity in the primal space.
    integer :: csimtype

    ! Absolute level where the simulation is executed.
    integer :: ilevel
  
    ! A t_fbsimPreconditioner structure that defines the preconditioner
    type(t_fbsimPreconditioner) :: rpreconditioner

    ! Pointer to the underlying space-time matrix that defines the timesteping etc.
    type(t_ccoptSpaceTimeMatrix), pointer :: p_rmatrix

    ! Level-info structure of the level where the simulation is executed.
    type(t_problem_lvl), pointer :: p_rlevelInfo => null()
    
    ! Pointer to the problem structure
    type(t_problem), pointer :: p_rproblem

    ! Parameters for the nonlinear iteration during a pure forwad simulation.
    type(t_fbsimNonlinearIteration) :: rnonlinearIteration

  end type

!</typeblock>

!</types>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine fbsim_getNLsimParameters (rparamList,sname,rnonlinearIteration)
  
!<description>
  ! Initialises the parameters of the nonlinear solver based on the settings
  ! in the parameter list rparamList.
!</description>

!<input>
  ! Parameter list with the parameters configuring the nonlinear solver
  type(t_parlist), intent(in) :: rparamList

  ! Name of the section in the parameter list containing the parameters
  ! of the nonlinear solver.
  character(LEN=*), intent(in) :: sname
!</input>

!<output>
  ! Nonlinar iteration structure saving data for the callback routines.
  ! Is filled with data.
  type(t_fbsimNonlinearIteration), intent(out) :: rnonlinearIteration
!</output>

!</subroutine>

    ! local variables
    type(t_parlstSection), pointer :: p_rsection
    character(len=SYS_STRLEN) :: sstring,snewton

    call parlst_getvalue_int (rparamList, sname, &
                              'ioutputLevel', rnonlinearIteration%ioutputLevel, 0)
    
    ! Clear auxiliary variables for the nonlinear iteration
    rnonlinearIteration%DresidualInit = 0.0_DP
    rnonlinearIteration%DresidualOld  = 0.0_DP
    
    call parlst_querysection(rparamList, sname, p_rsection) 

    if (.not. associated(p_rsection)) then
      call output_line ('Cannot create nonlinear solver; no section '''//&
          trim(sname)//'''!', &
          OU_CLASS_ERROR,OU_MODE_STD,'fbsim_getNLsimParameters')
      call sys_halt()
    end if

    ! Get stopping criteria of the nonlinear iteration
    call parlst_getvalue_double (p_rsection, 'depsD', &
                                 rnonlinearIteration%DepsNL(1), 0.1_DP)

    call parlst_getvalue_double (p_rsection, 'depsDiv', &
                                 rnonlinearIteration%DepsNL(2), 0.1_DP)

    call parlst_getvalue_double (p_rsection, 'depsUR', &
                                 rnonlinearIteration%DepsNL(3), 0.1_DP)

    call parlst_getvalue_double (p_rsection, 'depsPR', &
                                 rnonlinearIteration%DepsNL(4), 0.1_DP)

    call parlst_getvalue_double (p_rsection, 'dDampingD', &
                                 rnonlinearIteration%DepsNL(5), 0.1_DP)

    call parlst_getvalue_int (p_rsection, 'nminIterations', &
        rnonlinearIteration%nminIterations, rnonlinearIteration%nminIterations)

    call parlst_getvalue_int (p_rsection, 'nmaxIterations', &
        rnonlinearIteration%nmaxIterations, rnonlinearIteration%nmaxIterations)

    call parlst_getvalue_int (p_rsection, 'ioutputLevel', &
        rnonlinearIteration%ioutputLevel, rnonlinearIteration%ioutputLevel)

    ! We write out the data of the nonlinear solver to the benchmark
    ! log file as well.
    rnonlinearIteration%coutputMode = OU_MODE_STD+OU_MODE_BENCHLOG

    ! Get information about the iteration.
      
    ! At first, ask the parameters in the INI/DAT file which type of 
    ! preconditioner is to be used. The data in the preconditioner structure
    ! is to be initialised appropriately!
    call parlst_getvalue_int (rparamList, sname, &
        'ctypeIteration', rnonlinearIteration%ctypeIteration, 1)

    ! We have even the extended, dynamic Newton as preconditioner.
    ! Put the parameters for the extended Newton from the DAT file
    ! into the Adaptive-Newton configuration block.
    
    call parlst_getvalue_string (rparamList, sname, &
                                'spreconditionerAdaptiveNewton', sstring, '')
    snewton = ''
    if (sstring .ne. '') read (sstring,*) snewton
    if (snewton .ne. '') then
      ! Initialise the parameters of the adaptive Newton
      call parlst_getvalue_int (rparamList, snewton, &
          'nminFixPointIterations', rnonlinearIteration% &
          radaptiveNewton%nminFixPointIterations, 0)

      call parlst_getvalue_int (rparamList, snewton, &
          'nmaxFixPointIterations', rnonlinearIteration% &
          radaptiveNewton%nmaxFixPointIterations, 999)

      call parlst_getvalue_double (rparamList, snewton, &
          'depsAbsNewton', rnonlinearIteration% &
          radaptiveNewton%depsAbsNewton, 1E-5_DP)

      call parlst_getvalue_double (rparamList, snewton, &
          'depsRelNewton', rnonlinearIteration% &
          radaptiveNewton%depsRelNewton, 1E99_DP)

      call parlst_getvalue_int (rparamList, snewton, &
          'cinexactNewton', rnonlinearIteration% &
          radaptiveNewton%cinexactNewton, 1)

      call parlst_getvalue_double (rparamList, snewton, &
          'dinexactNewtonEpsRel', rnonlinearIteration% &
          radaptiveNewton%dinexactNewtonEpsRel, 1.0E-2_DP)

      call parlst_getvalue_double (rparamList, snewton, &
          'dinexactNewtonExponent', rnonlinearIteration% &
          radaptiveNewton%dinexactNewtonExponent, 2.0_DP)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fbsim_init (rproblem, nlmin, nlmax, csimtype, rsimsolver)
  
!<description>
  ! Initialises a rsimsolver structure based on the general problem structure
  ! for performing forward/backward sweeps.
!</description>
  
!<input>
  ! The general problem structure.
  type(t_problem), intent(inout), target :: rproblem

  ! Minimum level allowed to be used by the preconditioner.
  integer, intent(in) :: nlmin

  ! Absolute level where the simulation is executed.
  integer, intent(in) :: nlmax

  ! Type of simulation solver.
  ! =0: Navier-Stokes forward solver on the primal solution vector.
  !     The dual solution influences the simulation.
  ! =1: Oseen forward solver on the primal solution vector.
  ! =2: Oseen backward solver simulation on the dual solution vector.
  integer, intent(in) :: csimtype
!</input>

!<output>
  ! A t_simSolver node with parameters for the iteration.
  type(t_simSolver), intent(out) :: rsimsolver 
!</output>

!</subroutine>

    integer :: cspace

    ! Take the information from the problem structure.
    rsimsolver%ilevel = nlmax
    rsimsolver%p_rlevelInfo => rproblem%RlevelInfo(nlmax)
    rsimsolver%p_rproblem => rproblem
    
    rsimsolver%csimtype = csimtype
    
    select case (csimtype)
    case (0)
      cspace = CCSPACE_PRIMAL
    case (1)
      cspace = CCSPACE_PRIMAL
    case (2)
      cspace = CCSPACE_DUAL
    end select
    
    ! Initialise the preconditioner
    call fbsim_initPreconditioner (rproblem, rproblem%rparamlist, "CC-LINEARSOLVER", &
        nlmin, nlmax, cspace, rsimsolver%rpreconditioner)
        
    ! Get the parameters of the nonlinear solver in case we need them
    call fbsim_getNLsimParameters (rproblem%rparamList,"CC-NONLINEAR",&
        rsimsolver%rnonlinearIteration)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fbsim_done (rsimsolver)
  
!<description>
  ! Cleans up a rsimsolver structure.
!</description>
  
!<inputoutput>
  ! A t_simSolver node with parameters for the iteration.
  type(t_simSolver), intent(inout) :: rsimsolver 
!</inputoutput>

!</subroutine>

    call fbsim_donePreconditioner (rsimsolver%rpreconditioner)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fbsim_setMatrix (rsimsolver,rmatrix)
  
!<description>
  ! Assigns a space-time matrix to the simulator. This matrix defines
  ! the forward and backward problem.
!</description>
  
!<input>
  ! Space-time matrix structure to use for the iteration.
  type(t_ccoptSpaceTimeMatrix), intent(in), target :: rmatrix
!</input>

!<inputoutput>
  ! A t_simSolver node with parameters for the iteration.
  type(t_simSolver), intent(inout) :: rsimsolver 
!</inputoutput>

!</subroutine>

    ! Take the information from the problem structure.
    rsimsolver%p_rmatrix => rmatrix

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fbsim_initPreconditioner (rproblem, rparamlist, ssection, &
      nlmin, nlmax, cspace, rpreconditioner)
  
!<description>
  ! Initialises the preconditioner for each timestep based on the parameters 
  ! in rparamlist.
!</description>
  
!<input>
  ! The general problem structure.
  type(t_problem), intent(inout), target :: rproblem

  ! Parameter list with the solver parameters.
  type(t_parlist), intent(in) :: rparamlist
  
  ! Name of the section configuring the linear solver.
  character(len=*), intent(in) :: ssection
  
  ! Minimum available refinement level.
  integer, intent(in) :: nlmin
  
  ! Maximum refinement level; corresponds to the level where to do the 
  ! preconditioning.
  integer, intent(in) :: nlmax

  ! Solution space, the preconditioner is applied to.
  ! CCSPACE_PRIMAL: BC for the primal space.
  ! CCSPACE_DUAL: BC for the dual space.
  ! CCSPACE_PRIMALDUAL: BC for primal and dual space.
  integer, intent(in) :: cspace
!</input>

!<inputoutput>
  ! Configuration block of the preconditioner.
  type(t_fbsimPreconditioner), intent(out) :: rpreconditioner 
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_parlstSection), pointer :: p_rsection
    integer :: nlevels, ilev, nsm
    
    integer :: isolverType,ismootherType,icoarseGridSolverType
    character(LEN=SYS_STRLEN) :: sstring,ssolverSection,ssmootherSection
    character(LEN=SYS_STRLEN) :: scoarseGridSolverSection,spreconditionerSection
    type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo
    type(t_linsolNode), pointer :: p_rpreconditioner, p_rsmoother
    type(t_linsolNode), pointer :: p_rsolverNode
    type(t_nonlinearSpatialMatrix) :: rnonlinearSpatialMatrix
    
    ! -------------------------------------------------------------------------
    ! Part 1:
    ! Basic initialisation of all parameters and level-independent stuff.
    ! Reading of parameters from the DAT file.
    ! -------------------------------------------------------------------------

    ! Fetch level information where the preconditioner works.
    rpreconditioner%NLMIN = nlmin
    rpreconditioner%NLMAX = nlmax
    rpreconditioner%p_RlevelInfo => rproblem%RlevelInfo
    
    rpreconditioner%cspace = cspace
    
    ! Currently the only preconditioner available is a linear solver.
    rpreconditioner%ctypePreconditioner = 1

    ! Set up a filter that modifies the block vectors/matrix
    ! according to boundary conditions.
    !
    ! Initialise the first filter of the filter chain as boundary
    ! implementation filter for defect vectors:
    rpreconditioner%RfilterChain(1)%ifilterType = &
        FILTER_DISCBCDEFREAL

    ! The second filter filters for boundary conditions of fictitious boundary
    ! components
    rpreconditioner%RfilterChain(2)%ifilterType = &
        FILTER_DISCBCDEFFICT
    
    ! The last element is by default DO-NOTHING but can be changed to an L2_0-filter
    ! if necessary.
    rpreconditioner%RfilterChain(3)%ifilterType = FILTER_DONOTHING
    rpreconditioner%bneedPressureDiagonalBlock = .false.
    rpreconditioner%bhasNeumann = .true.
    
    ! Initialise a standard interlevel projection structure for every level
    call mlprj_initPrjHierarchy(rpreconditioner%rprjHierarchy,rpreconditioner%NLMIN,rpreconditioner%NLMAX)
    do ilev=rpreconditioner%NLMIN,rpreconditioner%NLMAX
    
      call mlprj_initPrjHierarchyLevel(rpreconditioner%rprjHierarchy,ilev,&
          rpreconditioner%p_RlevelInfo(ilev)%rdiscretisation)
    
      ! Initialise the projection structure with data from the INI/DAT
      ! files. This allows to configure prolongation/restriction.
      if (ilev .gt. rpreconditioner%NLMIN) then
        call fbsim_getProlRest (&
          rpreconditioner%rprjHierarchy%p_Rprojection(ilev-rpreconditioner%NLMIN+1), &
          rparamList, 'CC-PROLREST')
      end if
          
    end do
    call mlprj_commitPrjHierarchy(rpreconditioner%rprjHierarchy)
    
    ! Prepare boundary condition structures
    allocate(rpreconditioner%p_RdiscreteBC(rpreconditioner%NLMIN:rpreconditioner%NLMAX))
    allocate(rpreconditioner%p_RdiscreteFBC(rpreconditioner%NLMIN:rpreconditioner%NLMAX))
    do ilev=rpreconditioner%NLMIN,rpreconditioner%NLMAX
      call bcasm_initDiscreteBC(rpreconditioner%p_RdiscreteBC(ilev))
      call bcasm_initDiscreteFBC(rpreconditioner%p_RdiscreteFBC(ilev))
    end do
    
    ! Check that there is a section called ssolverName - otherwise we
    ! cannot create anything!
    
    call parlst_querysection(rparamList, ssection, p_rsection) 
    
    if (.not. associated(p_rsection)) then
      call output_line ('Cannot create linear solver; no section '''//trim(ssection)//&
                        '''!', OU_CLASS_ERROR,OU_MODE_STD,'fbsim_initPreconditioner')
      call sys_halt()
    end if
    
    ! Get the parameters that configure the solver type
    
    call parlst_getvalue_int (p_rsection, 'isolverType', isolverType, 1)
    call parlst_getvalue_int (p_rsection, 'ismootherType', ismootherType, 3)
    call parlst_getvalue_int (p_rsection, 'icoarseGridSolverType', &
        icoarseGridSolverType, 1)
        
    rpreconditioner%isolverType = isolverType
    rpreconditioner%ismootherType = ismootherType
    rpreconditioner%icoarseGridSolverType = icoarseGridSolverType

    call parlst_getvalue_string (p_rsection, 'ssolverSection', sstring,'')
    read (sstring,*) ssolverSection
    call parlst_getvalue_string (p_rsection, 'ssmootherSection', sstring,'')
    read (sstring,*) ssmootherSection
    call parlst_getvalue_string (p_rsection, 'scoarseGridSolverSection', sstring,'')
    read (sstring,*) scoarseGridSolverSection
    
    ! Which type of solver do we have?
    
    select case (isolverType)
    
    case (0)
    
      ! This is the UMFPACK solver. Very easy to initialise. No parameters at all.
      call linsol_initUMFPACK4 (p_rsolverNode)
    
    case (1)
    
      ! Multigrid solver. This is a little bit harder.
      !
      ! In a first step, initialise the main solver node for all our levels.
      nlevels = rpreconditioner%NLMAX - rpreconditioner%NLMIN + 1
      
      call linsol_initMultigrid2 (p_rsolverNode,nlevels,&
          rpreconditioner%RfilterChain)
      
      ! Manually trim the coarse grid correction in Multigrid to multiply the 
      ! pressure equation with -1. This (un)symmetrises the operator and gives
      ! much better convergence rates.
      call cgcor_release(p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection)
      call cgcor_init(p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection,3)
      p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%p_DequationWeights(3) &
          = -1.0_DP

      ! Init standard solver parameters and extended multigrid parameters
      ! from the DAT file.
      call linsolinit_initParams (p_rsolverNode,rparamList,ssolverSection,&
          LINSOL_ALG_UNDEFINED)
      call linsolinit_initParams (p_rsolverNode,rparamList,ssolverSection,&
          LINSOL_ALG_MULTIGRID2)
      
      ! Ok, now we have to initialise all levels. First, we create a coarse
      ! grid solver and configure it.
      call linsol_getMultigrid2Level (p_rsolverNode,1,p_rlevelInfo)
      
      select case (icoarseGridSolverType)
      case (0)
        ! UMFPACK coarse grid solver. Easy.
        call linsol_initUMFPACK4 (p_rlevelInfo%p_rcoarseGridSolver)
        
      case (1)
        ! Defect correction with diagonal VANKA preconditioning.
        !
        ! Create VANKA and initialise it with the parameters from the DAT file.
        call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVST)
        
        call parlst_getvalue_string (rparamList, scoarseGridSolverSection, &
            'spreconditionerSection', sstring, '')
        read (sstring,*) spreconditionerSection
        call linsolinit_initParams (p_rpreconditioner,rparamList,&
            spreconditionerSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rpreconditioner,rparamList,&
            spreconditionerSection,p_rpreconditioner%calgorithm)
        
        ! Create the defect correction solver, attach VANKA as preconditioner.
        call linsol_initDefCorr (p_rlevelInfo%p_rcoarseGridSolver,p_rpreconditioner,&
            rpreconditioner%RfilterChain)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,rparamList,&
            scoarseGridSolverSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,rparamList,&
            scoarseGridSolverSection,p_rpreconditioner%calgorithm)
        
        ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
        rpreconditioner%bneedVirtTransposedDonCoarse = .true.

      case (2)
        ! Defect correction with full VANKA preconditioning.
        !
        ! Create VANKA and initialise it with the parameters from the DAT file.
        call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVST)
        
        call parlst_getvalue_string (rparamList, scoarseGridSolverSection, &
            'spreconditionerSection', sstring, '')
        read (sstring,*) spreconditionerSection
        call linsolinit_initParams (p_rpreconditioner,rparamList,&
            spreconditionerSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rpreconditioner,rparamList,&
            spreconditionerSection,p_rpreconditioner%calgorithm)
        
        ! Create the defect correction solver, attach VANKA as preconditioner.
        call linsol_initDefCorr (p_rlevelInfo%p_rcoarseGridSolver,p_rpreconditioner,&
            rpreconditioner%RfilterChain)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,rparamList,&
            scoarseGridSolverSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,rparamList,&
            scoarseGridSolverSection,p_rpreconditioner%calgorithm)
        
        ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
        rpreconditioner%bneedVirtTransposedDonCoarse = .true.

      case (3)
        ! BiCGStab with diagonal VANKA preconditioning.
        !
        ! Create VANKA and initialise it with the parameters from the DAT file.
        call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DNAVST)
        
        call parlst_getvalue_string (rparamList, scoarseGridSolverSection, &
           'spreconditionerSection', sstring, '')
        read (sstring,*) spreconditionerSection
        call linsolinit_initParams (p_rpreconditioner,rparamList,&
            spreconditionerSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rpreconditioner,rparamList,&
            spreconditionerSection,p_rpreconditioner%calgorithm)
        
        ! Create the defect correction solver, attach VANKA as preconditioner.
        call linsol_initBiCGStab (p_rlevelInfo%p_rcoarseGridSolver,p_rpreconditioner,&
            rpreconditioner%RfilterChain)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,rparamList,&
            scoarseGridSolverSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,rparamList,&
            scoarseGridSolverSection,p_rlevelInfo%p_rcoarseGridSolver%calgorithm)
        
        ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
        rpreconditioner%bneedVirtTransposedDonCoarse = .true.

      case (4)
        ! BiCGStab with full VANKA preconditioning.
        !
        ! Create VANKA and initialise it with the parameters from the DAT file.
        call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVST)
        
        call parlst_getvalue_string (rparamList, scoarseGridSolverSection, &
           'spreconditionerSection', sstring, '')
        read (sstring,*) spreconditionerSection
        call linsolinit_initParams (p_rpreconditioner,rparamList,&
            spreconditionerSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rpreconditioner,rparamList,&
            spreconditionerSection,p_rpreconditioner%calgorithm)
        
        ! Create the defect correction solver, attach VANKA as preconditioner.
        call linsol_initBiCGStab (p_rlevelInfo%p_rcoarseGridSolver,p_rpreconditioner,&
            rpreconditioner%RfilterChain)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,rparamList,&
            scoarseGridSolverSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,rparamList,&
            scoarseGridSolverSection,p_rlevelInfo%p_rcoarseGridSolver%calgorithm)

        ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
        rpreconditioner%bneedVirtTransposedDonCoarse = .true.

      case (5)
        ! BiCGStab with full VANKA preconditioning.
        !
        ! Create VANKA and initialise it with the parameters from the DAT file.
        call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_GENERAL)
        
        call parlst_getvalue_string (rparamList, scoarseGridSolverSection, &
           'spreconditionerSection', sstring, '')
        read (sstring,*) spreconditionerSection
        call linsolinit_initParams (p_rpreconditioner,rparamList,&
            spreconditionerSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rpreconditioner,rparamList,&
            spreconditionerSection,p_rpreconditioner%calgorithm)
        
        ! Create the defect correction solver, attach VANKA as preconditioner.
        call linsol_initBiCGStab (p_rlevelInfo%p_rcoarseGridSolver,p_rpreconditioner,&
            rpreconditioner%RfilterChain)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,rparamList,&
            scoarseGridSolverSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,rparamList,&
            scoarseGridSolverSection,p_rlevelInfo%p_rcoarseGridSolver%calgorithm)

      case (6)
        ! BiCGStab with diagonal VANKA preconditioning, new implementation
        ! for general elements
        !
        ! Create VANKA and initialise it with the parameters from the DAT file.
        call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DNAVST)
        
        call parlst_getvalue_string (rparamList, scoarseGridSolverSection, &
           'spreconditionerSection', sstring, '')
        read (sstring,*) spreconditionerSection
        call linsolinit_initParams (p_rpreconditioner,rparamList,&
            spreconditionerSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rpreconditioner,rparamList,&
            spreconditionerSection,p_rpreconditioner%calgorithm)
        
        ! Create the defect correction solver, attach VANKA as preconditioner.
        call linsol_initBiCGStab (p_rlevelInfo%p_rcoarseGridSolver,p_rpreconditioner,&
            rpreconditioner%RfilterChain)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,rparamList,&
            scoarseGridSolverSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,rparamList,&
            scoarseGridSolverSection,p_rlevelInfo%p_rcoarseGridSolver%calgorithm)

        ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
        rpreconditioner%bneedVirtTransposedDonCoarse = .true.

      case default
      
        call output_line ('Unknown coarse grid solver.', &
            OU_CLASS_ERROR,OU_MODE_STD,'cc_initLinearSolver')
        call sys_halt()
          
      end select
      
      ! Save the reference to the coarse grid solver.
      rpreconditioner%p_rcoarseGridsolverNode => p_rlevelInfo%p_rcoarseGridSolver
      
      ! Now after the coarse grid solver is done, we turn to the smoothers
      ! on all levels. Their initialisation is similar to the coarse grid
      ! solver. Note that we use the same smoother on all levels, for 
      ! presmoothing as well as for postsmoothing.
      
      do ilev = 2,nlevels

        ! Initialise the smoothers.
        select case (ismootherType)
        
        case (0:9)

          nullify(p_rsmoother)
        
          ! This is some kind of VANKA smoother. Initialise the correct one.
          select case (ismootherType)
          case (0)
            call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_GENERAL)
          case (1)
            call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_GENERALDIRECT)
          case (2)
            call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVST)

            ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
            rpreconditioner%bneedVirtTransposedD = .true.

          case (3)
            call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTDIRECT)

            ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
            rpreconditioner%bneedVirtTransposedD = .true.

          case (4)
            call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVST)

            ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
            rpreconditioner%bneedVirtTransposedD = .true.

          case (5)
            call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTDIRECT)

            ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
            rpreconditioner%bneedVirtTransposedD = .true.

          case (6)
            call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVST)
            call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,&
                rpreconditioner%RfilterChain)

            ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
            rpreconditioner%bneedVirtTransposedD = .true.

          case (7)
            call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DNAVST)
            call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,&
                rpreconditioner%RfilterChain)

            ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
            rpreconditioner%bneedVirtTransposedD = .true.

          case (8)
            call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DNAVST)
            call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,&
                rpreconditioner%RfilterChain)

            ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
            rpreconditioner%bneedVirtTransposedD = .true.

          case (9)
            call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DNAVST)

            ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
            rpreconditioner%bneedVirtTransposedD = .true.

          end select
          
          ! Initialise the parameters -- if there are any.
          call linsolinit_initParams (p_rsmoother,rparamList,&
              ssmootherSection,LINSOL_ALG_UNDEFINED)
          call linsolinit_initParams (p_rsmoother,rparamList,&
              ssmootherSection,p_rsmoother%calgorithm)
          
          ! Convert to a smoother with a defined number of smoothing steps.
          call parlst_getvalue_int (rparamList, ssmootherSection, &
                    'nsmoothingSteps', nsm, 4)
          call linsol_convertToSmoother (p_rsmoother,nsm)
          
          ! Put the smoother into the level info structure as presmoother
          ! and postsmoother
          call linsol_getMultigrid2Level (p_rsolverNode,ilev,p_rlevelInfo)
          p_rlevelInfo%p_rpresmoother => p_rsmoother
          p_rlevelInfo%p_rpostsmoother => p_rsmoother
          
          ! Set up the interlevel projection structure for the projection from/to
          ! the lower level.
          call linsol_initProjMultigrid2Level(p_rlevelInfo,&
              rpreconditioner%rprjHierarchy%p_Rprojection(ilev))
          
        case default
        
          call output_line ('Unknown smoother.', &
              OU_CLASS_ERROR,OU_MODE_STD,'cc_initLinearSolver')
          call sys_halt()
          
        end select
      
      end do

      ! Get information about adaptive matrix generation from INI/DAT files
      call parlst_getvalue_int (rparamList, 'CC-DISCRETISATION', &
          'iAdaptiveMatrix', rpreconditioner%iadaptiveMatrices, 0)
                                
      call parlst_getvalue_double(rparamList, 'CC-DISCRETISATION', &
          'dAdMatThreshold', rpreconditioner%dAdMatThreshold, 20.0_DP)

    case (2:3)
    
      ! VANKA smoother: 1 step defect correction with nmaxIterations steps VANKA.
      ! ismootherType defines the type of smoother to use.
      select case (ismootherType)
      
      case (0:9)

        nullify(p_rsmoother)
      
        ! This is some kind of VANKA smoother. Initialise the correct one.
        select case (ismootherType)
        case (0)
          call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_GENERAL)
        case (1)
          call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_GENERALDIRECT)
        case (2)
          call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVST)

          ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
          rpreconditioner%bneedVirtTransposedD = .true.

        case (3)
          call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTDIRECT)

          ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
          rpreconditioner%bneedVirtTransposedD = .true.

        case (4)
          call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVST)

          ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
          rpreconditioner%bneedVirtTransposedD = .true.

        case (5)
          call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTDIRECT)

          ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
          rpreconditioner%bneedVirtTransposedD = .true.

        case (6)
          call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVST)
          call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,&
              rpreconditioner%RfilterChain)

          ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
          rpreconditioner%bneedVirtTransposedD = .true.

        case (7)
          call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVST)
          call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,&
              rpreconditioner%RfilterChain)

          ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
          rpreconditioner%bneedVirtTransposedD = .true.

        case (8)
          call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVST)
          call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,&
              rpreconditioner%RfilterChain)

        case (9)
          call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVST)

        end select
        
        ! Initialise the parameters -- if there are any.
        call linsolinit_initParams (p_rsmoother,rparamList,&
            ssmootherSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rsmoother,rparamList,&
            ssmootherSection,p_rsmoother%calgorithm)
        
        ! Convert to a smoother with a defined number of smoothing steps.
        call parlst_getvalue_int (rparamList, ssmootherSection, &
                  'nsmoothingSteps', nsm, 4)
        call linsol_convertToSmoother (p_rsmoother,nsm)
        
      case default
      
        call output_line ('Unknown smoother.', &
            OU_CLASS_ERROR,OU_MODE_STD,'cc_initLinearSolver')
        call sys_halt()
        
      end select

      ! Init defect correction, 1-step... with that smoother
      if (isolverType .eq. 2) then
        call linsol_initDefCorr (p_rsolverNode,p_rsmoother,rpreconditioner%RfilterChain)
      else
        call linsol_initBiCGStab (p_rsolverNode,p_rsmoother,rpreconditioner%RfilterChain)
      end if
      call linsolinit_initParams (p_rsolverNode,rparamList,&
          ssolverSection,LINSOL_ALG_UNDEFINED)
      call linsolinit_initParams (p_rsolverNode,rparamList,&
          ssolverSection,p_rsolverNode%calgorithm)
    
    end select    

    ! Put the final solver node to the preconditioner structure.
    rpreconditioner%p_rsolverNode => p_rsolverNode
    
    ! -------------------------------------------------------------------------
    ! Part 2:
    ! Matrix allocation & basic setup of the preconditioner.
    ! -------------------------------------------------------------------------
    
    ! Prepare preconditioner matrices on each level, attach them top the
    ! solver and do a structure-initialisation of the solver.
    !
    ! Allocate the matrix array.
    allocate(rpreconditioner%p_RmatrixPrecondFullSpace(&
        rpreconditioner%NLMIN:rpreconditioner%NLMAX))

    ! Get a dummy structure for a full matrix.
    call cc_getFullMatrixDummy(rproblem%rphysicsPrimal,rnonlinearSpatialMatrix)
    
    ! Deactivate the offdiagonal submatrices, we don't need them.
    ! Saves some memory.
    call cc_disableSubmatrix (rnonlinearSpatialMatrix,1,2)
    call cc_disableSubmatrix (rnonlinearSpatialMatrix,2,1)
    
    ! Allocate memory on each level
    do ilev = rpreconditioner%NLMIN,rpreconditioner%NLMAX
      call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rproblem,&
          rpreconditioner%p_RlevelInfo(ilev)%rdiscretisation,&
          rpreconditioner%p_RlevelInfo(ilev)%rstaticInfo)
    
      call cc_assembleMatrix (CCMASM_ALLOCMEM,CCMASM_MTP_AUTOMATIC,&
          rpreconditioner%p_RmatrixPrecondFullSpace(ilev),&
          rnonlinearSpatialMatrix)
    end do

    ! Allocate the actual preconditioner matrices. These will be created
    ! later as submatrices of p_RmatrixPrecondFullSpace.
    allocate(rpreconditioner%p_RmatrixPrecond(&
        rpreconditioner%NLMIN:rpreconditioner%NLMAX))
        
    ! -------------------------------------------------------------------------
    ! Part 3:
    ! Allocation of temp vectors
    ! -------------------------------------------------------------------------
    
    allocate(rpreconditioner%p_RtempVec(3,rpreconditioner%NLMIN:rpreconditioner%NLMAX))
    
    select case (rpreconditioner%cspace)
    case (CCSPACE_PRIMAL,CCSPACE_DUAL)
      ! Create temp vectors in the size of the primal space
      do ilev = rpreconditioner%NLMIN,rpreconditioner%NLMAX
        call lsysbl_createVectorBlock(&
            rpreconditioner%p_RlevelInfo(ilev)%rdiscretisationPrimal,&
            rpreconditioner%p_RtempVec(1,ilev),.false.)
        call lsysbl_createVectorBlock(&
            rpreconditioner%p_RlevelInfo(ilev)%rdiscretisationPrimal,&
            rpreconditioner%p_RtempVec(2,ilev),.false.)
        call lsysbl_createVectorBlock(&
            rpreconditioner%p_RlevelInfo(ilev)%rdiscretisationPrimal,&
            rpreconditioner%p_RtempVec(3,ilev),.false.)
      end do

    case (CCSPACE_PRIMALDUAL)
      ! Create temp vectors in the size of the full space
      do ilev = rpreconditioner%NLMIN,rpreconditioner%NLMAX
        call lsysbl_createVectorBlock(&
            rpreconditioner%p_RlevelInfo(ilev)%rdiscretisation,&
            rpreconditioner%p_RtempVec(1,ilev),.false.)
        call lsysbl_createVectorBlock(&
            rpreconditioner%p_RlevelInfo(ilev)%rdiscretisation,&
            rpreconditioner%p_RtempVec(2,ilev),.false.)
        call lsysbl_createVectorBlock(&
            rpreconditioner%p_RlevelInfo(ilev)%rdiscretisation,&
            rpreconditioner%p_RtempVec(3,ilev),.false.)
      end do
    end select
        
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fbsim_getProlRest (rprojection, rparamList, sname)
  
!<description>
  ! Initialises an existing interlevel projection structure rprojection
  ! with parameters from the INI/DAT files. sname is the section in the
  ! parameter list containing parameters about prolongation restriction.
!</description>

!<input>
  ! Parameter list that contains the parameters from the INI/DAT file(s).
  type(t_parlist), intent(IN) :: rparamList
  
  ! Name of the section in the parameter list containing the parameters
  ! of the prolongation/restriction.
  character(LEN=*), intent(IN) :: sname
!</input>

!<output>
  ! An interlevel projection block structure containing an initial
  ! configuration of prolongation/restriction. The structure is modified
  ! according to the parameters in the INI/DAT file(s).
  type(t_interlevelProjectionBlock), intent(INOUT) :: rprojection
!</output>

!</subroutine>

    ! local variables
    type(t_parlstSection), pointer :: p_rsection
    integer :: i1
    real(DP) :: d1

    ! Check that there is a section called sname - otherwise we
    ! cannot create anything!
    
    call parlst_querysection(rparamList, sname, p_rsection) 

    if (.not. associated(p_rsection)) then
      ! We use the default configuration; stop here.
      return
    end if
    
    ! Now take a look which parameters appear in that section.

    ! Prolongation/restriction order for velocity components
    call parlst_getvalue_int (p_rsection,'iinterpolationOrderVel',i1,-1)
    
    if (i1 .ne. -1) then
      ! Initialise order of prolongation/restriction for velocity components
      rprojection%RscalarProjection(:,1:NDIM2D)%iprolongationOrder  = i1
      rprojection%RscalarProjection(:,1:NDIM2D)%irestrictionOrder   = i1
      rprojection%RscalarProjection(:,1:NDIM2D)%iinterpolationOrder = i1
    end if

    ! Prolongation/restriction order for pressure
    call parlst_getvalue_int (p_rsection,'iinterpolationOrderPress',i1,-1)
    
    if (i1 .ne. -1) then
      ! Initialise order of prolongation/restriction for pressure components
      rprojection%RscalarProjection(:,NDIM2D+1)%iprolongationOrder  = i1
      rprojection%RscalarProjection(:,NDIM2D+1)%irestrictionOrder   = i1
      rprojection%RscalarProjection(:,NDIM2D+1)%iinterpolationOrder = i1
    end if
    
    ! Prolongation/restriction variant for velocity components
    ! in case of Q1~ discretisation
    call parlst_getvalue_int (p_rsection,'iinterpolationVariantVel',i1,0)
    
    if (i1 .ne. -1) then
      rprojection%RscalarProjection(:,1:NDIM2D)%iprolVariant  = i1
      rprojection%RscalarProjection(:,1:NDIM2D)%irestVariant  = i1
    end if
    
    ! Aspect-ratio indicator in case of Q1~ discretisation
    ! with extended prolongation/restriction
    call parlst_getvalue_int (p_rsection,'iintARIndicatorEX3YVel',i1,1)
    
    if (i1 .ne. 1) then
      rprojection%RscalarProjection(:,1:NDIM2D)%iprolARIndicatorEX3Y  = i1
      rprojection%RscalarProjection(:,1:NDIM2D)%irestARIndicatorEX3Y  = i1
    end if

    ! Aspect-ratio bound for switching to constant prolongation/restriction
    ! in case of Q1~ discretisation with extended prolongation/restriction
    call parlst_getvalue_double (p_rsection,'dintARboundEX3YVel',d1,20.0_DP)
    
    if (d1 .ne. 20.0_DP) then
      rprojection%RscalarProjection(:,1:NDIM2D)%dprolARboundEX3Y  = d1
      rprojection%RscalarProjection(:,1:NDIM2D)%drestARboundEX3Y  = d1
    end if

  end subroutine
    
  ! ***************************************************************************

!<subroutine>

  subroutine fbsim_updateDiscreteBCprec (rproblem,rpreconditioner,dtime)
  
!<description>
  ! Updates the discrete boundary conditions in the preconditioner
!</description>
  
!<input>
  ! The general problem structure.
  type(t_problem), intent(inout), target :: rproblem

  ! Current simulation time.
  real(dp), intent(in) :: dtime
!</input>
  
!<inputoutput>
  ! Configuration block of the preconditioner.
  type(t_fbsimPreconditioner), intent(inout) :: rpreconditioner 
!</inputoutput>

!</subroutine>

    ! local variables
    logical :: bneumann
    integer :: ilev

    ! Initialise the collection for the assembly.
    call cc_initCollectForAssembly (rproblem,dtime,rproblem%rcollection)

    ! Clear the BC`s and reassemble on all levels.
    do ilev = rpreconditioner%NLMIN,rpreconditioner%NLMAX
      call bcasm_clearDiscreteBC(rpreconditioner%p_RdiscreteBC(ilev))
      call bcasm_clearDiscreteFBC(rpreconditioner%p_RdiscreteFBC(ilev))
      
      call cc_assembleBDconditions (rproblem,dtime,&
          rpreconditioner%p_RlevelInfo(ilev)%rdiscretisationPrimal,&
          rpreconditioner%cspace,rpreconditioner%p_RdiscreteBC(ilev),rproblem%rcollection,bneumann)
      call cc_assembleFBDconditions (rproblem,dtime,&
          rpreconditioner%p_RlevelInfo(ilev)%rdiscretisationPrimal,&
          rpreconditioner%cspace,rpreconditioner%p_RdiscreteFBC(ilev),rproblem%rcollection)
    end do

    ! Clean up the collection (as we are done with the assembly, that's it.
    call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)

    ! Do we have Neumann boundary?
    ! The Neumann flag on the maximum level decides upon that.
    ! This may actually change from level to level but we simplify here
    ! and use the filter only depending on the max. level.
    rpreconditioner%bhasNeumann = bneumann
    if (.not. bneumann) then
      ! Pure Dirichlet problem -- Neumann boundary for the pressure.
      ! Filter the pressure to avoid indefiniteness.
      rpreconditioner%RfilterChain(3)%ifilterType = FILTER_TOL20
      rpreconditioner%RfilterChain(3)%itoL20component = 3
      rpreconditioner%bneedPressureDiagonalBlock = .true.
    else
      rpreconditioner%RfilterChain(3)%ifilterType = FILTER_DONOTHING
      rpreconditioner%bneedPressureDiagonalBlock = .false.
    end if
    

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fbsim_donePreconditioner (rpreconditioner)
  
!<description>
  ! Cleans up the preconditioner.
!</description>
  
!<inputoutput>
  ! Configuration block of the preconditioner.
  type(t_fbsimPreconditioner), intent(inout) :: rpreconditioner 
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ilev
    
    ! Cancel if not initialised.
    if (rpreconditioner%ctypePreconditioner .eq. 0) return

    ! Release the linear solver
    call linsol_releaseSolver (rpreconditioner%p_rsolverNode)

    ! Clean up data about the projection
    call mlprj_releasePrjHierarchy(rpreconditioner%rprjHierarchy)
    
    ! Release matrices on each level
    do ilev = rpreconditioner%NLMAX,rpreconditioner%NLMIN,-1
      call lsysbl_releaseMatrix (&
          rpreconditioner%p_RmatrixPrecondFullSpace(ilev))
      call lsysbl_releaseMatrix (&
          rpreconditioner%p_RmatrixPrecond(ilev))
    end do
    deallocate(rpreconditioner%p_RmatrixPrecondFullSpace)
    deallocate(rpreconditioner%p_RmatrixPrecond)

    ! Release boundary condition structures
    do ilev = rpreconditioner%NLMAX,rpreconditioner%NLMIN,-1
      call bcasm_releaseDiscreteFBC(rpreconditioner%p_RdiscreteFBC(ilev))
      call bcasm_releaseDiscreteBC(rpreconditioner%p_RdiscreteBC(ilev))
    end do
    deallocate(rpreconditioner%p_RdiscreteFBC)
    deallocate(rpreconditioner%p_RdiscreteBC)
    
    ! Release temp vectors
    do ilev = rpreconditioner%NLMAX,rpreconditioner%NLMIN,-1
      call lsysbl_releaseVector(rpreconditioner%p_RtempVec(1,ilev))
      call lsysbl_releaseVector(rpreconditioner%p_RtempVec(2,ilev))
      call lsysbl_releaseVector(rpreconditioner%p_RtempVec(3,ilev))
    end do
    
    rpreconditioner%ctypePreconditioner = 0

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fbsim_preparePrecMatrixAssembly (rpreconditioner,ilev,&
      rnonlinearSpatialMatrixTemplate,rnonlinearSpatialMatrix)

!<description>
  ! Prepares a rnonlinearSpatialMatrix structure for the assembly according
  ! to a preconditioner on level ilev. rpreconditioner specifies a couple of preconditioner
  ! flags that configure the shape of the system matrix that the preconditioner
  ! needs. rnonlinearSpatialMatrixTemplate is a template structure defining the
  ! weights for all levels. The routine generates a new structure
  ! rnonlinearSpatialMatrix using the weights from rnonlinearSpatialMatrixTemplate
  ! and incorporating special assembly specific flags.
  !
  ! fbsim_initNonlinMatrix must have been called prior to this routine to
  ! initialise the basic matrix. fbsim_preparePrecondMatrixAssembly will then
  ! add assembly-specific parameters of the preconditioner.
!</description>

!<input>
  ! Current assembly level.
  integer, intent(in) :: ilev
  
  ! Structure with assembly-specific parameters of the preconditioner.
  type(t_fbsimPreconditioner), intent(in) :: rpreconditioner

  ! Template nonlinear matrix structure.
  type(t_nonlinearSpatialMatrix), intent(in) :: rnonlinearSpatialMatrixTemplate
!</input>

!<inputoutput>
  ! New nonlinear matrix structure based on rnonlinearSpatialMatrixTemplate.
  ! Assembly-specific parameters are added.
  type(t_nonlinearSpatialMatrix), intent(out) :: rnonlinearSpatialMatrix
!</inputoutput>
              
!</subroutine>

    ! Copy the template structure.
    rnonlinearSpatialMatrix = rnonlinearSpatialMatrixTemplate
    
    ! Parameters for adaptive matrices for Q1~ with anisotropic elements
    rnonlinearSpatialMatrix%iadaptiveMatrices = rpreconditioner%iadaptiveMatrices
    rnonlinearSpatialMatrix%dadmatthreshold = rpreconditioner%dadmatthreshold
    
    ! Pointers for that level.
    rnonlinearSpatialMatrix%p_rdiscretisation => &
        rpreconditioner%p_RlevelInfo(ilev)%rdiscretisation
    rnonlinearSpatialMatrix%p_rstaticInfo => &
        rpreconditioner%p_RlevelInfo(ilev)%rstaticInfo
    
    ! Depending on the level, we have to set up information about
    ! transposing B-matrices.
    if (ilev .eq. rpreconditioner%nlmin) then
      rnonlinearSpatialMatrix%bvirtualTransposedD = rpreconditioner%bneedVirtTransposedDonCoarse
    else
      rnonlinearSpatialMatrix%bvirtualTransposedD = rpreconditioner%bneedVirtTransposedD
    end if
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fbsim_assemblePrecMatrices (rpreconditioner,rnonlinearSpatialMatrix,&
      rx1,rx2,rx3,bincorporateBC)

!<description>
  ! Assembles on every level a matrix for the preconditioner in each timestep.
  ! The output is written to the preconditioner matrices in rpreconditioner.
!</description>

!<input>
  ! Weights defining the core equation.
  type(t_nonlinearSpatialMatrix), intent(IN) :: rnonlinearSpatialMatrix

  ! Current iteration vector of the 'previous' timestep. May be undefined
  ! if there is no previous timestep. 
  type(t_vectorBlock), intent(IN), target :: rx1

  ! Current iteration vector. 
  type(t_vectorBlock), intent(IN), target :: rx2

  ! Current iteration vector of the 'next' timestep. May be undefined
  ! if there is no previous timestep. 
  type(t_vectorBlock), intent(IN), target :: rx3
  
  ! Defines if boundary conditions are implemented into the matrix.
  ! If this is .true., the boundary conditions for the current time
  ! must be available in rpreconditioner!
  logical :: bincorporateBC
!</input>

!<inputoutput>
  ! Spatial preconditioner structure where the preconditioner matrices
  ! should be updated.
  ! If bincorporateBC=.true., the structure must provide the boundary conditions!
  ! p_RmatrixPrecondFullSpace and p_RmatrixPrecond are assembled, boundary
  ! conditions attached.
  type(t_fbsimPreconditioner), intent(inout), target :: rpreconditioner
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ilev
    type(t_matrixBlock), pointer :: p_rmatrix,p_rmatrixFine
    type(t_vectorBlock), pointer :: p_rvectorFine1,p_rvectorFine2,p_rvectorFine3
    type(t_vectorBlock), pointer :: p_rvectorCoarse1,p_rvectorCoarse2,p_rvectorCoarse3
    type(t_nonlinearSpatialMatrix) :: rlocalNonlSpatialMatrix
    integer, dimension(1), parameter :: Irows = (/1/)

    ! DEBUG!!!
    real(dp), dimension(:), pointer :: p_Da,p_vec,p_def
    
    ! DEBUG!!!
    !call lsysbl_getbase_double (rd,p_def)
    !call lsysbl_getbase_double (rx,p_vec)

    ! On all levels, we have to set up the nonlinear system matrix,
    ! so that the linear solver can be applied to it.
    
    nullify(p_rmatrix)

    do ilev=rpreconditioner%NLMAX,rpreconditioner%NLMIN,-1
    
      ! Get the matrix on the current level.
      ! Shift the previous matrix to the pointer of the fine grid matrix.
      p_rmatrixFine => p_rmatrix
      p_rmatrix => rpreconditioner%p_RmatrixPrecondFullSpace(ilev)
      
      ! DEBUG!!!
      call lsyssc_getbase_double (p_rmatrix%RmatrixBlock(1,1),p_Da)
    
      ! On the highest level, we use rx as solution to build the nonlinear
      ! matrix. On lower levels, we have to create a solution
      ! on that level from a fine-grid solution before we can use
      ! it to build the matrix!
      if (ilev .eq. rpreconditioner%NLMAX) then
      
        p_rvectorCoarse1 => rx1
        p_rvectorCoarse2 => rx2
        p_rvectorCoarse3 => rx3
        
      else
        ! We have to discretise a level hierarchy and are on a level < NLMAX.

        ! Get the temporary vector on level i. Will receive the solution
        ! vector on that level. 
        p_rvectorCoarse1 => rpreconditioner%p_rtempVec(1,ilev)
        p_rvectorCoarse2 => rpreconditioner%p_rtempVec(2,ilev)
        p_rvectorCoarse3 => rpreconditioner%p_rtempVec(3,ilev)
        
        ! Get the solution vector on level i+1. This is either the temporary
        ! vector on that level, or the solution vector on the maximum level.
        if (ilev .lt. rpreconditioner%NLMAX-1) then
          p_rvectorFine1 => rpreconditioner%p_rtempVec(1,ilev+1)
          p_rvectorFine2 => rpreconditioner%p_rtempVec(2,ilev+1)
          p_rvectorFine3 => rpreconditioner%p_rtempVec(3,ilev+1)
        else
          p_rvectorFine1 => rx1
          p_rvectorFine2 => rx2
          p_rvectorFine3 => rx3
        end if

        ! Interpolate the solution from the finer grid to the coarser grid.
        ! The interpolation is configured in the interlevel projection
        ! structure we got from the collection.
        call mlprj_performInterpolationHier (rpreconditioner%rprjHierarchy,ilev+1,&
            p_rvectorCoarse1,p_rvectorFine1)
        call mlprj_performInterpolationHier (rpreconditioner%rprjHierarchy,ilev+1,&
            p_rvectorCoarse2,p_rvectorFine2)
        call mlprj_performInterpolationHier (rpreconditioner%rprjHierarchy,ilev+1,&
            p_rvectorCoarse3,p_rvectorFine3)

        ! Apply the filter chain to the temp vector.
        ! This implements the boundary conditions that are attached to it.
        ! NOTE: Deactivated for standard CC2D compatibility -- and because
        ! it has to be checked whether the correct boundary conditions
        ! are attached to that vector!
        ! CALL filter_applyFilterChainVec (p_rvectorCoarse, p_RfilterChain)

      end if

      ! Generate a local nonlinear matrix structure for the current level.
      call fbsim_preparePrecMatrixAssembly (rpreconditioner,ilev,&
          rnonlinearSpatialMatrix,rlocalNonlSpatialMatrix)

      ! Assemble the matrix.
      ! If we are on a lower level, we can specify a 'fine-grid' matrix.
      if (ilev .eq. rpreconditioner%NLMAX) then
        call cc_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
            p_rmatrix,rlocalNonlSpatialMatrix,&
            p_rvectorCoarse1,p_rvectorCoarse2,p_rvectorCoarse3)
      else
        call cc_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
            p_rmatrix,rlocalNonlSpatialMatrix,&
            p_rvectorCoarse1,p_rvectorCoarse2,p_rvectorCoarse3,&
            p_rmatrixFine)
      end if

    end do
    
    if (rpreconditioner%bneedPressureDiagonalBlock) then
      
      ! The 3,3-matrix must exist! This is ensured by the initialisation routine.
      !
      ! We have a pure Dirichlet problem. This may give us some difficulties
      ! in the case, the preconditioner uses a direct solver (UMFPACK).
      ! In this case, we have to include a unit vector to the pressure
      ! matrix to make the problem definite!
      if (rpreconditioner%isolverType .eq. 0) then
        p_rmatrix => rpreconditioner%p_RmatrixPrecondFullSpace(rpreconditioner%NLMAX)
        
        ! Include a unit vector to the matrix part of the pressure in
        ! the primal equation -- as long as there is not a full identity
        ! matrix in the pressure matrix (what would be the case for 
        ! the initial condition).
        if (rlocalNonlSpatialMatrix%Dkappa(1,1) .eq. 0.0_DP) then
          ! Switch the pressure matrix on and clear it; we don't know what is inside.
          p_rmatrix%RmatrixBlock(3,3)%dscaleFactor = 1.0_DP
          call lsyssc_clearMatrix (p_rmatrix%RmatrixBlock(3,3))
          call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,1),Irows)
          call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,2),Irows)
          call mmod_replaceLinesByUnit(p_rmatrix%RmatrixBlock(3,3),Irows)
          call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,4),Irows)
          call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,5),Irows)
          call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,6),Irows)
        end if

        ! Also in the dual equation, as the BC type coincides
        if (rlocalNonlSpatialMatrix%Dkappa(2,2) .eq. 0.0_DP) then
          ! Switch the pressure matrix on and clear it; we don't know what is inside.
          p_rmatrix%RmatrixBlock(6,6)%dscaleFactor = 1.0_DP
          call lsyssc_clearMatrix (p_rmatrix%RmatrixBlock(6,6))
          call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,1),Irows)
          call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,2),Irows)
          call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,3),Irows)
          call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,4),Irows)
          call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,5),Irows)
          call mmod_replaceLinesByUnit(p_rmatrix%RmatrixBlock(6,6),Irows)
        end if
        
      end if
      
      if (rpreconditioner%isolverType .eq. 1) then
      
        ! If we have a MG solver, We also check the coarse grid solver for 
        ! the same thing!
        ! What we don't check is the smoother, thus we assume that smoothers
        ! are always solvers that allow the applicance of a filter chain.
        if (rpreconditioner%icoarseGridSolverType .eq. 0) then
          p_rmatrix => rpreconditioner%p_RmatrixPrecondFullSpace(rpreconditioner%NLMIN)
          
          ! Include a unit vector to the matrix part of the pressure in
          ! the primal equation -- as long as there is not a full identity
          ! matrix in the pressure matrix (what would be the case for 
          ! the initial condition).
          if (rlocalNonlSpatialMatrix%Dkappa(1,1) .eq. 0.0_DP) then
            ! Switch the pressure matrix on and clear it; we don't know what is inside.
            p_rmatrix%RmatrixBlock(3,3)%dscaleFactor = 1.0_DP
            call lsyssc_clearMatrix (p_rmatrix%RmatrixBlock(3,3))
            call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,1),Irows)
            call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,2),Irows)
            call mmod_replaceLinesByUnit(p_rmatrix%RmatrixBlock(3,3),Irows)
            call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,4),Irows)
            call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,5),Irows)
            call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,6),Irows)
          end if

          ! Also in the dual equation, as the BC type coincides
          if (rlocalNonlSpatialMatrix%Dkappa(2,2) .eq. 0.0_DP) then
            ! Switch the pressure matrix on and clear it; we don't know what is inside.
            p_rmatrix%RmatrixBlock(6,6)%dscaleFactor = 1.0_DP
            call lsyssc_clearMatrix (p_rmatrix%RmatrixBlock(6,6))
            call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,1),Irows)
            call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,2),Irows)
            call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,3),Irows)
            call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,4),Irows)
            call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,5),Irows)
            call mmod_replaceLinesByUnit(p_rmatrix%RmatrixBlock(6,6),Irows)
          end if
          
        end if
        
      end if
        
    end if        
    
    ! Extract the correct submatrices for the preconiditioner from the full matrix.
    ! The extracted submatrix identifies either the diagonal of the primal
    ! space, the dual space or the full matrix.
    do ilev=rpreconditioner%NLMIN,rpreconditioner%NLMAX
    
      select case (rpreconditioner%cspace)
      
      case (CCSPACE_PRIMAL)
      
        ! Get the primal sub-operators
        call lsysbl_deriveSubmatrix (&
            rpreconditioner%p_RmatrixPrecondFullSpace(ilev),&
            rpreconditioner%p_RmatrixPrecond(ilev),&
            LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE,1,3)
            
        call lsysbl_assignDiscrDirectMat (&
            rpreconditioner%p_RmatrixPrecond(ilev),&
            rpreconditioner%p_RlevelInfo(ilev)%rdiscretisationPrimal)
            
      case (CCSPACE_DUAL)

        ! Get the dual sub-operators
        call lsysbl_deriveSubmatrix (&
            rpreconditioner%p_RmatrixPrecondFullSpace(ilev),&
            rpreconditioner%p_RmatrixPrecond(ilev),&
            LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE,4,6)
            
        call lsysbl_assignDiscrDirectMat (&
            rpreconditioner%p_RmatrixPrecond(ilev),&
            rpreconditioner%p_RlevelInfo(ilev)%rdiscretisationPrimal)
            
      case (CCSPACE_PRIMALDUAL)

        ! Take the full matrix.
        call lsysbl_duplicateMatrix (&
            rpreconditioner%p_RmatrixPrecondFullSpace(ilev),&
            rpreconditioner%p_RmatrixPrecond(ilev),&
            LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)

        call lsysbl_assignDiscrDirectMat (&
            rpreconditioner%p_RmatrixPrecond(ilev),&
            rpreconditioner%p_RlevelInfo(ilev)%rdiscretisation)
      
      end select
      
      ! Attach boundary conditions.
      call lsysbl_assignDiscreteBC (&
          rpreconditioner%p_RmatrixPrecond(ilev),&
          rpreconditioner%p_RdiscreteBC(ilev))

      call lsysbl_assignDiscreteFBC (&
          rpreconditioner%p_RmatrixPrecond(ilev),&
          rpreconditioner%p_RdiscreteFBC(ilev))
      
      ! Boundary conditions
      ! ---------------------------------------------------
      ! We implement the BC`s in the preconditioner matrix.
      if (bincorporateBC) then

        ! Call the matrix filter for the boundary conditions to include the BC's
        ! into the matrix.
        call matfil_discreteBC (rpreconditioner%p_RmatrixPrecond(ilev),&
          rpreconditioner%p_RdiscreteBC(ilev))
        call matfil_discreteFBC (rpreconditioner%p_RmatrixPrecond(ilev),&
          rpreconditioner%p_RdiscreteFBC(ilev))
        
        ! 'Nonlinear' boundary conditions like slip boundary conditions
        ! are not implemented with a filter chain into a matrix.
        ! Call the appropriate matrix filter of 'nonlinear' boundary
        ! conditions manually:
        call matfil_discreteNLSlipBC (rpreconditioner%p_RmatrixPrecond(ilev),&
          .true.,rpreconditioner%p_RdiscreteBC(ilev))
        
        ! DEBUG!!!
        !CALL matio_writeBlockMatrixHR (p_rmatrix, 'matrix',&
        !                              .TRUE., 0, 'matrix.txt','(E20.10)')
      end if

    end do
    
  end subroutine

  ! ***************************************************************************

  !<subroutine>

    subroutine fbsim_precondSpaceDefect (rpreconditioner,rd,bsuccess)
  
    use linearsystemblock
    use collection
    
  !<description>
    ! Applies the spatial preconditioner to a defect vector rd.
  !</description>

  !<inputoutput>
    ! Spatial preconditioner structure that defines all parameters how to perform
    ! preconditioning.
    type(t_fbsimPreconditioner), intent(inout) :: rpreconditioner

    ! Defect vector b-A(rx)x. This must be replaced by J^{-1} rd by a preconditioner.
    type(t_vectorBlock), intent(INOUT)            :: rd

    ! If the preconditioning was a success. Is normally automatically set to
    ! TRUE. If there is an error in the preconditioner, this flag can be
    ! set to FALSE. In this case, the nonlinear solver breaks down with
    ! the error flag set to 'preconditioner broke down'.
    logical, intent(INOUT)                        :: bsuccess
  !</inputoutput>
  
  !</subroutine>

      ! local variables
      integer :: ierror
      type(t_linsolNode), pointer :: p_rsolverNode 

      ! DEBUG!!!
      real(dp), dimension(:), pointer :: p_def
      call lsysbl_getbase_double (rd,p_def)

      select case (rpreconditioner%ctypePreconditioner)
      case (1)
        ! Preconditioning with a linear solver.
        !
        ! Get the solver node
      
        p_rsolverNode => rpreconditioner%p_rsolverNode

        ! DEBUG!!!
        !CALL matio_writeBlockMatrixHR (Rmatrices(rpreconditioner%NLMIN), 'matrix',&
        !                               .TRUE., 0, 'matrixstat.txt','(E10.2)')
        
        call linsol_setMatrices(p_rsolverNode,&
            rpreconditioner%p_RmatrixPrecond(rpreconditioner%NLMIN:rpreconditioner%NLMAX))
        
        ! DEBUG!!!
        !DO i=rnonlinearIteration%NLMIN,rnonlinearIteration%NLMAX
        !  CALL storage_getbase_double (Rmatrices(i)% &
        !      RmatrixBlock(4,1)%h_Da,p_Ddata)
        !END DO
        
        ! Initialise data of the solver. This in fact performs a numeric
        ! factorisation of the matrices in UMFPACK-like solvers.
        call linsol_initStructure (rpreconditioner%p_rsolverNode,ierror)
        call linsol_initData (p_rsolverNode, ierror)
        if (ierror .ne. LINSOL_ERR_NOERROR) then
          print *,'linsol_initData failed!'
          call sys_halt()
        end if
        
        ! Solve the system. As we want to solve Ax=b with
        ! b being the real RHS and x being the real solution vector,
        ! we use linsol_solveAdaptively. If b is a defect
        ! RHS and x a defect update to be added to a solution vector,
        ! we would have to use linsol_precondDefect instead.
        call linsol_precondDefect (p_rsolverNode,rd)

        ! Release the numeric factorisation of the matrix.
        ! We don't release the symbolic factorisation, as we can use them
        ! for the next iteration.
        call linsol_doneData (p_rsolverNode)
        call linsol_doneStructure (p_rsolverNode)
        
        ! Did the preconditioner work?
        bsuccess = p_rsolverNode%iresult .eq. 0
        
        if (bsuccess) then
          ! Filter the final defect.
          call vecfil_discreteBCdef (rd,rpreconditioner%p_RdiscreteBC(rpreconditioner%NLMAX))
          call vecfil_discreteFBCdef (rd,rpreconditioner%p_RdiscreteFBC(rpreconditioner%NLMAX))
        end if
        
        if (p_rsolverNode%dfinalDefect .gt. p_rsolverNode%dinitialDefect*0.99_DP) then
          ! Ignore the correction, it cannot be good enough!
          call output_line (&
            'Space-Preconditioner: Warning. Solution ignored for missing accuracy.')
            
          call lsysbl_clearVector (rd)
        end if
        
      end select
      
    end subroutine

! ***************************************************************************

!<subroutine>

  subroutine fbsim_getDefectNorm (rvector,rrhs,rdefect,Dresiduals)
  
!<description>
  ! Calculates a couple of norms from a given solution, defect and RHS vector
  ! of the nonlinear system. This can be used to check convergence criteria
  ! etc.
!</description>

!<input>
  ! The solution vector which is modified later during the nonlinear iteration.
  type(t_vectorBlock), intent(in) :: rvector

  ! The right-hand-side vector to use in the equation
  type(t_vectorBlock), intent(in) :: rrhs
  
  ! A defect vector calculated with rvector and rrhs
  type(t_vectorBlock), intent(in) :: rdefect
!</input>

!<output>
  ! An array receiving different defect norms calculated by the above vectors.
  ! Dresiduals(1) = RESU   = ||defect_u|| / ||rhs|| = velocity residual
  ! Dresiduals(2) = RESDIV = ||p|| / ||u||          = divergence residual
  real(DP), dimension(:), intent(out) :: Dresiduals
!</output>

!</subroutine>

    ! local variables
    real(DP) :: dresF,DresTmp(2),dnormU
    integer, dimension(2) :: Cnorms

    !-----------------------------------------------------------------------
    !     Compute the relative l2-norms  RESU,RESDIV
    !-----------------------------------------------------------------------

    Cnorms(:) = LINALG_NORMEUCLID

    ! RESF := max ( ||F1||_E , ||F2||_E )

    call lsysbl_vectorNormBlock (rrhs,Cnorms,DresTmp)
    dresF = max(DresTmp(1),DresTmp(2))
    if (dresF .lt. 1.0E-8_DP) dresF = 1.0_DP

    !               || (D1,D2) ||_E
    ! RESU = -----------------------------
    !        max ( ||F1||_E , ||F2||_E )

    call lsysbl_vectorNormBlock (rdefect,Cnorms,DresTmp)
    Dresiduals(1) = sqrt(DresTmp(1)**2+DresTmp(2)**2)/dresF

    ! DNORMU = || (U1,U2) ||_l2 

    call lsysbl_vectorNormBlock (rvector,Cnorms,DresTmp)
    dnormU = sqrt(DresTmp(1)**2+DresTmp(2)**2)
    if (dnormU .lt. 1.0E-8_DP) dnormU = 1.0_DP

    !             || DP ||_E
    ! RESDIV = ----------------
    !          || (U1,U2) ||_E

    Dresiduals(2) = &
        lsyssc_vectorNorm (rdefect%RvectorBlock(3),LINALG_NORMEUCLID) / dnormU
        
  end subroutine

  ! ***************************************************************************

  subroutine fbsim_resNormCheck (rnonlinearIteration,&
      ite,rx,rb,rd,bconvergence,bdivergence)

  use linearsystemblock
  use collection
  
!<description>
  ! Residual norm calculation & printing routine.
  ! This routine is called each time the norm of the residuum was calculated.
  ! It has to check the current residuum for convergence and/or divergence
  ! and can print the residuum to screen.
!</description>

!<inputoutput>
  ! Reference to the nonlinear iteration structure that configures the
  ! main nonlinear equation. Intermediate data is changed during the iteration.
  type(t_fbsimNonlinearIteration), intent(inout)   :: rnonlinearIteration

  ! Number of current iteration. Is set to 0 when the callback routine
  ! is called the first time. In this situation, rd describes the initial
  ! defect vector.
  integer, intent(in)                           :: ite

  ! Current iteration vector
  type(t_vectorBlock), intent(in), target       :: rx

  ! Right hand side vector of the equation.
  type(t_vectorBlock), intent(in), target       :: rb

  ! Defect vector b-A(x)x.
  type(t_vectorBlock), intent(in), target       :: rd
!</inputoutput>

!<output>
  ! Is set to TRUE if the residuum rd
  ! is within a desired tolerance, so that the solver should treat
  ! the iteration as 'converged'.
  logical, intent(out)                        :: bconvergence

  ! Is TRUE if the residuum rd
  ! is out of a desired tolerance, so that the solver should treat
  ! the iteration as 'diverged'.
  logical, intent(out)                        :: bdivergence
!</output>

    ! local variables
    real(DP), dimension(3) :: Dresiduals
    real(DP) :: dresOld,drhoNL,ddelP,ddelU,dtmp,dresU,dresDIV,dres,dresINIT
    real(DP) :: depsD,depsDiv,depsUR,depsPR,depsRES
    integer, dimension(3) :: Cnorms

    ! Calculate norms of the solution/defect vector
    call fbsim_getDefectNorm (rx,rb,rd,Dresiduals)
    
    Dresiduals(3) = sqrt(Dresiduals(1)**2 + Dresiduals(2)**2)

    ! In the first iteration (initial defect), print the norm of the defect
    ! and save the norm of the initial residuum to the structure
    if (ite .eq. 1) then
    
      call output_separator (OU_SEP_MINUS,coutputMode=rnonlinearIteration%coutputMode)
      call output_line (' IT  RELU     RELP     DEF-U    DEF-DIV'// &
                        '  DEF-TOT  RHONL    OMEGNL   RHOMG',&
                        coutputMode=rnonlinearIteration%coutputMode)
      call output_separator (OU_SEP_MINUS,coutputMode=rnonlinearIteration%coutputMode)     
      call output_line ('  0                   '// &
          trim(sys_sdEP(Dresiduals(1),9,2))//&
          trim(sys_sdEP(Dresiduals(2),9,2))//&
          trim(sys_sdEP(Dresiduals(3),9,2)),coutputMode=rnonlinearIteration%coutputMode)
      call output_separator (OU_SEP_MINUS,coutputMode=rnonlinearIteration%coutputMode)     

      rnonlinearIteration%DresidualInit (1:2) = Dresiduals(1:2)
      rnonlinearIteration%DresidualOld (1:2) = Dresiduals(1:2)

      rnonlinearIteration%dinitialDefectTotal = sqrt(Dresiduals(1)**2 + &
                                                     Dresiduals(2)**2)
      rnonlinearIteration%dfinalDefectTotal = rnonlinearIteration%dinitialDefectTotal

      bconvergence = .false.
      bdivergence = .false.
    
    else
      ! In the other iterations, calculate the relative change and
      ! test the convergence criteria.
    
      ! Old defect:
      dresOld = sqrt(rnonlinearIteration%DresidualOld(1)**2 + &
                      rnonlinearIteration%DresidualOld(2)**2)
      
      ! Initial defect
      dresINIT = sqrt(rnonlinearIteration%DresidualInit(1)**2 + &
                      rnonlinearIteration%DresidualInit(2)**2)
                      
      ! dresInit=0 may hardly occur -- except when we expect 'no flow'.
      ! But to prevent a check against "something<=0" in this case below,
      ! set dresInit to something <> 0.
      if (dresINIT .eq. 0.0_DP) dresINIT = 1.0_DP

      ! Replace the 'old' residual by the current one
      rnonlinearIteration%DresidualOld(1:2) = Dresiduals(1:2)

      ! Nonlinear convergence rate
      drhoNL = (Dresiduals(3)/dresINIT) ** (1.0_DP/real(ite,DP))
      
      ! Calculate norms of the solution/defect vector, calculated above
      dresU   = Dresiduals(1)
      dresDIV = Dresiduals(2)
      dres    = sqrt(dresU**2 + dresDIV**2)
      
      ! Calculate relative maximum changes 
      ! This simply calculates some postprocessing values of the relative
      ! change in the solution.
      !
      ! Maximum norm of solution vector:
      !
      !              || YP ||_max       || Pnew - Pold ||_max
      !   DELP := ------------------- = ---------------------
      !              || P ||_max           || Pnew ||_max
      !
      !
      ! Relative change of solution vector:
      !
      !            || (Y1,Y2) ||_max    || Unew - Uold ||_max
      !   DELU := ------------------- = --------------------- 
      !           || (KU1,KU2) ||_max       || Unew ||_max
      !
      ! The norms || YP ||_max, || Yi ||_max are saved in the nonlinear
      ! iteration structure from the last preconditioning!
      ! The other MAX-norms have to be calculated from U...
      
      Cnorms(:) = LINALG_NORMMAX
      call lsysbl_vectorNormBlock (rx,Cnorms,Dresiduals)

      dtmp = max(Dresiduals(1),Dresiduals(2))
      if (dtmp .lt. 1.0E-8_DP) dtmp = 1.0_DP
      ddelU = max(rnonlinearIteration%DresidualCorr(1),&
                  rnonlinearIteration%DresidualCorr(2))/dtmp
      
      dtmp = Dresiduals(3)
      if (dtmp .lt. 1.0E-8_DP) dtmp = 1.0_DP
      ddelP = rnonlinearIteration%DresidualCorr(3)/dtmp
      
      ! Check if the nonlinear iteration can prematurely terminate.
      !        
      ! Get the stopping criteria from the parameters.
      ! Use the DepsNL data according to the initialisation above.
      depsD   = rnonlinearIteration%DepsNL(1)
      depsDiv = rnonlinearIteration%DepsNL(2)
      depsUR  = rnonlinearIteration%DepsNL(3)
      depsPR  = rnonlinearIteration%DepsNL(4)
      depsRES = rnonlinearIteration%DepsNL(5)*dresINIT ! -> ddampingD
      
      ! All residual information calculated.
      ! Check for divergence; use a 'NOT' for better NaN handling.
      bdivergence = .not. (dres/dresINIT .lt. 1E5)
      
      ! Check for convergence
      if((ddelU .le. depsUR).and.(ddelP .le. depsPR)   .and. &
          (dresU .le. depsD) .and.(dresDiv .le. depsDiv).and. &
          (dres .le. depsRES)) then
        bconvergence = .true.
      else
        bconvergence = .false.
      end if

      if ((ddelU .lt. SYS_EPSREAL_DP*1E2_DP) .and. &
          (ddelP .lt. SYS_EPSREAL_DP*1E2_DP)) then
        ! We are hard on machine exactness, so stop the iteraton
        bconvergence =.true.
      end if
      
      ! Print residual information
      call output_line ( &
          trim(sys_si(ite,3))//' '//&
          trim(sys_sdEP(ddelU,9,2))// &
          trim(sys_sdEP(ddelP,9,2))// &
          trim(sys_sdEP(dresU,9,2))// &
          trim(sys_sdEP(dresDIV,9,2))// &
          trim(sys_sdEP(dres,9,2))// &
          trim(sys_sdEP(drhoNL,9,2))// &
          trim(sys_sdEP(rnonlinearIteration%domegaNL,9,2))// &
          trim(sys_sdEP(rnonlinearIteration%drhoLinearSolver,9,2)), &
          coutputMode=rnonlinearIteration%coutputMode)
      
    end if
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fbsim_simulate (rsimsolver, rsolvector, rrhsvector, &
      ifirststep, ilaststep, roseenSolution, bsuccess)
  
!<description>
  ! Performs a forward simulation through the primal solution or a
  ! backward sweep through the dual solution. The forward sweep can
  ! be performed linearly (Oseen equation) or nonlinearly ((Navier-)Stokes).
!</description>
  
!<input>
  ! Solver configuration of the simulation solver.
  type(t_simSolver), intent(inout), target :: rsimsolver

  ! A space-time RHS vector for all timesteps.
  type(t_spacetimevector), intent(inout) :: rrhsvector
  
  ! OPTIONAL: Evaluation point of the nonlinearity in the primal space.
  ! This is not used in a pure forward simulation.
  type(t_spacetimevector), intent(inout), optional :: roseenSolution
  
  ! First timestep to simulate
  integer, intent(in) :: ifirststep

  ! Last timestep to simulate
  integer, intent(in) :: ilaststep
!</input>

!<inputoutput>
  ! An initial space-time solution vector. This is modified according to the iteration.
  type(t_spacetimevector), intent(inout) :: rsolvector
  
  ! TRUE if the iteration was successful.
  logical, intent(out) :: bsuccess
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: istep
    real(dp) :: dtstep,dtime
    type(t_vectorBlock) :: roseensol1, roseensol2, roseensol3
    type(t_vectorBlock) :: rprevsol, rcurrentsol, rprevrhs, rcurrentrhs
    type(t_vectorBlock) :: rdefect
    type(t_vectorBlock) :: rdefectPrimal,rdefectDual
    type(t_problem_lvl), pointer :: p_rlevelInfo 
    real(dp) :: dtheta
    type(t_nonlinearSpatialMatrix) :: rnonlinearSpatialMatrix
    type(t_ccoptSpaceTimeMatrix), pointer :: p_rspaceTimeMatrix
    
    ! nonlinear solver
    integer :: ite
    real(dp) :: dres,dresInit,dtempdef
    logical :: bconvergence,bdivergence
    type(t_ccoptSpaceTimeMatrix), pointer :: rspaceTimeMatrixPrecond
    type(t_ccDynamicNewtonControl), pointer :: p_radaptiveNewton
    
    ! DEBUG!!!
    real(dp), dimension(:), pointer :: p_Dsol,p_Drhs,p_Doseen1,p_Doseen2,p_Doseen3
    real(dp), dimension(:), pointer :: p_Ddefect,p_DdefP,p_DdefD
    real(dp), dimension(:), pointer :: p_DprevSol,p_DprevRhs
    
    ! Fetch some necessary information from the structures.
    p_rspaceTimeMatrix => rsimsolver%p_rmatrix
    dtheta = rsimsolver%p_rproblem%rtimedependence%dtimeStepTheta
    dtstep = p_rspaceTimeMatrix%p_rspaceTimeDiscr%rtimeDiscr%dtstep
    
    ! Get the level information structure of the current level
    p_rlevelInfo => rsimsolver%p_rlevelInfo
    
    ! Create some temp vectors
    call lsysbl_createVectorBlock(p_rlevelInfo%rdiscretisation,rprevsol)
    call lsysbl_createVectorBlock(p_rlevelInfo%rdiscretisation,rcurrentsol)
    call lsysbl_createVectorBlock(p_rlevelInfo%rdiscretisation,roseensol1)
    call lsysbl_createVectorBlock(p_rlevelInfo%rdiscretisation,roseensol2)
    call lsysbl_createVectorBlock(p_rlevelInfo%rdiscretisation,roseensol3)
    call lsysbl_createVectorBlock(p_rlevelInfo%rdiscretisation,rprevrhs)
    call lsysbl_createVectorBlock(p_rlevelInfo%rdiscretisation,rcurrentrhs)

    call lsysbl_createVectorBlock(p_rlevelInfo%rdiscretisation,rdefect)
    
    ! DEBUG!!!
    call lsysbl_getbase_double (roseensol1,p_Doseen1)
    call lsysbl_getbase_double (roseensol2,p_Doseen2)
    call lsysbl_getbase_double (roseensol3,p_Doseen3)
    call lsysbl_getbase_double (rcurrentsol,p_Dsol)
    call lsysbl_getbase_double (rcurrentrhs,p_Drhs)
    call lsysbl_getbase_double (rprevrhs,p_DprevRhs)
    call lsysbl_getbase_double (rprevsol,p_DprevSol)
    
    ! Create a temp matrix.
    
    ! Create subvectors for the primal/dual part of the defect.
    call lsysbl_deriveSubvector(rdefect,rdefectPrimal,1,3,.true.)
    call lsysbl_deriveSubvector(rdefect,rdefectDual,4,6,.true.)
    call lsysbl_assignDiscrDirectVec (rdefectPrimal,p_rlevelInfo%rdiscretisationPrimal)
    call lsysbl_assignDiscrDirectVec (rdefectDual,p_rlevelInfo%rdiscretisationPrimal)

    ! DEBUG!!!
    call lsysbl_getbase_double (rdefect,p_Ddefect)
    call lsysbl_getbase_double (rdefectPrimal,p_DdefP)
    call lsysbl_getbase_double (rdefectDual,p_DdefD)
    
    ! What to do?
    select case (rsimsolver%csimtype)
    case (0)
      ! (Navier-)Stokes forward simulation. This iteration type is nonlinear 
      ! in each timestep.

      ! Assign BC`s to the defect vector.
      call lsysbl_assignDiscreteBC (rdefectPrimal,&
          rsimSolver%rpreconditioner%p_RdiscreteBC(rsimSolver%ilevel))
      call lsysbl_assignDiscreteFBC (rdefectPrimal,&
          rsimSolver%rpreconditioner%p_RdiscreteFBC(rsimSolver%ilevel))
          
      ! Derive a "preconditioner" space-time matrix from the basic space-time
      ! matrix. This may or may not include the Newton part - but by default it
      ! does not.
      rspaceTimeMatrixPrecond = p_rspaceTimeMatrix
      rspaceTimeMatrixPrecond%cmatrixType = 0

      ! Loop through the timesteps.
      do istep = ifirststep,ilaststep
      
        ! Current time step?
        dtime = &
            p_rspaceTimeMatrix%p_rspaceTimeDiscr%rtimeDiscr%dtimeInit + (istep-1) * dtstep

        ! Update the boundary conditions to the current time
        call fbsim_updateDiscreteBCprec (rsimsolver%p_rproblem,rsimsolver%rpreconditioner,dtime)

        ! Initialise the previous and current timestep. 
        if (istep .eq. ifirststep) then

          ! Initial condition is given. Read it.
          call sptivec_getTimestepData (rrhsvector, ifirststep, rcurrentrhs)
          call sptivec_getTimestepData (rsolvector, ifirststep, rcurrentsol)
          
          ! The "Oseen" solutions coincide with the current solution.
          ! The "next" Oseen solution is always zero as we do not have
          ! a "next" solution in a pure forward simulation.
          call lsysbl_clearVector (roseensol1)
          call lsysbl_copyVector (rcurrentsol,roseensol2)
          call lsysbl_clearVector (roseensol3)

        else  
                
          ! Solution.
          call lsysbl_copyVector (rcurrentsol,rprevsol)
          call lsysbl_copyVector (rcurrentsol,roseensol1)
          call sptivec_getTimestepData (rsolvector, istep, rcurrentsol)
          call lsysbl_copyVector (rcurrentsol,roseensol2)

          ! RHS.
          call lsysbl_copyVector (rcurrentrhs,rprevrhs)
          call sptivec_getTimestepData (rrhsvector, istep, rcurrentrhs)
          
        end if
        
        ! Nonlinear loop in the timestep
        ! ------------------------------
        ! Start the nonlinear solver
        
        do ite = 1,rsimsolver%rnonlinearIteration%nmaxIterations
        
          ! 1.) Create the nonlinear defect
          ! -------------------------------
        
          ! Get the RHS.
          call lsysbl_copyVector (rcurrentrhs,rdefect)
          
          ! If we are not in the first step...
          if (istep .ne. ifirststep) then
            ! ... subtract the offdiagonal, corresponding to the previous timestep.
            call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rsimsolver%p_rproblem,&
                p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
                p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
            call cc_setupMatrixWeights (rsimsolver%p_rproblem,p_rspaceTimeMatrix,dtheta,&
              istep-1,-1,rnonlinearSpatialMatrix)
            
            call cc_assembleDefect (rnonlinearSpatialMatrix,rprevsol,&
              rdefect,1.0_DP,roseensol1,roseensol2,roseensol3)
          end if

          ! Subtract A*primal/dual solution from the rhs to create the primal defect.
          call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rsimsolver%p_rproblem,&
              p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
              p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
          call cc_setupMatrixWeights (rsimsolver%p_rproblem,p_rspaceTimeMatrix,dtheta,&
            istep-1,0,rnonlinearSpatialMatrix)
          call cc_disableSubmatrix (rnonlinearSpatialMatrix,2,1)
          call cc_disableSubmatrix (rnonlinearSpatialMatrix,2,2)
          
          call cc_assembleDefect (rnonlinearSpatialMatrix,rcurrentsol,&
            rdefect,1.0_DP,roseensol1,roseensol2,roseensol3)

          ! Implement the boundary conditions        
          call vecfil_discreteBCdef(rdefectPrimal)
          
          ! 2.) Check for convergence
          ! -------------------------
          
          ! Calculate the nonlinear defect.
          ! In the first iteration, the result is also the initial defect.
          call fbsim_resNormCheck (rsimSolver%rnonlinearIteration,&
              ite,rcurrentsol,rcurrentrhs,rdefect,bconvergence,bdivergence)
          
          ! Perform at least nminIterations iterations
          if (ite .lt. rsimSolver%rnonlinearIteration%nminIterations) bconvergence = .false.
        
          ! Check for convergence
          if (bconvergence) then
            bsuccess = .true.
            exit
          end if
        
          ! Check for divergence
          if (bdivergence) then
            if (rsimSolver%rnonlinearIteration%ioutputLevel .ge. 0) then
              call output_line ('NLSOL: Iteration '//&
                  trim(sys_siL(ite,10))//' canceled, divergence detected!',&
                  coutputMode=rsimSolver%rnonlinearIteration%coutputMode)
            end if
            bsuccess = .false.
            exit
          end if
          
          ! 3.) Preconditioning to the defect
          ! ---------------------------------
        
          ! Based on the preconditioner space-time matrix, initialise the weights
          ! for the preconditioner and assemble the subblock of this preconditioner.
          ! Switch off all matrices except the primal diagonal block since that is
          ! the only matrix we need here.

          call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rsimsolver%p_rproblem,&
              p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
              p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
              
          ! Newton or standard matrix?
          select case (rsimsolver%rnonlinearIteration%ctypeIteration)
          case (0)
            ! Standard iteration.
            rspaceTimeMatrixPrecond%cmatrixType = 0
            
          case (1)
            ! Newton iteration
            rspaceTimeMatrixPrecond%cmatrixType = 1
            
          case (2)
            ! Adaptive Newton iteration. That's a bit more involving.
            ! We start being a standard iteration.
            rspaceTimeMatrixPrecond%cmatrixType = 0
            
            p_radaptiveNewton => rsimsolver%rnonlinearIteration%radaptiveNewton
            
            if (ite .gt. p_radaptiveNewton%nmaxFixPointIterations) then
              ! Force Newton to be used.
              rspaceTimeMatrixPrecond%cmatrixType = 1
            else
              if (ite .gt. p_radaptiveNewton%nminFixPointIterations) then
                ! In this case, the residuum of the last iterate decides on 
                ! whether to use Newton or not.
                dresInit = rsimsolver%rnonlinearIteration%dinitialDefectTotal
                dres = rsimsolver%rnonlinearIteration%dfinalDefectTotal
                if ((dres .lt. p_radaptiveNewton%depsAbsNewton) .and. &
                    (dres .lt. p_radaptiveNewton%depsRelNewton * dresInit)) then
                  rspaceTimeMatrixPrecond%cmatrixType = 1
                end if
              end if
              
              ! Otherwise: Use fixpoint iteration...
              
              ! Do we have to apply the inexact Newton?
              if (p_radaptiveNewton%cinexactNewton .ne. 0) then
              
                ! Adaptive stopping criterion / inexact Newton active.
                !
                ! Determine the stopping criterion for the linear solver.
                ! This is an adaptive stopping criterion depending on the current
                ! defect. In detail, for the inexact Newton we have
                !
                !   |b-Ax_{i+1}|         ( |b-Ax_i| ) exp             ( |b-Ax_i| )
                !   ------------ = min { ( -------- )     , depsrel * ( -------- ) }
                !     |b-Ax_0|           ( |b-Ax_0| )                 ( |b-Ax_0| )
                !
                ! see e.g. [Michael Hinze, Habilitation, p. 51]
                !
                ! If Newton is not active, we taje the formula
                !
                !   |b-Ax_{i+1}|             ( |b-Ax_i| )
                !   ------------ = depsrel * ( -------- ) 
                !     |b-Ax_0|               ( |b-Ax_0| )
                !
                ! to always gain depsrel.
                ! Switch off the relative stopping criterion in the linear solver:
                
                rsimSolver%rpreconditioner%p_rsolverNode%istoppingCriterion = 0
                
                ! Just for safetyness, gain at least one digit.
                rsimSolver%rpreconditioner%p_rsolverNode%depsRel = 1.0E-1_DP
                
                ! Calculate the new absolute stopping criterion:
                dresInit = rsimsolver%rnonlinearIteration%dinitialDefectTotal
                dres = rsimsolver%rnonlinearIteration%dfinalDefectTotal
                
                dtempdef = dres / dresInit
                
                if (rspaceTimeMatrixPrecond%cmatrixType .eq. 1) then
                  rsimSolver%rpreconditioner%p_rsolverNode%depsAbs = &
                      MIN(dtempDef**p_radaptiveNewton%dinexactNewtonExponent, &
                          p_radaptiveNewton%dinexactNewtonEpsRel*dtempdef) * dresInit
                else      
                  rsimSolver%rpreconditioner%p_rsolverNode%depsAbs = &
                      p_radaptiveNewton%dinexactNewtonEpsRel*dtempdef*dresInit
                end if
                
                ! If we have a multigrid solver, we also have to take care for
                ! the coarse grid solver!
                if (associated(rsimSolver%rpreconditioner%p_rcoarseGridSolverNode)) then
                  ! For the coarse grid solver, we choose the same stopping criterion.
                  ! But just for safetyness, the coarse grid solver should gain at least
                  ! one digit!
                  rsimSolver%rpreconditioner%p_rcoarseGridSolverNode%istoppingCriterion = 0
                  rsimSolver%rpreconditioner%p_rcoarseGridSolverNode%depsRel = 1.0E-1_DP
                  rsimSolver%rpreconditioner%p_rcoarseGridSolverNode%depsAbs = &
                      rsimSolver%rpreconditioner%p_rsolverNode%depsAbs
                end if
                
              end if
              
            end if
            
          end select
              
          call cc_setupMatrixWeights (rsimsolver%p_rproblem,rspaceTimeMatrixPrecond,dtheta,&
            istep-1,0,rnonlinearSpatialMatrix)

          call cc_disableSubmatrix (rnonlinearSpatialMatrix,2,1)
          call cc_disableSubmatrix (rnonlinearSpatialMatrix,2,2)
          call cc_disableSubmatrix (rnonlinearSpatialMatrix,1,2)
          
          ! Assemble the preconditioner matrices on all levels.
          call fbsim_assemblePrecMatrices (rsimsolver%rpreconditioner,rnonlinearSpatialMatrix,&
              roseensol1,roseensol2,roseensol3,.true.)
              
          ! Do preconditioning to the defect
          call fbsim_precondSpaceDefect (rsimsolver%rpreconditioner,rdefectPrimal,bsuccess)
          
          ! 4.) Update of the solution
          ! --------------------------
          
          if (bsuccess) then
            ! Combine to get the new solution vector.
            call lsyssc_vectorLinearComb(&
                rdefectPrimal%RvectorBlock(1),rcurrentsol%RvectorBlock(1),1.0_DP,1.0_DP)
            call lsyssc_vectorLinearComb(&
                rdefectPrimal%RvectorBlock(2),rcurrentsol%RvectorBlock(2),1.0_DP,1.0_DP)
            call lsyssc_vectorLinearComb(&
                rdefectPrimal%RvectorBlock(3),rcurrentsol%RvectorBlock(3),1.0_DP,1.0_DP)
                
            ! This is also the new Oseen solution
            call lsysbl_copyVector (rcurrentsol,roseenSol2)
          else
            exit
          end if

        end do ! ite

        ! Nonlinear loop finished. Now save the solution or stop incase of an error.

        if (bsuccess) then
          ! Save the solution.
          call sptivec_setTimestepData (rsolvector, istep, rcurrentsol)
        else  
          exit
        end if
      
      end do ! istep
      
    case (1)
    
      ! Oseen forward simulation in the primal space.
      
      ! Assign BC`s to the defect vector.
      call lsysbl_assignDiscreteBC (rdefectPrimal,&
          rsimSolver%rpreconditioner%p_RdiscreteBC(rsimSolver%ilevel))
      call lsysbl_assignDiscreteFBC (rdefectPrimal,&
          rsimSolver%rpreconditioner%p_RdiscreteFBC(rsimSolver%ilevel))
      
      ! Loop through the timesteps.
      do istep = ifirststep,ilaststep
      
        ! Current time step?
        dtime = &
            p_rspaceTimeMatrix%p_rspaceTimeDiscr%rtimeDiscr%dtimeInit + (istep-1) * dtstep

        ! Update the boundary conditions to the current time
        call fbsim_updateDiscreteBCprec (rsimsolver%p_rproblem,rsimsolver%rpreconditioner,dtime)

        if (istep .eq. ifirststep) then
        
          ! Initial condition is given.
          call sptivec_getTimestepData (rrhsvector, ifirststep, rcurrentrhs)
          call sptivec_getTimestepData (rsolvector, ifirststep, rcurrentsol)
          call lsysbl_clearVector (roseensol1)
          call sptivec_getTimestepData (roseenSolution, ifirststep, roseensol2)
          call sptivec_getTimestepData (roseenSolution, ifirststep+1, roseensol3)

        else

          ! Initialise the previous and current timestep. 
          ! Solution.
          call lsysbl_copyVector (rcurrentsol,rprevsol)
          call sptivec_getTimestepData (rsolvector, istep, rcurrentsol)

          ! RHS.
          call lsysbl_copyVector (rcurrentrhs,rprevrhs)
          call sptivec_getTimestepData (rrhsvector, istep, rcurrentrhs)
          
          ! Get the Oseen solutions for the assembly of the nonlinearity -- 
          ! if there is any.
          call lsysbl_copyVector (roseensol2,roseensol1)
          call lsysbl_copyVector (roseensol3,roseensol2)
          if (istep .lt. roseenSolution%NEQtime) then
            call sptivec_getTimestepData (roseenSolution, istep+1, roseensol3)
          end if
        
        end if
        
        ! Create the RHS.
        call lsysbl_copyVector (rcurrentrhs,rdefect)

        ! If we are not in the first step...
        if (istep .ne. ifirststep) then
          ! ... subtract the offdiagonal, corresponding to the previous timestep.
          call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rsimsolver%p_rproblem,&
              p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
              p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
          call cc_setupMatrixWeights (rsimsolver%p_rproblem,p_rspaceTimeMatrix,dtheta,&
            istep-1,-1,rnonlinearSpatialMatrix)
          
          call cc_assembleDefect (rnonlinearSpatialMatrix,rprevsol,&
            rdefect,1.0_DP,roseensol1,roseensol2,roseensol3)
        end if

        ! Subtract A*primal/dual solution from the rhs to create the primal defect.
        call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rsimsolver%p_rproblem,&
            p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
            p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
        call cc_setupMatrixWeights (rsimsolver%p_rproblem,p_rspaceTimeMatrix,dtheta,&
          istep-1,0,rnonlinearSpatialMatrix)
        call cc_disableSubmatrix (rnonlinearSpatialMatrix,2,1)
        call cc_disableSubmatrix (rnonlinearSpatialMatrix,2,2)
        
        call cc_assembleDefect (rnonlinearSpatialMatrix,rcurrentsol,&
          rdefect,1.0_DP,roseensol1,roseensol2,roseensol3)

        ! Implement the boundary conditions        
        call vecfil_discreteBCdef(rdefectPrimal)
        
        ! Disable everything except the primal matrix and
        ! use the remaining weights to create the matrix for the preconditioner
        ! on every level.
        call cc_disableSubmatrix (rnonlinearSpatialMatrix,1,2)
        
        ! Assemble the preconditioner matrices on all levels.
        call fbsim_assemblePrecMatrices (rsimsolver%rpreconditioner,rnonlinearSpatialMatrix,&
            roseensol1,roseensol2,roseensol3,.true.)
            
        ! Do preconditioning to the defect
        call fbsim_precondSpaceDefect (rsimsolver%rpreconditioner,rdefectPrimal,bsuccess)
        
        if (bsuccess) then
          ! Combine to get the new solution vector.
          call lsyssc_vectorLinearComb(&
              rdefectPrimal%RvectorBlock(1),rcurrentsol%RvectorBlock(1),1.0_DP,1.0_DP)
          call lsyssc_vectorLinearComb(&
              rdefectPrimal%RvectorBlock(2),rcurrentsol%RvectorBlock(2),1.0_DP,1.0_DP)
          call lsyssc_vectorLinearComb(&
              rdefectPrimal%RvectorBlock(3),rcurrentsol%RvectorBlock(3),1.0_DP,1.0_DP)
          call sptivec_setTimestepData (rsolvector, istep, rcurrentsol)
        else
          exit
        end if
      
      end do ! istep
    
    case (2)
    
      ! Oseen backward simulation.
      
      ! Assign BC`s to the defect vector.
      call lsysbl_assignDiscreteBC (rdefectDual,&
          rsimSolver%rpreconditioner%p_RdiscreteBC(rsimSolver%ilevel))
      call lsysbl_assignDiscreteFBC (rdefectDual,&
          rsimSolver%rpreconditioner%p_RdiscreteFBC(rsimSolver%ilevel))

      do istep = ilaststep,ifirststep,-1
      
        ! Current time step?
        dtime = &
            p_rspaceTimeMatrix%p_rspaceTimeDiscr%rtimeDiscr%dtimeInit + (istep-1) * dtstep

        ! Update the boundary conditions to the current time
        call fbsim_updateDiscreteBCprec (rsimsolver%p_rproblem,rsimsolver%rpreconditioner,dtime)
        
        if (istep .eq. ilastStep) then
        
          ! Terminal condition is given.
          call sptivec_getTimestepData (rrhsvector, ilaststep, rcurrentrhs)
          call sptivec_getTimestepData (rsolvector, ilaststep, rcurrentsol)
          call sptivec_getTimestepData (roseenSolution, ilaststep, roseensol2)
          call sptivec_getTimestepData (roseenSolution, ilaststep-1, roseensol1)

        else

          ! Initialise the previous and current timestep. 
          ! Solution.
          call lsysbl_copyVector (rcurrentsol,rprevsol)
          call sptivec_getTimestepData (rsolvector, istep, rcurrentsol)

          ! RHS.
          call lsysbl_copyVector (rcurrentrhs,rprevrhs)
          call sptivec_getTimestepData (rrhsvector, istep, rcurrentrhs)
          
          ! Get the Oseen solutions for the assembly of the nonlinearity -- 
          ! if there is any.
          call lsysbl_copyVector (roseensol2,roseensol3)
          call lsysbl_copyVector (roseensol1,roseensol2)
          if (istep .ge. 2) then
            call sptivec_getTimestepData (roseenSolution, istep-1, roseensol1)
          end if
          
        end if
        
        ! Create the RHS.
        call lsysbl_copyVector (rcurrentrhs,rdefect)

        ! If we are not in the last step...
        if (istep .ne. ilaststep) then
          ! ... subtract the offdiagonal, corresponding to the next timestep.
          call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rsimsolver%p_rproblem,&
              p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
              p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
          call cc_setupMatrixWeights (rsimsolver%p_rproblem,p_rspaceTimeMatrix,dtheta,&
            istep-1,1,rnonlinearSpatialMatrix)
          
          call cc_assembleDefect (rnonlinearSpatialMatrix,rprevsol,&
            rdefect,1.0_DP,roseensol1,roseensol2,roseensol3)
        end if

        ! Subtract A*primal/dual solution from the rhs to create the primal defect.
        call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rsimsolver%p_rproblem,&
            p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
            p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
        call cc_setupMatrixWeights (rsimsolver%p_rproblem,p_rspaceTimeMatrix,dtheta,&
          istep-1,0,rnonlinearSpatialMatrix)
        call cc_disableSubmatrix (rnonlinearSpatialMatrix,1,1)
        call cc_disableSubmatrix (rnonlinearSpatialMatrix,1,2)
        
        call cc_assembleDefect (rnonlinearSpatialMatrix,rcurrentsol,&
          rdefect,1.0_DP,roseensol1,roseensol2,roseensol3)
        
        ! Implement the boundary conditions        
        call vecfil_discreteBCdef(rdefectDual)

        ! Disable everything except the dual matrix and
        ! use the remaining weights to create the matrix for the preconditioner
        ! on every level.
        call cc_disableSubmatrix (rnonlinearSpatialMatrix,2,1)
        
        ! Assemble the preconditioner matrices on all levels.
        call fbsim_assemblePrecMatrices (rsimsolver%rpreconditioner,rnonlinearSpatialMatrix,&
            roseensol1,roseensol2,roseensol3,.true.)
            
        ! Do preconditioning to the defect
        call fbsim_precondSpaceDefect (rsimsolver%rpreconditioner,rdefectDual,bsuccess)
        
        if (bsuccess) then
          ! Combine to get the new solution vector.
          call lsyssc_vectorLinearComb(&
              rdefectDual%RvectorBlock(1),rcurrentsol%RvectorBlock(4),1.0_DP,1.0_DP)
          call lsyssc_vectorLinearComb(&
              rdefectDual%RvectorBlock(2),rcurrentsol%RvectorBlock(5),1.0_DP,1.0_DP)
          call lsyssc_vectorLinearComb(&
              rdefectDual%RvectorBlock(3),rcurrentsol%RvectorBlock(6),1.0_DP,1.0_DP)
          call sptivec_setTimestepData (rsolvector, istep, rcurrentsol)
        else
          exit
        end if
      
      end do ! istep
      
    end select
    
    ! Release the temp vectors
    call lsysbl_releaseVector (roseensol3)
    call lsysbl_releaseVector (roseensol2)
    call lsysbl_releaseVector (roseensol1)
    call lsysbl_releaseVector (rcurrentsol)
    call lsysbl_releaseVector (rprevsol)

  end subroutine

end module
