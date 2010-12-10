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
!#
!# Auxiliary routines:
!#
!# 1.) fbsim_initPreconditioner
!#     -> Initialise a spatial preconditioner
!#
!# 2.) fbsim_updateDiscreteBCprec
!#     -> Update boundary conditions in a spatial preconditioner
!#
!# 3.) fbsim_donePreconditioner
!#     -> Releasea spatial preconditioner
!#
!# 4.) fbsim_assemblePrecMatrices
!#     -> Assemble matrices in a spatial preconditioner
!#
!# 5.) fbsim_getNLsimParameters
!#     -> Reads parameters for the nonlinear solver from a DAT file structure
!#
!# 6.) fbsim_precondSpaceDefect
!#     -> Preconditioning of a nonlinear defect vector
!#
!# 7.) fbsim_getDefectNorm
!#     -> Calculate a nonlinear defect norm
!#
!# 8.) fbsim_resNormCheck
!#     -> Check if a nonlinear iteration converged
!# </purpose>
!##############################################################################

module forwardbackwardsimulation

    !filter chain?
    !solvers appropriate for both, primal and primal/dual problems?
    !if yes -> replaces spacetimepreconditioner in the spacetimelinearsolver!

  use fsystem
  use genoutput
  use basicgeometry
  use element
  use cubature
  use triangulation
  use linearsystemscalar
  use linearsystemblock
  use filtersupport
  use multilevelprojection
  use linearsolver
  use nonlinearsolver
  use paramlist
  use discretebc
  use discretefbc
  use linearalgebra
  use bcassembly
  use spatialdiscretisation
  use coarsegridcorrection
  use linearsolverautoinitialise
  use vectorfilters
  use matrixfilters
  use multilevelprojection
  use matrixmodification
  use ucd
  use spdiscprojection
  use collection
  use statistics
  
  use fespacehierarchybase
  use fespacehierarchy

  ! Include main structures; used for matrix assembly!
  use constantsoptc
  use assemblytemplates
  use assemblytemplatesoptc
  use structuresoptc
  
  use structuresoptflow
  
  use timediscretisation
  use spatialbcdef
  use spacematvecassembly
  use spacetimevectors
  use spacetimelinearsystem
  
  use user_callback

  implicit none
  
  private

!<constants>

!<constantblock description="ID's that define the type of the solver.">

  ! Nonlinear forward simulation.
  integer, parameter, public :: FBSIM_SOLVER_NLFORWARD = 0

  ! Linear forward simulation with prescribed nonlinearity.
  integer, parameter, public :: FBSIM_SOLVER_LINFORWARD = 1
  
  ! Linear backward simulation with prescribed nonlinearity.
  integer, parameter, public :: FBSIM_SOLVER_LINBACKWARD = 2

!</constantblock>

!</constants>

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
  
    ! Minimum number of usul fix point iteration before to switch to
    ! preconfitioning with the Newton matrix. (IFIXMIN)

    integer :: nminFixPointIterations = 0

    ! Maximum number of usul fix point iteration before to switch to
    ! preconfitioning with the Newton matrix. (IFIXMAX)

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
  type t_fbsimNonlinearIteration
  
    ! Norm of initial residuum for the complete solution vector.
    real(DP) :: dinitialDefectTotal = 0.0_DP

    ! Norm of final residuum for the complete solution vector.
    real(DP) :: dfinalDefectTotal = 0.0_DP

    ! Norm of initial residuum for the complete solution vector
    ! before the last preconditioning.
    real(DP) :: dlastInitPrecDef = 0.0_DP

    ! Norm of initial residuum for the complete solution vector
    ! fter the last preconditioning.
    real(DP) :: dlastPrecDef = 0.0_DP

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
    ! =-1: undefined
    ! =0: simple linear solver
    ! =1: nonlinear defect correction solver.
    ! =2: Newton solver
    ! =3: adaptive Newton with parameters in radaptiveNewton
    integer :: ctypeIteration = -1
    
    ! Output mode. Used for printing messages.
    ! =OU_MODE_STD: Print messages to the terminal and probably to a log 
    ! file (if a log file is opened).
    integer(I32) :: coutputmode = OU_MODE_STD

    ! Parameters for the adaptive Newton iteration.
    type(t_ccDynamicNewtonControl) :: radaptiveNewton
    
    ! <!-- ---------- -->
    ! <!-- STATISTICS -->
    ! <!-- ---------- -->
    
    ! Total number of nonlinear iterations
    integer :: nnonlinearIterations = 0
    
    ! Total number of linear iterations
    integer :: nlinearIterations = 0
    
    ! Total time for the solver
    real(DP) :: dtimeTotal = 0.0_DP
    
    ! Total time for nonlinear iteration
    real(DP) :: dtimeNonlinearSolver = 0.0_DP
    
    ! Total time for linear solver
    real(DP) :: dtimeLinearSolver = 0.0_DP
    
    ! Total time for matrix assembly
    real(DP) :: dtimeMatrixAssembly = 0.0_DP

    ! Total time for defect calculation
    real(DP) :: dtimeDefectCalculation = 0.0_DP

    ! Total time for matrix assembly
    real(DP) :: dtimePostprocessing = 0.0_DP

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
    integer :: nlmin = 0
    
    ! Maximum refinement level
    integer :: nlmax = 0
    
    ! A pointer to a hierarchy of space assembly template structures.
    type(t_staticSpaceAsmHierarchy), pointer :: p_rspaceAsmTemplHier => null()

    ! A pointer to a hierarchy of space assembly template structures.
    type(t_staticSpaceAsmHierarchyOptC), pointer :: p_rspaceAsmTemplHierOptC => null()
    
    ! A solver node that accepts parameters for the linear solver    
    type(t_linsolNode), pointer :: p_rsolverNode => null()

    ! If the linear solver contains a multigrid subsolver, this is a reference
    ! to the coarse grid solver.
    type(t_linsolNode), pointer :: p_rcoarseGridSolverNode => null()

    ! A filter chain that is used for implementing boundary conditions or other
    ! things when invoking the linear solver.
    type(t_filterChain), dimension(6) :: RfilterChain

    ! A hierarchy of space levels for the space that is used.
    type(t_feHierarchy), pointer :: p_rfeHierarchy => null()
    
    ! A time dscretisation structure that defines the weights in time.
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr => null()

    ! Pointer to interlevel projection hierarchy.
    type(t_interlevelProjectionHier), pointer :: p_rprjHierarchy => null()
    
    ! An array of matrices used for preconditioning in space.
    ! These matrices are 6x6 matrices for both, primal and dual space.
    ! The actual preconditioner matrices for either the primal or the
    ! dual space can be found in p_RmatrixPrecond and are shares submatrices
    ! of these matrices.
    type(t_matrixBlock), dimension(:), pointer :: p_RmatrixPrecondFullSpace

    ! An array of matrices used for preconditioning in space,
    ! only for primal or dual space.
    type(t_matrixBlock), dimension(:), pointer :: p_RmatrixPrecond
    
    ! Boundary conditions to use.
    type(t_optcBDC), pointer :: p_rboundaryConditions => null()
    
    ! Discrete boundary conditions on all levels corresponding to the
    ! discretisation of p_RmatrixPrecond.
    type(t_discreteBC), dimension(:), pointer :: p_RdiscreteBC => null()
    
    ! Discrete fictitious BC`s on all levels corresponding to the
    ! discretisation of p_RmatrixPrecond.
    type(t_discreteFBC), dimension(:), pointer :: p_RdiscreteFBC => null()
    
    ! An array of 3x#levels temp vectors.
    ! Note: The temp vectors on the maximum level are created but
    ! actually not used. The application can use them to save the
    ! evaluation point of the nonlinearity there!
    type(t_vectorBlock), dimension(:,:), pointer :: p_RtempVec
    
    ! Assembly data on all levels.
    type(t_spatialMatrixDiscrData), dimension(:), pointer :: p_RassemblyData
    
    ! Whether there are Neumann boundary components in the BC`s.
    logical :: bhasNeumann
    
    !<!-- Statistical data -->
    
    ! Returns the number of iterations, a linear subsolver needed 
    ! for preconditioning
    integer :: nlinearIterations = 0
    
    ! Time needed for solving linear systems.
    real(DP) :: dtimeLinearSolver = 0.0_DP
    
  end type

!</typeblock>

!<typeblock>
  ! Configures the postprocessing during the iteration.
  type t_fbsimPostprocessing

    ! <!-- Input parameters -->
  
    ! Type of output file to generate from solutions.
    ! 0=disabled
    ! 1=GMV
    ! 2=AVS
    ! 3=Paraview (VTK)
    ! 4=Matlab
    integer :: ioutputUCD = 0
    
    ! Filename for UCD output.
    ! A timestep index '.0001','.0002',... is appended to this.
    character(len=SYS_STRLEN) :: sfilenameUCD = "./gmv/u"

    ! Write solution vector to a sequence of files.
    ! =0: don't write
    ! =1: write out, use formatted output.
    ! =2: write out, use unformatted output.

    integer :: cwriteSolution = 0

    ! Name/path of file for writing solution

    character(len=SYS_STRLEN) :: swriteSolutionFilename = ""
  
    ! <!-- the following parameters are automatically maintained during a simulation -->
    
    ! Space that is available in rsolution. One of the CCSPACE_xxxx constants.
    integer :: cspace = CCSPACE_PRIMAL
    
    ! Underlying space discretisation
    type(t_blockDiscretisation), pointer :: p_rspaceDiscr

    ! Underlying time discretisation
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr
    
    ! Discrete boundary conditions
    type(t_discreteBC) :: rdiscreteBC
    
    ! Discrete fictitious boundary conditions
    type(t_discreteFBC) :: rdiscreteFBC
  
    ! Boundary conditions to use.
    type(t_optcBDC), pointer :: p_rboundaryConditions => null()
    
    ! User defined parameter list; if not NULL(), this is added
    ! as comments to postprocessing files.
    type(t_parlist), pointer :: p_rparlist => null()

  end type
!</typeblock>

  public :: t_fbsimPostprocessing

!<typeblock>

  ! Configuration block for the simulation solver.
  type t_simSolver
  
    ! Type of simulation solver. An FBSIM_SOLVER_xxxx constant.
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
    integer :: csimtype = FBSIM_SOLVER_NLFORWARD

    ! Assembly data on the level where to solve.
    type(t_spatialMatrixDiscrData) :: rdiscrData

    ! Absolute level where the simulation is executed.
    integer :: ilevel
    
    ! Output level of the solver.
    ! =-1: no output
    ! =0: no output, only errors
    ! =1: basic output
    ! =2: standard output
    integer :: ioutputLevel = 2
  
    ! A t_fbsimPreconditioner structure that defines the preconditioner
    type(t_fbsimPreconditioner) :: rpreconditioner

    ! A t_fbsimPreconditioner structure that defines an alternative preconditioner
    type(t_fbsimPreconditioner) :: rpreconditionerAlternative

    ! Pointer to the underlying space-time matrix that defines the timesteping etc.
    type(t_ccoptSpaceTimeMatrix), pointer :: p_rmatrix => null()

    ! Level-info structure of the level where the simulation is executed.
    type(t_staticSpaceAsmTemplates), pointer :: p_rspaceAsmTemplHier => null()
    
    ! Boundary conditions to use.
    type(t_optcBDC), pointer :: p_rboundaryConditions => null()
    
    ! Global program settings
    type(t_settings_optflow), pointer :: p_rsettings

    ! Parameters for the nonlinear iteration during a pure forward simulation.
    type(t_fbsimNonlinearIteration) :: rnonlinearIteration

    ! Postprocessing settings
    type(t_fbsimPostprocessing) :: rpostprocessing

    ! User defined parameter list; if not NULL(), this is added
    ! as comments to postprocessing files.
    type(t_parlist), pointer :: p_rparlist => null()
    
    ! <!-- ---------- -->
    ! <!-- STATISTICS -->
    ! <!-- ---------- -->
    
    ! Total number of nonlinear iterations
    integer :: nnonlinearIterations = 0
    
    ! Total number of linear iterations
    integer :: nlinearIterations = 0
    
    ! Total time for the solver
    real(DP) :: dtimeTotal = 0.0_DP
    
    ! Total time for nonlinear iteration
    real(DP) :: dtimeNonlinearSolver = 0.0_DP
    
    ! Total time for linear solver
    real(DP) :: dtimeLinearSolver = 0.0_DP
    
    ! Total time for matrix assembly
    real(DP) :: dtimeMatrixAssembly = 0.0_DP

    ! Total time for defect calculation
    real(DP) :: dtimeDefectCalculation = 0.0_DP

    ! Total time for matrix assembly
    real(DP) :: dtimePostprocessing = 0.0_DP

  end type

!</typeblock>

!</types>

  public :: t_ccDynamicNewtonControl
  public :: t_fbsimNonlinearIteration
  public :: t_fbsimPreconditioner
  public :: t_simSolver
  
  public :: fbsim_getNLsimParamsLinSol
  public :: fbsim_getNLsimParameters
  public :: fbsim_setMatrix

  public :: fbsim_init
  public :: fbsim_simulate
  public :: fbsim_done
  
  public :: fbsim_initPreconditioner
  public :: fbsim_updateDiscreteBCprec
  public :: fbsim_precondSpaceDefect
  public :: fbsim_donePreconditioner
  
  public :: fbsim_assemblePrecMatrices
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine fbsim_getNLsimParamsLinSol (rnonlinearIteration)
  
!<description>
  ! Initialises the parameters in rnonlinearIteration with dummy parameters
  ! to reflect a simple linear solver without nonlinear iteration.
  ! The min/max number of iterations are set to 1, stopping criteria
  ! are switched off.
  !
  ! The resulting structure can be used for forward/backward simulations
  ! which do not involve a nonlinear iteration in each timestep.
  ! For simulations that need a nonlinear iteration in each timestep,
  ! the user can call fbsim_getNLsimParameters to initialise the structure
  ! based on a parameter list for a real nonlinear iteration.
!</description>

!<output>
  ! Nonlinar iteration structure saving data for the callback routines.
  ! Is filled with data.
  type(t_fbsimNonlinearIteration), intent(out) :: rnonlinearIteration
!</output>

!</subroutine>

    ! Clear auxiliary variables for the nonlinear iteration
    rnonlinearIteration%DresidualInit = 0.0_DP
    rnonlinearIteration%DresidualOld  = 0.0_DP
    
    ! Set number of iterations to 1.
    rnonlinearIteration%nminIterations = 1
    rnonlinearIteration%nmaxIterations = 1

    ! We write out the data of the nonlinear solver to the benchmark
    ! log file as well.
    rnonlinearIteration%coutputMode = OU_MODE_STD+OU_MODE_BENCHLOG
    
    ! The solver is a simple linear solver without a nonlinear iteration
    rnonlinearIteration%ctypeIteration = 0

  end subroutine

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
        'spreconditionerAdaptiveNewton', sstring, "",bdequote=.true.)
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

  subroutine fbsim_init (rsettings, rparlist, ssection, nlmin, nlmax, csimtype, rsimsolver,&
      rphysics)
  
!<description>
  ! Initialises a rsimsolver structure based on the general problem structure
  ! for performing forward/backward sweeps.
!</description>
  
!<input>
  ! The structure of the main solver
  type(t_settings_optflow), intent(inout), target :: rsettings
  
  ! Param,eter list containing the parameters for the solver.
  type(t_parlist), intent(in) :: rparlist
  
  ! Name of the section containing the parameters of the solver in
  ! each timestep. This section identifies either a linear or a nonlinear
  ! solver, depending on csimtype.
  character(len=*), intent(in) :: ssection

  ! Minimum space level allowed to be used by the preconditioner.
  ! Relative to the global space-hierarchy in rsettings.
  integer, intent(in) :: nlmin

  ! Level where the simulation is carried out.
  ! Relative to the global space-hierarchy in rsettings.
  integer, intent(in) :: nlmax

  ! Type of simulation solver.
  ! =0: Navier-Stokes forward solver on the primal solution vector.
  !     The dual solution influences the simulation.
  !     ssection must identify a parameter block of a nonlinear solver.
  ! =1: Oseen forward solver on the primal solution vector.
  !     ssection must identify a parameter block of a linear solver.
  ! =2: Oseen backward solver simulation on the dual solution vector.
  !     ssection must identify a parameter block of a linear solver.
  integer, intent(in) :: csimtype
  
  ! OPTIONAL: Alternative physics definition to use.
  ! If not present, the standard global physics settings are used.
  type(t_settings_physics), intent(in), optional :: rphysics
!</input>

!<output>
  ! A t_simSolver node with parameters for the iteration.
  type(t_simSolver), intent(out) :: rsimsolver 
!</output>

!</subroutine>

    integer :: cspace
    character(len=SYS_STRLEN) :: ssectionLinSol,ssectionLinSol2

    ! Reset statistrics
    rsimsolver%nnonlinearIterations = 0
    rsimsolver%nlinearIterations = 0
    rsimsolver%dtimeNonlinearSolver = 0.0_DP
    rsimsolver%dtimeLinearSolver = 0.0_DP
    rsimsolver%dtimeMatrixAssembly = 0.0_DP

    ! Take the information from the problem structure.
    rsimsolver%ilevel = nlmax
    rsimsolver%p_rspaceAsmTemplHier => rsettings%rspaceAsmHierarchy%p_RasmTemplList(nlmax)
    rsimsolver%p_rsettings => rsettings
    
    rsimsolver%csimtype = csimtype
    
    select case (csimtype)
    case (FBSIM_SOLVER_NLFORWARD)
      cspace = CCSPACE_PRIMAL
      
      ! Get the parameters of the nonlinear solver 
      call fbsim_getNLsimParameters (rparlist,ssection, &  !"CC-NONLINEAR",&
          rsimsolver%rnonlinearIteration)

      ! Get the section defining the linear solver.
      call parlst_getvalue_string (rparlist, ssection,"slinearSolver", &
          ssectionLinSol, "CC-LINEARSOLVER",bdequote=.true.)

      ! Get the section defining the alternative linear solver.
      call parlst_getvalue_string (rparlist, ssection,"slinearSolverAlternative", &
          ssectionLinSol2, "",bdequote=.true.)

    case (FBSIM_SOLVER_LINFORWARD)
      cspace = CCSPACE_PRIMAL
      
      ! Initialise the nonlinear solver structure as dummy.
      call fbsim_getNLsimParamsLinSol(rsimsolver%rnonlinearIteration)
      
      ssectionLinSol = ssection
      
    case (FBSIM_SOLVER_LINBACKWARD)
      cspace = CCSPACE_DUAL

      ! Initialise the nonlinear solver structure as dummy.
      call fbsim_getNLsimParamsLinSol(rsimsolver%rnonlinearIteration)

      ssectionLinSol = ssection
      
    end select
    
    ! Initialise the preconditioner.
    ! Do not initialise the time discretisation; this will be initialised
    ! during setmatix!
    call fbsim_initPreconditioner (rsimsolver%rpreconditioner, rsettings, &
        nlmin, nlmax, cspace, rparamlist=rparlist, ssection=ssectionLinSol)
        
    ! Probably initialise the alternative preconditioner.
    if (ssectionLinSol2 .ne. "") then
      call fbsim_initPreconditioner (rsimsolver%rpreconditionerAlternative, rsettings, &
          nlmin, nlmax, cspace, rparamlist=rparlist, ssection=ssectionLinSol2)
    end if
        
    ! Get the assembly data on our level that allows us to create nonlinear
    ! matrices.
    call smva_getDiscrData (rsettings, nlmax, rsimsolver%rdiscrData, rphysics)
        
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

    ! Release the preconditioner.
    call fbsim_donePreconditioner (rsimsolver%rpreconditioner)
    
    if (rsimsolver%rpreconditionerAlternative%ctypePreconditioner .ne. 0) then
      ! Release also the alternative preconditioner.
      call fbsim_donePreconditioner (rsimsolver%rpreconditionerAlternative)
    end if

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
    
    ! Remember a pointer to the underlying time discretisation.
    rsimSolver%rpreconditioner%p_rtimeDiscr => rmatrix%rdiscrData%p_rtimeDiscr

    if (rsimsolver%rpreconditionerAlternative%ctypePreconditioner .ne. 0) then
      ! THe same for the alternative preconditioner
      rsimSolver%rpreconditionerAlternative%p_rtimeDiscr => rmatrix%rdiscrData%p_rtimeDiscr
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fbsim_getMatrixAssemblyFlags (rpreconditioner,ilev,rassemblyFlags)

!<description>
  ! Create a standard matrix assembly structure based on a preconditioner structure.
!</description>

!<input>
  ! Preconditioner structure.
  type(t_fbsimPreconditioner), intent(in) :: rpreconditioner
  
  ! Current space level.
  integer, intent(in) :: ilev
!</input>

!<inputoutput>
  ! Matrix assembly structure to initialise.
  type(t_matrixAssemblyFlags), intent(out) :: rassemblyFlags
!</inputoutput>

!</subroutine>

    ! local variables

    ! Prepare an assembly flags structure for the assembly.
    rassemblyFlags%iadaptiveMatrices = rpreconditioner%iadaptiveMatrices
    rassemblyFlags%dadmatthreshold = rpreconditioner%dadmatthreshold
    if (ilev .eq. rpreconditioner%nlmin) then
      rassemblyFlags%bvirtualTransposedD = rpreconditioner%bneedVirtTransposedDonCoarse
    else
      rassemblyFlags%bvirtualTransposedD = rpreconditioner%bneedVirtTransposedD
    end if

  end subroutine
      
  ! ***************************************************************************

!<subroutine>

  subroutine fbsim_initPreconditioner (rpreconditioner, rsettings, &
      nlmin, nlmax, cspace, rtimeDiscr, rboundaryConditions, rparamlist, ssection)
  
!<description>
  ! Initialises the preconditioner for each timestep based on the parameters 
  ! in rparamlist.
!</description>
  
!<input>
  ! The general problem structure.
  type(t_settings_optflow), intent(inout), target :: rsettings

  ! Minimum available space refinement level.
  integer, intent(in) :: nlmin
  
  ! Maximum space refinement level; corresponds to the level where to do the 
  ! preconditioning.
  integer, intent(in) :: nlmax
  
  ! Solution space, the preconditioner is applied to.
  ! CCSPACE_PRIMAL: Primal space.
  ! CCSPACE_DUAL: Dual space.
  ! CCSPACE_PRIMALDUAL: Primal + dual space.
  integer, intent(in) :: cspace

  ! OPTIONAL: Time discretisation structure.
  ! If not specified, the caller must specify it later before
  ! the simulation!
  type(t_timeDiscretisation), intent(in), optional, target :: rtimeDiscr

  ! OPTIONAL: Boundary conditions
  ! If not specified, the caller must specify it later before
  ! the simulation!
  type(t_optcBDC), intent(in), optional, target :: rboundaryConditions

  ! OPTIONAL: Parameter list with the solver parameters.
  ! Of not specified, no solver is created. The structure rpreconditioner 
  ! can then only be used to create matrices but not to solve systems.
  type(t_parlist), intent(in), optional :: rparamlist
  
  ! OPTIONAL: Name of the section configuring the linear solver.
  ! Of not specified, no solver is created. The structure rpreconditioner 
  ! can then only be used to create matrices but not to solve systems.
  character(len=*), intent(in), optional :: ssection
  
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
    type(t_spatialMatrixNonlinearData) :: rnonlinearity
    type(t_matrixAssemblyFlags) :: rassemblyFlags
    
    ! -------------------------------------------------------------------------
    ! Part 1:
    ! Basic initialisation of all parameters and level-independent stuff.
    ! Reading of parameters from the DAT file.
    ! -------------------------------------------------------------------------

    ! Fetch level information where the preconditioner works.
    rpreconditioner%NLMIN = nlmin
    rpreconditioner%NLMAX = nlmax
    rpreconditioner%p_rspaceAsmTemplHier => rsettings%rspaceAsmHierarchy
    rpreconditioner%p_rspaceAsmTemplHierOptC => rsettings%rspaceAsmHierarchyOptC
    
    rpreconditioner%cspace = cspace
    
    ! Currently the only preconditioner available is a linear solver.
    rpreconditioner%ctypePreconditioner = 1

    ! Set up a filter that modifies the block vectors/matrix
    ! according to boundary conditions.
    select case (cspace)
    case (CCSPACE_PRIMAL, CCSPACE_DUAL)
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
      
    case (CCSPACE_PRIMALDUAL)
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

      ! Initialise the first filter of the filter chain as boundary
      ! implementation filter for defect vectors:
      rpreconditioner%RfilterChain(4)%ifilterType = &
          FILTER_DISCBCDEFREAL

      ! The second filter filters for boundary conditions of fictitious boundary
      ! components
      rpreconditioner%RfilterChain(5)%ifilterType = &
          FILTER_DISCBCDEFFICT
      
      ! The last element is by default DO-NOTHING but can be changed to an L2_0-filter
      ! if necessary.
      rpreconditioner%RfilterChain(6)%ifilterType = FILTER_DONOTHING
    end select
    rpreconditioner%bneedPressureDiagonalBlock = .false.
    rpreconditioner%bhasNeumann = .true.
    
    ! Depending on the space, initialise the projection hierarchy with primal
    ! and/or dual space parameters.
    select case (rpreconditioner%cspace)
    
      case (CCSPACE_PRIMAL,CCSPACE_DUAL)
        ! Only primal or dual space.
        rpreconditioner%p_rfeHierarchy => rsettings%rfeHierPrimal
            
        ! Get the interlevel projection structure for space prolongation/restriction
        rpreconditioner%p_rprjHierarchy => rsettings%rprjHierSpacePrimal

      case (CCSPACE_PRIMALDUAL) 
        ! Full primal/dual space
        rpreconditioner%p_rfeHierarchy => rsettings%rfeHierPrimalDual

        ! Get the interlevel projection structure for space prolongation/restriction
        rpreconditioner%p_rprjHierarchy => rsettings%rprjHierSpacePrimalDual
          
    end select
    
    ! Initialise the time discretisation if present.
    if (present(rtimeDiscr)) then
      rpreconditioner%p_rtimeDiscr => rtimeDiscr
    else
      nullify(rpreconditioner%p_rtimeDiscr)
    end if
    
    ! The same for the boundary conditions
    if (present(rboundaryConditions)) then
      rpreconditioner%p_rboundaryConditions => rboundaryConditions
    else
      nullify(rpreconditioner%p_rboundaryConditions)
    end if
    
    ! Get assembly data for nonlinear matrices on all levels.
    allocate(rpreconditioner%p_RassemblyData(rpreconditioner%NLMIN:rpreconditioner%NLMAX))
    do ilev=rpreconditioner%NLMIN,rpreconditioner%NLMAX
      call smva_getDiscrData (rsettings, ilev, &
          rpreconditioner%p_RassemblyData(ilev))
    end do
        
    ! Prepare boundary condition structures
    allocate(rpreconditioner%p_RdiscreteBC(rpreconditioner%NLMIN:rpreconditioner%NLMAX))
    allocate(rpreconditioner%p_RdiscreteFBC(rpreconditioner%NLMIN:rpreconditioner%NLMAX))
    do ilev=rpreconditioner%NLMIN,rpreconditioner%NLMAX
      call bcasm_initDiscreteBC(rpreconditioner%p_RdiscreteBC(ilev))
      call bcasm_initDiscreteFBC(rpreconditioner%p_RdiscreteFBC(ilev))
    end do
    
    if (present(rparamList) .and. present(ssection)) then
    
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

      call parlst_getvalue_string (p_rsection, 'ssolverSection', ssolverSection,&
          "",bdequote=.true.)
      call parlst_getvalue_string (p_rsection, 'ssmootherSection', ssmootherSection,&
          "",bdequote=.true.)
      call parlst_getvalue_string (p_rsection, 'scoarseGridSolverSection', &
          scoarseGridSolverSection,"",bdequote=.true.)
      
      ! Which type of solver do we have?
      
      select case (isolverType)
      
      case (0)
      
        ! This is the UMFPACK solver. Very easy to initialise. No parameters at all.
        call linsol_initUMFPACK4 (p_rsolverNode)
        !p_rsolverNode%p_rsubnodeUmfpack4%imatrixDebugOutput = 1
      
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

        select case (cspace)
        case (CCSPACE_PRIMAL, CCSPACE_DUAL)
          call cgcor_init(p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection,3)
          p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%p_DequationWeights(3) &
              = -1.0_DP
        case (CCSPACE_PRIMALDUAL)
          call cgcor_init(p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection,6)
          p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%p_DequationWeights(3) &
              = -1.0_DP
          p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%p_DequationWeights(6) &
              = -1.0_DP
        end select

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
          select case (cspace)
          case (CCSPACE_PRIMAL, CCSPACE_DUAL)
            call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVST)
            rpreconditioner%bneedVirtTransposedDonCoarse = .true.
          case (CCSPACE_PRIMALDUAL)
            call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG2)
            rpreconditioner%bneedVirtTransposedDonCoarse = .false.
          end select
          
          call parlst_getvalue_string (rparamList, scoarseGridSolverSection, &
              'spreconditionerSection', sstring, "",bdequote=.true.)
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
          !rpreconditioner%bneedVirtTransposedDonCoarse = .true.
          
        case (2)
          ! Defect correction with full VANKA preconditioning.
          !
          ! Create VANKA and initialise it with the parameters from the DAT file.
          select case (cspace)
          case (CCSPACE_PRIMAL, CCSPACE_DUAL)
            call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVST)
          case (CCSPACE_PRIMALDUAL)
            call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOC)
          end select
          
          call parlst_getvalue_string (rparamList, scoarseGridSolverSection, &
              'spreconditionerSection', sstring, "", bdequote=.true.)
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
          select case (cspace)
          case (CCSPACE_PRIMAL, CCSPACE_DUAL)
            call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DNAVST)
            rpreconditioner%bneedVirtTransposedDonCoarse = .true.
          case (CCSPACE_PRIMALDUAL)
            call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG2)
            rpreconditioner%bneedVirtTransposedDonCoarse = .false.
          end select
          
          call parlst_getvalue_string (rparamList, scoarseGridSolverSection, &
            'spreconditionerSection', sstring, "", bdequote=.true.)
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
          !rpreconditioner%bneedVirtTransposedDonCoarse = .true.

        case (4)
          ! BiCGStab with full VANKA preconditioning.
          !
          ! Create VANKA and initialise it with the parameters from the DAT file.
          select case (cspace)
          case (CCSPACE_PRIMAL, CCSPACE_DUAL)
            call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVST)
          case (CCSPACE_PRIMALDUAL)
            call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOC)
          end select
          
          call parlst_getvalue_string (rparamList, scoarseGridSolverSection, &
            'spreconditionerSection', sstring, "", bdequote=.true.)
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
          select case (cspace)
          case (CCSPACE_PRIMAL, CCSPACE_DUAL)
            call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_GENERAL)
          case (CCSPACE_PRIMALDUAL)
            call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_GENERAL)
          end select
          
          call parlst_getvalue_string (rparamList, scoarseGridSolverSection, &
            'spreconditionerSection', sstring, "", bdequote=.true.)
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
          select case (cspace)
          case (CCSPACE_PRIMAL, CCSPACE_DUAL)
            call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DNAVST)
          case (CCSPACE_PRIMALDUAL)
            call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG2)
          end select
          
          call parlst_getvalue_string (rparamList, scoarseGridSolverSection, &
            'spreconditionerSection', sstring, "", bdequote=.true.)
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
          
          case (0:10)

            nullify(p_rsmoother)
          
            ! This is some kind of VANKA smoother. Initialise the correct one.
            select case (ismootherType)
            case (0)
              call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_GENERAL)
            case (1)
              call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_GENERALDIRECT)
              
            case (2)
              select case (cspace)
              case (CCSPACE_PRIMAL, CCSPACE_DUAL)
                call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DNAVST)
                rpreconditioner%bneedVirtTransposedD = .true.
              case (CCSPACE_PRIMALDUAL)
                call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG2)
                rpreconditioner%bneedVirtTransposedD = .false.
              end select

              ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
              !rpreconditioner%bneedVirtTransposedD = .true.

            case (3)
              select case (cspace)
              case (CCSPACE_PRIMAL, CCSPACE_DUAL)
                call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DNAVSTDIRECT)
              case (CCSPACE_PRIMALDUAL)
                call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAGDIR)
              end select

              ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
              rpreconditioner%bneedVirtTransposedD = .true.

            case (4)
              select case (cspace)
              case (CCSPACE_PRIMAL, CCSPACE_DUAL)
                call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVST)
              case (CCSPACE_PRIMALDUAL)
                call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTOC)
              end select

              ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
              rpreconditioner%bneedVirtTransposedD = .true.

            case (5)
              select case (cspace)
              case (CCSPACE_PRIMAL, CCSPACE_DUAL)
                call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTDIRECT)
              case (CCSPACE_PRIMALDUAL)
                call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIRECT)
              end select

              ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
              rpreconditioner%bneedVirtTransposedD = .true.

            case (6)
              select case (cspace)
              case (CCSPACE_PRIMAL, CCSPACE_DUAL)
                call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVST)
              case (CCSPACE_PRIMALDUAL)
                call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOC)
              end select
              call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,&
                  rpreconditioner%RfilterChain)

              ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
              rpreconditioner%bneedVirtTransposedD = .true.

            case (7)
              select case (cspace)
              case (CCSPACE_PRIMAL, CCSPACE_DUAL)
                call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DNAVST)
                rpreconditioner%bneedVirtTransposedD = .true.
              case (CCSPACE_PRIMALDUAL)
                call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG2)
                rpreconditioner%bneedVirtTransposedD = .false.
              end select
              call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,&
                  rpreconditioner%RfilterChain)

              ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
              !rpreconditioner%bneedVirtTransposedD = .true.

            case (8)
              select case (cspace)
              case (CCSPACE_PRIMAL, CCSPACE_DUAL)
                call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DNAVST)
                rpreconditioner%bneedVirtTransposedD = .true.
              case (CCSPACE_PRIMALDUAL)
                call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG2)
                rpreconditioner%bneedVirtTransposedD = .false.
              end select
              call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,&
                  rpreconditioner%RfilterChain)

              ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
              !rpreconditioner%bneedVirtTransposedD = .true.

            case (9)
              select case (cspace)
              case (CCSPACE_PRIMAL, CCSPACE_DUAL)
                call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DNAVST)
                rpreconditioner%bneedVirtTransposedD = .true.
              case (CCSPACE_PRIMALDUAL)
                call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG2)
                rpreconditioner%bneedVirtTransposedD = .false.
              end select

              ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
              !rpreconditioner%bneedVirtTransposedD = .true.

            case (10)
              select case (cspace)
              case (CCSPACE_PRIMAL, CCSPACE_DUAL)
                call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_GENERAL)
                rpreconditioner%bneedVirtTransposedD = .false.
              case (CCSPACE_PRIMALDUAL)
                call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_GENERAL)
                rpreconditioner%bneedVirtTransposedD = .false.
              end select

              call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,&
                  rpreconditioner%RfilterChain)

              ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
              !rpreconditioner%bneedVirtTransposedD = .true.

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
                rpreconditioner%p_rprjHierarchy%p_Rprojection(ilev-nlmin+1))
            
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
            select case (cspace)
            case (CCSPACE_PRIMAL, CCSPACE_DUAL)
              call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVST)
            case (CCSPACE_PRIMALDUAL)
              call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTOC)
            end select

            ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
            rpreconditioner%bneedVirtTransposedD = .true.

          case (3)
            select case (cspace)
            case (CCSPACE_PRIMAL, CCSPACE_DUAL)
              call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTDIRECT)
            case (CCSPACE_PRIMALDUAL)
              call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIRECT)
            end select

            ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
            rpreconditioner%bneedVirtTransposedD = .true.

          case (4)
            select case (cspace)
            case (CCSPACE_PRIMAL, CCSPACE_DUAL)
              call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVST)
            case (CCSPACE_PRIMALDUAL)
              call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG)
            end select

            ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
            rpreconditioner%bneedVirtTransposedD = .true.

          case (5)
            select case (cspace)
            case (CCSPACE_PRIMAL, CCSPACE_DUAL)
              call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTDIRECT)
            case (CCSPACE_PRIMALDUAL)
              call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAGDIR)
            end select

            ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
            rpreconditioner%bneedVirtTransposedD = .true.

          case (6)
            select case (cspace)
            case (CCSPACE_PRIMAL, CCSPACE_DUAL)
              call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVST)
            case (CCSPACE_PRIMALDUAL)
              call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOC)
            end select
            
            call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,&
                rpreconditioner%RfilterChain)

            ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
            rpreconditioner%bneedVirtTransposedD = .true.

          case (7)
            select case (cspace)
            case (CCSPACE_PRIMAL, CCSPACE_DUAL)
              call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVST)
            case (CCSPACE_PRIMALDUAL)
              call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG)
            end select
            
            call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,&
                rpreconditioner%RfilterChain)

            ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
            rpreconditioner%bneedVirtTransposedD = .true.

          case (8)
            select case (cspace)
            case (CCSPACE_PRIMAL, CCSPACE_DUAL)
              call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVST)
            case (CCSPACE_PRIMALDUAL)
              call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG2)
            end select
            
            call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,&
                rpreconditioner%RfilterChain)

          case (9)
            select case (cspace)
            case (CCSPACE_PRIMAL, CCSPACE_DUAL)
              call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVST)
            case (CCSPACE_PRIMALDUAL)
              call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG2)
            end select

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
      
    else
    
      ! Otherwise: no linear solver.
      nullify(rpreconditioner%p_rsolverNode)
    
    end if
    
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

    ! Allocate memory on each level
    do ilev = rpreconditioner%NLMIN,rpreconditioner%NLMAX
      ! Get the matrix assembly flags.
      call fbsim_getMatrixAssemblyFlags (rpreconditioner,ilev,rassemblyFlags)
    
      ! Initialise the matrix. Specify a 'dummy' nonlinearity (empty rnonlinearity)
      ! because there is no nonlinearity during memory allocation.
      call smva_initNonlinMatrix (rnonlinearSpatialMatrix,&
          rpreconditioner%p_RassemblyData(ilev),rnonlinearity)

      ! Get a dummy structure for a full matrix.
      call stlin_getFullMatrixDummy(rsettings%rphysicsPrimal,rnonlinearSpatialMatrix)
      
        ! Probably deactivate the offdiagonal submatrices, we don't need them.
        ! Saves some memory.
      select case (cspace)
      case (CCSPACE_PRIMAL,CCSPACE_DUAL)
        call stlin_disableSubmatrix (rnonlinearSpatialMatrix,1,2)
        call stlin_disableSubmatrix (rnonlinearSpatialMatrix,2,1)
      end select
      
      ! Allocate memory.
      call smva_assembleMatrix (CCMASM_ALLOCMEM,CCMASM_MTP_AUTOMATIC,rassemblyFlags,&
          rnonlinearSpatialMatrix,rpreconditioner%p_RmatrixPrecondFullSpace(ilev))
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
    
    ! Create temp vectors in the size of the FE space
    do ilev = rpreconditioner%NLMIN,rpreconditioner%NLMAX
      call lsysbl_createVectorBlock(&
          rpreconditioner%p_rfeHierarchy%p_rfeSpaces(ilev)%p_rdiscretisation,&
          rpreconditioner%p_RtempVec(1,ilev),.false.)
      call lsysbl_createVectorBlock(&
          rpreconditioner%p_rfeHierarchy%p_rfeSpaces(ilev)%p_rdiscretisation,&
          rpreconditioner%p_RtempVec(2,ilev),.false.)
      call lsysbl_createVectorBlock(&
          rpreconditioner%p_rfeHierarchy%p_rfeSpaces(ilev)%p_rdiscretisation,&
          rpreconditioner%p_RtempVec(3,ilev),.false.)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fbsim_updateDiscreteBCprec (rglobalData,rpreconditioner,dtimePrimal,dtimeDual)
  
!<description>
  ! Updates the discrete boundary conditions in the preconditioner
!</description>
  
!<input>
  ! Global settings for callback routines.
  type(t_globalData), intent(inout), target :: rglobalData

  ! Current simulation time for the primal equation.
  real(dp), intent(in) :: dtimePrimal

  ! Current simulation time for the dual equation.
  real(dp), intent(in) :: dtimeDual
!</input>
  
!<inputoutput>
  ! Configuration block of the preconditioner.
  type(t_fbsimPreconditioner), intent(inout) :: rpreconditioner 
!</inputoutput>

!</subroutine>

    ! local variables
    logical :: bneumann
    integer :: ilev

    ! Clear the BC`s and reassemble on all levels.
    do ilev = rpreconditioner%NLMIN,rpreconditioner%NLMAX
      call bcasm_clearDiscreteBC(rpreconditioner%p_RdiscreteBC(ilev))
      call bcasm_clearDiscreteFBC(rpreconditioner%p_RdiscreteFBC(ilev))
      
      call sbc_assembleBDconditions (rpreconditioner%p_rboundaryConditions,dtimePrimal,dtimeDual,&
          rpreconditioner%p_rfeHierarchy%p_rfeSpaces(ilev)%p_rdiscretisation,&
          rpreconditioner%p_rtimeDiscr,&
          rpreconditioner%cspace,rpreconditioner%p_RdiscreteBC(ilev),&
          rglobalData,bneumann)
      call sbc_assembleFBDconditions (dtimePrimal,&
          rpreconditioner%p_rfeHierarchy%p_rfeSpaces(ilev)%p_rdiscretisation,&
          rpreconditioner%p_rtimeDiscr,&
          rpreconditioner%cspace,rpreconditioner%p_RdiscreteFBC(ilev),&
          rglobalData)
    end do

    ! Do we have Neumann boundary?
    ! The Neumann flag on the maximum level decides upon that.
    ! This may actually change from level to level but we simplify here
    ! and use the filter only depending on the max. level.
    rpreconditioner%bhasNeumann = bneumann
    if (.not. bneumann) then
      ! Pure Dirichlet problem -- Neumann boundary for the pressure.
      ! Filter the pressure to avoid indefiniteness.
      rpreconditioner%bneedPressureDiagonalBlock = .true.
      rpreconditioner%RfilterChain(3)%ifilterType = FILTER_DONOTHING
      rpreconditioner%RfilterChain(6)%ifilterType = FILTER_DONOTHING
      select case (rpreconditioner%cspace)
      case (CCSPACE_PRIMAL)
        rpreconditioner%RfilterChain(3)%ifilterType = FILTER_TOL20
        rpreconditioner%RfilterChain(3)%itoL20component = 3
      case (CCSPACE_DUAL)
        rpreconditioner%RfilterChain(6)%ifilterType = FILTER_TOL20
        rpreconditioner%RfilterChain(6)%itoL20component = 6
      case (CCSPACE_PRIMALDUAL)
        rpreconditioner%RfilterChain(3)%ifilterType = FILTER_TOL20
        rpreconditioner%RfilterChain(3)%itoL20component = 3
        rpreconditioner%RfilterChain(6)%ifilterType = FILTER_TOL20
        rpreconditioner%RfilterChain(6)%itoL20component = 6
      end select
    else
      rpreconditioner%RfilterChain(3)%ifilterType = FILTER_DONOTHING
      rpreconditioner%RfilterChain(6)%ifilterType = FILTER_DONOTHING
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
    if (associated(rpreconditioner%p_rsolverNode)) then
      call linsol_releaseSolver (rpreconditioner%p_rsolverNode)
    end if

    ! Clean up data about the projection
    nullify(rpreconditioner%p_rprjHierarchy)
    
    ! Release assembly data.
    deallocate(rpreconditioner%p_RassemblyData)
    
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
    deallocate(rpreconditioner%p_RtempVec)
    
    rpreconditioner%ctypePreconditioner = 0

  end subroutine

!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine fbsim_preparePrecMatrixAssembly (rpreconditioner,ilev,&
!      rnonlinearSpatialMatrixTemplate,rnonlinearSpatialMatrix)
!
!!<description>
!  ! Prepares a rnonlinearSpatialMatrix structure for the assembly according
!  ! to a preconditioner on level ilev. rpreconditioner specifies a couple of preconditioner
!  ! flags that configure the shape of the system matrix that the preconditioner
!  ! needs. rnonlinearSpatialMatrixTemplate is a template structure defining the
!  ! weights for all levels. The routine generates a new structure
!  ! rnonlinearSpatialMatrix using the weights from rnonlinearSpatialMatrixTemplate
!  ! and incorporating special assembly specific flags.
!  !
!  ! fbsim_initNonlinMatrix must have been called prior to this routine to
!  ! initialise the basic matrix. fbsim_preparePrecondMatrixAssembly will then
!  ! add assembly-specific parameters of the preconditioner.
!!</description>
!
!!<input>
!  ! Current assembly level.
!  integer, intent(in) :: ilev
!  
!  ! Structure with assembly-specific parameters of the preconditioner.
!  type(t_fbsimPreconditioner), intent(in) :: rpreconditioner
!
!  ! Template nonlinear matrix structure.
!  type(t_nonlinearSpatialMatrix), intent(in) :: rnonlinearSpatialMatrixTemplate
!!</input>
!
!!<inputoutput>
!  ! New nonlinear matrix structure based on rnonlinearSpatialMatrixTemplate.
!  ! Assembly-specific parameters are added.
!  type(t_nonlinearSpatialMatrix), intent(out) :: rnonlinearSpatialMatrix
!!</inputoutput>
!              
!!</subroutine>
!
!    ! Copy the template structure.
!    rnonlinearSpatialMatrix = rnonlinearSpatialMatrixTemplate
!    
!    ! Parameters for adaptive matrices for Q1~ with anisotropic elements
!    rnonlinearSpatialMatrix%iadaptiveMatrices = rpreconditioner%iadaptiveMatrices
!    rnonlinearSpatialMatrix%dadmatthreshold = rpreconditioner%dadmatthreshold
!    
!    ! Pointers for that level.
!    rnonlinearSpatialMatrix%p_rdiscretisation => &
!        rpreconditioner%p_rspaceAsmTempl(ilev)%rdiscretisation
!    rnonlinearSpatialMatrix%p_rstaticInfo => &
!        rpreconditioner%p_rspaceAsmTempl(ilev)%rstaticInfo
!    
!    ! Depending on the level, we have to set up information about
!    ! transposing B-matrices.
!    if (ilev .eq. rpreconditioner%nlmin) then
!      rnonlinearSpatialMatrix%bvirtualTransposedD = rpreconditioner%bneedVirtTransposedDonCoarse
!    else
!      rnonlinearSpatialMatrix%bvirtualTransposedD = rpreconditioner%bneedVirtTransposedD
!    end if
!    
!  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fbsim_assemblePrecMatrices (rpreconditioner,ieqTime,ioffdiag,&
      rspaceTimeMatrix,rnonlinearData,bincorporateBC)

!<description>
  ! Assembles on every level a matrix for the preconditioner in each timestep.
  ! The output is written to the preconditioner matrices in rpreconditioner.
!</description>

!<input>
  ! Data that specifies the nonlinearity on the maximum level of the preconditioner.
  type(t_spatialMatrixNonlinearData), intent(in), target :: rnonlinearData

  ! Id of the iterate, relative to the space-time matrix.
  integer, intent(in) :: ieqTime
  
  ! Determines the index of the offdiagonal in the space-time matrix which should
  ! be assembled. =0 identifies the diagonal, =-1 the lower subdiagonal,
  ! =1 the upper subdiagonal.
  integer, intent(in) :: ioffdiag
  
  ! Space-time matrix that configures the weights.
  type(t_ccoptSpaceTimeMatrix), intent(in) :: rspaceTimeMatrix

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
    type(t_nonlinearSpatialMatrix) :: rnonlinearSpatialMatrix
    type(t_matrixAssemblyFlags) :: rassemblyFlags
    type(t_spatialMatrixNonlinearData), target :: rlocalNonlinearity
    integer, dimension(1), parameter :: Irows = (/1/)

    ! DEBUG!!!
    !real(dp), dimension(:), pointer :: p_Da,p_vec,p_def
    
    ! DEBUG!!!
    !call lsysbl_getbase_double (rd,p_def)
    !call lsysbl_getbase_double (rx,p_vec)
    
    ! On all levels, we have to set up the nonlinear system matrix,
    ! so that the linear solver can be applied to it.
    
    nullify(p_rmatrix)

    do ilev=rpreconditioner%NLMAX,rpreconditioner%NLMIN,-1
    
      ! Prepare an assembly flags structure for the assembly.
      call fbsim_getMatrixAssemblyFlags (rpreconditioner,ilev,rassemblyFlags)

      ! Get the matrix on the current level.
      ! Shift the previous matrix to the pointer of the fine grid matrix.
      p_rmatrixFine => p_rmatrix
      p_rmatrix => rpreconditioner%p_RmatrixPrecondFullSpace(ilev)
      
      ! DEBUG!!!
      !call lsyssc_getbase_double (p_rmatrix%RmatrixBlock(1,1),p_Da)
    
      ! On the highest level, we use rx as solution to build the nonlinear
      ! matrix. On lower levels, we have to create a solution
      ! on that level from a fine-grid solution before we can use
      ! it to build the matrix!
      if (ilev .eq. rpreconditioner%NLMAX) then
      
        p_rvectorCoarse1 => rnonlinearData%p_rvector1
        p_rvectorCoarse2 => rnonlinearData%p_rvector2
        p_rvectorCoarse3 => rnonlinearData%p_rvector3
        
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
          p_rvectorFine1 => rnonlinearData%p_rvector1
          p_rvectorFine2 => rnonlinearData%p_rvector2
          p_rvectorFine3 => rnonlinearData%p_rvector3
        end if

        ! Interpolate the solution from the finer grid to the coarser grid.
        ! The interpolation is configured in the interlevel projection
        ! structure we got from the collection.
        call mlprj_performInterpolationHier (rpreconditioner%p_rprjHierarchy,ilev+1,&
            p_rvectorCoarse1,p_rvectorFine1)
        call mlprj_performInterpolationHier (rpreconditioner%p_rprjHierarchy,ilev+1,&
            p_rvectorCoarse2,p_rvectorFine2)
        call mlprj_performInterpolationHier (rpreconditioner%p_rprjHierarchy,ilev+1,&
            p_rvectorCoarse3,p_rvectorFine3)

        ! Apply the filter chain to the temp vector.
        ! This implements the boundary conditions that are attached to it.
        ! NOTE: Deactivated for standard CC2D compatibility -- and because
        ! it has to be checked whether the correct boundary conditions
        ! are attached to that vector!
        ! CALL filter_applyFilterChainVec (p_rvectorCoarse, p_RfilterChain)

      end if
      
      ! Prepare a local nonlinearity structure for the matrix assembly.
      call smva_initNonlinearData (rlocalNonlinearity,&
          p_rvectorCoarse1,p_rvectorCoarse2,p_rvectorCoarse3)

      ! Create a nonlinear matrix on the current level.
      call smva_initNonlinMatrix (rnonlinearSpatialMatrix,&
          rpreconditioner%p_RassemblyData(ilev),rlocalNonlinearity)
    
      ! Initialise the weights according to the diagonal block of the global
      ! space-time matrix. Disable all matrix blocks which do not belong to
      ! the space we are currently processing.
      call stlin_setupMatrixWeights (rspaceTimeMatrix,&
          ieqTime,ioffdiag,rnonlinearSpatialMatrix)
        
      select case (rpreconditioner%cspace)
      case (CCSPACE_PRIMAL)
        call stlin_disableSubmatrix (rnonlinearSpatialMatrix,2,1)
        call stlin_disableSubmatrix (rnonlinearSpatialMatrix,2,2)
        call stlin_disableSubmatrix (rnonlinearSpatialMatrix,1,2)
      case (CCSPACE_DUAL)
        call stlin_disableSubmatrix (rnonlinearSpatialMatrix,1,1)
        call stlin_disableSubmatrix (rnonlinearSpatialMatrix,1,2)
        call stlin_disableSubmatrix (rnonlinearSpatialMatrix,2,1)
      end select
      !<- letzte änderung!
    
      ! Assemble the matrix.
      ! If we are on a lower level, we can specify a 'fine-grid' matrix.
      if (ilev .eq. rpreconditioner%NLMAX) then
        call smva_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,rassemblyFlags,&
            rnonlinearSpatialMatrix,p_rmatrix)
      else
        call smva_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,rassemblyFlags,&
            rnonlinearSpatialMatrix,p_rmatrix,p_rmatrixFine)
      end if

    end do
    
    if (rpreconditioner%bneedPressureDiagonalBlock .and. bincorporateBC) then
      
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
        if (rnonlinearSpatialMatrix%Dkappa(1,1) .eq. 0.0_DP) then
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
        if (rnonlinearSpatialMatrix%Dkappa(2,2) .eq. 0.0_DP) then
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
          if (rnonlinearSpatialMatrix%Dkappa(1,1) .eq. 0.0_DP) then
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
          if (rnonlinearSpatialMatrix%Dkappa(2,2) .eq. 0.0_DP) then
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
            
      case (CCSPACE_DUAL)

        ! Get the dual sub-operators
        call lsysbl_deriveSubmatrix (&
            rpreconditioner%p_RmatrixPrecondFullSpace(ilev),&
            rpreconditioner%p_RmatrixPrecond(ilev),&
            LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE,4,6)
            
      case (CCSPACE_PRIMALDUAL)

        ! Take the full matrix.
        call lsysbl_duplicateMatrix (&
            rpreconditioner%p_RmatrixPrecondFullSpace(ilev),&
            rpreconditioner%p_RmatrixPrecond(ilev),&
            LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)

      end select

      call lsysbl_assignDiscrDirectMat (&
          rpreconditioner%p_RmatrixPrecond(ilev),&
          rpreconditioner%p_rfeHierarchy%p_rfeSpaces(ilev)%p_rdiscretisation)
      
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

  subroutine fbsim_precondSpaceDefect (rpreconditioner,rd,bsuccess,rnonlinearIteration)

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
  type(t_vectorBlock), intent(inout)            :: rd

  ! If the preconditioning was a success. Is normally automatically set to
  ! TRUE. If there is an error in the preconditioner, this flag can be
  ! set to FALSE. In this case, the nonlinear solver breaks down with
  ! the error flag set to 'preconditioner broke down'.
  logical, intent(inout) :: bsuccess

  ! OPTIONAL: Reference to the nonlinear solver structure.
  ! If present, convergence rates and needed time are written to here.
  type(t_fbsimNonlinearIteration), intent(inout), optional :: rnonlinearIteration
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ierror
    type(t_linsolNode), pointer :: p_rsolverNode
    type(t_timer) :: rtimer

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
      
      ! Initialise data of the solver. This in fact performs a numeric
      ! factorisation of the matrices in UMFPACK-like solvers.
      call linsol_initStructure (rpreconditioner%p_rsolverNode,ierror)
      call linsol_initData (p_rsolverNode, ierror)
      if (ierror .ne. LINSOL_ERR_NOERROR) then
        print *,'linsol_initData failed!'
        call sys_halt()
      end if
      
      ! Init statistics
      call stat_clearTimer (rtimer)
      
      ! Solve the system. As we want to solve Ax=b with
      ! b being the real RHS and x being the real solution vector,
      ! we use linsol_solveAdaptively. If b is a defect
      ! RHS and x a defect update to be added to a solution vector,
      ! we would have to use linsol_precondDefect instead.
      call stat_startTimer (rtimer)
      call linsol_precondDefect (p_rsolverNode,rd)
      call stat_stopTimer (rtimer)
      
      ! Return statistical data in the preconditioner structure.
      rpreconditioner%nlinearIterations = p_rsolverNode%iiterations
      rpreconditioner%dtimeLinearSolver = rtimer%delapsedReal

      ! Release the numeric factorisation of the matrix.
      ! We don't release the symbolic factorisation, as we can use them
      ! for the next iteration.
      call linsol_doneData (p_rsolverNode)
      call linsol_doneStructure (p_rsolverNode)
      
      ! Remember convergence rate for output
      if (present(rnonlinearIteration)) then
        rnonlinearIteration%drhoLinearSolver = p_rsolverNode%dconvergenceRate
        rnonlinearIteration%dlastInitPrecDef = p_rsolverNode%dinitialDefect
        rnonlinearIteration%dlastPrecDef = p_rsolverNode%dfinalDefect

        rnonlinearIteration%dtimeLinearSolver = rtimer%delapsedReal
        rnonlinearIteration%nlinearIterations = p_rsolverNode%iiterations
      end if

      ! Did the preconditioner work?
      bsuccess = p_rsolverNode%iresult .eq. 0
      
      if (bsuccess) then
        ! Filter the final defect.
        call vecfil_discreteBCdef (rd,rpreconditioner%p_RdiscreteBC(rpreconditioner%NLMAX))
        call vecfil_discreteFBCdef (rd,rpreconditioner%p_RdiscreteFBC(rpreconditioner%NLMAX))
        if (.not. rpreconditioner%bhasNeumann) then
          call vecfil_normaliseToL20Sca (rd%RvectorBlock(3))
          if (rpreconditioner%cspace .eq. CCSPACE_PRIMALDUAL) then
            call vecfil_normaliseToL20Sca (rd%RvectorBlock(6))
          end if
        end if
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
    if (ite .eq. 0) then
    
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

      rnonlinearIteration%dinitialDefectTotal = sqrt(Dresiduals(1)**2 + Dresiduals(2)**2)
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

      if ((ddelU .lt. SYS_EPSREAL*1E2_DP) .and. &
          (ddelP .lt. SYS_EPSREAL*1E2_DP)) then
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

  subroutine fbsim_getNonlinearDefect (cspace,rspaceTimeMatrix,ieqTime,&
      rdiscrData,rnonlinearData,rsolution,rrhs,rdefect)
  
!<description>
  ! Calculates the nonlinear defect in a timestep.
  !
  ! Note: The boundary conditions are not implemented into the defect.
!</description>
  
!<input>
  ! Underlying space. One of the CCSPACE_xxxx constants.
  integer, intent(in) :: cspace

  ! Underlying Space-time matrix.
  type(t_ccoptSpaceTimeMatrix), intent(in) :: rspaceTimeMatrix
  
  ! Number of the time solution whose nonlinear defect should be calculated.
  integer, intent(in) :: ieqTime

  ! Assembly data 
  type(t_spatialMatrixDiscrData), intent(in) :: rdiscrData
  
  ! Nonlinearity
  type(t_spatialMatrixNonlinearData), intent(in) :: rnonlinearData

  ! Solution vector
  type(t_vectorBlock), intent(in) :: rsolution

  ! Right hand side vector
  type(t_vectorBlock), intent(in) :: rrhs
!</input>

!<inputoutput>
  ! Defect vector
  type(t_vectorBlock), intent(inout) :: rdefect
!</inputoutput>

!</subroutine>

  ! Nonlinear matrix
  type(t_nonlinearSpatialMatrix) :: rnonlinearSpatialMatrix

    ! Prepare subtraction of A*primal/dual solution from the rhs to 
    ! create the primal defect.
    call smva_initNonlinMatrix (rnonlinearSpatialMatrix,rdiscrData,rnonlinearData)
    call stlin_setupMatrixWeights (rspaceTimeMatrix,&
        ieqTime,0,rnonlinearSpatialMatrix)

    ! Get the rhs, disable not used submatrices
    select case (cspace)
    case (CCSPACE_PRIMAL)
      ! Only primal space, set the dual defect to zero
      call lsyssc_copyVector (rrhs%RvectorBlock(1),rdefect%RvectorBlock(1))
      call lsyssc_copyVector (rrhs%RvectorBlock(2),rdefect%RvectorBlock(2))
      call lsyssc_copyVector (rrhs%RvectorBlock(3),rdefect%RvectorBlock(3))
      call lsyssc_clearVector (rdefect%RvectorBlock(4))
      call lsyssc_clearVector (rdefect%RvectorBlock(5))
      call lsyssc_clearVector (rdefect%RvectorBlock(6))

      call stlin_disableSubmatrix (rnonlinearSpatialMatrix,2,1)
      call stlin_disableSubmatrix (rnonlinearSpatialMatrix,2,2)
      call stlin_disableSubmatrix (rnonlinearSpatialMatrix,1,2)
      
    case (CCSPACE_DUAL)
      ! Only dual space, set the primal defect to zero
      call lsyssc_clearVector (rdefect%RvectorBlock(1))
      call lsyssc_clearVector (rdefect%RvectorBlock(2))
      call lsyssc_clearVector (rdefect%RvectorBlock(3))
      call lsyssc_copyVector (rrhs%RvectorBlock(4),rdefect%RvectorBlock(4))
      call lsyssc_copyVector (rrhs%RvectorBlock(5),rdefect%RvectorBlock(5))
      call lsyssc_copyVector (rrhs%RvectorBlock(6),rdefect%RvectorBlock(6))

      call stlin_disableSubmatrix (rnonlinearSpatialMatrix,1,1)
      call stlin_disableSubmatrix (rnonlinearSpatialMatrix,1,2)
      call stlin_disableSubmatrix (rnonlinearSpatialMatrix,2,1)
    case default
      call lsysbl_copyVector (rrhs,rdefect)
    end select
    
    call smva_assembleDefect (rnonlinearSpatialMatrix,rsolution,&
      rdefect,1.0_DP)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fbsim_initpostprocessing (rpostproc,cspace,rboundaryConditions,rspaceTimeMatrix)
  
!<description>
  ! Initialises the postprocessing
!</description>

!<input>
  ! Space that is available in rsolution. One of the CCSPACE_xxxx constants.
  integer, intent(in) :: cspace
  
  ! Boundary conditions to use.
  type(t_optcBDC), intent(in), target  :: rboundaryConditions

  ! Reference to the space-time matrix defining the space and time discretisation
  ! on the maximum level.
  type(t_ccoptSpaceTimeMatrix), intent(in) :: rspaceTimeMatrix
!</input>
  
!<inputoutput>
  ! Parameters about the postprocessing.
  type(t_fbsimPostprocessing), intent(inout) :: rpostproc
!</inputoutput>
  
!</subroutine>

    ! Fetch data
    if (rpostproc%cspace .eq. CCSPACE_UNDEFINED) then
      ! Allow a caller to specify the type of postprocessing, do not overwrite
      ! previous data.
      rpostproc%cspace = cspace
    end if
    rpostproc%p_rspaceDiscr => rspaceTimeMatrix%rdiscrData%p_rspaceDiscr
    rpostproc%p_rtimeDiscr => rspaceTimeMatrix%rdiscrData%p_rtimeDiscr
    rpostproc%p_rboundaryConditions => rboundaryConditions
    call bcasm_initDiscreteBC(rpostproc%rdiscreteBC)
    call bcasm_initDiscreteFBC(rpostproc%rdiscreteFBC)

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine fbsim_donepostprocessing (rpostproc)
  
!<description>
  ! Clean up the postprocessing structure
!</description>
  
!<inputoutput>
  ! Parameters about the postprocessing.
  type(t_fbsimPostprocessing), intent(inout) :: rpostproc
!</inputoutput>
  
!</subroutine>

    ! local variables
    rpostproc%cspace = CCSPACE_UNDEFINED
    nullify(rpostproc%p_rspaceDiscr)
    nullify(rpostproc%p_rtimeDiscr)
    nullify(rpostproc%p_rboundaryConditions)
    call bcasm_releaseDiscreteBC(rpostproc%rdiscreteBC)
    call bcasm_releaseDiscreteFBC(rpostproc%rdiscreteFBC)

  end subroutine
  
!******************************************************************************

!<subroutine>

  subroutine fbsim_writeUCD (rpostprocessing,rvector,istep,dtimePrimal,dtimeDual,rsettings)

!<description>
  ! Writes an UCD postprocessing file as configured in the DAT file.
  ! (-> GMV, AVS, Paraview,...)
!</description>
  
!<input>
  ! Solution vector.
  type(t_vectorBlock), intent(in) :: rvector
  
  ! Id of the timestep.
  integer, intent(in) :: istep
  
  ! Simulation time of the primal equation
  ! Must be set to 0 in stationary simulations.
  real(DP), intent(in) :: dtimePrimal

  ! Simulation time of the dual equation
  ! Must be set to 0 in stationary simulations.
  real(DP), intent(in) :: dtimeDual

  ! The structure of the main solver
  type(t_settings_optflow), intent(inout), target :: rsettings
!</input>

!<inputoutput>
  ! Postprocessing structure. Must have been initialised prior
  ! to calling this routine.
  ! The time stamp of the last written out GMV is updated.
  type(t_fbsimPostprocessing), intent(inout) :: rpostprocessing  
!</inputoutput>

!</subroutine>

    ! We need some more variables for postprocessing - i.e. writing
    ! a GMV file.
    real(DP), dimension(:), pointer :: p_Ddata,p_Ddata2

    ! A pointer to the triangulation.
    type(t_triangulation), pointer :: p_rtriangulation
    
    ! A vector accepting Q1 data
    type(t_vectorBlock) :: rprjVector
    
    ! A discretisation structure for Q1
    type(t_blockDiscretisation) :: rprjDiscretisation
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    
    integer(I32) :: ieltype
    
    character(SYS_STRLEN) :: sfile
    
    if (rpostprocessing%ioutputUCD .eq. 0) return

    ! The solution vector is probably not in the way, GMV likes it!
    ! GMV for example does not understand Q1~ vectors!
    ! Therefore, we first have to convert the vector to a form that
    ! GMV understands.
    ! GMV understands only Q1/P1 solutions! So the task is now to
    ! create a Q1/P1 solution from rvector and write that out.
    !
    ! For this purpose, first create a 'derived' simple discretisation
    ! structure based on Q1/P1 by copying the main guiding block 
    ! discretisation structure and modifying the discretisation
    ! structures of the two velocity subvectors:
    
    call spdiscr_duplicateBlockDiscr(rvector%p_rblockDiscr,rprjDiscretisation)
    
    call spdiscr_deriveDiscr_triquad (&
                 rvector%p_rblockDiscr%RspatialDiscr(1), &
                 EL_P1, EL_Q1, CUB_TRZ_T, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscr(1))

    call spdiscr_deriveDiscr_triquad (&
                 rvector%p_rblockDiscr%RspatialDiscr(2), &
                 EL_P1, EL_Q1, CUB_TRZ_T, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscr(2))

    call spdiscr_deriveDiscr_triquad (&
                 rvector%p_rblockDiscr%RspatialDiscr(3), &
                 EL_P0, EL_Q0, CUB_TRZ_T, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscr(3))

    call spdiscr_deriveDiscr_triquad (&
                 rvector%p_rblockDiscr%RspatialDiscr(4), &
                 EL_P1, EL_Q1, CUB_TRZ_T, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscr(4))

    call spdiscr_deriveDiscr_triquad (&
                 rvector%p_rblockDiscr%RspatialDiscr(5), &
                 EL_P1, EL_Q1, CUB_TRZ_T, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscr(5))

    call spdiscr_deriveDiscr_triquad (&
                 rvector%p_rblockDiscr%RspatialDiscr(6), &
                 EL_P0, EL_Q0, CUB_TRZ_T, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscr(6))
                 
    ! The pressure discretisation substructure stays the old.
    !
    ! Now set up a new solution vector based on this discretisation,
    ! allocate memory.
    call lsysbl_createVectorBlock (rprjDiscretisation,rprjVector,.false.)
    
    ! Then take our original solution vector and convert it according to the
    ! new discretisation:
    call spdp_projectSolution (rvector,rprjVector)

    ! Discretise the boundary conditions according to the Q1/Q1/Q0 
    ! discretisation for implementing them into a solution vector.
    call sbc_assembleBDconditions (rpostprocessing%p_rboundaryConditions,dtimePrimal,dtimeDual,&
        rprjDiscretisation,rpostprocessing%p_rtimeDiscr,&
        rpostprocessing%cspace,rpostprocessing%rdiscreteBC,rsettings%rglobalData)
    call sbc_assembleFBDconditions (dtimePrimal,&
        rprjDiscretisation,rpostprocessing%p_rtimeDiscr,&
        rpostprocessing%cspace,rpostprocessing%rdiscreteFBC,rsettings%rglobalData)
    
    ! Filter the solution vector to implement discrete BC`s.
    call vecfil_discreteBCsol (rprjVector,rpostprocessing%rdiscreteBC)
    call vecfil_discreteFBCsol (rprjVector,rpostprocessing%rdiscreteFBC)
    
    ! Create the actual filename
    sfile = trim(adjustl(rpostprocessing%sfilenameUCD))//'.'//sys_si0(istep,5)
                                 
    ! Now we have a Q1/Q1/Q0 solution in rprjVector -- on the level NLMAX.
    ! The next step is to project it down to level ilevelUCD.
    ! Due to the fact that solutions are usually 2-level-ordered,
    ! this can be shortened by taking only the first NVT vertices
    ! of the solution vector!
    
    ! From the attached discretisation, get the underlying triangulation
    ! of that level
    p_rtriangulation => rpostprocessing%p_rspaceDiscr%p_rtriangulation
    
    ! Start UCD export to GMV file:
    call output_lbrk ()
    call output_line ('Writing GMV file: '//sfile)
    
    select case (rpostprocessing%ioutputUCD)
    case (1)
      call ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,sfile)

    case (2)
      call ucd_startAVS (rexport,UCD_FLAG_STANDARD,p_rtriangulation,sfile)
          
    case (3)
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,p_rtriangulation,sfile)
          
    case default
      call output_line ('Invalid UCD ooutput type.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'fbsim_writeUCD')
      stop
    end select
        
    ! Set the simulation time.
    call ucd_setSimulationTime (rexport,dtimePrimal)
    
    ! Write the configuration of the application as comment block
    ! to the output file.
    call ucd_addCommentLine (rexport,'Configuration:')
    call ucd_addCommentLine (rexport,'---------------')
    call ucd_addParameterList (rexport,rsettings%p_rparlist)
    call ucd_addCommentLine (rexport,'---------------')

    ! Get the velocity field
    select case (rpostprocessing%cspace)
    case (CCSPACE_PRIMAL)
      call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
      call lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
      
      ! Write the velocity field
      call ucd_addVarVertBasedVec (rexport,'velocity_p',&
          p_Ddata(1:p_rtriangulation%NVT),p_Ddata2(1:p_rtriangulation%NVT))
      
      ! Write out cell based or node based pressure.
      call lsyssc_getbase_double (rprjVector%RvectorBlock(3),p_Ddata)
      ieltype = rvector%p_rblockDiscr%RspatialDiscr(3)% &
                RelementDistr(1)%celement
                
      if ((elem_getPrimaryElement(ieltype) .eq. EL_Q1) .or. &
          ((elem_getPrimaryElement(ieltype) .eq. EL_P1))) then
        call ucd_addVariableVertexBased (rexport,'pressure_p',UCD_VAR_STANDARD, &
            p_Ddata(1:p_rtriangulation%NVT))
      else
        call ucd_addVariableElementBased (rexport,'pressure_p',UCD_VAR_STANDARD, &
            p_Ddata(1:p_rtriangulation%NEL))
      end if
     
    case (CCSPACE_DUAL)
      call lsyssc_getbase_double (rprjVector%RvectorBlock(4),p_Ddata)
      call lsyssc_getbase_double (rprjVector%RvectorBlock(5),p_Ddata2)
      
      ! Write the velocity field
      call ucd_addVarVertBasedVec (rexport,'velocity_d',&
          p_Ddata(1:p_rtriangulation%NVT),p_Ddata2(1:p_rtriangulation%NVT))
      
      ! Write cell based or node based pressure.
      call lsyssc_getbase_double (rprjVector%RvectorBlock(6),p_Ddata)
      ieltype = rvector%p_rblockDiscr%RspatialDiscr(6)% &
                RelementDistr(1)%celement
                
      if ((elem_getPrimaryElement(ieltype) .eq. EL_Q1) .or. &
          ((elem_getPrimaryElement(ieltype) .eq. EL_P1))) then
        call ucd_addVariableVertexBased (rexport,'pressure_d',UCD_VAR_STANDARD, &
            p_Ddata(1:p_rtriangulation%NVT))
      else
        call ucd_addVariableElementBased (rexport,'pressure_d',UCD_VAR_STANDARD, &
            p_Ddata(1:p_rtriangulation%NEL))
      end if

    case default
      ! Primal solution
      call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
      call lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
      
      ! Write the velocity field
      call ucd_addVarVertBasedVec (rexport,'velocity_p',&
          p_Ddata(1:p_rtriangulation%NVT),p_Ddata2(1:p_rtriangulation%NVT))
      
      ! Write pressure
      call lsyssc_getbase_double (rprjVector%RvectorBlock(3),p_Ddata)
      call ucd_addVariableElementBased (rexport,'pressure',UCD_VAR_STANDARD, &
          p_Ddata(1:p_rtriangulation%NEL))

      ! If we have a simple Q1 discretisation in the pressure, write it out as it is
      if (rvector%p_rblockDiscr%RspatialDiscr(3)%ccomplexity .eq. SPDISC_UNIFORM) then
        ieltype = rvector%p_rblockDiscr%RspatialDiscr(3)% &
                  RelementDistr(1)%celement
                  
        if (elem_getPrimaryElement(ieltype) .eq. EL_Q1) then
          call lsyssc_getbase_double (rvector%RvectorBlock(3),p_Ddata)
          call ucd_addVariableVertexBased (rexport,'pressure',UCD_VAR_STANDARD, &
              p_Ddata(1:p_rtriangulation%NEL))
        end if
      end if

      ! Dual solution
      call lsyssc_getbase_double (rprjVector%RvectorBlock(4),p_Ddata)
      call lsyssc_getbase_double (rprjVector%RvectorBlock(5),p_Ddata2)
      
      ! Write the velocity field
      call ucd_addVarVertBasedVec (rexport,'velocity_d',&
          p_Ddata(1:p_rtriangulation%NVT),p_Ddata2(1:p_rtriangulation%NVT))
      
      ! Write pressure
      call lsyssc_getbase_double (rprjVector%RvectorBlock(6),p_Ddata)
      call ucd_addVariableElementBased (rexport,'pressure_d',UCD_VAR_STANDARD, &
          p_Ddata(1:p_rtriangulation%NEL))

      ! If we have a simple Q1 discretisation in the pressure, write it out as it is
      if (rvector%p_rblockDiscr%RspatialDiscr(6)%ccomplexity .eq. SPDISC_UNIFORM) then
        ieltype = rvector%p_rblockDiscr%RspatialDiscr(6)% &
                  RelementDistr(1)%celement
                  
        if (elem_getPrimaryElement(ieltype) .eq. EL_Q1) then
          call lsyssc_getbase_double (rvector%RvectorBlock(6),p_Ddata)
          call ucd_addVariableVertexBased (rexport,'pressure_d',UCD_VAR_STANDARD, &
              p_Ddata(1:p_rtriangulation%NEL))
        end if
      end if

    end select    

    ! Write the file to disc, that is it.
    call ucd_write (rexport)
    call ucd_release (rexport)
    
    ! Release the auxiliary vector
    call lsysbl_releaseVector (rprjVector)
    
    ! Release the discretisation structure.
    call spdiscr_releaseBlockDiscr (rprjDiscretisation)
    
  end subroutine

!******************************************************************************

!<subroutine>

  subroutine fbsim_writeRAW (rpostprocessing,rvector,istep,rsettings)

!<description>
  ! Writes out the solution as raw data file to disc.
!</description>
  
!<input>
  ! Solution vector.
  type(t_vectorBlock), intent(in) :: rvector
  
  ! Id of the timestep.
  integer, intent(in) :: istep
  
  ! The structure of the main solver
  type(t_settings_optflow), intent(inout), target :: rsettings
!</input>

!<inputoutput>
  ! Postprocessing structure. Must have been initialised prior
  ! to calling this routine.
  ! The time stamp of the last written out GMV is updated.
  type(t_fbsimPostprocessing), intent(inout) :: rpostprocessing  
!</inputoutput>

!</subroutine>

    ! local variables
    character(len=SYS_STRLEN) :: sfilename
    
    ! Create a filename.
    if (rpostprocessing%swriteSolutionFilename .eq. "") return
    sfilename = trim(rpostprocessing%swriteSolutionFilename)//sys_si0L(istep,5)
    
    ! Write out the file.
    select case (rpostprocessing%cwriteSolution) 
    case (1)
      ! Formatted
      call vecio_writeBlockVectorHR (rvector, "vector", .true.,&
          0, sfilename, "(E20.10)")
    case (2)
      ! Unformatted
      call vecio_writeBlockVectorHR (rvector, "vector", .true.,&
          0, sfilename)
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fbsim_postprocessing (rpostproc,rsolution,istep,dtimePrimal,dtimeDual,rsettings)
  
!<description>
  ! Performs postprocessing to a solution vector.
!</description>
  
!<input>
  ! Solution vector for postprocessing
  type(t_vectorBlock), intent(in) :: rsolution

  ! Id of the timestep.
  integer, intent(in) :: istep
  
  ! Simulation time for the primal equation
  real(DP), intent(in) :: dtimePrimal

  ! Simulation time for the dual equation
  real(DP), intent(in) :: dtimeDual

  ! Reference to the problem structure
  type(t_settings_optflow), intent(inout) :: rsettings
!</input>

!<inputoutput>
  ! Parameters about the postprocessing.
  type(t_fbsimPostprocessing), intent(inout) :: rpostproc
!</inputoutput>
  
!</subroutine>

    ! Output to a visualisation file?
    if (rpostproc%ioutputUCD .ne. 0) then
      call fbsim_writeUCD (rpostproc,rsolution,istep-1,dtimePrimal,dtimeDual,rsettings)
    end if
    
    if (rpostproc%cwriteSolution .ne. 0) then
      call fbsim_writeRAW (rpostproc,rsolution,istep-1,rsettings)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fbsim_simulate (rsimsolver, rsolvector, rrhsvector, &
      rboundaryConditions, ifirstinterval, ilastinterval, roseenSolution, &
      domega,bsuccess)
  
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

  ! Boundary conditions to use.
  type(t_optcBDC), intent(in), target :: rboundaryConditions
  
  ! First time interval to simulate. 
  integer, intent(in) :: ifirstinterval

  ! Last time interval to simulate
  ! ifirstinterval=1 and ilastinterval=0 simulates the initial condition.
  integer, intent(in) :: ilastinterval
  
  ! Damping parameter. Standard = 1.0.
  real(dp), intent(in) :: domega
!</input>

!<inputoutput>
  ! An initial space-time solution vector. This is modified according to the iteration.
  ! the first subvector must contain a proper initial condition.
  type(t_spacetimevector), intent(inout) :: rsolvector
  
  ! TRUE if the iteration was successful.
  logical, intent(out),optional :: bsuccess
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iiterate,NEQtime
    logical :: blocalsuccess
    real(dp) :: dtimePrimal,dtimeDual,dweightold,dweightnew,dtstep
    type(t_vectorBlock) :: rrhs
    type(t_vectorBlock) :: rdefect
    type(t_vectorBlock) :: rdefectPrimal,rdefectDual,rdefectBackup
    type(t_staticSpaceAsmTemplates), pointer :: p_rspaceAsmTempl 
    type(t_nonlinearSpatialMatrix) :: rnonlinearSpatialMatrix
    type(t_ccoptSpaceTimeMatrix), pointer :: p_rspaceTimeMatrix
    type(t_spatialMatrixNonlinearData) :: rnonlinearData
    
    ! nonlinear solver
    integer :: ite
    real(dp) :: dres,dresInit,dtempdef
    logical :: bconvergence,bdivergence
    integer, dimension(3) :: Cnorms
    type(t_ccoptSpaceTimeMatrix) :: rspaceTimeMatrixPrecond
    type(t_ccDynamicNewtonControl), pointer :: p_radaptiveNewton
    
    ! Statistical data
    type(t_timer) :: rtimerNonlinearSolver,rtimerMatrixAssembly
    type(t_timer) :: rtimerTotal,rtimerPostProc,rtimerDefectCalc
    
    ! Vector pool for quick access
    type(t_spaceTimeVectorAccess) :: rsolAccess, roseenAccess, rrhsAccess
    type(t_vectorBlock), pointer :: p_rprevsol, p_rcurrentsol, p_rprevrhs, p_rcurrentrhs
    type(t_vectorBlock), pointer :: p_roseensol1, p_roseensol2, p_roseensol3
    
    ! DEBUG!!!
    real(dp), dimension(:), pointer :: p_Dsol,p_Drhs,p_Doseen1,p_Doseen2,p_Doseen3
    real(dp), dimension(:), pointer :: p_Dactrhs
    real(dp), dimension(:), pointer :: p_Ddefect,p_DdefP,p_DdefD
    real(dp), dimension(:), pointer :: p_DprevSol,p_DprevRhs
    
    ! Fetch some necessary information from the structures.
    p_rspaceTimeMatrix => rsimsolver%p_rmatrix
    
    ! Get the level information structure of the current level
    p_rspaceAsmTempl => p_rspaceTimeMatrix%rdiscrData%p_rstaticSpaceAsmTempl
    
    ! Put the BDC's to the solver structure as well as to the preconditioner.
    rsimsolver%p_rboundaryConditions => rboundaryConditions
    rsimSolver%rpreconditioner%p_rboundaryConditions => rboundaryConditions
    
    if (rsimSolver%rpreconditionerAlternative%ctypePreconditioner .ne. 0) then
      ! Also to the alternative preconditioner.
      rsimSolver%rpreconditionerAlternative%p_rboundaryConditions => rboundaryConditions
    end if
    
    ! Create vector pools for the space-time vectors to have quicker access.
    call sptivec_createAccessPool(rrhsvector,rrhsAccess,4)
    call sptivec_createAccessPool(rsolvector,rsolAccess,4)
    if (present(roseenSolution)) then
      call sptivec_createAccessPool(roseenSolution,roseenAccess,4)
    end if
    
    call lsysbl_createVectorBlock(rsimsolver%rdiscrData%p_rdiscrPrimalDual,rdefect)
    
    ! DEBUG!!!
!    call lsysbl_getbase_double (roseensol1,p_Doseen1)
!    call lsysbl_getbase_double (roseensol2,p_Doseen2)
!    call lsysbl_getbase_double (roseensol3,p_Doseen3)
!    call lsysbl_getbase_double (rcurrentsol,p_Dsol)
!    call lsysbl_getbase_double (rcurrentrhs,p_Drhs)
!    call lsysbl_getbase_double (rprevrhs,p_DprevRhs)
!    call lsysbl_getbase_double (rprevsol,p_DprevSol)
    
    ! Create a temp matrix.
    
    ! Create subvectors for the primal/dual part of the defect.
    call lsysbl_deriveSubvector(rdefect,rdefectPrimal,1,3,.true.)
    call lsysbl_deriveSubvector(rdefect,rdefectDual,4,6,.true.)
    call lsysbl_assignDiscrDirectVec (rdefectPrimal,rsimsolver%rdiscrData%p_rdiscrPrimal)
    call lsysbl_assignDiscrDirectVec (rdefectDual,rsimsolver%rdiscrData%p_rdiscrPrimal)

    ! create a backup vector in case there is an alternative preconditioner.
    call lsysbl_createVectorBlock(rsimsolver%rdiscrData%p_rdiscrPrimalDual,rdefectBackup)

    ! DEBUG!!!
    call lsysbl_getbase_double (rdefect,p_Ddefect)
    call lsysbl_getbase_double (rdefectPrimal,p_DdefP)
    call lsysbl_getbase_double (rdefectDual,p_DdefD)

    blocalsuccess = .false.

    ! Assign BC`s to the vectors.
    call lsysbl_assignDiscreteBC (rdefectPrimal,&
        rsimSolver%rpreconditioner%p_RdiscreteBC(rsimSolver%ilevel))
    call lsysbl_assignDiscreteFBC (rdefectPrimal,&
        rsimSolver%rpreconditioner%p_RdiscreteFBC(rsimSolver%ilevel))

    call sptivec_bindDiscreteBCtoBuffer(rsolAccess,&
        rsimSolver%rpreconditioner%p_RdiscreteBC(rsimSolver%ilevel),&
        rsimSolver%rpreconditioner%p_RdiscreteFBC(rsimSolver%ilevel))

    call sptivec_bindDiscreteBCtoBuffer (rrhsAccess,&
        rsimSolver%rpreconditioner%p_RdiscreteBC(rsimSolver%ilevel),&
        rsimSolver%rpreconditioner%p_RdiscreteFBC(rsimSolver%ilevel))

    ! Prepare the postprocessing
    call fbsim_initpostprocessing (rsimsolver%rpostprocessing,&
        rsimsolver%rpreconditioner%cspace,rsimsolver%p_rboundaryConditions,&
        rsimsolver%p_rmatrix)
    
    ! Initialise timers & statistics
    call stat_clearTimer(rtimerNonlinearSolver)
    call stat_clearTimer(rtimerMatrixAssembly)
    call stat_clearTimer(rtimerPostproc)
    call stat_clearTimer(rtimerDefectCalc)
    call stat_clearTimer(rtimerTotal)
    
    rsimsolver%nlinearIterations = 0
    rsimsolver%dtimeLinearSolver = 0.0_DP
    
    call stat_startTimer(rtimerTotal)
    
    ! What to do?
    select case (rsimsolver%csimtype)
    case (FBSIM_SOLVER_NLFORWARD)
      ! (Navier-)Stokes forward simulation. This iteration type is nonlinear 
      ! in each timestep.

      ! Create a RHS vector.
      call lsysbl_createVectorBlock(rsimsolver%rdiscrData%p_rdiscrPrimalDual,rrhs)
      call lsysbl_assignDiscreteBC (rrhs,&
          rsimSolver%rpreconditioner%p_RdiscreteBC(rsimSolver%ilevel))
      call lsysbl_assignDiscreteFBC (rrhs,&
          rsimSolver%rpreconditioner%p_RdiscreteFBC(rsimSolver%ilevel))

      ! DEBUG!!!
      call lsysbl_getbase_double (rrhs,p_Dactrhs)

      ! Derive a "preconditioner" space-time matrix from the basic space-time
      ! matrix. This may or may not include the Newton part - but by default it
      ! does not.
      rspaceTimeMatrixPrecond = p_rspaceTimeMatrix
      rspaceTimeMatrixPrecond%cmatrixType = 0
      
      ! Focus on the initial solution.
      call tdiscr_getTimestep(p_rspaceTimeMatrix%rdiscrData%p_rtimeDiscr,&
          ifirstinterval-1,dtimePrimal,dtstep)
      dtimeDual = dtimePrimal - (1.0_DP-p_rspaceTimeMatrix%rdiscrData%p_rtimeDiscr%dtheta)*dtstep

      call output_separator (OU_SEP_MINUS,&
          coutputMode=rsimSolver%rnonlinearIteration%coutputMode)
      call output_line ('Time-Iterate '//trim(sys_siL(ifirstinterval,6))// &
          ', Time = '//trim(sys_sdL(dtimePrimal,5)), &
          coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
          
      call fbsim_updateDiscreteBCprec (rsimsolver%p_rsettings%rglobalData,&
          rsimsolver%rpreconditioner,dtimePrimal,dtimeDual)

      ! Get the start vector      
      call sptivec_getVectorFromPool(rsolAccess, ifirstinterval, p_rcurrentsol)
     
      ! Initial postprocessing
      call stat_startTimer (rtimerPostProc)
      call fbsim_postprocessing (rsimSolver%rpostprocessing,p_rcurrentsol,ifirstinterval,&
          dtimePrimal,dtimeDual,rsimsolver%p_rsettings)
      call stat_stopTimer (rtimerPostProc)
          
      ! Loop through the timesteps. Ignore the initial solution, this is the
      ! initial condition.
      do iiterate = ifirstinterval+1,ilastinterval+1
      
        ! Timestep size? Current time step?
        call tdiscr_getTimestep(p_rspaceTimeMatrix%rdiscrData%p_rtimeDiscr,iiterate-1,&
            dtimePrimal,dtstep)
        dtimeDual = dtimePrimal - (1.0_DP-p_rspaceTimeMatrix%rdiscrData%p_rtimeDiscr%dtheta)*dtstep

        if (rsimSolver%ioutputLevel .ge. 1) then
          call output_separator (OU_SEP_MINUS,&
              coutputMode=rsimSolver%rnonlinearIteration%coutputMode)
          call output_line ('Time-Iterate '//trim(sys_siL(iiterate,6))// &
              ', Time = '//trim(sys_sdL(dtimePrimal,5))// &
              ', Stepsize = '//trim(sys_sdL(dtstep,5)),&
              coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
        end if

        ! Update the boundary conditions to the current time
        call fbsim_updateDiscreteBCprec (rsimsolver%p_rsettings%rglobalData,&
            rsimsolver%rpreconditioner,dtimePrimal,dtimeDual)

        call stat_startTimer(rtimerDefectCalc)
        if (iiterate .gt. ifirstinterval) then
          
          ! Get the solution from the access pool.
          call sptivec_getVectorFromPool(rsolAccess, iiterate-1, p_rprevsol)
          call sptivec_getVectorFromPool(rsolAccess, iiterate, p_rcurrentsol)
          call vecfil_discreteBCsol(p_rcurrentsol)
          p_roseensol1 => p_rprevsol
          p_roseensol2 => p_rcurrentsol
          p_roseensol3 => p_rcurrentsol

          ! Create a nonlinear-data structure that specifies the evaluation point
          ! of nonlinear matrices. the evaluation point is given by the three vectors
          ! roseensolX.
          call smva_initNonlinearData (rnonlinearData,p_roseensol1,p_roseensol2,p_roseensol3)

          ! RHS.
          call sptivec_getVectorFromPool(rrhsAccess, iiterate-1, p_rprevrhs)
          call sptivec_getVectorFromPool(rrhsAccess, iiterate, p_rcurrentrhs)
        
          ! Calculate the RHS according to the timestep scheme.
          call tdiscr_getTimestepWeights(p_rspaceTimeMatrix%rdiscrData%p_rtimeDiscr,&
              iiterate-1,dweightold,dweightnew)
          call lsysbl_copyVector (p_rcurrentrhs,rrhs)
          call lsysbl_vectorLinearComb (p_rprevrhs,rrhs,dweightold,dweightnew)

          !! If we are not in the first step...
          ! ... subtract the offdiagonal, corresponding to the previous timestep.
          call smva_initNonlinMatrix (rnonlinearSpatialMatrix,&
              rsimsolver%rdiscrData,rnonlinearData)
          call stlin_setupMatrixWeights (p_rspaceTimeMatrix,&
              iiterate,-1,rnonlinearSpatialMatrix)
          
          call smva_assembleDefect (rnonlinearSpatialMatrix,p_rprevsol,rrhs,1.0_DP)
        else
          ! Initial condition is given. Read it.
          call sptivec_getVectorFromPool(rrhsAccess,ifirstinterval, p_rcurrentrhs)
          call sptivec_getVectorFromPool(rsolAccess,ifirstinterval, p_rcurrentsol)
          call vecfil_discreteBCsol(p_rcurrentsol)
          
          ! All oseen solutions point to the current solution vector.
          p_roseensol1 => p_rcurrentsol
          p_roseensol2 => p_rcurrentsol
          p_roseensol3 => p_rcurrentsol

          ! Create a nonlinear-data structure that specifies the evaluation point
          ! of nonlinear matrices. the evaluation point is given by the three vectors
          ! roseensolX.
          call smva_initNonlinearData (rnonlinearData,p_roseensol1,p_roseensol2,p_roseensol3)

          ! Current RHS
          call lsysbl_copyVector (p_rcurrentrhs,rrhs)
        end if
        call stat_stopTimer(rtimerDefectCalc)
        
        ! Filter the current RHS in rdefect as well as the solution vector
        ! to implement boundary conditions.
        call vecfil_discreteBCrhs(rrhs)
        call vecfil_discreteBCsol(p_rcurrentsol)

        ! Start with ite=0, the 0th iteration; this represents the start point
        ite = 0

        ! Create the initial nonlinear defect
        call stat_startTimer(rtimerDefectCalc)
        call fbsim_getNonlinearDefect (rsimSolver%rpreconditioner%cspace,p_rspaceTimeMatrix,&
            iiterate,rsimsolver%rdiscrData,rnonlinearData,p_rcurrentsol,rrhs,rdefect)
        call stat_stopTimer(rtimerDefectCalc)

        ! Implement the boundary conditions        
        call vecfil_discreteBCdef(rdefectPrimal)

        ! Initial check for convergence/divergence        
        call fbsim_resNormCheck (rsimSolver%rnonlinearIteration,&
            ite,p_rcurrentsol,p_rcurrentrhs,rdefect,bconvergence,bdivergence)

        ! Perform at least nminIterations iterations
        if (ite .lt. rsimSolver%rnonlinearIteration%nminIterations) bconvergence = .false.

        ! Nonlinear loop in the timestep
        ! ------------------------------
        
        if ((.not. bconvergence) .and. (.not. bdivergence)) then
        
          ! Start the nonlinear solver
          call stat_startTimer (rtimerNonlinearSolver)
          do ite = 1,rsimsolver%rnonlinearIteration%nmaxIterations

            ! Preconditioning to the defect
            ! ---------------------------------
            ! Create a space-time matrix that suits our simulation.
            !        
            ! Newton or standard matrix?
            select case (rsimsolver%rnonlinearIteration%ctypeIteration)
            case (0,1)
              ! Standard iteration.
              rspaceTimeMatrixPrecond%cmatrixType = 0
              
            case (2)
              ! Newton iteration
              rspaceTimeMatrixPrecond%cmatrixType = 1
              
            case (3)
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
                
            ! Assemble the preconditioner matrices on all levels.
            call stat_startTimer (rtimerMatrixAssembly)
            call fbsim_assemblePrecMatrices (rsimsolver%rpreconditioner,&
                iiterate,0,rspaceTimeMatrixPrecond,rnonlinearData,.true.)
            call stat_stopTimer (rtimerMatrixAssembly)
            
            if (rsimsolver%rpreconditionerAlternative%ctypePreconditioner .ne. 0) then
              ! There is an alternative preconditioner, so make a backup of our
              ! defect in case we have to apply that.
              call lsysbl_copyVector (rdefect,rdefectBackup)
            end if
                
            ! Do preconditioning to the defect
            call fbsim_precondSpaceDefect (rsimsolver%rpreconditioner,&
                rdefectPrimal,blocalsuccess,rsimsolver%rnonlinearIteration)
            
            ! Gather statistics
            rsimsolver%dtimeLinearSolver = rsimsolver%dtimeLinearSolver + &
                rsimsolver%rnonlinearIteration%dtimeLinearSolver
            rsimsolver%nlinearIterations = rsimsolver%nlinearIterations + &
                rsimsolver%rnonlinearIteration%nlinearIterations
            
            ! If bsuccess=false, the preconditioner had an error.
            if (.not. blocalsuccess) then
              if (rsimSolver%rnonlinearIteration%ioutputLevel .ge. 0) then
                call output_line ('fbsim_simulate: Iteration '//&
                    trim(sys_siL(ite,10))//' canceled as the preconditioner went down!',&
                    coutputMode=rsimSolver%rnonlinearIteration%coutputMode)
              end if
              if (rsimsolver%rpreconditionerAlternative%ctypePreconditioner .eq. 0) then
                blocalsuccess = .false.
                exit
              else
                if (rsimSolver%rnonlinearIteration%ioutputLevel .ge. 0) then
                  call output_line ('fbsim_simulate: Trying alternative preconditioner...',&
                      coutputMode=rsimSolver%rnonlinearIteration%coutputMode)
                end if
                
                ! Restore the backup of the defect.
                call lsysbl_copyVector (rdefectBackup,rdefect)
                
                ! Prepeare the alternative preconditioner.
                call fbsim_updateDiscreteBCprec (rsimsolver%p_rsettings%rglobalData,&
                    rsimsolver%rpreconditionerAlternative,dtimePrimal,dtimeDual)

                ! Assemble the preconditioner matrices on all levels.
                call stat_startTimer (rtimerMatrixAssembly)
                call fbsim_assemblePrecMatrices (rsimsolver%rpreconditionerAlternative,&
                    iiterate,0,rspaceTimeMatrixPrecond,rnonlinearData,.true.)
                call stat_stopTimer (rtimerMatrixAssembly)

                ! Do preconditioning to the defect
                call fbsim_precondSpaceDefect (rsimsolver%rpreconditionerAlternative,&
                    rdefectPrimal,blocalsuccess,rsimsolver%rnonlinearIteration)
                
                ! Gather statistics
                rsimsolver%dtimeLinearSolver = rsimsolver%dtimeLinearSolver + &
                    rsimsolver%rnonlinearIteration%dtimeLinearSolver
                rsimsolver%nlinearIterations = rsimsolver%nlinearIterations + &
                    rsimsolver%rnonlinearIteration%nlinearIterations

                if (.not. blocalsuccess) then
                  if (rsimSolver%rnonlinearIteration%ioutputLevel .ge. 0) then
                    call output_line ('fbsim_simulate: Iteration '//&
                        trim(sys_siL(ite,10))//' canceled as the alternative preconditioner went down!',&
                        coutputMode=rsimSolver%rnonlinearIteration%coutputMode)
                  end if
                  blocalsuccess = .false.
                  exit
                end if

              end if
            end if

            ! Update of the solution
            ! --------------------------
            
            if (blocalsuccess) then
              ! We do not damp here.
              rsimSolver%rnonlinearIteration%domegaNL = 1.0_DP
            
              ! Combine to get the new solution vector.
              call lsyssc_vectorLinearComb(&
                  rdefectPrimal%RvectorBlock(1),p_rcurrentsol%RvectorBlock(1),1.0_DP,1.0_DP)
              call lsyssc_vectorLinearComb(&
                  rdefectPrimal%RvectorBlock(2),p_rcurrentsol%RvectorBlock(2),1.0_DP,1.0_DP)
              call lsyssc_vectorLinearComb(&
                  rdefectPrimal%RvectorBlock(3),p_rcurrentsol%RvectorBlock(3),1.0_DP,1.0_DP)
                  
            else
              exit
            end if
            
            ! Calculate the max-norm of the correction vector.
            ! This is used for the stopping criterium in resNormCheck!
            Cnorms(:) = LINALG_NORMMAX
            call lsysbl_vectorNormBlock(rdefectPrimal,Cnorms,&
                rsimsolver%rnonlinearIteration%DresidualCorr)

            ! Create the new nonlinear defect
            call stat_startTimer(rtimerDefectCalc)
            call fbsim_getNonlinearDefect (rsimSolver%rpreconditioner%cspace,p_rspaceTimeMatrix,&
                iiterate,rsimsolver%rdiscrData,rnonlinearData,p_rcurrentsol,rrhs,rdefect)
            call stat_stopTimer(rtimerDefectCalc)

            ! Implement the boundary conditions        
            call vecfil_discreteBCdef(rdefectPrimal)

            ! Check for convergence/divergence        
            call fbsim_resNormCheck (rsimSolver%rnonlinearIteration,&
                ite,p_rcurrentsol,p_rcurrentrhs,rdefect,bconvergence,bdivergence)
          
            ! Perform at least nminIterations iterations
            if (ite .lt. rsimSolver%rnonlinearIteration%nminIterations) bconvergence = .false.
          
            ! Check for convergence
            if (bconvergence) then
              blocalsuccess = .true.
              exit
            end if
          
            ! Check for divergence
            if (bdivergence) then
              if (rsimSolver%rnonlinearIteration%ioutputLevel .ge. 0) then
                call output_line ('fbsim_simulate: Iteration '//&
                    trim(sys_siL(ite,10))//' canceled, divergence detected!',&
                    coutputMode=rsimSolver%rnonlinearIteration%coutputMode)
              end if
              blocalsuccess = .false.
              exit
            end if
            
          end do ! ite
          
          call stat_stopTimer (rtimerNonlinearSolver)
          
        end if

        ! ite may be maxite+1 in case the loop passed through.
        rsimSolver%nnonlinearIterations = rsimSolver%nnonlinearIterations + &
            max(ite,rsimsolver%rnonlinearIteration%nmaxIterations)
        
        ! Nonlinear loop finished. Now save the solution or stop incase of an error.

        if (blocalsuccess) then
          ! Save the solution.
          ! Write the solution back from the pool into the space-time vector.
          ! p_rcurrentsol is connected to the pool vector with index iiterate.
          call sptivec_commitVecInPool (rsolAccess,iiterate)

          ! Postprocessing
          call stat_startTimer (rtimerPostProc)
          call fbsim_postprocessing (rsimSolver%rpostprocessing,p_rcurrentsol,iiterate,&
              dtimePrimal,dtimeDual,rsimsolver%p_rsettings)
          call stat_stopTimer (rtimerPostProc)
        else  
          exit
        end if
      
      end do ! iiterate
  
      ! Release the temporary RHS.    
      call lsysbl_releaseVector (rrhs)
      
    case (FBSIM_SOLVER_LINFORWARD)
    
      ! Oseen forward simulation in the primal space.
      NEQtime = ilastinterval+1
      
      ! Assign BC`s to the defect vector.
      call lsysbl_assignDiscreteBC (rdefectPrimal,&
          rsimSolver%rpreconditioner%p_RdiscreteBC(rsimSolver%ilevel))
      call lsysbl_assignDiscreteFBC (rdefectPrimal,&
          rsimSolver%rpreconditioner%p_RdiscreteFBC(rsimSolver%ilevel))
      
      ! Loop through the timesteps.
      do iiterate = ifirstinterval,NEQtime
      
        ! Current time step?
        call tdiscr_getTimestep(p_rspaceTimeMatrix%rdiscrData%p_rtimeDiscr,iiterate-1,&
            dtimePrimal,dtstep)
        dtimeDual = dtimePrimal - (1.0_DP-p_rspaceTimeMatrix%rdiscrData%p_rtimeDiscr%dtheta)*dtstep

        if (rsimSolver%ioutputLevel .ge. 1) then
          call output_line ("fbsim_simulate: Forward Iteration "//&
              trim(sys_siL(iiterate,10))//" of ["//&
              trim(sys_siL(ifirstinterval,10))//".."//&
              trim(sys_siL(NEQtime,10))//"], Time="//&
              trim(sys_sdL(dtimePrimal,10)),&
              coutputMode=rsimSolver%rnonlinearIteration%coutputMode)
        end if

        ! Update the boundary conditions to the current time
        call fbsim_updateDiscreteBCprec (rsimsolver%p_rsettings%rglobalData,&
            rsimsolver%rpreconditioner,dtimePrimal,dtimeDual)

        call stat_startTimer(rtimerDefectCalc)
        if (iiterate .gt. ifirstinterval) then

          ! Get the solution from the access pool.
          call sptivec_getVectorFromPool(rsolAccess, iiterate-1, p_rprevsol)
          call sptivec_getVectorFromPool(rsolAccess, iiterate, p_rcurrentsol)

          ! RHS.
          call sptivec_getVectorFromPool(rrhsAccess, iiterate-1, p_rprevrhs)
          call sptivec_getVectorFromPool(rrhsAccess, iiterate, p_rcurrentrhs)
          
          ! Get the Oseen solutions for the assembly of the nonlinearity -- 
          ! if there is any.
          call sptivec_getVectorFromPool(roseenAccess, iiterate-1, p_roseensol1)
          call sptivec_getVectorFromPool(roseenAccess, iiterate, p_roseensol2)
          if (iiterate .lt. NEQtime) then
            call sptivec_getVectorFromPool(roseenAccess, iiterate+1, p_roseensol3)
          end if
          
          ! Create a nonlinear-data structure that specifies the evaluation point
          ! of nonlinear matrices. the evaluation point is given by the three vectors
          ! roseensolX.
          call smva_initNonlinearData (rnonlinearData,p_roseensol1,p_roseensol2,p_roseensol3)
          
          ! Calculate the defect. For that purpose, start with the
          ! defect. Remember that the time stepping scheme is already
          ! incorporated to our RHS, so we do not have to take care of
          ! any timestepping weights concerning the RHS!
          call lsysbl_copyVector (p_rcurrentrhs,rdefect)
        
          ! ... subtract the offdiagonal, corresponding to the previous timestep.
          call smva_initNonlinMatrix (rnonlinearSpatialMatrix,&
              rsimsolver%rdiscrData,rnonlinearData)
          call stlin_setupMatrixWeights (p_rspaceTimeMatrix,&
              iiterate,-1,rnonlinearSpatialMatrix)
          call stlin_disableSubmatrix (rnonlinearSpatialMatrix,2,1)
          call stlin_disableSubmatrix (rnonlinearSpatialMatrix,2,2)
          
          call smva_assembleDefect (rnonlinearSpatialMatrix,p_rprevsol,&
              rdefect,1.0_DP)
            
        else
        
          ! Initial condition is given.
          call sptivec_getVectorFromPool(rrhsAccess, iiterate, p_rcurrentrhs)
          call sptivec_getVectorFromPool(rsolAccess, iiterate, p_rcurrentsol)

          call sptivec_getVectorFromPool(roseenAccess, iiterate, p_roseensol2)
          call sptivec_getVectorFromPool(roseenAccess, iiterate+1, p_roseensol3)
          p_roseensol1 => p_roseensol2 ! dummy
          
          ! Create a nonlinear-data structure that specifies the evaluation point
          ! of nonlinear matrices. the evaluation point is given by the three vectors
          ! roseensolX.
          call smva_initNonlinearData (rnonlinearData,p_roseensol1,p_roseensol2,p_roseensol3)
          
          ! There is no previous timestep. Load the RHS into the defect vector
          call lsysbl_copyVector (p_rcurrentrhs,rdefect)

        end if

        ! Subtract A*primal/dual solution from the rhs to create the primal defect.
        call smva_initNonlinMatrix (rnonlinearSpatialMatrix,&
            rsimsolver%rdiscrData,rnonlinearData)
        call stlin_setupMatrixWeights (p_rspaceTimeMatrix,&
            iiterate,0,rnonlinearSpatialMatrix)
        call stlin_disableSubmatrix (rnonlinearSpatialMatrix,2,1)
        call stlin_disableSubmatrix (rnonlinearSpatialMatrix,2,2)
        
        call smva_assembleDefect (rnonlinearSpatialMatrix,p_rcurrentsol,&
          rdefect,1.0_DP)

        ! Implement the boundary conditions        
        call vecfil_discreteBCdef(rdefectPrimal)
        call stat_stopTimer(rtimerDefectCalc)
        
        ! Assemble the preconditioner matrices on all levels.
        call stat_startTimer (rtimerMatrixAssembly)
        call fbsim_assemblePrecMatrices (rsimsolver%rpreconditioner,&
            iiterate,0,p_rspaceTimeMatrix,rnonlinearData,.true.)
        call stat_stopTimer (rtimerMatrixAssembly)
            
        ! Do preconditioning to the defect
        call fbsim_precondSpaceDefect (rsimsolver%rpreconditioner,&
            rdefectPrimal,blocalsuccess,rsimsolver%rnonlinearIteration)
            
        ! Gather statistics
        rsimsolver%dtimeLinearSolver = rsimsolver%dtimeLinearSolver + &
            rsimsolver%rnonlinearIteration%dtimeLinearSolver
        rsimsolver%nlinearIterations = rsimsolver%nlinearIterations + &
            rsimsolver%rnonlinearIteration%nlinearIterations

        if (rsimSolver%ioutputLevel .ge. 2) then
          call output_line ("fbsim_simulate: Forward Iteration "//&
              trim(sys_siL(iiterate,10))//". Residual damped from "//&
              trim(sys_sdEL(rsimsolver%rnonlinearIteration%dlastInitPrecDef,3))&
              //" to "//&
              trim(sys_sdEL(rsimsolver%rnonlinearIteration%dlastPrecDef,3)),&
              coutputMode=rsimSolver%rnonlinearIteration%coutputMode)
        end if
        
        if (blocalsuccess) then
          ! Combine to get the new solution vector.
          call lsyssc_vectorLinearComb(&
              rdefectPrimal%RvectorBlock(1),p_rcurrentsol%RvectorBlock(1),domega,1.0_DP)
          call lsyssc_vectorLinearComb(&
              rdefectPrimal%RvectorBlock(2),p_rcurrentsol%RvectorBlock(2),domega,1.0_DP)
          call lsyssc_vectorLinearComb(&
              rdefectPrimal%RvectorBlock(3),p_rcurrentsol%RvectorBlock(3),domega,1.0_DP)

          ! Write the solution back from the pool into the space-time vector.
          ! p_rcurrentsol is connected to the pool vector with index iiterate.
          call sptivec_commitVecInPool (rsolAccess,iiterate)

          ! Postprocessing          
          call stat_startTimer (rtimerPostProc)
          call fbsim_postprocessing (rsimSolver%rpostprocessing,p_rcurrentsol,iiterate,&
              dtimePrimal,dtimeDual,rsimsolver%p_rsettings)
          call stat_stopTimer (rtimerPostProc)

        else
          !exit
        end if
      
      end do ! iiterate
    
    case (FBSIM_SOLVER_LINBACKWARD)
    
      ! Oseen backward simulation.
      NEQtime = ilastinterval+1
      
      ! Assign BC`s to the defect vector.
      call lsysbl_assignDiscreteBC (rdefectDual,&
          rsimSolver%rpreconditioner%p_RdiscreteBC(rsimSolver%ilevel))
      call lsysbl_assignDiscreteFBC (rdefectDual,&
          rsimSolver%rpreconditioner%p_RdiscreteFBC(rsimSolver%ilevel))

      do iiterate = NEQtime,ifirstinterval,-1
      
        ! Current time step?
        call tdiscr_getTimestep(p_rspaceTimeMatrix%rdiscrData%p_rtimeDiscr,iiterate,&
            dtstep=dtstep,dtimestart=dtimePrimal)
        dtimeDual = dtimePrimal - (1.0_DP-p_rspaceTimeMatrix%rdiscrData%p_rtimeDiscr%dtheta)*dtstep

        if (rsimSolver%ioutputLevel .ge. 1) then
          call output_line ("fbsim_simulate: Backward Iteration "//&
              trim(sys_siL(iiterate,10))//" of ["//&
              trim(sys_siL(ifirstinterval,10))//".."//&
              trim(sys_siL(NEQtime,10))//"], Time="//&
              trim(sys_sdL(dtimePrimal,10)),&
              coutputMode=rsimSolver%rnonlinearIteration%coutputMode)
        end if

        ! Update the boundary conditions to the current time
        call fbsim_updateDiscreteBCprec (rsimsolver%p_rsettings%rglobalData,&
            rsimsolver%rpreconditioner,dtimePrimal,dtimeDual)
        
        call stat_startTimer(rtimerDefectCalc)
        if (iiterate .lt. NEQtime) then

          ! Initialise the previous and current timestep. 
          ! Solution.
          call sptivec_getVectorFromPool(rsolAccess, iiterate+1, p_rprevsol)
          call sptivec_getVectorFromPool(rsolAccess, iiterate, p_rcurrentsol)

          ! RHS.
          call sptivec_getVectorFromPool(rrhsAccess, iiterate+1, p_rprevrhs)
          call sptivec_getVectorFromPool(rrhsAccess, iiterate, p_rcurrentrhs)
          
          ! Get the Oseen solutions for the assembly of the nonlinearity -- 
          ! if there is any.
          call sptivec_getVectorFromPool(roseenAccess, iiterate+1, p_roseensol3)
          call sptivec_getVectorFromPool(roseenAccess, iiterate, p_roseensol2)
          if (iiterate .gt. ifirstinterval) then
            call sptivec_getVectorFromPool(roseenAccess, iiterate-1, p_roseensol1)
          end if

          ! Create a nonlinear-data structure that specifies the evaluation point
          ! of nonlinear matrices. the evaluation point is given by the three vectors
          ! roseensolX.
          call smva_initNonlinearData (rnonlinearData,p_roseensol1,p_roseensol2,p_roseensol3)

          ! Calculate the defect. For that purpose, start with the
          ! defect. Remember that the time stepping scheme is already
          ! incorporated to our RHS, so we do not have to take care of
          ! any timestepping weights concerning the RHS!
          call lsysbl_copyVector (p_rcurrentrhs,rdefect)

          ! ... subtract the offdiagonal, corresponding to the next timestep.
          call smva_initNonlinMatrix (rnonlinearSpatialMatrix,&
              rsimsolver%rdiscrData,rnonlinearData)
          call stlin_setupMatrixWeights (p_rspaceTimeMatrix,&
              iiterate,1,rnonlinearSpatialMatrix)
          call stlin_disableSubmatrix (rnonlinearSpatialMatrix,1,1)
          call stlin_disableSubmatrix (rnonlinearSpatialMatrix,1,2)
          
          call smva_assembleDefect (rnonlinearSpatialMatrix,p_rprevsol,&
              rdefect,1.0_DP)
            
        else
            
          ! Terminal condition is given.
          call sptivec_getVectorFromPool(rrhsAccess, NEQtime, p_rcurrentrhs)
          call sptivec_getVectorFromPool(rsolAccess, NEQtime, p_rcurrentsol)

          call sptivec_getVectorFromPool(roseenAccess, NEQtime, p_roseensol2)
          call sptivec_getVectorFromPool(roseenAccess, NEQtime-1, p_roseensol1)
          p_roseensol3 => p_roseensol2 ! dummy

          ! Create a nonlinear-data structure that specifies the evaluation point
          ! of nonlinear matrices. the evaluation point is given by the three vectors
          ! roseensolX.
          call smva_initNonlinearData (rnonlinearData,p_roseensol1,p_roseensol2,p_roseensol3)

          ! There is next timestep. Load the RHS into the defect vector
          call lsysbl_copyVector (p_rcurrentrhs,rdefect)

        end if

        ! Subtract A*primal/dual solution from the rhs to create the primal defect.
        call smva_initNonlinMatrix (rnonlinearSpatialMatrix,&
            rsimsolver%rdiscrData,rnonlinearData)
        call stlin_setupMatrixWeights (p_rspaceTimeMatrix,&
            iiterate,0,rnonlinearSpatialMatrix)
        call stlin_disableSubmatrix (rnonlinearSpatialMatrix,1,1)
        call stlin_disableSubmatrix (rnonlinearSpatialMatrix,1,2)
        
        call smva_assembleDefect (rnonlinearSpatialMatrix,p_rcurrentsol,&
            rdefect,1.0_DP)
        call stat_stopTimer(rtimerDefectCalc)
        
        ! Implement the boundary conditions        
        call vecfil_discreteBCdef(rdefectDual)

        ! Assemble the preconditioner matrices on all levels.
        call stat_startTimer (rtimerMatrixAssembly)
        call fbsim_assemblePrecMatrices (rsimsolver%rpreconditioner,&
            iiterate,0,p_rspaceTimeMatrix,rnonlinearData,.true.)
        call stat_stopTimer (rtimerMatrixAssembly)
            
        ! Do preconditioning to the defect
        call fbsim_precondSpaceDefect (rsimsolver%rpreconditioner,&
            rdefectDual,blocalsuccess,rsimsolver%rnonlinearIteration)
        
        ! Gather statistics
        rsimsolver%dtimeLinearSolver = rsimsolver%dtimeLinearSolver + &
            rsimsolver%rnonlinearIteration%dtimeLinearSolver
        rsimsolver%nlinearIterations = rsimsolver%nlinearIterations + &
            rsimsolver%rnonlinearIteration%nlinearIterations

        if (rsimSolver%ioutputLevel .ge. 2) then
          call output_line ("fbsim_simulate: Backward Iteration "//&
              trim(sys_siL(iiterate,10))//". Residual damped from "//&
              trim(sys_sdEL(rsimsolver%rnonlinearIteration%dlastInitPrecDef,3))&
              //" to "//&
              trim(sys_sdEL(rsimsolver%rnonlinearIteration%dlastPrecDef,3)),&
              coutputMode=rsimSolver%rnonlinearIteration%coutputMode)
        end if

        if (blocalsuccess) then
          ! Combine to get the new solution vector.
          call lsyssc_vectorLinearComb(&
              rdefectDual%RvectorBlock(1),p_rcurrentsol%RvectorBlock(4),domega,1.0_DP)
          call lsyssc_vectorLinearComb(&
              rdefectDual%RvectorBlock(2),p_rcurrentsol%RvectorBlock(5),domega,1.0_DP)
          call lsyssc_vectorLinearComb(&
              rdefectDual%RvectorBlock(3),p_rcurrentsol%RvectorBlock(6),domega,1.0_DP)

          ! Write the solution back from the pool into the space-time vector.
          ! p_rcurrentsol is connected to the pool vector with index iiterate.
          call sptivec_commitVecInPool (rsolAccess,iiterate)

          ! Postprocessing          
          call stat_startTimer (rtimerPostProc)
          call fbsim_postprocessing (rsimSolver%rpostprocessing,p_rcurrentsol,iiterate,&
              dtimePrimal,dtimeDual,rsimsolver%p_rsettings)
          call stat_stopTimer (rtimerPostProc)
          
        else
          !exit
        end if
      
      end do ! iiterate
      
    end select
    
    call stat_stopTimer(rtimerTotal)
    
    ! Gather statistics
    rsimSolver%dtimeNonlinearSolver = rtimerNonlinearSolver%delapsedReal
    rsimSolver%dtimeMatrixAssembly = rtimerMatrixAssembly%delapsedReal
    rsimSolver%dtimeTotal = rtimerTotal%delapsedReal
    rsimSolver%dtimePostprocessing = rtimerPostProc%delapsedReal
    rsimSolver%dtimeDefectCalculation = rtimerDefectCalc%delapsedReal

    if (present(bsuccess)) then
      bsuccess = blocalsuccess
    end if
    
    ! Release postprocessing stuff.
    call fbsim_donepostprocessing (rsimsolver%rpostprocessing)
    
    ! Release the temp vectors
    call lsysbl_releaseVector (rdefect)
    call lsysbl_releaseVector (rdefectBackup)

    if (present(roseenSolution)) then
      call sptivec_releaseAccessPool(roseenAccess)
    end if
    call sptivec_releaseAccessPool(rsolAccess)
    call sptivec_releaseAccessPool(rrhsAccess)

  end subroutine

end module
