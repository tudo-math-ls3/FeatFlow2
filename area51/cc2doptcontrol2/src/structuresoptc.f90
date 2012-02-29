!##############################################################################
!# ****************************************************************************
!# <name> structuresoptc </name>
!# ****************************************************************************
!#
!# <purpose>
!# Underlying structures of the space-time optimal control solver.
!# </purpose>
!##############################################################################

module structuresoptc

  use fsystem
  use storage
  use boundary
  use triangulation
  use paramlist
  use discretebc
  use discretefbc
  use fparser
  
  use collection
  use paramlist
  use linearsystemblock

  use analyticsolution
  
  use timediscretisation
  
  implicit none
  
  private
  
!<types>

!<typeblock>

  ! Type block encapsuling constraints, only in space.
  type t_optcconstraintsSpace
  
    ! Type of constraints to apply to the control u.
    ! =0: No constraints.
    ! =1: Constant constraints on u active:
    !     dumin1 <= u_1 <= dumax1, dumin2 <= u_2 <= dumax2
    !     Implementation by DOF.
    ! =2: Constant constraints on u active:
    !     dumin1 <= u_1 <= dumax1, dumin2 <= u_2 <= dumax2
    !     Implemented by cubature point.
    ! =3: Constant constraints on u active: 
    !     dumin1 <= u_1 <= dumax1, dumin2 <= u_2 <= dumax2.
    !     Implemented by DOF. Newton matrix is approximative.
    ! =4: Constant constraints on u active: 
    !     dumin1 <= u_1 <= dumax1, dumin2 <= u_2 <= dumax2.
    !     Implemented by cubature point (more exact). Adaptive integration.
    integer :: ccontrolConstraints = 0

    ! Type of definition of the constraints if ccontrolConstraints <> 0.
    ! =0: constants specified in dumin1/2, dumax1/2.
    ! =1: analytical functions defined in p_rumin1/2, p_rumax1/2.
    integer :: ccontrolConstraintsType = 0

    ! Constraints on u_1
    real(DP) :: dumin1 = -1.0E10
    real(DP) :: dumax1 = 1.0E10
    
    ! Constraints in u_2
    real(DP) :: dumin2 = -1.0E10
    real(DP) :: dumax2 = 1.0E10

    ! Analytical constraints for u_1
    type(t_anSolution), pointer :: p_rumin1 => null()
    type(t_anSolution), pointer :: p_rumax1 => null()
    
    ! Analytical constraints for u_2
    type(t_anSolution), pointer :: p_rumin2 => null()
    type(t_anSolution), pointer :: p_rumax2 => null()
    
    ! Discrete constraints for u_1
    type(t_vectorBlock), pointer :: p_rvectorumin => null()
    type(t_vectorBlock), pointer :: p_rvectorumax => null()


    ! Type of constraints to apply to the state y.
    ! =0: No constraints.
    ! =1: Constant constraints on u active:
    !     dymin1 <= y_1 <= dymax1, dymin2 <= y_2 <= dymax2
    !     Implementation by DOF.
    ! =2: Constant constraints on y active:
    !     dymin1 <= y_1 <= dymax1, dymin2 <= y_2 <= dymax2
    !     Implemented by cubature point.
    integer :: cstateConstraints = 0

    ! Type of definition of the constraints if ccontrolConstraints <> 0.
    ! =0: constants specified in dumin1/2, dumax1/2.
    ! =1: analytical functions defined in p_rumin1/2, p_rumax1/2.
    integer :: cstateConstraintsType = 0

    ! Regularisation parameter for the Moreau-Yosida regularisation
    real(DP) :: dstateConstrReg = 1.0_DP

    ! Constraints on u_1
    real(DP) :: dymin1 = -1.0E10
    real(DP) :: dymax1 = 1.0E10
    
    ! Constraints in u_2
    real(DP) :: dymin2 = -1.0E10
    real(DP) :: dymax2 = 1.0E10

    ! Analytical constraints for u_1
    type(t_anSolution), pointer :: p_rymin1 => null()
    type(t_anSolution), pointer :: p_rymax1 => null()
    
    ! Analytical constraints for u_2
    type(t_anSolution), pointer :: p_rymin2 => null()
    type(t_anSolution), pointer :: p_rymax2 => null()
    
    ! Discrete constraints for u_1
    type(t_vectorBlock), pointer :: p_rvectorymin => null()
    type(t_vectorBlock), pointer :: p_rvectorymax => null()
    
    
    ! Time used for the evaluation of the analytical constraints
    ! p_rumin1/2, p_rumax1/2.
    real(DP) :: dconstrainsTime = 0.0_DP
    
  end type
  
!</types>

  public :: t_optcconstraintsSpace

!<types>

  ! Type block encapsuling constraints, in space/time.
  type t_optcconstraintsSpaceTime
  
    ! Type of constraints to apply to the control u.
    ! =0: No constraints.
    ! =1: Constant constraints on u active:
    !     dumin1 <= u_1 <= dumax1, dumin2 <= u_2 <= dumax2
    !     Implementation by DOF.
    ! =2: Constant constraints on u active:
    !     dumin1 <= u_1 <= dumax1, dumin2 <= u_2 <= dumax2
    !     Implemented by cubature point.
    ! =3: Constant constraints on u active: 
    !     dumin1 <= u_1 <= dumax1, dumin2 <= u_2 <= dumax2.
    !     Implemented by DOF. Newton matrix is approximative.
    ! =4: Constant constraints on u active: 
    !     dumin1 <= u_1 <= dumax1, dumin2 <= u_2 <= dumax2.
    !     Implemented by cubature point (more exact). Adaptive integration.
    integer :: ccontrolConstraints = 0

    ! Type of definition of the constraints if ccontrolConstraints <> 0.
    ! =0: constants specified in dumin1/2, dumax1/2.
    ! =1: analytical functions defined in p_rumin1/2, p_rumax1/2.
    integer :: ccontrolConstraintsType = 0

    ! Constraints on u_1
    real(DP) :: dumin1 = -1.0E10
    real(DP) :: dumax1 = 1.0E10

    ! Analytical constraints for u_1.
    type(t_anSolution), pointer :: p_rumin1 => null()
    type(t_anSolution), pointer :: p_rumax1 => null()

    ! Constraints in u_2
    real(DP) :: dumin2 = -1.0E10
    real(DP) :: dumax2 = 1.0E10
  
    ! Analytical constraints for u_2
    type(t_anSolution), pointer :: p_rumin2 => null()
    type(t_anSolution), pointer :: p_rumax2 => null()
  
    ! Type of constraints to apply to the state y.
    ! =0: No constraints.
    ! =1: Constant constraints on u active:
    !     dymin1 <= y_1 <= dymax1, dymin2 <= y_2 <= dymax2
    !     Implementation by DOF.
    ! =2: Constant constraints on y active:
    !     dymin1 <= y_1 <= dymax1, dymin2 <= y_2 <= dymax2
    !     Implemented by cubature point.
    integer :: cstateConstraints = 0

    ! Type of definition of the constraints if ccontrolConstraints <> 0.
    ! =0: constants specified in dumin1/2, dumax1/2.
    ! =1: analytical functions defined in p_rumin1/2, p_rumax1/2.
    integer :: cstateConstraintsType = 0

    ! Regularisation parameter for the Moreau-Yosida regularisation
    real(DP) :: dstateConstrReg = 1.0_DP

    ! Constraints on u_1
    real(DP) :: dymin1 = -1.0E10
    real(DP) :: dymax1 = 1.0E10
    
    ! Constraints in u_2
    real(DP) :: dymin2 = -1.0E10
    real(DP) :: dymax2 = 1.0E10

    ! Analytical constraints for u_1
    type(t_anSolution), pointer :: p_rymin1 => null()
    type(t_anSolution), pointer :: p_rymax1 => null()
    
    ! Analytical constraints for u_2
    type(t_anSolution), pointer :: p_rymin2 => null()
    type(t_anSolution), pointer :: p_rymax2 => null()
    
  end type

!</typeblock>

  public :: t_optcconstraintsSpaceTime

!<typeblock>

  ! This block saves parameters of the optimal control problem
  type t_settings_optcontrol
  
    ! $\alpha$ parameter of the optimal control functional
    real(DP) :: dalphaC = 1.0_DP
    
    ! $\alpha$ parameter of the optimal control functional via boundary control
    real(DP) :: dbetaC = 1.0_DP
    
    ! Penalty parameter for the dirichlet boundary control
    real(DP) :: ddirichletBCPenalty = 100.0_DP
    
    ! $\gamma$ parameter of the nonstationary optimal control functional
    real(DP) :: dgammaC = 0.0_DP
  
    ! Formulation of the Space-time problem.
    ! =0: usual formulation as specified in the DFG applicance
    ! =1: Formulation for the generation of reference results from papers
    ! The two formulations differ in a "-"-sign in front of the dual velocity.
    integer :: ispaceTimeFormulation = 0
  
    ! Whether to treat the convection explicitly or implicitly.
    ! =0: Treat the convection implicitely.
    ! =1: Treat the convection explicitely. (This is the formulation of the paper
    !     of Baerwolff and Hinze!)
    integer :: iconvectionExplicit = 0

    ! Type of scaling of the global system.
    ! =0: no scaling.
    ! =1: System scaled by Delta(t) to normalise diagonal blocks.
    integer :: csystemScaling = 0

    !<!-- --------------- -->
    !<!-- TARGET FUNCTION -->
    !<!-- --------------- -->

    ! Analytic solution defining the target function
    type(t_anSolution) :: rtargetFunction
    
    !<!-- ----------- -->
    !<!-- CONSTRAINTS -->
    !<!-- ----------- -->
  
    ! Possible constraints in the problem.
    type(t_optcconstraintsSpaceTime) :: rconstraints

    !<!-- ---------------- -->
    !<!-- OBSERVATION AREA -->
    !<!-- ---------------- -->
    
    ! Observation area. If this points to NULL(), the whole
    ! domain is observed. Otherwise, this is an array with four entries
    ! in the form
    !   p_DobservationArea = (x1 y1 x2 y2)
    ! specifying the area to be observed.
    real(DP), dimension(:), pointer :: p_DobservationArea => null()
    
    
  end type

!</typeblock>

  public :: t_settings_optcontrol

!<typeblock>

  ! This type block encapsules all physical constants and configuration
  ! parameters for the primal equation. This includes e.g. the type of the equation,
  ! viscosity parameter etc.
  type t_settings_physics
  
    ! Type of problem.
    ! =0: Navier-Stokes 2D.
    ! =1: Stokes 2D.
    integer :: cequation
    
    ! Type of subproblem of the main problem. Depending on iequationType.
    ! If iequationType=0 or =1:
    ! =0: (Navier-)Stokes with gradient tensor
    ! =1: (Navier-)Stokes with deformation tensor
    integer :: isubEquation
    
    ! Model for the viscosity.
    ! =0: Constant viscosity.
    ! =1: Power law: nu = nu_0 * z^(dviscoexponent/2 - 1), nu_0 = 1/RE, z=||D(u)||^2+dviscoEps
    ! =2: Bingham fluid: nu = nu_0 + dviscoyield / sqrt(|D(u)||^2+dviscoEps^2), nu_0 = 1/RE
    integer :: cviscoModel
        
    ! Exponent parameter for the viscosity model
    real(DP) :: dviscoexponent

    ! Epsilon regularisation for the viscosity model
    real(DP) :: dviscoEps
    
    ! Yield stress for Bingham fluid
    real(DP) :: dviscoYield

    ! Viscosity parameter nu = 1/Re if cviscoModel=0 (constant viscosity).
    ! Otherwise not used.
    real(DP) :: dnuConst
      
  end type

!</typeblock>

  public :: t_settings_physics

!<typeblock>

  ! Type block that specifies settings about the refinement of a mesh.
  type t_settings_refinement
  
    ! Type of refinement.
    ! =0: 2-level refinement
    integer :: crefType = 0

    ! Number of pre-refinements with 2-level refinement to calculate the
    ! coarse mesh in rmeshHierarchy from rtriaCoarse.
    integer :: npreref = 0
    
    ! Total number of levels.
    integer :: nlevels = 0
    
    ! Output level that defines which information to print
    ! during the refinement.
    ! =-1: no information.
    ! = 0: only errors/warnings.
    ! = 1: basic information.
    ! = 2: standard information.
    integer :: coutputlevel = 2
    
  end type

!</typeblock>

  public :: t_settings_refinement

!<typeblock>

  ! This type block encapsules the settings for the stabilisation of the convection.
  type t_settings_stabil
  
    ! Type of stabilization of the convective term.
    ! 0=Streamline Diffusion
    ! 1=upwind
    ! 2=unified edge oriented jump stabilisation
    ! 3=Streamline Diffusion (new implementation)
    ! 4=unified edge oriented jump stabilisation (new implementation)
    !   iUpwind1 holds for the primal equation, iUpwind2 for the dual one.
    ! 5=UEO stabilisation with precomputed matrix and optimised matrix stencils
    integer :: cupwind = 0
    
    ! Relaxation parameter for upwind/Streamline Diffusion/Jump stabilisation.
    ! Standard values: Streamline diffusion=1.0, Upwind=0.1, Jump stabil=0.01.
    ! dUpsam1 holds for the primal equation, dUpsam2 for the dual one.
    real(DP) :: dupsam = 0.0_DP
    
    ! Defines whether or not the EOJ stabilisation is applied to the boundary.
    ! =0: Stabilisation is not applied to the boundary.
    ! =1: Stabilisation is applied to all bonudary edges (default).
    integer :: ceojStabilOnBoundary = 1
    
    ! Definbes whether the convection operator of the dual equation is included
    ! into the solution.
    ! =0: no, =1: yes
    integer :: cconvectionOnBoundaryDefect = 1

    ! Definbes whether the convection operator of the dual equation is set up
    ! for the preconditioners on the boundary.
    ! =0: no, =1: yes
    integer :: cconvectionOnBoundaryMatrix = 1

  end type

!</typeblock>

  public :: t_settings_stabil

!<typeblock>

  ! This type encapsules all debug flags that may be of use somewhere
  ! in the program.
  type t_optcDebugFlags
  
    ! If the following constant is set from 1.0 to 0.0, the primal system is
    ! decoupled from the dual system!
    real(dp) :: dprimalDualCoupling = 1.0

    ! If the following constant is set from 1.0 to 0.0, the dual system is
    ! decoupled from the primal system!
    real(dp) :: ddualPrimalCoupling = 1.0

    ! If the following parameter is set from 1.0 to 0.0, the terminal
    ! condition between the primal and dual equation is decoupled, i.e.
    ! the dual equation gets independent from the primal one.
    real(dp) :: dterminalCondDecoupled = 1.0

    ! If the following parameter is set from 1.0 to 0.0, the time coupling
    ! is disabled, resulting in a stationary simulation in every timestep.
    real(dp) :: dtimeCoupling = 1.0
    
    ! Timestep-scheme for the pressure.
    ! =0: Use the same timestep scheme as for the velocity.
    ! =1: Calculate the pressure fully implicitely (standard)
    ! No effect if Implicit-Euler is used.
    ! Difference can only be seen for higher order time-discretisation.
    integer :: ipressureFullyImplicit = 1
  
    ! Tell the space-time UMFPACK (if present in the lienar solver) to write
    ! out the global matrix.
    integer :: cumfpackWriteMatrix = 0

    ! Filename of the file in which UMFPACK should write the global matrix
    character(len=SYS_STRLEN) :: sumfpackMatrixFilename = "./matrix.txt"
  
    ! Additional weight for the convective operator. May be used to
    ! switch off the convection.
    ! Standard value = 1.0
    real(DP) :: dweightConvection = 1.0_DP

    ! Additional weight for the term -y grad(lambda) in the dual equation.
    ! Standard value = 1.0.
    real(DP) :: dweightDualConvection = 1.0_DP
    
    ! Additional weight for the term grad(y)^T lambda in the dual equation.
    ! Standard value = 1.0.
    real(DP) :: dweightDualNewtonT = 1.0_DP
    
    ! Additional weight for the term (y*n)lambda in the natural boundary condition
    ! of the dual equation.
    ! Standard value = 1.0
    real(DP) :: dweightNaturalBdcDual = 1.0_DP

    ! Modification to the discrete RHS.
    ! =0: No modification (standard).
    ! =1: Disturb the primal velocity RHS in all DOF's with a random value.
    ! =2: Disturb the primal velocity RHS in all DOF's except for the boundary DOF's
    !     with a random value.
    ! =3: Disturb the dual velocity RHS in all DOF's with a random value.
    ! =4: Disturb the dual velocity RHS in all DOF's except for the boundary DOF's
    !     with a random value.
    ! =5: Disturb the velocity RHS in all DOF's with a random value.
    ! =6: Disturb the velocity RHS in all DOF's except for the boundary DOF's
    !     with a random value.
    integer :: crhsmodification = 0

    ! Maximum error to be introduced to the RHS if crhsmodification=1/2.
    real(DP) :: drhsrandomMax = 1E-13
    
  end type

!</typeblock>

  public :: t_optcDebugFlags

!<typeblock>

  ! Type encapsuling the boundary conditions.
  type t_optcBDC
  
    ! Physics of the problem
    type(t_settings_physics), pointer :: p_rphysics => null()

    ! Name of the section in rparamList describing the boundary conditions
    character(len=SYS_STRLEN) :: ssectionBdExpressions = ""
    character(len=SYS_STRLEN) :: ssectionBdConditions = ""

    ! A param list that saves the boundary conditions.
    type(t_parlist), pointer :: p_rparamListBDC => null()

  end type

!</typeblock>

  public :: t_optcBDC

!<typeblock>

  ! A collection of global data which is passen to the callback routines
  ! in user_callback.f90
  type t_globalData
    ! An application specific parameter list.
    type(t_parlist), pointer :: p_rparlist => null()
    
    ! Pointer to the coarse time discretisation.
    type(t_timeDiscretisation), pointer :: p_rtimeCoarse => null()
    
    ! Reference to the physics parameters
    type(t_settings_physics), pointer :: p_rphysics => null()

    ! Reference to the optimal control parameters.
    type(t_settings_optcontrol), pointer :: p_rsettingsOptControl => null()
    
    ! Pointer to the right hand side.
    type(t_anSolution), pointer :: p_rrhs => null()

    ! Pointer to the target function.
    type(t_anSolution), pointer :: p_rtargetFunction => null()
  end type

!</typeblock>

  public :: t_globalData


!!<typeblock>
!
!  ! Application-specific type block for the Nav.St. problem
!  type t_problem
!
!    ! Output level during the initialisation phase.
!    integer :: MSHOW_Initialisation
!
!    ! Output level of the application.
!    integer :: MT_OutputLevel
!
!    ! Physics of the problem
!    type(t_settings_physics) :: rphysicsPrimal
!
!    ! A parameter block for everything that controls the optimal control
!    ! problem.
!    type(t_settings_optcontrol) :: roptcontrol
!
!    ! An object for saving the domain:
!    type(t_boundary) :: rboundary
!
!    ! A collection object that saves structural data and some
!    ! problem-dependent information which is e.g. passed to
!    ! callback routines.
!    type(t_collection)                    :: rcollection
!
!    ! A param list that saves all parameters from the DAT/INI file(s).
!    type(t_parlist)                       :: rparamList
!
!    ! A t_problem_oneshot structure that is only valid during the
!    ! execution of the command parser.
!    type(t_problem_oneshot)               :: rdataOneshot
!
!  end type
!
!!</typeblock>

!</types>

end module
