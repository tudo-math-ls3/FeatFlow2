!##############################################################################
!# ****************************************************************************
!# <name> structuresoptc </name>
!# ****************************************************************************
!#
!# <purpose>
!# Underlying structures of the space-time optimal control solver.
!#
!# Routines in this module:
!#
!# 1.) soptc_initParOptControl
!#     -> Reads the parameters for the optimal control from the parameter list
!#        and stores them in the structure.
!#
!# 2.) soptc_doneParOptControl
!#     -> Releases data allocated in soptc_initParOptControl.
!#
!# 3.) soptc_initStepPathFollower
!#     -> Adapts the parameters in the optimal control structure
!#        according to path following information.
!# </purpose>
!##############################################################################

module structuresoptc

  use fsystem
  use storage
  use genoutput
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
    
    ! Discrete constraints for u_1, u_2
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
    
    ! Discrete constraints for u_1, u_2
    type(t_vectorBlock), pointer :: p_rvectorymin => null()
    type(t_vectorBlock), pointer :: p_rvectorymax => null()
    
    
    ! Time used for the evaluation of the analytical control constraints
    ! p_rumin1/2, p_rumax1/2.
    real(DP) :: dstateConstrTime = 0.0_DP

    ! Time used for the evaluation of the analytical state constraints
    ! p_rymin1/2, p_rymax1/2.
    real(DP) :: dcontrolConstrTime = 0.0_DP
    
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

  ! This structure collects a set of path following data.
  ! During the nonlinear iteration, this data is automatically
  ! adapted from nonlinear to nonlinear step.
  type t_optcpathfollowingdata
    
    ! List of $\alpha$ parameter of the optimal control functional parameters. 
    ! If this is not =>NULL, a path following algorithm will be applied. 
    ! Every nonlinear step, the next value is taken.
    real(DP), dimension(:), pointer :: p_DalphaC => null()
    
    ! List of $\beta$ parameter of the optimal control functional parameters. 
    ! If this is not =>NULL, a path following algorithm will be applied. 
    ! Every nonlinear step, the next value is taken.
    real(DP), dimension(:), pointer :: p_DbetaC => null()
    
    ! List of $\gamma$ parameter of the optimal control functional parameters. 
    ! If this is not =>NULL, a path following algorithm will be applied. 
    ! Every nonlinear step, the next value is taken.
    real(DP), dimension(:), pointer :: p_DgammaC => null()
  
  end type

!</typeblock>

  public :: t_optcpathfollowingdata

!<typeblock>

  ! This block saves parameters of the optimal control problem
  type t_settings_optcontrol
  
    ! $\alpha$ parameter of the optimal control functional
    real(DP) :: dalphaC = 1.0_DP
    
    ! $\beta$ parameter of the optimal control functional via boundary control
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

    !<!-- ---------------- -->
    !<!-- OBSERVATION AREA -->
    !<!-- ---------------- -->

    ! Contains data for path following strategies
    type(t_optcpathfollowingdata) :: rpathfollowingdata
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
  public :: soptc_initParOptControl
  public :: soptc_doneParOptControl
  public :: soptc_initStepPathFollower

contains

  ! ***************************************************************************

!<subroutine>

  subroutine sys_getParameter (sstring,stoken,istart)
  
!<description>
  ! Gets a parameter out of a list of parameters.
!</description>

!<input>
  ! A string with parameters, e.g. "0.1 0.2 0.3",...
  character(len=*), intent(in) :: sstring
  
  ! String that receives the parameter.
  character(len=*), intent(out) :: stoken
!</input>

!<inputoutput>
  ! On input: Must be set to =1 for the first call.
  !   Position from where to search for the next parameter
  ! On output:
  !   =0, if this was the last parameter
  !   Otherwise: position in sstring where the next parameter starts.
  integer, intent(inout) :: istart
!</inputoutput>

!</subroutine>

    integer :: i,slen
    
    ! If istart=0, do not do anything. There is nothing in the string
    stoken = ""
    if (istart .le. 0) return
    
    ! If the string length is =0, finish.
    slen = len_trim(sstring)
    if (slen .eq. 0) return
    
    ! Skip whitespaces.
    do while (istart .le. slen)
      if (sstring(istart:istart) .eq. " ") then
        istart = istart + 1
      else
        exit
      end if
    end do
    
    ! End of the string?
    i = istart
    do while (i .le. slen)
      if (sstring(i:i) .ne. " ") then
        i = i + 1
      else
        exit
      end if
    end do
    
    ! Copy
    if (istart .le. slen) then
      stoken = sstring(istart:i)
    end if
    
    ! Put istart behind the parameter or set it to 0.
    if (i .ge. slen) then
      ! No more parameters
      istart = 0
    else
      ! Put istart behind the parameter -- on the whitespace which
      ! is skipped in the next loop.
      istart = i+1
    end if
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sys_countParameters (sstring,ntokens)
  
!<description>
  ! Returns the number of parameters in sstring
!</description>

!<input>
  ! A string with parameters, e.g. "0.1 0.2 0.3",...
  character(len=*), intent(in) :: sstring
!</input>

!<output>
  ! Number of parameters in sstring
  integer, intent(out) :: ntokens
!</output>

!</subroutine>

    integer :: slen, istart
    
    ! If the string length is =0, finish.
    ntokens = 0
    
    slen = len_trim(sstring)
    if (slen .eq. 0) return
    
    ! Find all substrings  
    istart = 1
    
    do
      ! Cancel if nothing left.
      if (istart .gt. slen) exit

      ! Skip whitespaces.
      do while (istart .le. slen)
        if (sstring(istart:istart) .eq. " ") then
          istart = istart + 1
        else
          exit
        end if
      end do
      
      ! Cancel if we reached the string end
      if (istart .gt. slen) exit
      
      ! End of the string?
      do while (istart .le. slen)
        if (sstring(istart:istart) .ne. " ") then
          istart = istart + 1
        else
          exit
        end if
      end do
      
      ! One more
      ntokens = ntokens + 1
      
    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine soptc_initParOptControl (rparlist,ssectionOptC,roptcontrol)
  
!<description>
  ! Reads the parameters for the optimal control from the parameter list
  ! and stores them in the structure.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section where the parameters of the optimal control can be found.
  character(len=*), intent(in) :: ssectionOptC
!</input>

!<inputoutput>
  ! Structure defining the parameters for the optimal control.
  type(t_settings_optcontrol), intent(inout) :: roptcontrol
!</inputoutput>

!</subroutine>
    
    character(len=SYS_STRLEN) :: sstring,stoken
    integer :: ntokens,istart,itoken

    ! Alpha/Gamma parameters.
    !
    ! These parameters may have two types: either it is a single
    ! parameter or a list of parameters. In the second case,
    ! this defines a path following strategy, i.e., the parameters
    ! change according to the list in every nonlinear iteration.
    call parlst_getvalue_string (rparlist,ssectionOptC,&
        "dalphaC",sstring,"-1.0",bdequote=.true.)

    ! Is this a list?
    call sys_countParameters (sstring,ntokens)
    
    if (ntokens .eq. 1) then
    
      read(sstring,*) roptcontrol%dalphaC
      
    else if (ntokens .gt. 1) then
      ! Read the list
      roptcontrol%dalphaC = 0.0_DP
      allocate(roptcontrol%rpathfollowingdata%p_DalphaC(ntokens))
      
      ! Loop through all tokens and save them
      istart = 1
      do itoken = 1,ntokens
        call sys_getParameter (sstring,stoken,istart)
        read(stoken,*) roptcontrol%rpathfollowingdata%p_DalphaC(itoken)
      end do
    end if
        
    call parlst_getvalue_string (rparlist,ssectionOptC,&
        "dbetaC",sstring,"-1.0",bdequote=.true.)

    ! Is this a list?
    call sys_countParameters (sstring,ntokens)

    if (ntokens .eq. 1) then
    
      read(sstring,*) roptcontrol%dbetaC
      
    else if (ntokens .gt. 1) then
      ! Read the list
      roptcontrol%dbetaC = 0.0_DP
      allocate(roptcontrol%rpathfollowingdata%p_DbetaC(ntokens))
      
      ! Loop through all tokens and save them
      istart = 1
      do itoken = 1,ntokens
        call sys_getParameter (sstring,stoken,istart)
        read(stoken,*) roptcontrol%rpathfollowingdata%p_DbetaC(itoken)
      end do
    end if

    call parlst_getvalue_string (rparlist,ssectionOptC,&
        "dgammaC",sstring,"0.0",bdequote=.true.)

    ! Is this a list?
    call sys_countParameters (sstring,ntokens)

    if (ntokens .eq. 1) then
    
      read(sstring,*) roptcontrol%dgammaC

    else if (ntokens .gt. 1) then
      ! Read the list
      roptcontrol%dgammaC = 0.0_DP
      allocate(roptcontrol%rpathfollowingdata%p_DgammaC(ntokens))
      
      ! Loop through all tokens and save them
      istart = 1
      do itoken = 1,ntokens
        call sys_getParameter (sstring,stoken,istart)
        read(stoken,*) roptcontrol%rpathfollowingdata%p_DgammaC(itoken)
      end do
    end if

            
    call parlst_getvalue_double (rparlist,ssectionOptC,&
        "ddirichletBCPenalty",roptcontrol%ddirichletBCPenalty,100.0_DP)
        
    ! Type of the formulation
    call parlst_getvalue_int (rparlist,ssectionOptC,&
        "ispaceTimeFormulation",roptcontrol%ispaceTimeFormulation,0)
    
    call parlst_getvalue_int (rparlist,ssectionOptC,&
        "iconvectionExplicit",roptcontrol%iconvectionExplicit,0)

    call parlst_getvalue_int (rparlist,ssectionOptC,&
        "csystemScaling",roptcontrol%csystemScaling,0)
    
    ! Parameters defining the control constraints
    call parlst_getvalue_int (rparlist,ssectionOptC,&
        "ccontrolConstraints",roptcontrol%rconstraints%ccontrolConstraints,0)

    roptcontrol%rconstraints%ccontrolConstraintsType = 0

    call parlst_getvalue_double (rparlist,ssectionOptC,&
        "dumin1",roptcontrol%rconstraints%dumin1,-1.0E10_DP)

    call parlst_getvalue_double (rparlist,ssectionOptC,&
        "dumax1",roptcontrol%rconstraints%dumax1,1.0E10_DP)

    call parlst_getvalue_double (rparlist,ssectionOptC,&
        "dumin2",roptcontrol%rconstraints%dumin2,-1.0E10_DP)

    call parlst_getvalue_double (rparlist,ssectionOptC,&
        "dumax2",roptcontrol%rconstraints%dumax2,1.0E10_DP)

    ! Parameters defining the state constraints
    call parlst_getvalue_int (rparlist,ssectionOptC,&
        "cstateConstraints",roptcontrol%rconstraints%cstateConstraints,0)

    roptcontrol%rconstraints%cstateConstraintsType = 0

    call parlst_getvalue_double (rparlist,ssectionOptC,&
        "dymin1",roptcontrol%rconstraints%dymin1,-1.0E10_DP)

    call parlst_getvalue_double (rparlist,ssectionOptC,&
        "dymax1",roptcontrol%rconstraints%dymax1,1.0E10_DP)

    call parlst_getvalue_double (rparlist,ssectionOptC,&
        "dymin2",roptcontrol%rconstraints%dymin2,-1.0E10_DP)

    call parlst_getvalue_double (rparlist,ssectionOptC,&
        "dymax2",roptcontrol%rconstraints%dymax2,1.0E10_DP)

    call parlst_getvalue_double (rparlist,ssectionOptC,&
        "dstateConstrReg",roptcontrol%rconstraints%dstateConstrReg,1.0_DP)
        
    ! Observation area
    call parlst_getvalue_string (rparlist,ssectionOptC,&
        "DobservationArea",sstring,"")
        
    if (sstring .ne. "") then
      ! Read the observation area. This is a box, the parameters
      ! have the format "x1, y1, x2, y2".
      allocate (roptcontrol%p_DobservationArea(4))
      read (sstring,*) &
          roptcontrol%p_DobservationArea(1), &
          roptcontrol%p_DobservationArea(2), &
          roptcontrol%p_DobservationArea(3), &
          roptcontrol%p_DobservationArea(4)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine soptc_doneParOptControl (roptcontrol)
  
!<description>
  ! Cleans up information in the optimal control structure.
!</description>

!<inputoutput>
  ! Structure defining the parameters for the optimal control.
  type(t_settings_optcontrol), intent(inout) :: roptcontrol
!</inputoutput>

!</subroutine>

    ! Release the observation area
    if (associated(roptcontrol%p_DobservationArea)) then
      deallocate(roptcontrol%p_DobservationArea)
    end if
    
    ! Release path following data
    if (associated(roptcontrol%rpathfollowingdata%p_DalphaC)) then
      deallocate(roptcontrol%rpathfollowingdata%p_DalphaC)
    end if

    if (associated(roptcontrol%rpathfollowingdata%p_DbetaC)) then
      deallocate(roptcontrol%rpathfollowingdata%p_DbetaC)
    end if

    if (associated(roptcontrol%rpathfollowingdata%p_DgammaC)) then
      deallocate(roptcontrol%rpathfollowingdata%p_DgammaC)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine soptc_initStepPathFollower (roptcontrol,istep,bprint)
  
!<description>
  ! Adapts the parameters in the optimal control structure
  ! according to path following information.
  ! The step istep defines the current setting of the parameters
  ! which are adapted via path following methods.
!</description>

!<input>
  ! Current step in the path following strategy
  integer, intent(in) :: istep
  
  ! Whether or not to print information about the path following parameters.
  logical, intent(in) :: bprint
!</input>

!<inputoutput>
  ! Structure defining the parameters for the optimal control.
  ! The parameters in this structure which are influenced
  ! via the path following data are modified according to
  ! step istep.
  type(t_settings_optcontrol), intent(inout) :: roptcontrol
!</inputoutput>

!</subroutine>

    ! Change the data in roptcontrol.
    ! istep defines the array position of the element to take.
    if (associated(roptcontrol%rpathfollowingdata%p_DalphaC)) then
      roptControl%dalphaC = roptControl%rpathfollowingdata%p_DalphaC(&
              min(size(roptcontrol%rpathfollowingdata%p_DalphaC),istep))
      if (bprint) then
        call output_line ("Path-Follower: Parameter ALPHA = "//&
            trim(sys_sdEL(roptcontrol%dalphaC,10)))
      end if
    end if

    if (associated(roptcontrol%rpathfollowingdata%p_DbetaC)) then
      roptcontrol%dbetaC = roptControl%rpathfollowingdata%p_DbetaC(&
              min(size(roptcontrol%rpathfollowingdata%p_DbetaC),istep))
      if (bprint) then
        call output_line ("Path-Follower: Parameter BETA  = "//&
            trim(sys_sdEL(roptcontrol%dbetaC,10)))
      end if
    end if

    if (associated(roptcontrol%rpathfollowingdata%p_DgammaC)) then
      roptcontrol%dgammaC = roptControl%rpathfollowingdata%p_DgammaC(&
              min(size(roptcontrol%rpathfollowingdata%p_DgammaC),istep))
      if (bprint) then
        call output_line ("Path-Follower: Parameter GAMMA = "//&
            trim(sys_sdEL(roptcontrol%dgammaC,10)))
      end if
    end if

  end subroutine

end module
