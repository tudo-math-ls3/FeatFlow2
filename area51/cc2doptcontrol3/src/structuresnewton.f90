!##############################################################################
!# ****************************************************************************
!# <name> structuresnewton </name>
!# ****************************************************************************
!#
!# <purpose>
!# General structures for Newton algorithms -- in space as well as in
!# space-time.
!# </purpose>
!##############################################################################

module structuresnewton

  use fsystem
  use genoutput
  use paramlist
  use linearalgebra
  use iterationcontrol
  
  implicit none
  
  private
  
!<constants>

!<constantblock description="Constants defining the partial Newton">

  ! Full Newton
  integer, parameter, public :: NEWTN_PN_FULLNEWTON = 0 

  ! Partial Newton
  integer, parameter, public :: NEWTN_PN_PARTIALNEWTON = 1

  ! Partial Newton in the dual equation, full Newton in the primal equation
  integer, parameter, public :: NEWTN_PN_PARTIALNEWTONDUAL = 2
  
!</constantblock>

!</constants>
  
!<type>

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
  
    ! Type of partial Newton to apply.
    ! =0: No partial Newton.
    ! =1: Use partial Newton until nmaxFixedPointIterations is reached.
    !     For primal/dual problems: Use partial Newton in the primal and 
    !     dual equation.
    ! =2: For primal/dual problems: Use full Newton in the primal and 
    !     partial Newton in the dual equation until nmaxFixedPointIterations 
    !     is reached.
    integer :: cpartialNewton = NEWTN_PN_FULLNEWTON

    ! Minimum number of partial Newton iterations before switching to
    ! the full Newton
    integer :: nminPartialNewtonIterations = 0

    ! Maximum number of partial Newton iterations before switching to
    ! the full Newton
    integer :: nmaxPartialNewtonIterations = 0

    ! Norm of absolute residuum before applying the full Newton.
    ! Newton is only applied
    ! if   ||absolute residuum|| < dtolAbsNewton
    ! and  ||relative residuum|| < dtolRelNewton.
    ! Otherwise, the usual fix point iteration is used.
    ! Value = 0.0: Disable check.
    real(DP) :: dtolAbsPartialNewton = 0.0_DP

    ! Norm of relative residuum before applying Newton.
    ! Newton is only applied
    ! if   ||absolute residuum|| < dtolAbsNewton
    ! and  ||relative residuum|| < dtolRelNewton.
    ! Otherwise, the usual fix point iteration is used.
    ! Standard value = 1E-1.
    ! Value = 0.0: Disable check.
    real(DP) :: dtolRelPartialNewton = 1.0E-1_DP
  
    ! Whether to use the inexact Newton iteration or not.
    ! The inexact Newton controls the stopping criterion of the linear
    ! solver according to the nonlinear residual.
    integer :: cinexactNewton = 1

    ! Stopping criterion for the linear solver in the inexact Newton iteration.
    ! Controls the minimum number of digits to gain in the linear solver
    ! per Newton iteration. 
    real(dp) :: dinexactNewtonTolRel = 1.0E-2_DP

    ! Exponent to control the stopping criterion for the linear solver in
    ! an inexact Newton iteration. =2 result in quadratic convergence,
    ! =1.5 in superlinear convergence. 
    real(dp) :: dinexactNewtonExponent = 2.0_DP

    ! Lower bound for the absolute residual. Subproblems are not solved
    ! more exact than this.
    real(DP) :: dinexactNewtonTolAbs = 1E-15_DP
    
    ! Number of smoothing steps after a Newton step.
    integer :: nsmoothingSteps = 0

    ! Step length strategy.
    ! =0: No step length control.
    ! =1: Damp with domega for the first nsteplengthSteps steps, then switch to full Newton.
    integer :: cstepLengthStrategy = 0

    ! Damping parameter
    real(DP) :: dstepLengthOmega = 1.0_DP

    ! Number of steps
    integer :: nstepLengthSteps = 0

  end type
  
!</typeblock>

  public :: t_ccDynamicNewtonControl

!<typeblock>
  
  ! Contains the parameters of the Newton iteration.
  type t_newtonParameters
  
    ! <!-- -------------- -->
    ! <!-- STATUS, OUTPUT -->
    ! <!-- -------------- -->

    ! OUTPUT: Result
    ! The result of the solution process.
    ! =-1: convergence criterion not reached.
    ! =0: success.
    ! =1: iteration broke down, diverging.
    ! =2: iteration broke down, preconditioner did not work.
    ! =3: error in the parameters.
    integer :: iresult = 0

    ! Output level of the solver.
    ! =-1: no output
    ! =0: no output, only errors
    ! =1: basic output
    ! =2: standard output
    integer :: ioutputLevel = 2

    ! Output mode. Used for printing messages.
    ! =OU_MODE_STD: Print messages to the terminal and probably to a log
    ! file (if a log file is opened).
    integer(I32) :: coutputmode = OU_MODE_STD
  
    ! <!-- ----------------------------------------- -->
    ! <!-- STOPPING CRITERIA, DAMPING PARAMETERS,... --> 
    ! <!-- ----------------------------------------- -->

    ! Type of the iteration.
    ! =-1: undefined
    ! =0: simple linear solver
    ! =1: nonlinear loop for a semilinear equation
    ! =2: Newton solver
    ! =3: adaptive Newton with parameters in radaptiveNewton
    integer :: ctypeIteration = -1
    
    ! Type of stopping criterion to use for standard convergence test. One of the
    ! NLSOL_STOP_xxxx constants.
    ! Note: This parameter is only evaluated in the stanard convergence test.
    ! If the caller of the nonlinear solver specifies a callback routine fcb_resNormCheck
    ! for checking the convergence, that callback routine must implement its own
    ! logic to handle relative and absolute convrgence criteria!
    ! integer :: cstoppingCriterion = NLSOL_STOP_STANDARD

    ! Type of norm to use in the residual checking of the total vector
    ! (cf. linearalgebra.f90).
    ! =0: euclidian norm, =1: l1-norm, =2: l2-norm, =3: MAX-norm
    ! integer :: iresNormTotal = 2

    ! RHS-vector is treated as zero if max(defect) < drhsZero
    real(DP) :: drhsZero = 1E-90_DP

    ! How to check residuals.
    ! =0: Check in Euclidian norm (like in old CC).
    ! =2: L2-norm
    integer :: iresNorm = LINALG_NORML2

    ! Damping parameter for the Newton iteration
    real(DP) :: domega = 1.0_DP
  
    ! <!-- ----------------------------- -->
    ! <!-- SUBSOLVERS AND OTHER SETTINGS -->
    ! <!-- ----------------------------- -->

    ! Parameters for the dynamic and inexact Newton algorithm
    type(t_ccDynamicNewtonControl) :: radaptiveNewton

  end type
  
!</typeblock>
  
  public :: t_newtonParameters

!<typeblock>
  
  ! Standard parameters for preconditioners in a Newton iteration
  type t_newtonPrecParameters
  
    ! <!-- -------------- -->
    ! <!-- STATUS, OUTPUT -->
    ! <!-- -------------- -->

    ! OUTPUT: Result
    ! The result of the solution process.
    ! =-1: convergence criterion not reached.
    ! =0: success.
    ! =1: iteration broke down, diverging.
    ! =2: iteration broke down, preconditioner did not work.
    ! =3: error in the parameters.
    integer :: iresult = 0

    ! Output level of the solver.
    ! =-1: no output
    ! =0: no output, only errors
    ! =1: basic output
    ! =2: standard output
    integer :: ioutputLevel = 2

    ! Output mode. Used for printing messages.
    ! =OU_MODE_STD: Print messages to the terminal and probably to a log
    ! file (if a log file is opened).
    integer(I32) :: coutputmode = OU_MODE_STD
  
    ! How to check residuals.
    ! =0: Check in Euclidian norm (like in old CC).
    ! =2: L2-norm
    integer :: iresNorm = LINALG_NORML2

    ! RHS-vector is treated as zero if max(defect) < drhsZero
    real(DP) :: drhsZero = 1E-90_DP

    ! Damping parameter
    real(DP) :: domega = 1.0_DP

  end type
  
!</typeblock>
  
  public :: t_newtonPrecParameters

!</types>

  ! Basic initialisation of the Newton solver
  public :: newtonit_initBasicParams
  
  ! Basic initialisation of the Newton solver
  public :: newtonlin_initBasicParams
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_initBasicParams (rsolver,riter,ssection,rparamList)
  
!<description>
  ! Initialises basic solver parameters according to a parameter list.
!</description>
  
!<input>
  ! Parameter list with the parameters configuring the nonlinear solver
  type(t_parlist), intent(in) :: rparamList

  ! Name of the section in the parameter list containing the parameters
  ! of the nonlinear solver.
  character(LEN=*), intent(in) :: ssection
!</input>

!<output>
  ! Solver structure receiving the parameters
  type(t_newtonParameters), intent(inout) :: rsolver
  
  ! Structure receiving iteration control parameters
  type(t_iterationControl), intent(inout) :: riter
!</output>

!</subroutine>

    ! local variables
    type(t_parlstSection), pointer :: p_rsection
    character(len=SYS_STRLEN) :: sstring,snewton

    call parlst_querysection(rparamList, ssection, p_rsection)

    if (.not. associated(p_rsection)) then
      call output_line ("Cannot create nonlinear solver; no section """//&
          trim(ssection)//"""!", &
          OU_CLASS_ERROR,OU_MODE_STD,"newtonit_initBasicParams")
      call sys_halt()
    end if

    ! Read iteration control parameters
    call itc_getParamsFromParlist(riter,ssection,rparamList)
    
    ! Read additional parameters
    call parlst_getvalue_double (p_rsection, "domega", &
                                 rsolver%domega, rsolver%domega)

    call parlst_getvalue_int (p_rsection, "ioutputLevel", &
        rsolver%ioutputLevel, rsolver%ioutputLevel)

    call parlst_getvalue_int (p_rsection, "iresnorm", &
        rsolver%iresnorm, rsolver%iresnorm)

    ! Get information about the iteration.
      
    ! At first, ask the parameters in the INI/DAT file which type of
    ! preconditioner is to be used. The data in the preconditioner structure
    ! is to be initialised appropriately!
    call parlst_getvalue_int (p_rsection, "ctypeIteration", &
        rsolver%ctypeIteration, rsolver%ctypeIteration)

    ! We have even the extended, dynamic Newton as preconditioner.
    ! Put the parameters for the extended Newton from the DAT file
    ! into the Adaptive-Newton configuration block.
    
    call parlst_getvalue_string (rparamList, ssection, &
        "ssectionAdaptiveNewton", sstring, "", bdequote=.true.)
    snewton = ""
    if (sstring .ne. "") read (sstring,*) snewton
    if (snewton .ne. "") then
      ! Initialise the parameters of the adaptive Newton
      call parlst_getvalue_int (rparamList, snewton, &
          "cpartialNewton", &
          rsolver%radaptiveNewton%cpartialNewton, &
          rsolver%radaptiveNewton%cpartialNewton)

      call parlst_getvalue_int (rparamList, snewton, &
          "nminPartialNewtonIterations", &
          rsolver%radaptiveNewton%nminPartialNewtonIterations, &
          rsolver%radaptiveNewton%nminPartialNewtonIterations)

      call parlst_getvalue_int (rparamList, snewton, &
          "nmaxPartialNewtonIterations", &
          rsolver%radaptiveNewton%nmaxPartialNewtonIterations, &
          rsolver%radaptiveNewton%nmaxPartialNewtonIterations)

      call parlst_getvalue_double (rparamList, snewton, &
          "dtolAbsPartialNewton", &
          rsolver%radaptiveNewton%dtolAbsPartialNewton, &
          rsolver%radaptiveNewton%dtolAbsPartialNewton)

      call parlst_getvalue_double (rparamList, snewton, &
          "dtolRelPartialNewton", &
          rsolver%radaptiveNewton%dtolRelPartialNewton, &
          rsolver%radaptiveNewton%dtolRelPartialNewton)

      call parlst_getvalue_int (rparamList, snewton, &
          "cinexactNewton", &
          rsolver%radaptiveNewton%cinexactNewton, &
          rsolver%radaptiveNewton%cinexactNewton)

      call parlst_getvalue_double (rparamList, snewton, &
          "dinexactNewtonTolRel", &
          rsolver%radaptiveNewton%dinexactNewtonTolRel, &
          rsolver%radaptiveNewton%dinexactNewtonTolRel)

      call parlst_getvalue_double (rparamList, snewton, &
          "dinexactNewtonExponent", &
          rsolver%radaptiveNewton%dinexactNewtonExponent, &
          rsolver%radaptiveNewton%dinexactNewtonExponent)

      call parlst_getvalue_double (rparamList, snewton, &
          "dinexactNewtonTolAbs", &
          rsolver%radaptiveNewton%dinexactNewtonTolAbs, &
          rsolver%radaptiveNewton%dinexactNewtonTolAbs)

      call parlst_getvalue_int (rparamList, snewton, &
          "nsmoothingSteps", &
          rsolver%radaptiveNewton%nsmoothingSteps, &
          rsolver%radaptiveNewton%nsmoothingSteps)

      call parlst_getvalue_int (rparamList, snewton, &
          "cstepLengthStrategy", &
          rsolver%radaptiveNewton%cstepLengthStrategy, &
          rsolver%radaptiveNewton%cstepLengthStrategy)

      call parlst_getvalue_int (rparamList, snewton, &
          "nstepLengthSteps", &
          rsolver%radaptiveNewton%nstepLengthSteps, &
          rsolver%radaptiveNewton%nstepLengthSteps)

      call parlst_getvalue_double (rparamList, snewton, &
          "dstepLengthOmega", &
          rsolver%radaptiveNewton%dstepLengthOmega, &
          rsolver%radaptiveNewton%dstepLengthOmega)

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_initBasicParams (rsolver,riter,ssection,rparamList)
  
!<description>
  ! Initialises basic solver parameters according to a parameter list.
!</description>
  
!<input>
  ! Parameter list with the parameters configuring the nonlinear solver
  type(t_parlist), intent(in) :: rparamList

  ! Name of the section in the parameter list containing the parameters
  ! of the nonlinear solver.
  character(LEN=*), intent(in) :: ssection
!</input>

!<output>
  ! Solver structure receiving the parameters
  type(t_newtonPrecParameters), intent(inout) :: rsolver

  ! Structure receiving iteration control parameters
  type(t_iterationControl), intent(inout) :: riter
!</output>

!</subroutine>

    ! local variables
    type(t_parlstSection), pointer :: p_rsection

    call parlst_querysection(rparamList, ssection, p_rsection)

    if (.not. associated(p_rsection)) then
      call output_line ("Cannot create linear solver; no section """//&
          trim(ssection)//"""!", &
          OU_CLASS_ERROR,OU_MODE_STD,"newtonlin_initBasicParams")
      call sys_halt()
    end if

    ! Read iteration control parameters
    call itc_getParamsFromParlist(riter,ssection,rparamList)

    call parlst_getvalue_int (p_rsection, "ioutputLevel", &
        rsolver%ioutputLevel, rsolver%ioutputLevel)

    call parlst_getvalue_int (p_rsection, "iresnorm", &
        rsolver%iresnorm, rsolver%iresnorm)

    call parlst_getvalue_double (p_rsection, "domega", &
        rsolver%domega, rsolver%domega)

  end subroutine

end module
