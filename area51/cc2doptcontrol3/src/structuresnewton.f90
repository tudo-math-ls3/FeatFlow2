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
  
  implicit none
  
  private
  
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
  
    ! Minimum number of fix point iterations before switching to
    ! preconditioning with the Newton matrix.

    integer :: nminFixedPointIterations = 0

    ! Maximum number of fix point iterations before switching to
    ! preconditioning with the Newton matrix.

    integer :: nmaxFixedPointIterations = 0

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
    ! per Newton iteration. 
    
    real(dp) :: dinexactNewtonEpsRel = 1.0E-2_DP

    ! Exponent to control the stopping criterion for the linear solver in
    ! an inexact Newton iteration. =2 result in quadratic convergence,
    ! =1.5 in superlinear convergence. 
    
    real(dp) :: dinexactNewtonExponent = 2.0_DP

    ! Lower bound for the absolute residual. Subproblems are not solved
    ! more exact than this.
    real(DP) :: dinexactNewtonEpsAbs = 1E-15_DP

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

    ! Relative stopping criterion. Stop iteration if
    ! !!defect!! < EPSREL * !!initial defect!!.
    ! =0: ignore, use absolute stopping criterion; standard = 1E-5
    ! Remark: do not set depsAbs=depsRel=0!
    real(DP) :: depsRel = 1E-5_DP

    ! Absolute stopping criterion. Stop iteration if
    ! !!defect!! < EPSABS.
    ! =0: ignore, use relative stopping criterion; standard = 1E-5
    ! Remark: do not set depsAbs=depsRel=0!
    real(DP) :: depsAbs = 1E-5_DP

    ! Relative divergence criterion.  Treat iteration as
    ! diverged if
    !   !!defect!! >= DIVREL * !!initial defect!!
    ! A value of SYS_INFINITY_DP disables the relative divergence check.
    ! standard = 1E3
    real(DP) :: ddivRel = 1E6_DP

    ! Absolute divergence criterion. Treat iteration as
    ! diverged if
    !   !!defect!! >= DIVREL
    ! A value of SYS_INFINITY_DP disables the absolute divergence check.
    ! standard = SYS_INFINITY_DP
    real(DP) :: ddivAbs = SYS_INFINITY_DP

    ! Minimum number of iterations top perform
    integer :: nminIterations = 1

    ! Maximum number of iterations top perform
    integer :: nmaxIterations = 50
  
    ! Damping parameter for the Newton iteration
    real(DP) :: domega = 1.0_DP
  
    ! How to check residuals.
    ! =0: Check in Euclidian norm (like in old CC).
    ! =2: L2-norm
    integer :: iresNorm = LINALG_NORML2

    ! <!-- ----------------------------- -->
    ! <!-- SUBSOLVERS AND OTHER SETTINGS -->
    ! <!-- ----------------------------- -->

    ! Parameters for the dynamic and inexact Newton algorithm
    type(t_ccDynamicNewtonControl) :: radaptiveNewton

    ! <!-- ---------- -->
    ! <!-- STATISTICS -->
    ! <!-- ---------- -->
    
    ! Initial residual in the control space.
    real(DP) :: dresInit = 0.0_DP

    ! Final residual in the control space.
    real(DP) :: dresFinal = 0.0_DP

    ! Total number of nonlinear iterations
    integer :: niterations = 0
    
    ! Total number of linear iterations
    integer :: nlinearIterations = 0
    
!    ! Total time for the solver
!    real(DP) :: dtimeTotal = 0.0_DP
!    
!    ! Total time for nonlinear iteration
!    real(DP) :: dtimeNonlinearSolver = 0.0_DP
!    
!    ! Total time for linear space-time solver
!    real(DP) :: dtimeLinearSolver = 0.0_DP
!
!    ! Total time for linear solver in space
!    real(DP) :: dtimeLinearSolverTimestep = 0.0_DP
!    
!    ! Total time for matrix assembly
!    real(DP) :: dtimeMatrixAssembly = 0.0_DP
!
!    ! Total time for defect calculation
!    real(DP) :: dtimeDefectCalculation = 0.0_DP
!
!    ! Total time for matrix assembly
!    real(DP) :: dtimePostprocessing = 0.0_DP

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
  
    ! <!-- ----------------------------------------- -->
    ! <!-- STOPPING CRITERIA, DAMPING PARAMETERS,... --> 
    ! <!-- ----------------------------------------- -->

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

    ! Relative stopping criterion. Stop iteration if
    ! !!defect!! < EPSREL * !!initial defect!!.
    ! =0: ignore, use absolute stopping criterion; standard = 1E-5
    ! Remark: do not set depsAbs=depsRel=0!
    real(DP) :: depsRel = 1E-5_DP

    ! Absolute stopping criterion. Stop iteration if
    ! !!defect!! < EPSABS.
    ! =0: ignore, use relative stopping criterion; standard = 1E-5
    ! Remark: do not set depsAbs=depsRel=0!
    real(DP) :: depsAbs = 1E-5_DP

    ! Relative divergence criterion.  Treat iteration as
    ! diverged if
    !   !!defect!! >= DIVREL * !!initial defect!!
    ! A value of SYS_INFINITY_DP disables the relative divergence check.
    ! standard = 1E3
    real(DP) :: ddivRel = 1E6_DP

    ! Absolute divergence criterion. Treat iteration as
    ! diverged if
    !   !!defect!! >= DIVREL
    ! A value of SYS_INFINITY_DP disables the absolute divergence check.
    ! standard = SYS_INFINITY_DP
    real(DP) :: ddivAbs = SYS_INFINITY_DP

    ! Minimum number of iterations top perform
    integer :: nminIterations = 1

    ! Maximum number of iterations top perform
    integer :: nmaxIterations = 50
  
    ! Damping parameter for the Newton iteration
    real(DP) :: domega = 1.0_DP
  
    ! How to check residuals.
    ! =0: Check in Euclidian norm (like in old CC).
    ! =2: L2-norm
    integer :: iresNorm = LINALG_NORML2

    ! RHS-vector is treated as zero if max(defect) < drhsZero
    real(DP) :: drhsZero = 1E-90_DP

    ! <!-- ---------- -->
    ! <!-- STATISTICS -->
    ! <!-- ---------- -->
    
    ! Initial residual in the control space.
    real(DP) :: dresInit = 0.0_DP

    ! Final residual in the control space.
    real(DP) :: dresFinal = 0.0_DP

    ! Convergence rate
    real(DP) :: dconvergenceRate = 0.0_DP

    ! Total number of iterations
    integer :: niterations = 0
    
!    ! Total time for the solver
!    real(DP) :: dtimeTotal = 0.0_DP
!    
!    ! Total time for nonlinear iteration
!    real(DP) :: dtimeNonlinearSolver = 0.0_DP
!    
!    ! Total time for linear space-time solver
!    real(DP) :: dtimeLinearSolver = 0.0_DP
!
!    ! Total time for linear solver in space
!    real(DP) :: dtimeLinearSolverTimestep = 0.0_DP
!    
!    ! Total time for matrix assembly
!    real(DP) :: dtimeMatrixAssembly = 0.0_DP
!
!    ! Total time for defect calculation
!    real(DP) :: dtimeDefectCalculation = 0.0_DP
!
!    ! Total time for matrix assembly
!    real(DP) :: dtimePostprocessing = 0.0_DP

  end type
  
!</typeblock>
  
  public :: t_newtonPrecParameters

!</types>

  ! Basic initialisation of the Newton solver
  public :: newtonit_initBasicParams
  
  ! Standard convergence check
  public :: newtonit_checkConvergence
  
  ! Check for stopping the iteration
  public :: newtonit_checkIterationStop
  
  ! Standard divergence check
  public :: newtonit_checkDivergence
  
  ! Basic initialisation of the Newton solver
  public :: newtonlin_initBasicParams
  
  ! Standard convergence check
  public :: newtonlin_checkConvergence
  
  ! Check for stopping the iteration
  public :: newtonlin_checkIterationStop
  
  ! Standard divergence check
  public :: newtonlin_checkDivergence
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_initBasicParams (rsolver,ssection,rparamList)
  
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
!</output>

!</subroutine>

    ! local variables
    type(t_parlstSection), pointer :: p_rsection
    character(len=SYS_STRLEN) :: sstring,snewton

    call parlst_querysection(rparamList, ssection, p_rsection)

    if (.not. associated(p_rsection)) then
      call output_line ("Cannot create nonlinear solver; no section ''"//&
          trim(ssection)//"''!", &
          OU_CLASS_ERROR,OU_MODE_STD,"newtonit_initBasicParams")
      call sys_halt()
    end if

    ! Get stopping criteria of the nonlinear iteration
    call parlst_getvalue_double (p_rsection, "depsRel", &
                                 rsolver%depsRel, rsolver%depsRel)

    call parlst_getvalue_double (p_rsection, "depsAbs", &
                                 rsolver%depsAbs, rsolver%depsAbs)

    call parlst_getvalue_double (p_rsection, "ddivRel", &
                                 rsolver%ddivRel, rsolver%ddivRel)

    call parlst_getvalue_double (p_rsection, "ddivAbs", &
                                 rsolver%ddivAbs, rsolver%ddivAbs)

    call parlst_getvalue_double (p_rsection, "domega", &
                                 rsolver%domega, rsolver%domega)

    call parlst_getvalue_int (p_rsection, "nminIterations", &
        rsolver%nminIterations, rsolver%nminIterations)

    call parlst_getvalue_int (p_rsection, "nmaxIterations", &
        rsolver%nmaxIterations, rsolver%nmaxIterations)

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
          "nminFixedPointIterations", &
          rsolver%radaptiveNewton%nminFixedPointIterations, &
          rsolver%radaptiveNewton%nminFixedPointIterations)

      call parlst_getvalue_int (rparamList, snewton, &
          "nmaxFixedPointIterations", &
          rsolver%radaptiveNewton%nmaxFixedPointIterations, &
          rsolver%radaptiveNewton%nmaxFixedPointIterations)

      call parlst_getvalue_double (rparamList, snewton, &
          "depsAbsNewton", &
          rsolver%radaptiveNewton%depsAbsNewton, &
          rsolver%radaptiveNewton%depsAbsNewton)

      call parlst_getvalue_double (rparamList, snewton, &
          "depsRelNewton", &
          rsolver%radaptiveNewton%depsRelNewton, &
          rsolver%radaptiveNewton%depsRelNewton)

      call parlst_getvalue_int (rparamList, snewton, &
          "cinexactNewton", &
          rsolver%radaptiveNewton%cinexactNewton, &
          rsolver%radaptiveNewton%cinexactNewton)

      call parlst_getvalue_double (rparamList, snewton, &
          "dinexactNewtonEpsRel", &
          rsolver%radaptiveNewton%dinexactNewtonEpsRel, &
          rsolver%radaptiveNewton%dinexactNewtonEpsRel)

      call parlst_getvalue_double (rparamList, snewton, &
          "dinexactNewtonExponent", &
          rsolver%radaptiveNewton%dinexactNewtonExponent, &
          rsolver%radaptiveNewton%dinexactNewtonExponent)

      call parlst_getvalue_double (rparamList, snewton, &
          "dinexactNewtonEpsAbs", &
          rsolver%radaptiveNewton%dinexactNewtonEpsAbs, &
          rsolver%radaptiveNewton%dinexactNewtonEpsAbs)

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  logical function newtonit_checkConvergence (rsolver)
  
!<description>
  ! Standard convergence check.
!</description>
  
!<result>
  ! Returns TRUE if the iteration is converged.
!</result>
  
!<inputoutput>
  ! Solver structure receiving the parameters.
  ! If convergence is reached, iresult is set.
  type(t_newtonParameters), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>

    newtonit_checkConvergence = .false.
    
    ! Check the residual.
    if (rsolver%niterations .ge. rsolver%nminIterations) then
      ! Absolute residual
      if (rsolver%depsAbs .gt. 0.0_DP) then
        if (rsolver%dresFinal .le. rsolver%depsAbs) then
          rsolver%iresult = 0
          newtonit_checkConvergence = .true.
          return
        end if
      end if

      ! Relative residual
      if (rsolver%depsRel .gt. 0.0_DP) then
        if (rsolver%dresFinal .le. rsolver%depsRel * rsolver%dresInit) then
          rsolver%iresult = 0
          newtonit_checkConvergence = .true.
          return
        end if
      end if

    end if

  end function

  ! ***************************************************************************

!<subroutine>

  logical function newtonit_checkIterationStop (rsolver,iite)
  
!<description>
  ! Checks whether or not to stop the iteration.
!</description>
  
!<result>
  ! Returns TRUE if the iteration should be stopped.
!</result>
  
!<input>
  ! Current iteration counter
  integer, intent(in) :: iite
!</input>

!<inputoutput>
  ! Solver structure receiving the parameters.
  ! If convergence is reached, iresult is set.
  type(t_newtonParameters), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>

    newtonit_checkIterationStop = .false.
    
    if (rsolver%niterations .ge. rsolver%nmaxIterations) then
      ! Maximum number of iterations reached.
      rsolver%iresult = -1
      newtonit_checkIterationStop = .true.
      return
    end if

  end function

  ! ***************************************************************************

!<subroutine>

  logical function newtonit_checkDivergence (rsolver)
  
!<description>
  ! Standard convergence check.
!</description>
  
!<result>
  ! Returns TRUE if the iteration is converged.
!</result>
  
!<inputoutput>
  ! Solver structure receiving the parameters.
  ! If divergence is reached, iresult is set.
  type(t_newtonParameters), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>

    newtonit_checkDivergence = .false.

    ! Absolute residual
    if (.not. (rsolver%dresFinal .lt. rsolver%ddivAbs)) then
      rsolver%iresult = 1
      newtonit_checkDivergence = .true.
      return
    end if
    
    ! Relative residual
    if (.not. (rsolver%dresFinal .le. rsolver%ddivRel * rsolver%dresInit)) then
      rsolver%iresult = 1
      newtonit_checkDivergence = .true.
      return
    end if

  end function

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_initBasicParams (rsolver,ssection,rparamList)
  
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
!</output>

!</subroutine>

    ! local variables
    type(t_parlstSection), pointer :: p_rsection

    call parlst_querysection(rparamList, ssection, p_rsection)

    if (.not. associated(p_rsection)) then
      call output_line ("Cannot create linear solver; no section ''"//&
          trim(ssection)//"''!", &
          OU_CLASS_ERROR,OU_MODE_STD,"newtonlin_initBasicParams")
      call sys_halt()
    end if

    ! Get stopping criteria of the nonlinear iteration
    call parlst_getvalue_double (p_rsection, "depsRel", &
                                 rsolver%depsRel, rsolver%depsRel)

    call parlst_getvalue_double (p_rsection, "depsAbs", &
                                 rsolver%depsAbs, rsolver%depsAbs)

    call parlst_getvalue_double (p_rsection, "ddivRel", &
                                 rsolver%ddivRel, rsolver%ddivRel)

    call parlst_getvalue_double (p_rsection, "ddivAbs", &
                                 rsolver%ddivAbs, rsolver%ddivAbs)

    call parlst_getvalue_double (p_rsection, "domega", &
                                 rsolver%domega, rsolver%domega)

    call parlst_getvalue_int (p_rsection, "nminIterations", &
        rsolver%nminIterations, rsolver%nminIterations)

    call parlst_getvalue_int (p_rsection, "nmaxIterations", &
        rsolver%nmaxIterations, rsolver%nmaxIterations)

    call parlst_getvalue_int (p_rsection, "ioutputLevel", &
        rsolver%ioutputLevel, rsolver%ioutputLevel)

    call parlst_getvalue_int (p_rsection, "iresnorm", &
        rsolver%iresnorm, rsolver%iresnorm)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  logical function newtonlin_checkConvergence (rsolver)
  
!<description>
  ! Standard convergence check.
!</description>
  
!<result>
  ! Returns TRUE if the iteration is converged.
!</result>
  
!<inputoutput>
  ! Solver structure receiving the parameters.
  ! If convergence is reached, iresult is set.
  type(t_newtonPrecParameters), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>

    newtonlin_checkConvergence = .false.
    
    ! Check the residual.
    if (rsolver%niterations .ge. rsolver%nminIterations) then
      ! Absolute residual
      if (rsolver%depsAbs .gt. 0.0_DP) then
        if (rsolver%dresFinal .le. rsolver%depsAbs) then
          rsolver%iresult = 0
          newtonlin_checkConvergence = .true.
          return
        end if
      end if

      ! Relative residual
      if (rsolver%depsRel .gt. 0.0_DP) then
        if (rsolver%dresFinal .le. rsolver%depsRel * rsolver%dresInit) then
          rsolver%iresult = 0
          newtonlin_checkConvergence = .true.
          return
        end if
      end if
      
    end if

  end function

  ! ***************************************************************************

!<subroutine>

  logical function newtonlin_checkIterationStop (rsolver)
  
!<description>
  ! Checks whether or not to stop the iteration.
!</description>
  
!<result>
  ! Returns TRUE if the iteration should be stopped.
!</result>
  
!<inputoutput>
  ! Solver structure receiving the parameters.
  ! If convergence is reached, iresult is set.
  type(t_newtonPrecParameters), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>

    newtonlin_checkIterationStop = .false.
    
    if (rsolver%niterations .ge. rsolver%nmaxIterations) then
      ! Maximum number of iterations reached.
      rsolver%iresult = -1
      newtonlin_checkIterationStop = .true.
      return
    end if

  end function

  ! ***************************************************************************

!<subroutine>

  logical function newtonlin_checkDivergence (rsolver)
  
!<description>
  ! Standard convergence check.
!</description>
  
!<result>
  ! Returns TRUE if the iteration is converged.
!</result>
  
!<inputoutput>
  ! Solver structure receiving the parameters.
  ! If divergence is reached, iresult is set.
  type(t_newtonPrecParameters), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>

    newtonlin_checkDivergence = .false.

    ! Absolute residual
    if (rsolver%dresFinal .ge. rsolver%ddivAbs) then
      rsolver%iresult = 1
      newtonlin_checkDivergence = .true.
      return
    end if
    
    ! Relative residual
    if (rsolver%dresFinal .ge. rsolver%ddivRel * rsolver%dresInit) then
      rsolver%iresult = 1
      newtonlin_checkDivergence = .true.
      return
    end if

  end function

end module
