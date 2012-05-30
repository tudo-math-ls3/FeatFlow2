!##############################################################################
!# ****************************************************************************
!# <name> newtoniteration </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module realises a (semismooth) Newton iteration in the control 
!# space of the optimal control problem.
!# </purpose>
!##############################################################################

module newtoniteration

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
  use assemblytemplates
  
  use spacematvecassembly
  
  use kktsystemspaces
  use kktsystem
  
  use newtoniterationlinear
  
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
    ! =1: nonlinear defect correction solver.
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
  
    ! User defined parameter list; if not NULL(), this is added
    ! as comments to postprocessing files.
    type(t_parlist), pointer :: p_rparlist => null()

    ! <!-- ----------------------------- -->
    ! <!-- SUBSOLVERS AND OTHER SETTINGS -->
    ! <!-- ----------------------------- -->

    ! Parameters for the dynamic and inexact Newton algorithm
    type(t_ccDynamicNewtonControl) :: radaptiveNewton

    ! Parameters for the linear space-time subsolver.
    type(t_linsolParameters) :: rlinsolParam
    
    ! <!-- ---------- -->
    ! <!-- STATISTICS -->
    ! <!-- ---------- -->
    
    ! Initial residual in the control space.
    real(DP) :: dresInit = 0.0_DP

    ! Final residual in the control space.
    real(DP) :: dresFinal = 0.0_DP

    ! Total number of nonlinear iterations
    integer :: nnonlinearIterations = 0
    
    ! Total number of linear iterations
    integer :: nlinearIterations = 0
    
    ! Total time for the solver
    real(DP) :: dtimeTotal = 0.0_DP
    
    ! Total time for nonlinear iteration
    real(DP) :: dtimeNonlinearSolver = 0.0_DP
    
    ! Total time for linear space-time solver
    real(DP) :: dtimeLinearSolver = 0.0_DP

    ! Total time for linear solver in space
    real(DP) :: dtimeLinearSolverTimestep = 0.0_DP
    
    ! Total time for matrix assembly
    real(DP) :: dtimeMatrixAssembly = 0.0_DP

    ! Total time for defect calculation
    real(DP) :: dtimeDefectCalculation = 0.0_DP

    ! Total time for matrix assembly
    real(DP) :: dtimePostprocessing = 0.0_DP

  end type
  
!</typeblock>
  
  public :: t_newtonParameters

!</types>

  ! Apply a Newton iteration in the control space
  public :: newtonit_solve
  
  ! Basic initialisation of the Newton solver
  public :: newtonit_init
  
  ! Structural initialisation
  public :: newtonit_initStructure
  
  ! Final initialisation
  public :: newtonit_initData
  
  ! Cleanup of data
  public :: newtonit_doneData

  ! Cleanup of structures
  public :: newtonit_doneStructure

  ! Final cleanup
  public :: newtonit_done
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_initParams (rsolver,ssection,rparamList)
  
!<description>
  ! Initialises the solver parameters according to a parameter list.
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
          OU_CLASS_ERROR,OU_MODE_STD,"newtonit_initParams")
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

    call parlst_getvalue_int (p_rsection, 'ioutputLevel', &
        rsolver%ioutputLevel, rsolver%ioutputLevel)

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

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_getResidual (rkktsystem,rresidual,dres)
  
!<description>
  ! Calculates the basic (unprecondiotioned) search direction of the 
  ! Newton algorithm.
!</description>
  
!<inputoutput>
  ! Structure defining the KKT system. The control in this structure
  ! defines the current 'state' of the Newton algorithm.
  ! On output, the 
  type(t_kktsystem), intent(inout) :: rkktsystem
  
  ! On output, this structure receives a representation of the search
  ! direction / residual in the Newton iteration.
  type(t_controlSpace), intent(inout) :: rresidual
!</inputoutput>

!<output>
  ! L2-Norm of the residual
  real(DP), intent(out) :: dres
!</output>

!</subroutine>

    ! -------------------------------------------------------------
    ! Step 1: Solve the primal and dual system.
    ! -------------------------------------------------------------

    ! Solve the primal equation, update the primal solution.
    call kkt_solvePrimal (rkktsystem)
    
    ! Solve the dual equation, update the dual solution.
    call kkt_solveDual (rkktsystem)

    ! -------------------------------------------------------------
    ! Step 2: Calculate the search direction   
    ! -------------------------------------------------------------

    ! The search direction is just the residual in the control equation.
    call kkt_calcControlRes (rkktsystem,rresidual,dres)

  end subroutine

    ! ***************************************************************************

!<subroutine>

  subroutine newtonit_updateControl (rnewtonParam,rkktsystem,rcorrection)
  
!<description>
  ! Calculates the basic (unprecondiotioned) search direction of the 
  ! Newton algorithm.
!</description>  

!<input>
  ! Parameters for the Newton iteration.
  ! The output parameters are changed according to the iteration.
  type(t_newtonParameters), intent(in) :: rnewtonParam

  ! The preconditioned search direction
  type(t_controlSpace), intent(inout) :: rcorrection
!</input>
  
!<inputoutput>
  ! Structure defining the KKT system. The control in this structure
  ! is updated according to the search direction.
  type(t_kktsystem), intent(inout) :: rkktsystem
!</inputoutput>

!</subroutine>

    ! Currectly, this is just a linear combination of the control variables.
    !
    !    u_n+1  =  u_n + g_n
    !
    ! Later, a step length control can be added here.
    
    call kktsp_controlLinearComb (&
        rcorrection,rnewtonParam%domega,rkktsystem%p_rcontrol,1.0_DP)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_solve (rnewtonParam,rkktsystem)
  
!<inputoutput>
  ! Parameters for the Newton iteration.
  ! The output parameters are changed according to the iteration.
  type(t_newtonParameters), intent(inout) :: rnewtonParam

  ! Structure defining the KKT system.
  ! The solutions in this structure are taken as initial
  ! values. On exit, the structure contains improved solutions.
  type(t_kktsystem), intent(inout) :: rkktsystem
!</inputoutput>

!</subroutine>
   
    ! local variables
    type(t_kktsystemDirDeriv) :: rkktsystemDirDeriv
    type(t_controlSpace) :: rdescentDir
    
    ! Prepare a structure that encapsules the directional derivative.
    
    ! Apply the Newton iteration
    rnewtonParam%nnonlinearIterations = 0
    
    do while (.true.)
    
      ! The Newton iteration reads
      !
      !    u_n+1  =  u_n  -  [J''(u_n)]^-1  J'(u_n)
      !
      ! or in other words,
      !
      !    u_n+1  =  u_n  +  [J''(u_n)]^-1  d_n
      !
      ! with the residual   d_n = -J'(u_n)   specifying a 'search direction'.
      
      ! -------------------------------------------------------------
      ! Get the current residual / search direction
      ! -------------------------------------------------------------

      ! Compute the basic (unpreconditioned) search direction d_n.
      call newtonit_getResidual (rkktsystem,rdescentDir,rnewtonParam%dresFinal)

      if (rnewtonParam%nnonlinearIterations .eq. 1) then
        ! Remember the initial residual
        rnewtonParam%dresInit = rnewtonParam%dresFinal
      end if

      ! -------------------------------------------------------------
      ! Check for convergence
      ! -------------------------------------------------------------
      if (rnewtonParam%nnonlinearIterations .ge. rnewtonParam%nminIterations) then
        ! Check the residual.
        !
        ! Absolute residual
        if (rnewtonParam%dresFinal .le. rnewtonParam%depsAbs) then
          rnewtonParam%iresult = 0
          exit
        end if
        
        ! Relative residual
        if (rnewtonParam%dresFinal .le. rnewtonParam%depsRel * rnewtonParam%dresInit) then
          rnewtonParam%iresult = 0
          exit
        end if
      end if
      
      ! -------------------------------------------------------------
      ! Check for divergence
      ! -------------------------------------------------------------
      ! Absolute residual
      if (rnewtonParam%dresFinal .ge. rnewtonParam%ddivAbs) then
        rnewtonParam%iresult = 1
        exit
      end if
      
      ! Relative residual
      if (rnewtonParam%dresFinal .ge. rnewtonParam%ddivRel * rnewtonParam%dresInit) then
        rnewtonParam%iresult = 1
        exit
      end if
      
      ! -------------------------------------------------------------
      ! Check other stopping criteria
      ! -------------------------------------------------------------

      if (rnewtonParam%nnonlinearIterations .ge. rnewtonParam%nmaxIterations) then
        ! Maximum number of iterations reached.
        rnewtonParam%iresult = -1
        exit
      end if
      
      ! -------------------------------------------------------------
      ! Preconditioning with the Newton matrix
      ! -------------------------------------------------------------
      
      ! Actual Newton iteration. Apply the Newton preconditioner
      ! to get the Newton search direction:
      !
      !    J''(u_n) g_n  =  d_n
      !
      call newtonlin_precondNewton (rnewtonParam%rlinsolParam,&
          rkktsystemDirDeriv,rdescentDir)
      
      ! -------------------------------------------------------------
      ! Update of the solution
      ! -------------------------------------------------------------

      ! Update the control according to the search direction:
      !
      !    u_n+1  =  u_n  +  g_n
      !
      ! or to any configured step-length control rule.
      call newtonit_updateControl (rnewtonParam,rkktsystem,rdescentDir)
      
      ! -------------------------------------------------------------
      ! Proceed with the next iteration
      ! -------------------------------------------------------------
      ! Next iteration
      rnewtonParam%nnonlinearIterations = rnewtonParam%nnonlinearIterations + 1
    
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_init (rlinsolParam,rparlist,ssection)
  
!<description>
  ! Basic initialisation of the Newton solver.
!</description>

!<input>
  ! Parameter list that contains the data for the preconditioning.
  type(t_parlist), intent(in) :: rparlist
  
  ! Entry Section in the parameter list containing the data of the 
  ! preconditioner.
  character(len=*), intent(in) :: ssection
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_newtonParameters), intent(out) :: rlinsolParam
!</inputoutput>

!</subroutine>
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_initStructure (rlinsolParam)
  
!<description>
  ! Structural initialisation of the Newton solver.
!</description>

!<inputoutput>
  ! Structure to be initialised.
  type(t_newtonParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>
   
  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_initData (rlinsolParam)
  
!<description>
  ! Final preparation of the Newton solver.
!</description>

!<inputoutput>
  ! Structure to be initialised.
  type(t_newtonParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_doneData (rlinsolParam)
  
!<description>
  ! Cleanup of the data initalised in newtonit_initData.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_newtonParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_doneStructure (rlinsolParam)
  
!<description>
  ! Cleanup of the data initalised in newtonit_initStructure.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_newtonParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_done (rlinsolParam)
  
!<description>
  ! Clean up the Newton iteration.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_newtonParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>
   
  end subroutine

end module
