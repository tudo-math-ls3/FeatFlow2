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
  use statistics
  use iterationcontrol
  
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
  use structuresoptflow
  use structuresoperatorasm
  use structuresnewton
  use assemblytemplates
  
  use spacematvecassembly
  use spacelinearsolver
  use spacesolver
  
  use spacetimehierarchy
  use spacetimeinterlevelprj
  use kktsystemspaces
  use kktsystem
  use kktsystemhierarchy
  use postprocessing
  
  use newtoniterationlinear

  implicit none
  
  private

!<types>

!<typeblock>

  ! Linear solver statistics
  type t_newtonitSolverStat
  
    ! Number of iterations necessary for the solver
    integer :: niterations = 0
    
    ! Total time necessary for the solver
    type(t_timer) :: rtotalTime

    ! Time necessary for solving the nonlinear forward problem
    type(t_timer) :: rtimeForward

    ! Time necessary for solving the backward problem
    type(t_timer) :: rtimeBackward
    
    ! Time necessary for the creation of nonlinear defects
    type(t_timer) :: rtimeDefect

    ! Time necessary for the prolongation of the solution to all levels
    type(t_timer) :: rtimeProlRest
    
    ! Time necessary for postprocessing
    type(t_timer) :: rtimePostprocessing

    ! Statistics of the subsolver in space used for calculating the residual.
    type(t_spaceslSolverStat) :: rspaceslSolverStat

    ! Statistics of the linear subsolver
    type(t_newtonlinSolverStat) :: rnewtonlinSolverStat
    
  end type

!</typeblock>

  public :: t_newtonitSolverStat

!<typeblock>
  
  ! Parameters for the space-time Newton algorithm.
  type t_spacetimeNewton

    ! <!-- --------------------------------------- -->
    ! <!-- GENRERAL PARAMETERS AND SOLVER SETTINGS -->
    ! <!-- --------------------------------------- -->

    ! General newton parameters
    type(t_newtonParameters) :: rnewtonParams
    
    ! Iteration control parameters
    type(t_iterationControl) :: riter

    ! Parameters of the OptFlow solver
    type(t_settings_optflow), pointer :: p_rsettingsSolver => null()

    ! <!-- ----------------------------- -->
    ! <!-- SUBSOLVERS AND OTHER SETTINGS -->
    ! <!-- ----------------------------- -->

    ! Defines a policy how to generate the initial condition of a timestep.
    ! =0: Always take zero
    ! =1: Propagate the solution of the previous/next timestep to the
    !     current one. (Default)
    ! =2: Take the solution of the last space-time iteration
    integer :: cspatialInitCondPolicy = SPINITCOND_PREVTIMESTEP

    ! Whether to postprocess intermediate solutions.
    ! =1: Calculate functional values, errors,...
    ! =2: Write postprocessing files with unique filename.
    ! =3: Calculate functional values, errors,... and
    !     write postprocessing files with unique filename.
    integer :: cpostprocessIterates = 1

    ! Parameters for the linear space-time subsolver.
    type(t_linsolParameters) :: rlinsolParam
    
    ! KKT subsolver hierarchy which encapsules all subsolvers used
    ! by the KKT system solvers.
    type(t_kktSubsolverSet) :: rkktSubsolvers
    
    ! <!-- -------------- -->
    ! <!-- TEMPORARY DATA -->
    ! <!-- -------------- -->
    
    ! Underlying KKT system hierarchy that defines the shape of the solutions.
    type(t_kktsystemHierarchy), pointer :: p_rkktsystemHierarchy => null()
    
    ! Hierarchy of solutions (nonlinearities) for the nonlinear iteration
    type(t_kktsystemHierarchy), pointer :: p_rsolutionHierarchy => null()
        
    ! Hierarchy of directional derivatives. Calculated during the Newton iteration.
    type(t_kktsystemDirDerivHierarchy), pointer :: p_rdirDerivHierarchy => null()
    
    ! Descent direction
    type(t_controlSpace), pointer :: p_rdescentDir => null()
    
  end type

!</typeblock>

   public :: t_spacetimeNewton

!</types>
  
  ! Basic initialisation of the Newton solver
  public :: newtonit_init
  
  ! Structural initialisation
  public :: newtonit_initStructure
  
  ! Apply a Newton iteration in the control space
  public :: newtonit_solve
  
  ! Cleanup of structures
  public :: newtonit_doneStructure

  ! Final cleanup
  public :: newtonit_done
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_init (rsolver,rsettingsSolver,ssection,rparamList)
  
!<description>
  ! Initialises the solver parameters according to a parameter list.
!</description>
  
!<input>
  ! Parameters of the OptFlow solver
  type(t_settings_optflow), intent(in), target :: rsettingsSolver
  
  ! Parameter list with the parameters configuring the nonlinear solver
  type(t_parlist), intent(in) :: rparamList

  ! Name of the section in the parameter list containing the parameters
  ! of the nonlinear solver.
  character(LEN=*), intent(in) :: ssection
!</input>

!<output>
  ! Solver structure receiving the parameters
  type(t_spacetimeNewton), intent(inout) :: rsolver
!</output>

!</subroutine>

    ! local variables
    character(LEN=SYS_STRLEN) :: ssolverLin

    ! Remember the solver settings for later use
    rsolver%p_rsettingsSolver => rsettingsSolver

    ! Initialise basic parameters
    call newtonit_initBasicParams (rsolver%rnewtonParams,rsolver%riter,ssection,rparamList)

    ! Get the sections with the parameters for the nonlinear / linear
    ! solver in space
    call parlst_getvalue_int (rparamList, ssection, "cspatialInitCondPolicy", &
        rsolver%cspatialInitCondPolicy, rsolver%cspatialInitCondPolicy)

    call parlst_getvalue_int (rparamList, ssection, "cpostprocessIterates", &
        rsolver%cpostprocessIterates, rsolver%cpostprocessIterates)
        
    call parlst_getvalue_string (rparamList, ssection, &
        "ssectionLinSolverSpaceTime", ssolverLin, "SPACETIME-LINSOLVER",bdequote=.true.)

    ! Initialise the KKT system subsolvers.
    !
    ! The definition of the lineraised forward/backward equation depends upon
    ! whether we use the full Newton approach or not.
    select case (rsolver%rnewtonParams%ctypeIteration)

    ! --------------
    ! Partial Newton
    ! --------------
    case (1)
      call kkt_initSubsolvers (rsolver%rkktSubsolvers,&
          rsettingsSolver,ssection,rparamList,0)
    
    ! ----------------------------
    ! Full Newton, adaptive Newton
    ! ----------------------------
    case (2,3)
      call kkt_initSubsolvers (rsolver%rkktSubsolvers,&
          rsettingsSolver,ssection,rparamList,1)
          
    case default
      call output_line ("Invalid nonlinear iteration",&
          OU_CLASS_ERROR,OU_MODE_STD,"newtonit_init")
      call sys_halt()
      
    end select

    ! Initialise the linear subsolver
    call newtonlin_init (rsolver%rlinsolParam,rsettingsSolver,&
        rsolver%rkktSubsolvers,rparamList,ssolverLin)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_getResidual (rsolver,rkktsystem,rresidual,dres,iresnorm,rstatistics)
  
!<description>
  ! Calculates the basic (unprecondiotioned) search direction of the 
  ! Newton algorithm.
!</description>
  
!<input>
  ! Type of norm. A LINALG_NORMxxxx constant.
  integer, intent(in) :: iresnorm
!</input>

!<inputoutput>
  ! Parameters for the Newton iteration.
  type(t_spacetimeNewton), intent(inout) :: rsolver

  ! Structure defining the KKT system. The control in this structure
  ! defines the current "state" of the Newton algorithm.
  ! On output, the primal and dual solution are updated.
  type(t_kktsystem), intent(inout) :: rkktsystem
  
  ! On output, this structure receives a representation of the search
  ! direction / residual in the Newton iteration.
  type(t_controlSpace), intent(inout) :: rresidual
!</inputoutput>


!<output>
  ! Statistic structure. 
  type(t_newtonitSolverStat), intent(out) :: rstatistics

  ! L2-Norm of the residual
  real(DP), intent(out) :: dres
!</output>

!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rvector
    real(DP), dimension(:), pointer :: p_Ddata
    type(t_spaceslSolverStat) :: rlocalStat
    
    call stat_startTimer (rstatistics%rtotalTime)
    
    call sptivec_getVectorFromPool(rkktsystem%p_rcontrol%p_rvectorAccess,1,p_rvector)
    call lsysbl_getbase_double (p_rvector,p_Ddata)

    ! -------------------------------------------------------------
    ! Step 1: Solve the primal system
    ! -------------------------------------------------------------

    if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
      call output_line ("Nonlin. space-time Residual: Solving the primal equation")
    end if

    ! Solve the primal equation, update the primal solution.
    output_iautoOutputIndent = output_iautoOutputIndent + 2

    call kkt_solvePrimal (rkktsystem,&
        rsolver%cspatialInitCondPolicy,rsolver%rkktSubsolvers,rlocalStat)

    output_iautoOutputIndent = output_iautoOutputIndent - 2
    
    call spacesl_sumStatistics(rlocalStat,rstatistics%rspaceslSolverStat)

    if (rsolver%rnewtonParams%ioutputLevel .ge. 3) then
      call output_line ("Nonlin. space-time Residual: Time for solving      : "//&
          trim(sys_sdL(rlocalStat%rtotalTime%delapsedReal,10)))
      
      call output_line ("Nonlin. space-time Residual: Time for space-defects: "//&
          trim(sys_sdL(rlocalStat%rtimeDefect%delapsedReal,10)))
      call output_line ("Nonlin. space-time Residual:   Time for space-RHS  : "//&
          trim(sys_sdL(rlocalStat%rtimeRHS%delapsedReal,10)))
      call output_line ("Nonlin. space-time Residual: Time for mat. assembly: "//&
          trim(sys_sdL(rlocalStat%rtimeMatrixAssembly%delapsedReal,10)))
      call output_line ("Nonlin. space-time Residual: Time for factorisation: "//&
          trim(sys_sdL(rlocalStat%rlssSolverStat%rtimeSymbolicFactorisation%delapsedReal+&
                       rlocalStat%rlssSolverStat%rtimeNumericFactorisation%delapsedReal,10)))
      call output_line ("Nonlin. space-time Residual: Time for space-solver : "//&
          trim(sys_sdL(rlocalStat%rlssSolverStat%rtotalTime%delapsedReal,10)))
    end if
    
    ! Add time to the time of the forward equation
    call stat_addTimers (rlocalStat%rtotalTime,rstatistics%rtimeForward)

    ! -------------------------------------------------------------
    ! Step 2: Solve the dual system
    ! -------------------------------------------------------------

    if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
      call output_line ("Nonlin. space-time Residual: Solving the dual equation")
    end if

    ! Solve the dual equation, update the dual solution.
    output_iautoOutputIndent = output_iautoOutputIndent + 2

    call kkt_solveDual (rkktsystem,&
        rsolver%cspatialInitCondPolicy,rsolver%rkktSubsolvers,rlocalStat)

    output_iautoOutputIndent = output_iautoOutputIndent - 2
    
    call spacesl_sumStatistics(rlocalStat,rstatistics%rspaceslSolverStat)
    
    ! Add time to the time of the backward equation
    call stat_addTimers (rlocalStat%rtotalTime,rstatistics%rtimeBackward)

    if (rsolver%rnewtonParams%ioutputLevel .ge. 3) then
      call output_line ("Nonlin. space-time Residual: Time for solving      : "//&
          trim(sys_sdL(rlocalStat%rtotalTime%delapsedReal,10)))

      call output_line ("Nonlin. space-time Residual: Time for space-defects: "//&
          trim(sys_sdL(rlocalStat%rtimeDefect%delapsedReal,10)))
      call output_line ("Nonlin. space-time Residual:   Time for space-RHS  : "//&
          trim(sys_sdL(rlocalStat%rtimeRHS%delapsedReal,10)))
      call output_line ("Nonlin. space-time Residual: Time for mat. assembly: "//&
          trim(sys_sdL(rlocalStat%rtimeMatrixAssembly%delapsedReal,10)))
      call output_line ("Nonlin. space-time Residual: Time for factorisation: "//&
          trim(sys_sdL(rlocalStat%rlssSolverStat%rtimeSymbolicFactorisation%delapsedReal+&
                       rlocalStat%rlssSolverStat%rtimeNumericFactorisation%delapsedReal,10)))
      call output_line ("Nonlin. space-time Residual: Time for space-solver : "//&
          trim(sys_sdL(rlocalStat%rlssSolverStat%rtotalTime%delapsedReal,10)))
    end if

    ! -------------------------------------------------------------
    ! Step 3: Calculate the intermediate control and the search direction   
    ! -------------------------------------------------------------

    ! The search direction is just the residual in the control equation.
    call kkt_calcControlRes (rkktsystem,rresidual,dres,iresnorm,&
        rsolver%rkktSubsolvers,rlocalStat)
    call spacesl_sumStatistics(rlocalStat,rstatistics%rspaceslSolverStat)

    call stat_stopTimer (rstatistics%rtotalTime)

  end subroutine

    ! ***************************************************************************

!<subroutine>

  subroutine newtonit_updateControl (rsolver,rkktsystem,rcorrection)
  
!<description>
  ! Calculates the basic (unprecondiotioned) search direction of the 
  ! Newton algorithm.
!</description>  

!<input>
  ! Parameters for the Newton iteration.
  ! The output parameters are changed according to the iteration.
  type(t_spacetimeNewton), intent(in) :: rsolver

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
        rcorrection,rsolver%rnewtonParams%domega,rkktsystem%p_rcontrol,1.0_DP)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_adNewton_setEps (rsolver,radNewtonParams,dresInit,dresLastIte)
  
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
  real(DP), intent(in) :: dresInit
  
  ! Residual obtained in the last nonlinear iteration.
  ! In the first call, this should be set to dresInit.
  real(DP), intent(in) :: dresLastIte
!</input>

!<inputoutput>
  ! Parameters for the Newton iteration.
  ! The output parameters are changed according to the iteration.
  type(t_spacetimeNewton), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>

    real(DP) :: depsAbs,depsRel,ddigitsGained, ddigitsToGain
    
    if ((dresInit .eq. 0.0_DP) .and. (dresLastIte .eq. 0.0_DP)) then
      
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
      
      ddigitsGained = dresLastIte/dresInit
      
      ddigitsToGain = min(radNewtonParams%dinexactNewtonTolRel*ddigitsGained,&
          ddigitsGained ** radNewtonParams%dinexactNewtonExponent)
          
      depsRel = 0.0_DP
      depsAbs = max(dresInit * ddigitsToGain,radNewtonParams%dinexactNewtonTolAbs)
      
      ! Do not gain too much.
      depsAbs = max(depsAbs,&
          max(dresInit * rsolver%riter%dtolrel * radNewtonParams%dinexactNewtonTolRel,&
              rsolver%riter%dtolabs * radNewtonParams%dinexactNewtonTolRel))
      
      if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
        call output_line ("Adaptive Newton: New stopping criterion. ||res|| < "//&
            trim(sys_sdEL(depsAbs,10)))
      end if
    end if
    
    ! Initialise the linear subsolver(s).
    call newtonlin_adNewton_setEps (rsolver%rlinsolParam,depsAbs,depsRel)
    
    ! Initialise the nonlinear and linear solver in space which are used
    ! for the calvulation of the current residual.
    !call spaceslh_setEps (rsolver%p_rsolverHierPrimal,depsAbs,depsRel)
    !call spaceslh_setEps (rsolver%p_rsolverHierDual,depsAbs,depsRel)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_solve (rsolver,rsolution,rstatistics)
  
!<description>
  ! Applies a Newton iteration to solve the space-time system.
!</description>
  
!<inputoutput>
  ! Parameters for the Newton iteration.
  ! The output parameters are changed according to the iteration.
  type(t_spacetimeNewton), intent(inout) :: rsolver

  ! Structure defining the solutions of the KKT system.
  type(t_kktsystem), intent(inout), target :: rsolution
!</inputoutput>

!<output>
  ! Statistic structure which receives the statistics of the iteration.
  type(t_newtonitSolverStat), intent(out) :: rstatistics
!</output>

!</subroutine>
   
    ! local variables
    type(t_controlSpace), pointer :: p_rdescentDir
    type(t_kktsystem), pointer :: p_rsolution
    type(t_kktsystemDirDeriv), pointer :: p_rsolutionDirDeriv
    type(t_timer) :: rtotalTime
    type(t_newtonlinSolverStat) :: rstatisticsLinSol
    real(DP) :: delapsedReal, dres
    type(t_newtonitSolverStat) :: rlocalStat

    ! Measure the total computational time
    call stat_startTimer(rtotalTime)
    
    ! Initialise data for the nonlinear iteration
    call newtonit_initData (rsolver,rsolution)

    ! Get a pointer to the solution and directional derivative on the maximum level
    call kkth_getKKTsystem (rsolver%p_rsolutionHierarchy,0,p_rsolution)
    call kkth_getKKTsystemDirDeriv (rsolver%p_rdirDerivHierarchy,0,p_rsolutionDirDeriv)
    
    ! Temporary vector(s)
    p_rdescentDir => rsolver%p_rdescentDir
    
    ! Prepare a structure that encapsules the directional derivative.
    
    ! Apply the Newton iteration
    call itc_initIteration(rsolver%riter)
    
    !rsolver%p_rsettingsSolver%rsettingsOptControl%rconstraints%rconstraintsDistCtrl%cconstraints = 0
    
    do while (.true.)
    
      ! The Newton iteration reads
      !
      !    u_n+1  =  u_n  -  [J''(u_n)]^-1  J'(u_n)
      !
      ! or in other words,
      !
      !    u_n+1  =  u_n  +  [J''(u_n)]^-1  d_n
      !
      ! with the residual   d_n = -J'(u_n)   specifying a "search direction".
      
      ! -------------------------------------------------------------
      ! Get the current residual / search direction
      ! -------------------------------------------------------------

      if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
        call output_line ("Space-time Newton: Calculating the residual")
      end if

      ! Compute the basic (unpreconditioned) search direction d_n.
      output_iautoOutputIndent = output_iautoOutputIndent + 2
      
      call newtonit_getResidual (rsolver,p_rsolution,p_rdescentDir,&
          dres,rsolver%rnewtonParams%iresnorm,rlocalStat)
      
      output_iautoOutputIndent = output_iautoOutputIndent - 2
      
      call newtonit_sumStatistics(rlocalStat,rstatistics,.false.)
      call stat_addTimers (rlocalStat%rtotalTime,rstatistics%rtimeDefect)
      
      if (rsolver%riter%cstatus .eq. ITC_STATUS_UNDEFINED) then
        ! Remember the initial residual
        call itc_initResidual(rsolver%riter,dres)
      else
        ! Push the residual, increase the iteration counter
        call itc_pushResidual(rsolver%riter,dres)
      end if

      if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
        call output_line ("Space-time Newton: Iteration "// &
            trim(sys_siL(rsolver%riter%niterations,10))// &
            ", ||res(u)|| = "// &
            trim(sys_sdEL(dres,10)))
      end if

      ! -------------------------------------------------------------
      ! Check for convergence
      ! -------------------------------------------------------------
      if (rsolver%riter%cstatus .ne. ITC_STATUS_CONTINUE) exit
      
      ! -------------------------------------------------------------
      ! Postprocessing.
      ! It has to be after newtonit_getResidual as newtonit_getResidual
      ! computes the corresponding y/lambda to the control u.
      ! We do it after the residual checks, so the final solution
      ! is not postprocessed.
      ! -------------------------------------------------------------
      if (rsolver%cpostprocessIterates .ge. 0) then
        call output_separator (OU_SEP_MINUS)

        call stat_startTimer (rstatistics%rtimePostprocessing)
        
        call optcpp_postprocessSubstep (rsolver%p_rsettingsSolver%rpostproc,&
            rsolver%p_rsettingsSolver%rspaceTimeHierPrimal%p_rfeHierarchy%nlevels,&
            rsolver%p_rsettingsSolver%rspaceTimeHierPrimal%p_rtimeHierarchy%nlevels,&
            p_rsolution,rsolver%rkktSubsolvers,rsolver%p_rsettingsSolver,&
            rsolver%cpostprocessIterates,rsolver%riter%niterations)
        
        call stat_stopTimer (rstatistics%rtimePostprocessing)
        
        call output_separator (OU_SEP_MINUS)
      end if

      ! -------------------------------------------------------------
      ! Adaptive Newton for the next iteration?
      ! -------------------------------------------------------------
      if (rsolver%rnewtonParams%ctypeIteration .eq. 3) then
        call newtonit_adNewton_setEps (&
            rsolver,rsolver%rnewtonParams%radaptiveNewton,&
            rsolver%riter%dresInitial,rsolver%riter%dresFinal)
      end if
      
      ! -------------------------------------------------------------
      ! Partial Newton for the next iteration?
      ! -------------------------------------------------------------

      ! By default use partial Newton.
      if (rsolver%rnewtonParams%ctypeIteration .eq. 1) then
        
        rsolver%rlinsolParam%ceqnflags = SPACESLH_EQNF_NONEWTON
      
      else
      
        select case (rsolver%rnewtonParams%radaptiveNewton%cpartialNewton)
        case (NEWTN_PN_FULLNEWTON)
          rsolver%rlinsolParam%ceqnflags = SPACESLH_EQNF_DEFAULT
          call output_line (&
              "Partial Newton not used. Better use partial Newton (dual) in the first step!",&
              OU_CLASS_WARNING,OU_MODE_STD,"newtonit_solve")

        case (NEWTN_PN_PARTIALNEWTON)
          rsolver%rlinsolParam%ceqnflags = SPACESLH_EQNF_NONEWTON

        case (NEWTN_PN_PARTIALNEWTONDUAL)
          rsolver%rlinsolParam%ceqnflags = SPACESLH_EQNF_NONEWTONDUAL
        
        end select
        
        ! Should we switch to full Newton?
        if (rsolver%rnewtonParams%radaptiveNewton%cpartialNewton .ne. NEWTN_PN_FULLNEWTON) then
        
          if (rsolver%riter%niterations .gt. &
              rsolver%rnewtonParams%radaptiveNewton%nmaxPartialNewtonIterations) then
            ! Full Newton
            rsolver%rlinsolParam%ceqnflags = SPACESLH_EQNF_DEFAULT
            
          else if (rsolver%riter%niterations .ge. &
              rsolver%rnewtonParams%radaptiveNewton%nminPartialNewtonIterations) then
            ! We are allowed to switch to the full Newton if some conditions
            ! are met.
            if (rsolver%rnewtonParams%radaptiveNewton%dtolAbsPartialNewton .gt. 0.0_DP) then
              ! Check the absolute residual
              if (dres .le. rsolver%rnewtonParams%radaptiveNewton%dtolAbsPartialNewton) then
                rsolver%rlinsolParam%ceqnflags = SPACESLH_EQNF_DEFAULT
              end if
            end if
            if (rsolver%rnewtonParams%radaptiveNewton%dtolRelPartialNewton .gt. 0.0_DP) then
              ! Check the relative residual
              if (dres .le. rsolver%rnewtonParams%radaptiveNewton%dtolRelPartialNewton* &
                  rsolver%riter%dresInitial) then
                rsolver%rlinsolParam%ceqnflags = SPACESLH_EQNF_DEFAULT
              end if
            end if
              
          end if
          
        end if
      
      end if      
      
      ! Print out the choice
      if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
        select case (rsolver%rlinsolParam%ceqnflags)
        case (SPACESLH_EQNF_DEFAULT)
          call output_line ("Space-time Newton: Full Newton selected.")
        case (SPACESLH_EQNF_NONEWTON)
          call output_line ("Space-time Newton: Partial Newton for primal and dual equation selected.")
        case (SPACESLH_EQNF_NONEWTONDUAL)
          call output_line ("Space-time Newton: Partial Newton for dual equation selected.")
        end select
      end if

      ! -------------------------------------------------------------
      ! Preconditioning with the Newton matrix
      ! -------------------------------------------------------------
      
      if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
        call output_line ("Space-time Newton: Preconditioning")
      end if

      ! Initialise data arrays in the linear subsolver.
      call stat_startTimer(rstatistics%rtimeProlRest)
      call newtonlin_initNonlinearData (rsolver%rlinsolParam,rsolver%p_rsolutionHierarchy)
      call stat_stopTimer(rstatistics%rtimeProlRest)
      
      call newtonlin_initData (rsolver%rlinsolParam,rsolver%p_rsolutionHierarchy)

      ! Actual Newton iteration. Apply the Newton preconditioner
      ! to get the Newton search direction:
      !
      !    J''(u_n) g_n  =  d_n
      !
      ! The control on the maximum level of p_rdirDerivHierarchy
      ! (identified by p_rsolutionDirDeriv%p_rcontrolLin) receives the result.
      output_iautoOutputIndent = output_iautoOutputIndent + 2
      
      call newtonlin_precond (rsolver%rlinsolParam,&
          rsolver%p_rdirDerivHierarchy,p_rdescentDir,rstatisticsLinSol)
      
      output_iautoOutputIndent = output_iautoOutputIndent - 2
      
      ! Statistics
      call newtonlin_sumStatistics (&
          rstatisticsLinSol,rstatistics%rnewtonlinSolverStat,NLIN_STYPE_LINSOL)
      
      ! Clean up data in the linear subsolver
      call newtonlin_doneData (rsolver%rlinsolParam)

      ! -------------------------------------------------------------
      ! Update of the solution
      ! -------------------------------------------------------------

      ! Update the control according to the search direction:
      !
      !    u_n+1  =  u_n  +  g_n
      !
      ! or to any configured step-length control rule.
      call newtonit_updateControl (&
          rsolver,p_rsolution,p_rsolutionDirDeriv%p_rcontrolLin)

!    ! DEBUG!!!      
!    call kktsp_dualLinearComb (&
!        p_rsolutionDirDeriv%p_rdualSolLin,p_rsolution%p_rdualSol,1.0_DP,1.0_DP)
!    call kkt_calcControlRes (p_rsolution,p_rdescentDir,&
!        rsolver%rnewtonParams%dresFinal,rsolver%rnewtonParams%iresnorm)
!    !-> dres = u+du + 1/alpha (lambda+dlambda) =0 -> is ok!

      ! -------------------------------------------------------------
      ! Proceed with the next iteration
      ! -------------------------------------------------------------
      if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
        call output_separator (OU_SEP_MINUS)
        call output_line ("Space-time Newton: Time for the linear solver = "//&
            sys_sdL(rstatisticsLinSol%rtotalTime%delapsedReal,10))
        call stat_sampleTimer(rtotalTime,delapsedReal)
        call output_line ("Space-time Newton: Computation time           = "//&
            sys_sdL(delapsedReal,10))
      end if

      call output_separator (OU_SEP_MINUS)
      
      !rsolver%p_rsettingsSolver%rsettingsOptControl%rconstraints%rconstraintsDistCtrl%cconstraints = 1

    end do
    
    ! Release data
    call newtonit_doneData (rsolver)
    
    call stat_stopTimer(rtotalTime)
    call stat_addTimers(rtotalTime,rstatistics%rtotalTime)
    
    rstatistics%niterations = rstatistics%niterations + rsolver%riter%niterations
    
    ! Statistics
    if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
      call output_lbrk()
      call output_line ("Space-time Newton: Statistics.")
      call output_lbrk()
      call itc_printStatistics(rsolver%riter)
    end if
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_initStructure (rsolver,rkktsystemHierarchy,&
      rprjHierSpaceTimePrimal,rprjHierSpaceTimeDual,rprjHierSpaceTimeControl)
  
!<description>
  ! Structural initialisation of the Newton solver.
!</description>

!<input>
  ! Defines the basic hierarchy of the solutions of the KKT system.
  ! This can be a "template" structure, i.e., memory for the solutions
  ! in rkktsystemHierarchy does not have to be allocated.
  type(t_kktsystemHierarchy), intent(in), target :: rkktsystemHierarchy

  ! Projection hierarchy for the interlevel projection in space/time, primal space.
  type(t_sptiProjHierarchyBlock), intent(in), target :: rprjHierSpaceTimePrimal

  ! Projection hierarchy for the interlevel projection in space/time, dual space.
  type(t_sptiProjHierarchyBlock), intent(in), target :: rprjHierSpaceTimeDual

  ! Projection hierarchy for the interlevel projection in space/time, control space.
  type(t_sptiProjHierarchyBlock), intent(in), target :: rprjHierSpaceTimeControl
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_spacetimeNewton), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>

    ! Remember the structure of the solutions.
    rsolver%p_rkktsystemHierarchy => rkktsystemHierarchy
    
    ! Initialise the structures of the linear subsolver.
    call newtonlin_initStructure (rsolver%rlinsolParam,rkktsystemHierarchy,&
        rprjHierSpaceTimePrimal,rprjHierSpaceTimeDual,rprjHierSpaceTimeControl)
   
  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_initData (rsolver,rsolution)
  
!<description>
  ! Final preparation of the Newton solver.
!</description>

!<input>
  ! Structure defining the solutions of the KKT system.
  type(t_kktsystem), intent(in), target :: rsolution
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_spacetimeNewton), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>

    ! Allocate memory for the descent direction vector
    allocate(rsolver%p_rdescentDir)

    ! Create a vector holding the descent direction (highest level)
    call kkth_initControl (rsolver%p_rdescentDir,rsolver%p_rkktsystemHierarchy,0)
    
    ! Create temporary memory for the nonlinearity and the search direction.
    allocate(rsolver%p_rsolutionHierarchy)

    call kkth_initHierarchy (rsolver%p_rsolutionHierarchy,&
        rsolver%p_rkktsystemHierarchy%p_roperatorAsmHier,&
        rsolver%p_rkktsystemHierarchy%p_rspaceTimeHierPrimal,&
        rsolver%p_rkktsystemHierarchy%p_rspaceTimeHierDual,&
        rsolver%p_rkktsystemHierarchy%p_rspaceTimeHierControl,.true.,&
        rsolver%p_rsettingsSolver%roptcBDC,rsolution)

    ! create tempoprary memory for the search direction connected with
    ! the above solution hierarchy.
    allocate(rsolver%p_rdirDerivHierarchy)
    
    call kkth_initHierarchyDirDeriv (&
        rsolver%p_rdirDerivHierarchy,rsolver%p_rsolutionHierarchy)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_doneData (rsolver)
  
!<description>
  ! Cleanup of the data initalised in newtonit_initData.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_spacetimeNewton), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>

    ! Release the descent direction, memory for the directional derivative,...
    call kktsp_doneControlVector (rsolver%p_rdescentDir)
    deallocate(rsolver%p_rdescentDir)
    
    call kkth_doneHierarchyDirDeriv (rsolver%p_rdirDerivHierarchy)
    deallocate(rsolver%p_rdirDerivHierarchy)

    call kkth_doneHierarchy (rsolver%p_rsolutionHierarchy)
    deallocate(rsolver%p_rsolutionHierarchy)
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_doneStructure (rsolver)
  
!<description>
  ! Cleanup of the data initalised in newtonit_initStructure.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_spacetimeNewton), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>
   
    ! Release structures in the subsolver
    call newtonlin_doneStructure (rsolver%rlinsolParam)
   
    ! Detach the solution structure
    nullify(rsolver%p_rkktsystemHierarchy)
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_done (rsolver)
  
!<description>
  ! Clean up the Newton iteration.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_spacetimeNewton), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>

    ! Release the linear subsolver.
    call newtonlin_done (rsolver%rlinsolParam)
    
    ! Release the KKT subsolvers.
    call kkt_doneSubsolvers (rsolver%rkktSubsolvers)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_clearStatistics(rstatistics)
  
!<description>
  ! Resets a statistic structure.
!</description>

!<inputoutput>
  ! Structure to be reset.
  type(t_newtonitSolverStat), intent(inout) :: rstatistics
!</inputoutput>

!</subroutine>

    rstatistics%niterations = 0
    
    call stat_clearTimer(rstatistics%rtotalTime)
    call stat_clearTimer(rstatistics%rtimeForward)
    call stat_clearTimer(rstatistics%rtimeBackward)
    call stat_clearTimer(rstatistics%rtimeDefect)
    call stat_clearTimer(rstatistics%rtimePostprocessing)
    
    call spacesl_clearStatistics (rstatistics%rspaceslSolverStat)
    call newtonlin_clearStatistics (rstatistics%rnewtonlinSolverStat)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_sumStatistics(rstatistics1,rstatistics2,btotalTime)
  
!<description>
  ! Sums up the data of rstatistics1 to the data in rstatistics2.
!</description>

!<input>
  ! Source structure
  type(t_newtonitSolverStat), intent(in) :: rstatistics1

  ! OPTIONAL: Whether or not to sum up the total time.
  ! If not present, TRUE is assumed.
  logical, intent(in), optional :: btotalTime
!</input>

!<inputoutput>
  ! Destination structure.
  type(t_newtonitSolverStat), intent(inout) :: rstatistics2
!</inputoutput>

!</subroutine>

    rstatistics2%niterations = rstatistics2%niterations + rstatistics1%niterations
    
    if (.not. present(btotalTime)) then
      call stat_addTimers(rstatistics1%rtotalTime,         rstatistics2%rtotalTime)
    else if (btotalTime) then
      call stat_addTimers(rstatistics1%rtotalTime,         rstatistics2%rtotalTime)
    end if
    call stat_addTimers(rstatistics1%rtimeForward,       rstatistics2%rtimeForward)
    call stat_addTimers(rstatistics1%rtimeBackward,      rstatistics2%rtimeBackward)
    call stat_addTimers(rstatistics1%rtimeDefect,        rstatistics2%rtimeDefect)
    call stat_addTimers(rstatistics1%rtimeProlRest,      rstatistics2%rtimeProlRest)
    call stat_addTimers(rstatistics1%rtimePostprocessing,rstatistics2%rtimePostprocessing)
    
    call spacesl_sumStatistics (&
        rstatistics1%rspaceslSolverStat,rstatistics2%rspaceslSolverStat)

    call newtonlin_sumStatistics (&
        rstatistics1%rnewtonlinSolverStat,rstatistics2%rnewtonlinSolverStat,NLIN_STYPE_LINSOL)

  end subroutine

end module
