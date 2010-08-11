!##############################################################################
!# ****************************************************************************
!# <name> transport_application </name>
!# ****************************************************************************
!#
!# <purpose>
!# This application solves the generic time-dependent conservation law
!#
!#   $$\frac{\partial u}{\partial t}+\nabla\cdot{\bf f}(u)+r(u) = b$$
!#
!# for a scalar quantity $u$ in the domain $\Omega$ in 1D, 2D and 3D.
!#
!#
!# The flux function ${\bf f}(u)$ can be nonlinear, e.g. for Burgers`
!# equation or linear as it is the case for linear transport.
!#
!# In the linear case, the flux function is typically given by
!#
!#   $${\bf f}(u):={\bf v}u-D\nabla u$$
!#
!# whereby $\bf v$ denotes an externally defined velocity profile and
!# $D$ is the vsicosity tensor. It may be anisotropic
!#
!#   $$D=\left[\begin{array}{ccc}
!#       d_{11} & d_{12} & d_{13}\\
!#       d_{21} & d_{22} & d_{23}\\
!#       d_{31} & d_{32} & d_{33}
!#       \end{array}\right]$$
!#
!# or it can also reduce to the scalar coefficient
!#
!#   $$D=dI,\quad d_{ij}=0\,\forall j\ne i$$
!#
!#
!# In the nonlinear case, the flux function depends on the unknown
!# solution $u$ or its first and/or second derivative:
!#
!#   $${\bf f}={\bf f}(u,\nabla u, \Delta u)$$
!#
!#
!# Moreover, this application can handle space-time simulations
!# in 1D and 2D, whereby the flux function is defined such that
!#
!#   $$\nabla\cdot{\bf f}(u):=u_t+\frac{\partial f(u)}{\partial x}$$
!#
!# The above equation stands for a 1D/2D conservation law
!# which is solved in the 2D/3D space-time domain $\Omega=$.
!#
!#
!# The reactive term $r(u)$ is ignored at the moment. However, some
!# proper linearisation is required to preserve positivity.
!#
!#
!# The load vector $b$ is constant throughout the simulation.
!#
!#
!# The spatial discretisation is perform by means of the algebraic
!# flux correction (AFC) paradigm by Kuzmin, Moeller and Turek. In
!# particular, high-resolution finite element schemes of TVD- and
!# FCT-type are available. For the temporal discretisation, the
!# two-level theta-scheme is employed, whereby $\theta\in(0,1]$.
!#
!# Dynamic mesh adaptation is based on the red-green strategy, whereby
!# mixed triangulations are supported. Error estimation can either be
!# accomplished using gradient-recovery techniques or following goal-
!# oriented strategies.
!#
!#
!# The following routines are available:
!#
!# 1.) transp_app
!#     -> The main routine of the application called from the main
!#        program. The routine gets all required information from the
!#        parameter list which needs to be initialised and filled in
!#        the main program. It then works black-box, that is, it
!#        determines the solution algorithm to be used and performs
!#        the simulation. The user should only have to modify this
!#        routine if another solution algorithm is implemented.
!#
!# 2.) transp_initSolvers
!#     -> Initializes the solve structures from the parameter list.
!#
!# 3.) transp_initProblem
!#     -> Initializes the global problem structure based on the
!#        parameter settings given by the parameter list. This routine
!#        is quite universal, that is, it prepares the internal
!#        structure of the global problem and generates a linked
!#        list of problem levels used in the multigrid hierarchy.
!#
!# 4.) transp_initProblemLevel
!#     -> Initializes the individual problem level based on the
!#        parameter settings given by the parameter list.
!#        This routine is called repeatedly by the global
!#        initialisation routine transp_initAllProblemLevels.
!#
!# 5.) transp_initAllProblemLevels
!#     -> Initializes ALL problem levels attached to the global
!#        problem structure based on the parameter settings
!#        given by the parameter list.
!#
!# 6.) transp_initSolution
!#     -> Initializes the solution vector based on the parameter
!#        settings given by the parameter list
!#
!# 7.) transp_initRHS
!#     -> Initializes the right-hand side vector based on the
!#        parameter settings given by the application descriptor
!#
!# 8.) transp_initTargetFunc
!#     -> Initializes the target functional for the dual problem
!#
!# 9.) transp_outputSolution
!#     -> Outputs the solution vector to file in UCD format
!#
!# 10.) transp_outputStatistics
!#      -> Outputs the application statitics
!#
!# 11.) transp_estimateTargetFuncError
!#      -> Estimates the error in the quantity of interest
!#
!# 12.) transp_estimateRecoveryError
!#      -> Estimates the solution error using recovery techniques
!#
!# 13.) transp_adaptTriangulation
!#      -> Performs h-adaptation for the given triangulation
!#
!# 14.) transp_solveTransientPrimal
!#      -> Solves the primal formulation of the time-dependent
!#         convection-diffusion-reaction equation.
!#
!# 15.) transp_solveTransientPrimalDual
!#      -> Solves the primal and the dual formulation of the time-dependent
!#         convection-diffusion-reaction equation using pseudo time-stepping.
!#
!# 16.) transp_solvePseudoTransientPrimal
!#      -> Solves the primal formulation of the steady
!#         convection-diffusion-reaction equation using pseudo time-stepping.
!#
!# 17.) transp_solvePseudoTransientPrimalDual
!#      -> Solves the primal and the dual formulation of the steady
!#         convection-diffusion-reaction equation using pseudo time-stepping.
!#
!# 18.) transp_solveSteadyStatePrimal
!#      -> Solves the primal formulation of the steady
!#         convection-diffusion-reaction equation directly
!#
!# 19.) transp_solveSteadyStatePrimalDual
!#      -> Solves the primal and the dual formulation of the steady
!#         convection-diffusion-reaction equation directly
!#
!#
!# The following auxiliary routines are available:
!#
!# 1.) transp_parseCmdlArguments
!#     -> Parses the list of commandline arguments and overwrites
!#        parameter values from the parameter files
!#
!# 2.) transp_adjustParameterlist
!#     -> Adjust the parameter list depending on internal data
!#
!#
!# Frequently asked questions
!# --------------------------
!#
!# 1.) How to add a new sub-model, e.g. 1D-Burgers` equation in space time?
!#
!#
!# 2.) How to add new finite elements in spatial dimension xD?
!#
!#     -> This application builds on the group finite element formulation
!#        suggested by Fletcher. In essence, the solution vector $u$ and
!#        the flux function ${\bf f}(u)$ is approximated using the same
!#        basis functions. As a result, the matrices can be assembled
!#        once at the beginning of the simulation and each time the mesh
!#        has changed (e.g. due to h-adaptation) so that no expensive
!#        numerical integration needs to be performed in each step.
!#        If you want to add a new finite element, then you have to
!#        modify subroutines transp_initProblemLevel which generates
!#        the global matrix structure and builds all FE matrices.
!#
!# </purpose>
!##############################################################################

module transport_application

  use afcstabilisation
  use bilinearformevaluation
  use boundary
  use boundarycondaux
  use boundaryfilter
  use collection
  use cubature
  use derivatives
  use element
  use flagship_basic
  use fparser
  use fsystem
  use genoutput
  use graph
  use groupfemscalar
  use hadaptaux
  use hadaptivity
  use linearformevaluation
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use pprocerror
  use pprocgradients
  use pprocindicator
  use pprocsolution
  use problem
  use scalarpde
  use solveraux
  use spatialdiscretisation
  use statistics
  use stdoperators
  use storage
  use timestep
  use timestepaux
  use ucd

  use transport_basic
  use transport_callback
  use transport_callback1d
  use transport_callback2d
  use transport_callback3d

  implicit none

  private
  public :: transp_app
  public :: transp_adaptTriangulation
  public :: transp_adjustParameterlist
  public :: transp_estimateRecoveryError
  public :: transp_estimateTargetFuncError
  public :: transp_initAllProblemLevels
  public :: transp_initProblem
  public :: transp_initProblemLevel
  public :: transp_initRHS
  public :: transp_initSolution
  public :: transp_initSolvers
  public :: transp_initTargetFunc
  public :: transp_outputSolution
  public :: transp_outputStatistics
  public :: transp_solveTransientPrimal
  public :: transp_solveTransientPrimalDual
  public :: transp_solvePseudoTransientPrimal
  public :: transp_solvePseudoTransientPrimalDual
  public :: transp_solveSteadyStatePrimal
  public :: transp_solveSteadyStatePrimalDual

contains

  !*****************************************************************************

!<subroutine>

  subroutine transp_app(rparlist, ssectionName)

!<description>
    ! This is the main application for the convection-diffusion-reaction
    ! benchmark problem. It is a so-called driver routine which can be
    ! used to start a standalone convection-diffusion-reaction application.
!</description>

!<input>
    ! name of the top-most section of the application
    character(len=*), intent(in) :: ssectionName
!</input>

!<inputoutput>
    ! parameter list
    type(t_parlist), intent(inout) :: rparlist
!</inputoutput>
!</subroutine>

    !****************************************************************************
    ! Structures required for this application

    ! Global collection which is used to pass arguments to callback routines
    type(t_collection) :: rcollection

    ! Global function parser which is used to evaluate analytical functions
    type(t_fparser) :: rfparser

    ! Boundary condition structure for the primal problem
    type(t_boundaryCondition) :: rbdrCondPrimal

    ! Boundary condition structure for the dual problem
    type(t_boundaryCondition) :: rbdrCondDual

    ! Problem structure which holds all internal data (vectors/matrices)
    ! for the convection-diffusion-reaction application
    type(t_problem) :: rproblem

    ! Time-stepping structure
    type(t_timestep) :: rtimestep

    ! Global solver structure
    type(t_solver) :: rsolver

    ! Solution vector for the primal problem
    type(t_vectorBlock) :: rsolutionPrimal

    ! Solution vector for the dual problem
    type(t_vectorBlock) :: rsolutionDual

    ! Timer for the total solution process
    type(t_timer) :: rtimerTotal

    ! Timer for the solution process
    type(t_timer) :: rtimerSolution

    ! Timer for the adaptation process
    type(t_timer) :: rtimerAdaptation

    ! Timer for the error estimation process
    type(t_timer) :: rtimerErrorEstimation

    ! Timer for the triangulation process
    type(t_timer) :: rtimerTriangulation

    ! Timer for the assembly of constant coefficient matrices
    type(t_timer) :: rtimerAssemblyCoeff

    ! Timer for the assembly of system matrices
    type(t_timer) :: rtimerAssemblyMatrix

    ! Timer for the assembly of residual/right-hand side vectors
    type(t_timer) :: rtimerAssemblyVector

    ! Timer for pre- and post-processing
    type(t_timer) :: rtimerPrePostprocess

    ! Parameter file and section names
    character(LEN=SYS_STRLEN) :: sindatfileName
    character(LEN=SYS_STRLEN) :: sbdrcondName
    character(LEN=SYS_STRLEN) :: algorithm

    ! local variables
    integer :: systemMatrix, ndimension


    ! Start total time measurement
    call stat_startTimer(rtimerTotal)

    !---------------------------------------------------------------------------
    ! Pre-processing
    !---------------------------------------------------------------------------

    ! Start time measurement
    call stat_startTimer(rtimerPrePostprocess, STAT_TIMERSHORT)

    ! Overwrite configuration from command line arguments. After this
    ! subroutine has been called, the parameter list remains unchanged
    ! unless the used updates some parameter values interactively.
    call transp_parseCmdlArguments(rparlist)
    call transp_adjustParameterlist(rparlist, ssectionName)

    ! Initialize global collection structure
    call collct_init(rcollection)

    !  Attach the parameter list and the timers to the collection
    call collct_setvalue_parlst(rcollection,&
        'rparlist', rparlist, .true.)
    call collct_setvalue_timer(rcollection,&
        'rtimerSolution', rtimerSolution, .true.)
    call collct_setvalue_timer(rcollection,&
        'rtimerAdaptation', rtimerAdaptation, .true.)
    call collct_setvalue_timer(rcollection,&
        'rtimerErrorEstimation', rtimerErrorEstimation, .true.)
    call collct_setvalue_timer(rcollection,&
        'rtimerTriangulation', rtimerTriangulation, .true.)
    call collct_setvalue_timer(rcollection,&
        'rtimerAssemblyCoeff', rtimerAssemblyCoeff, .true.)
    call collct_setvalue_timer(rcollection,&
        'rtimerAssemblyMatrix', rtimerAssemblyMatrix, .true.)
    call collct_setvalue_timer(rcollection,&
        'rtimerAssemblyVector', rtimerAssemblyVector, .true.)
    call collct_setvalue_timer(rcollection,&
        'rtimerPrePostprocess', rtimerPrePostprocess, .true.)

    ! Create function parser
    call fparser_create(rfparser, 100)

    ! Read in all constants, predefined expressions
    ! and functions from the parameter file
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'indatfile', sindatfileName)
    call fparser_parseFileForKeyword(rfparser,&
        sindatfileName, 'defconst', FPAR_CONSTANT)
    call fparser_parseFileForKeyword(rfparser,&
        sindatfileName, 'defexpr', FPAR_EXPRESSION)
    call fparser_parseFileForKeyword(rfparser,&
        sindatfileName, 'deffunc', FPAR_FUNCTION)

    ! Attach the function parser to the collection
    call collct_setvalue_pars(rcollection, 'rfparser', rfparser, .true.)

    ! Initialize the solver structures
    call transp_initSolvers(rparlist, ssectionName, rtimestep, rsolver)

    ! Initialize the abstract problem structure
    call transp_initProblem(rparlist, ssectionName,&
        solver_getMinimumMultigridlevel(rsolver),&
        solver_getMaximumMultigridlevel(rsolver), rproblem)

    ! Initialize the individual problem levels
    call transp_initAllProblemLevels(rparlist, ssectionName,&
        rproblem, rcollection)

    ! Prepare internal data arrays of the solver structure
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'systemMatrix', systemMatrix)
    call flagship_updateSolverMatrix(rproblem%p_rproblemLevelMax,&
        rsolver, systemMatrix, SYSTEM_INTERLEAVEFORMAT, UPDMAT_ALL)
    call solver_updateStructure(rsolver)

    ! Stop time measurement for pre-processing
    call stat_stopTimer(rtimerPrePostprocess)


    !---------------------------------------------------------------------------
    ! Solution algorithm
    !---------------------------------------------------------------------------

    if (rtimestep%dfinalTime .gt. rtimestep%dinitialTime) then

      ! Get global configuration from parameter list
      call parlst_getvalue_string(rparlist,&
          ssectionName, 'algorithm', algorithm)
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'ndimension', ndimension)
      call parlst_getvalue_string(rparlist,&
          ssectionName, 'sprimalbdrcondname', sbdrcondName)

      ! The boundary conditions for the primal problem are required
      ! for all solution strategies. Hence, set the pointer and
      ! initialize boundary conditions from the parameter file.
      call bdrc_readBoundaryCondition(rbdrCondPrimal,&
          sindatfileName, '['//trim(sbdrcondName)//']',&
          ndimension, transp_getBdrCondExprNumber)

      ! What solution algorithm should be applied?
      if (trim(algorithm) .eq. 'transient_primal') then
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Solve the primal formulation for
        ! the time-dependent problem
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call transp_solveTransientPrimal(rparlist, ssectionName,&
            rbdrCondPrimal, rproblem, rtimestep, rsolver,&
            rsolutionPrimal, rcollection)

        call transp_outputSolution(rparlist, ssectionName,&
            rproblem%p_rproblemLevelMax,&
            rsolutionPrimal=rsolutionPrimal,&
            dtime=rtimestep%dTime)


      elseif (trim(algorithm) .eq. 'transient_primaldual') then
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Solve the primal and dual formulation for
        ! the time-dependent problem
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call parlst_getvalue_string(rparlist, ssectionName,&
            'sdualbdrcondname', sbdrcondName)
        call bdrc_readBoundaryCondition(rbdrCondDual,&
            sindatfileName, '['//trim(sbdrcondName)//']',&
            ndimension, transp_getBdrCondExprNumber)

        call transp_solveTransientPrimalDual(rparlist, ssectionName,&
            rbdrCondPrimal, rbdrCondDual, rproblem, rtimestep,&
            rsolver, rsolutionPrimal, rsolutionDual, rcollection)

        call transp_outputSolution(rparlist, ssectionName,&
            rproblem%p_rproblemLevelMax,&
            rsolutionPrimal, rsolutionDual, rtimestep%dTime)


      elseif (trim(algorithm) .eq. 'pseudotransient_primal') then
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Solve the primal formulation for
        ! the pseudo time-dependent problem
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call transp_solvePseudoTransientPrimal(rparlist, ssectionName,&
            rbdrCondPrimal, rproblem, rtimestep, rsolver,&
            rsolutionPrimal, rcollection)

        call transp_outputSolution(rparlist, ssectionName,&
            rproblem%p_rproblemLevelMax,&
            rsolutionPrimal=rsolutionPrimal,&
            dtime=rtimestep%dTime)


      elseif (trim(algorithm) .eq. 'pseudotransient_primaldual') then
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Solve the primal and dual formulation for
        ! the pseudo time-dependent problem
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call parlst_getvalue_string(rparlist, ssectionName,&
            'sdualbdrcondname', sbdrcondName)
        call bdrc_readBoundaryCondition(rbdrCondDual,&
            sindatfileName, '['//trim(sbdrcondName)//']',&
            ndimension, transp_getBdrCondExprNumber)

        call transp_solvePseudoTransientPrimalDual(rparlist, ssectionName,&
            rbdrCondPrimal, rbdrCondDual, rproblem, rtimestep,&
            rsolver, rsolutionPrimal, rsolutionDual, rcollection)

        call transp_outputSolution(rparlist, ssectionName,&
            rproblem%p_rproblemLevelMax,&
            rsolutionPrimal, rsolutionDual, rtimestep%dTime)


      elseif (trim(algorithm) .eq. 'stationary_primal') then
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Solve the primal formulation for
        ! the stationary problem
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call transp_solveSteadyStatePrimal(rparlist, ssectionName,&
            rbdrCondPrimal, rproblem, rtimestep, rsolver,&
            rsolutionPrimal, rcollection)

        call transp_outputSolution(rparlist, ssectionName,&
            rproblem%p_rproblemLevelMax,&
            rsolutionPrimal=rsolutionPrimal,&
            dtime=rtimestep%dTime)


      elseif (trim(algorithm) .eq. 'stationary_primaldual') then
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Solve the primal and dual formulation for
        ! the stationary problem
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call parlst_getvalue_string(rparlist, ssectionName,&
            'sdualbdrcondname', sbdrcondName)
        call bdrc_readBoundaryCondition(rbdrCondDual,&
            sindatfileName, '['//trim(sbdrcondName)//']',&
            ndimension, transp_getBdrCondExprNumber)

        call transp_solveSteadyStatePrimalDual(rparlist, ssectionName,&
            rbdrCondPrimal, rbdrCondDual, rproblem, rtimestep,&
            rsolver, rsolutionPrimal, rsolutionDual, rcollection)

        call transp_outputSolution(rparlist, ssectionName,&
            rproblem%p_rproblemLevelMax,&
            rsolutionPrimal, rsolutionDual, rtimestep%dTime)

      else
        call output_line(trim(algorithm)//' is not a valid solution algorithm!',&
            OU_CLASS_ERROR,OU_MODE_STD,'transp_app')
        call sys_halt()
      end if

    else

      ! Just output the computational mesh and exit
      call transp_outputSolution(rparlist, ssectionName,&
          rproblem%p_rproblemLevelMax)

    end if


    !---------------------------------------------------------------------------
    ! Post-processing
    !---------------------------------------------------------------------------

    ! Start time measurement for pre-processing
    call stat_startTimer(rtimerPrepostProcess, STAT_TIMERSHORT)

    ! Release time-stepping
    call tstep_infoTimestep(rtimestep, .false.)
    call tstep_releaseTimestep(rtimestep)

    ! Release solver
    call solver_infoSolver(rsolver, .false.)
    call solver_releaseSolver(rsolver)

    ! Release problem structure
    call problem_releaseProblem(rproblem)

    ! Release boundary conditions
    call bdrc_release(rbdrCondPrimal)
    call bdrc_release(rbdrCondDual)

    ! Release vectors
    call lsysbl_releaseVector(rsolutionPrimal)
    call lsysbl_releaseVector(rsolutionDual)

    ! Release function parser
    call fparser_release(rfparser)

    ! Stop time measurement for pre-processing
    call stat_stopTimer(rtimerPrePostprocess)

    ! Stop time measurement for total time measurement
    call stat_stopTimer(rtimerTotal)

    ! Output statistics
    call transp_outputStatistics(rtimerTotal, rcollection)

    ! Release collection
    call collct_done(rcollection)

  end subroutine transp_app

  !*****************************************************************************

!<subroutine>

  subroutine transp_initSolvers(rparlist, ssectionName,&
      rtimestep, rsolver)

!<description>
    ! This subroutine initializes the time-stepping structure
    ! and the top-level solver structure using the
    ! parameter settings defined in the parameter list
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName
!</input>

!<output>
    ! time-stepping structure
    type(t_timestep), intent(out) :: rtimestep

    ! solver struchture
    type(t_solver), intent(out) :: rsolver
!</output>
!</subroutine>

    ! section name for the top-level solver
    character(LEN=SYS_STRLEN) :: ssolverName

    ! section name for time-stepping scheme
    character(LEN=SYS_STRLEN) :: stimestepName

    ! local variables
    integer :: nlmin, nlmax


    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist, ssectionName, 'timestep', stimestepName)
    call parlst_getvalue_string(rparlist, ssectionName, 'solver',   ssolverName)

    ! Initialize time-stepping
    call tstep_createTimestep(rparlist, stimestepName, rtimestep)

    ! Initialize solver structure
    call solver_createSolver(rparlist, ssolverName, rsolver)
    if (rsolver%csolverType .eq. SV_FMG) then
      nlmin = rsolver%p_solverMultigrid%nlmin
      nlmax = rsolver%p_solverMultigrid%nlmax
      call solver_adjustHierarchy(rsolver, nlmin, nlmax)
    else
      call solver_adjustHierarchy(rsolver)
    end if
    call solver_updateStructure(rsolver)

  end subroutine transp_initSolvers

  !*****************************************************************************

!<subroutine>

  subroutine transp_initProblem(rparlist, ssectionName,&
      nlmin, nlmax, rproblem)

!<description>
    ! This subroutine initializes the abstract problem structure
    ! using the parameters settings defined in the parameter list
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName

    ! minimum/maximum problem level
    integer, intent(in) :: nlmin, nlmax
!</input>

!<output>
    ! problem structure
    type(t_problem), intent(out) :: rproblem
!</output>
!</subroutine>

    ! section names
    character(LEN=SYS_STRLEN) :: sconvection
    character(LEN=SYS_STRLEN) :: sdiffusion

    ! abstract problem descriptor
    type(t_problemDescriptor) :: rproblemDescriptor

    ! pointer to the problem level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! local variables
    integer :: convectionAFC
    integer :: diffusionAFC
    integer :: iconvToTria
    integer :: ivelocitytype


    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'trifile', rproblemDescriptor%trifile)
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'prmfile', rproblemDescriptor%prmfile, '')
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'ndimension', rproblemDescriptor%ndimension)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'ivelocitytype', ivelocitytype)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'iconvtotria', iconvToTria, 0)
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'diffusion', sdiffusion)
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'convection', sconvection)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'convectionAFC', convectionAFC)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'diffusionAFC', diffusionAFC)

    ! Set additional problem descriptor
    rproblemDescriptor%ndiscretisation = 1 ! one discretisation
    rproblemDescriptor%nafcstab        = 2 ! convective and diffusive stabilisation
    rproblemDescriptor%nlmin           = nlmin
    rproblemDescriptor%nlmax           = nlmax
    rproblemDescriptor%nmatrixScalar   = rproblemDescriptor%ndimension + 7
    rproblemDescriptor%nmatrixBlock    = 0
    rproblemDescriptor%nvectorScalar   = 0
    rproblemDescriptor%nvectorBlock    = merge(1,0,transp_hasVelocityVector(ivelocitytype)) ! velocity field

    ! Check if quadrilaterals should be converted to triangles
    if (iconvToTria .ne. 0) then
      rproblemDescriptor%iproblemSpec = rproblemDescriptor%iproblemSpec &
                                      + PROBDESC_MSPEC_CONVTRIANGLES
    end if

    ! Initialize problem structure
    call problem_initProblem(rproblemDescriptor, rproblem)

    ! Loop over all problem levels
    p_rproblemLevel => rproblem%p_rproblemLevelMax
    do while(associated(p_rproblemLevel))

      ! Initialize the stabilisation structure for convection (if required)
      if (convectionAFC > 0)&
          call afcstab_initFromParameterlist(rparlist, sconvection,&
          p_rproblemLevel%Rafcstab(convectionAFC))

      ! Initialize the stabilisation structure for diffusion (if required)
      if (diffusionAFC > 0)&
          call afcstab_initFromParameterlist(rparlist, sdiffusion,&
          p_rproblemLevel%Rafcstab(diffusionAFC))

      ! Switch to next coarser level
      p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
    end do

  end subroutine transp_initProblem

  !*****************************************************************************

!<subroutine>

  subroutine transp_initProblemLevel(rparlist, ssectionName,&
      rproblemLevel, rcollection)

!<description>
    ! This subroutine initielizes the individual problem level. It
    ! generates the discretisation, the template matrix and the
    ! coefficient matrices as duplicates of the template matrix.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout), target :: rproblemLevel

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</output>
!</subroutine>

    ! local variables
    integer :: templateMatrix
    integer :: systemMatrix
    integer :: jacobianMatrix
    integer :: transportMatrix
    integer :: consistentMassMatrix
    integer :: lumpedMassMatrix
    integer :: coeffMatrix_CX
    integer :: coeffMatrix_CY
    integer :: coeffMatrix_CZ
    integer :: coeffMatrix_S
    integer :: convectionAFC
    integer :: diffusionAFC
    integer :: discretisation
    integer :: ijacobianFormat
    integer :: celement
    integer :: imatrixFormat

    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    type(t_triangulation) , pointer :: p_rtriangulation
    type(t_boundary) , pointer :: p_rboundary
    type(t_fparser), pointer :: p_rfparser

    integer :: i,j,nsumcubRefBilForm,nsumcubRefLinForm,nsumcubRefEval

    ! Retrieve application specific parameters from the parameter list
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'templatematrix', templateMatrix)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'systemmatrix', systemMatrix)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'jacobianmatrix', jacobianMatrix)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'transportmatrix', transportMatrix)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'consistentmassmatrix', consistentMassMatrix)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'lumpedmassmatrix', lumpedMassMatrix)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_s', coeffMatrix_S)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_cx', coeffMatrix_CX)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_cy', coeffMatrix_CY)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_cz', coeffMatrix_CZ)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'convectionAFC', convectionAFC)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'diffusionAFC', diffusionAFC)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'discretisation', discretisation)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'celement', celement)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'nsumcubRefBilForm', nsumcubRefBilForm, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'nsumcubRefLinForm', nsumcubRefLinForm, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'nsumcubRefEval', nsumcubRefEval, 0)

    ! Set pointers to triangulation and boundary structure
    p_rtriangulation  => rproblemLevel%rtriangulation
    p_rboundary       => rproblemLevel%p_rproblem%rboundary


    ! Create discretisation structure
    if (discretisation > 0) then

      ! Initialize the discretisation structure
      p_rdiscretisation => rproblemLevel%Rdiscretisation(discretisation)
      if (p_rdiscretisation%ndimension .eq. 0) then
        call spdiscr_initBlockDiscr(p_rdiscretisation, 1, p_rtriangulation)
      end if

      ! Get spatial dimension
      select case(p_rdiscretisation%ndimension)
      case (NDIM1D)
        select case(celement)
        case (-1,1,11)
          ! P1=Q1 finite elements
          call spdiscr_initDiscr_simple(&
              p_rdiscretisation%RspatialDiscr(1),&
              EL_E001_1D, SPDISC_CUB_AUTOMATIC,&
              p_rtriangulation, p_rboundary)

        case(-2,2,12)
          ! P2=Q2 finite elements
          call spdiscr_initDiscr_simple(&
              p_rdiscretisation%RspatialDiscr(1),&
              EL_E002_1D, SPDISC_CUB_AUTOMATIC,&
              p_rtriangulation, p_rboundary)

        case DEFAULT
          call output_line('Unsupported element type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'transp_initProblemLevel')
          call sys_halt()
        end select

      case (NDIM2D)
        select case(celement)
        case (1)
          ! P1 finite elements
          call spdiscr_initDiscr_simple(&
              p_rdiscretisation%RspatialDiscr(1),&
              EL_E001, SPDISC_CUB_AUTOMATIC,&
              p_rtriangulation, p_rboundary)

        case (2)
          ! P2 finite elements
          call spdiscr_initDiscr_simple(&
              p_rdiscretisation%RspatialDiscr(1),&
              EL_E002, SPDISC_CUB_AUTOMATIC,&
              p_rtriangulation, p_rboundary)

        case (11)
          ! Q1 finite elements
          call spdiscr_initDiscr_simple(&
              p_rdiscretisation%RspatialDiscr(1),&
              EL_E011, SPDISC_CUB_AUTOMATIC,&
              p_rtriangulation, p_rboundary)

        case (12)
          ! Q2 finite elements
          call spdiscr_initDiscr_simple(&
              p_rdiscretisation%RspatialDiscr(1),&
              EL_E013, SPDISC_CUB_AUTOMATIC,&
              p_rtriangulation, p_rboundary)

        case (-1)
          ! mixed P1/Q1 finite elements
          call spdiscr_initDiscr_triquad(&
              p_rdiscretisation%RspatialDiscr(1), EL_E001, EL_E011,&
              SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC,&
              p_rtriangulation, p_rboundary)

        case (-2)
          ! mixed P2/Q2 finite elements
          call spdiscr_initDiscr_triquad(&
              p_rdiscretisation%RspatialDiscr(1), EL_E002, EL_E013,&
              SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC,&
              p_rtriangulation, p_rboundary)

        case DEFAULT
          call output_line('Unsupproted element type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'transp_initProblemLevel')
          call sys_halt()
        end select

      case (NDIM3D)
        select case(celement)
        case (1)
          ! P1 finite elements
          call spdiscr_initDiscr_simple(&
              p_rdiscretisation%RspatialDiscr(1),&
              EL_E001_3D, SPDISC_CUB_AUTOMATIC,&
              p_rtriangulation, p_rboundary)

        case (11)
          ! Q1 finite elements
          call spdiscr_initDiscr_simple(&
              p_rdiscretisation%RspatialDiscr(1),&
              EL_E010_3D, SPDISC_CUB_AUTOMATIC,&
              p_rtriangulation, p_rboundary)

        case DEFAULT
          call output_line('Unsupproted element type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'transp_initProblemLevel')
          call sys_halt()
        end select

      case DEFAULT
        call output_line('Invalid number of spatial dimensions',&
            OU_CLASS_ERROR,OU_MODE_STD,'transp_initProblemLevel')
        call sys_halt()
      end select

      ! Enforce using summed cubature formula (if any)
      do i = 1, p_rdiscretisation%ncomponents
        do j = 1, p_rdiscretisation%RspatialDiscr(i)%inumFESpaces
          p_rdiscretisation%RspatialDiscr(i)%RelementDistr(j)%ccubTypeBilForm =&
              cub_getSummedCubType(p_rdiscretisation%RspatialDiscr(i)&
              %RelementDistr(j)%ccubTypeBilForm, nsumcubRefBilForm)
          p_rdiscretisation%RspatialDiscr(i)%RelementDistr(j)%ccubTypeLinForm =&
              cub_getSummedCubType(p_rdiscretisation%RspatialDiscr(i)&
              %RelementDistr(j)%ccubTypeLinForm, nsumcubRefLinForm)
          p_rdiscretisation%RspatialDiscr(i)%RelementDistr(j)%ccubTypeEval =&
              cub_getSummedCubType(p_rdiscretisation%RspatialDiscr(i)&
              %RelementDistr(j)%ccubTypeEval, nsumcubRefEval)
        end do
      end do

    end if


    ! If the template matrix has no structure data then generate the
    ! finite element matrix sparsity structure based on the spatial
    ! descretisation and store it as the template matrix. Otherwise we
    ! assume that the template matrix has been generated externally.
    if (.not.lsyssc_hasMatrixStructure(rproblemLevel%Rmatrix(templateMatrix))) then
      call parlst_getvalue_int(rparlist, ssectionName, 'imatrixFormat', imatrixFormat)
      call bilf_createMatrixStructure(p_rdiscretisation&
          %RspatialDiscr(1), imatrixFormat, rproblemLevel&
          %Rmatrix(templateMatrix))

    end if


    ! Create system matrix as duplicate of the template matrix
    if (systemMatrix > 0) then
      if (lsyssc_isMatrixStructureShared(rproblemLevel%Rmatrix(systemMatrix),&
                                         rproblemLevel%Rmatrix(templateMatrix))) then
        call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(systemMatrix),&
            rproblemLevel%Rmatrix(templateMatrix), .false., .false., .true.)

      else
        call lsyssc_duplicateMatrix(&
            rproblemLevel%Rmatrix(templateMatrix),&
            rproblemLevel%Rmatrix(systemMatrix),&
            LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

      end if
    end if


    ! Create transport matrix as duplicate of the template matrix
    if (transportMatrix > 0) then
      if (lsyssc_isMatrixStructureShared(rproblemLevel%Rmatrix(transportMatrix),&
                                         rproblemLevel%Rmatrix(templateMatrix))) then
        call lsyssc_resizeMatrix(&
            rproblemLevel%Rmatrix(transportMatrix),&
            rproblemLevel%Rmatrix(templateMatrix),&
            .false., .false., .true.)

      else
        call lsyssc_duplicateMatrix(&
            rproblemLevel%Rmatrix(templateMatrix),&
            rproblemLevel%Rmatrix(transportMatrix),&
            LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

      end if
    end if


    ! Create Jacobian matrix. This is a little bit tricky. If the
    ! Jacobian matrix has the same sparsity pattern as the template
    ! matrix, we can just create the Jacobian matrix as a duplicate of
    ! the template matrix. If the Jacobian matrix has an extended
    ! sparsity pattern we must create it by using the template matrix
    if (jacobianMatrix > 0) then

      ! What format do we have for the Jacobian matrix?
      call parlst_getvalue_int(rparlist, ssectionName, 'ijacobianFormat', ijacobianFormat)

      if (lsyssc_hasMatrixStructure(rproblemLevel%Rmatrix(jacobianMatrix))) then
        if (ijacobianFormat .eq. 0) then
          call lsyssc_resizeMatrix(&
              rproblemLevel%Rmatrix(jacobianMatrix),&
              rproblemLevel%Rmatrix(templateMatrix),&
              .false., .false., .true.)
        else
          call afcstab_generateExtSparsity(&
              rproblemLevel%Rmatrix(templateMatrix),&
              rproblemLevel%Rmatrix(jacobianMatrix))
        end if
      else
        if (ijacobianFormat .eq. 0) then
          call lsyssc_duplicateMatrix(&
              rproblemLevel%Rmatrix(templateMatrix),&
              rproblemLevel%Rmatrix(jacobianMatrix),&
              LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

        else
          call afcstab_generateExtSparsity(&
              rproblemLevel%Rmatrix(templateMatrix),&
              rproblemLevel%Rmatrix(jacobianMatrix))

        end if
      end if
    end if


    ! Create consistent (and lumped) mass matrix as duplicate of the template matrix
    if (consistentMassMatrix > 0) then
      if (lsyssc_isMatrixStructureShared(rproblemLevel%Rmatrix(consistentMassMatrix),&
                                         rproblemLevel%Rmatrix(templateMatrix))) then
        call lsyssc_resizeMatrix(&
            rproblemLevel%Rmatrix(consistentMassMatrix),&
            rproblemLevel%Rmatrix(templateMatrix),&
            .false., .false., .true.)

      else
        call lsyssc_duplicateMatrix(&
            rproblemLevel%Rmatrix(templateMatrix),&
            rproblemLevel%Rmatrix(consistentMassMatrix),&
            LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

      end if
      call stdop_assembleSimpleMatrix(&
          rproblemLevel%Rmatrix(consistentMassMatrix), DER_FUNC, DER_FUNC)


      if (lumpedMassMatrix > 0) then
        call lsyssc_duplicateMatrix(&
            rproblemLevel%Rmatrix(consistentMassMatrix),&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)

        call lsyssc_lumpMatrixScalar(&
            rproblemLevel%Rmatrix(lumpedMassMatrix), LSYSSC_LUMP_DIAG)

      end if
    elseif (lumpedMassMatrix > 0) then
      call lsyssc_duplicateMatrix(&
          rproblemLevel%Rmatrix(templateMatrix),&
          rproblemLevel%Rmatrix(lumpedMassMatrix),&
          LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

      call stdop_assembleSimpleMatrix(&
          rproblemLevel%Rmatrix(lumpedMassMatrix), DER_FUNC, DER_FUNC)

      call lsyssc_lumpMatrixScalar(&
          rproblemLevel%Rmatrix(lumpedMassMatrix), LSYSSC_LUMP_DIAG)

    end if


    ! Create diffusion matrix as duplicate of the template matrix
    if (coeffMatrix_S > 0) then
      if (lsyssc_isMatrixStructureShared(rproblemLevel%Rmatrix(coeffMatrix_S),&
                                         rproblemLevel%Rmatrix(templateMatrix))) then
        call lsyssc_resizeMatrix(&
            rproblemLevel%Rmatrix(coeffMatrix_S),&
            rproblemLevel%Rmatrix(templateMatrix),&
            .false., .false., .true.)

      else
        call lsyssc_duplicateMatrix(&
            rproblemLevel%Rmatrix(templateMatrix),&
            rproblemLevel%Rmatrix(coeffMatrix_S),&
            LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

      end if

      ! Get function parser from collection
      p_rfparser => collct_getvalue_pars(rcollection, 'rfparser')

      select case(p_rtriangulation%ndim)
      case (NDIM1D)
        call initDiffusionMatrix1D(p_rfparser,&
            rproblemLevel%Rmatrix(coeffMatrix_S))
      case (NDIM2D)
        call initDiffusionMatrix2D(p_rfparser,&
            rproblemLevel%Rmatrix(coeffMatrix_S))
      case (NDIM3D)
        call initDiffusionMatrix3D(p_rfparser,&
            rproblemLevel%Rmatrix(coeffMatrix_S))
      case DEFAULT
        call lsyssc_releaseMatrix(rproblemLevel%Rmatrix(coeffMatrix_S))
      end select
    end if


    ! Create coefficient matrix (phi, dphi/dx) duplicate of the template matrix
    if (coeffMatrix_CX > 0) then
      if (lsyssc_isMatrixStructureShared(rproblemLevel%Rmatrix(coeffMatrix_CX),&
                                         rproblemLevel%Rmatrix(templateMatrix))) then
        call lsyssc_resizeMatrix(&
            rproblemLevel%Rmatrix(coeffMatrix_CX),&
            rproblemLevel%Rmatrix(templateMatrix),&
            .false., .false., .true.)

      else
        call lsyssc_duplicateMatrix(&
            rproblemLevel%Rmatrix(templateMatrix),&
            rproblemLevel%Rmatrix(coeffMatrix_CX),&
            LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

      end if
      call stdop_assembleSimpleMatrix(&
          rproblemLevel%Rmatrix(coeffMatrix_CX), DER_DERIV3D_X, DER_FUNC)

    end if


    ! Create coefficient matrix (phi, dphi/dy) duplicate of the template matrix
    if (coeffMatrix_CY > 0) then
      if (lsyssc_isMatrixStructureShared(rproblemLevel%Rmatrix(coeffMatrix_CY),&
                                         rproblemLevel%Rmatrix(templateMatrix))) then
        call lsyssc_resizeMatrix(&
            rproblemLevel%Rmatrix(coeffMatrix_CY),&
            rproblemLevel%Rmatrix(templateMatrix),&
            .false., .false., .true.)

      else
        call lsyssc_duplicateMatrix(&
            rproblemLevel%Rmatrix(templateMatrix),&
            rproblemLevel%Rmatrix(coeffMatrix_CY),&
            LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

      end if
      call stdop_assembleSimpleMatrix(&
          rproblemLevel%Rmatrix(coeffMatrix_CY), DER_DERIV3D_Y, DER_FUNC)

    end if


    ! Create coefficient matrix (phi, dphi/dz) duplicate of the template matrix
    if (coeffMatrix_CZ > 0) then
      if (lsyssc_isMatrixStructureShared(rproblemLevel%Rmatrix(coeffMatrix_CZ),&
                                         rproblemLevel%Rmatrix(templateMatrix))) then
        call lsyssc_resizeMatrix(&
            rproblemLevel%Rmatrix(coeffMatrix_CZ),&
            rproblemLevel%Rmatrix(templateMatrix),&
            .false., .false., .true.)

      else
        call lsyssc_duplicateMatrix(&
            rproblemLevel%Rmatrix(templateMatrix),&
            rproblemLevel%Rmatrix(coeffMatrix_CZ),&
            LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

      end if
      call stdop_assembleSimpleMatrix(&
          rproblemLevel%Rmatrix(coeffMatrix_CZ), DER_DERIV3D_Z, DER_FUNC)

    end if


    ! Resize stabilisation structure if necessary and remove the
    ! indicator for the subdiagonal edge structure. If they are
    ! needed, then they are re-generated on-the-fly.
    if (convectionAFC > 0) then
      if (rproblemLevel%Rafcstab(convectionAFC)%iSpec .eq. AFCSTAB_UNDEFINED) then
        call gfsc_initStabilisation(&
            rproblemLevel%Rmatrix(templateMatrix),&
            rproblemLevel%Rafcstab(convectionAFC),&
            p_rdiscretisation)

      else
        call afcstab_resizeStabilisation(&
            rproblemLevel%Rafcstab(convectionAFC),&
            rproblemLevel%Rmatrix(templateMatrix))

        rproblemLevel%Rafcstab(convectionAFC)%iSpec =&
            iand(rproblemLevel%Rafcstab(convectionAFC)%iSpec,&
            not(AFCSTAB_HAS_OFFDIAGONALEDGES))
      end if
    end if


    ! The same applies to the diffusive stabilisation structure
    if (diffusionAFC > 0) then
      if (rproblemLevel%Rafcstab(diffusionAFC)%iSpec .eq. AFCSTAB_UNDEFINED) then
        call gfsc_initStabilisation(&
            rproblemLevel%Rmatrix(templateMatrix),&
            rproblemLevel%Rafcstab(diffusionAFC),&
            p_rdiscretisation)

      else
        call afcstab_resizeStabilisation(&
            rproblemLevel%Rafcstab(diffusionAFC),&
            rproblemLevel%Rmatrix(templateMatrix))

        rproblemLevel%Rafcstab(diffusionAFC)%iSpec =&
            iand(rproblemLevel%Rafcstab(diffusionAFC)%iSpec,&
            not(AFCSTAB_HAS_OFFDIAGONALEDGES))
      end if
    end if

    ! Set update notification in problem level structure
    rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                     PROBLEV_MSPEC_UPDATE)

  contains

    ! Here, the initialisation routines for the diffusion matrix follow.

    !**************************************************************

    subroutine initDiffusionMatrix1D(rfparser, rmatrix)

      type(t_fparser), intent(in) :: rfparser
      type(t_matrixScalar), intent(inout) :: rmatrix

      ! local variables
      character(LEN=SYS_STRLEN) :: sdiffusionName
      real(DP), dimension(1) :: Dunity = (/1.0_DP/)
      real(DP) :: dalpha
      integer :: icomp,idiffusiontype

      ! Retrieve data from parameter list
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'idiffusiontype', idiffusiontype)

      select case(idiffusiontype)
      case (DIFFUSION_ISOTROPIC,&
            DIFFUSION_ANISOTROPIC)

        ! Retrieve name/number of expression describing the diffusion coefficient
        call parlst_getvalue_string(rparlist, ssectionName, 'sdiffusionname',&
                                    sdiffusionName, isubString=1)
        icomp = fparser_getFunctionNumber(rfparser, sdiffusionName)

        ! Evaluate the constant coefficient from the function parser
        call fparser_evalFunction(rfparser, icomp, Dunity, dalpha)

        ! Assemble the Laplace matrix multiplied by the negative value
        ! of the physical diffusion parameter alpha
        call stdop_assembleLaplaceMatrix(rmatrix, .true., -dalpha)

      case DEFAULT
        call lsyssc_clearMatrix(rmatrix)
      end select

    end subroutine initDiffusionMatrix1D

    !**************************************************************

    subroutine initDiffusionMatrix2D(rfparser, rmatrix)

      type(t_fparser), intent(in) :: rfparser
      type(t_matrixScalar), intent(inout) :: rmatrix

      ! local variables
      type(t_bilinearform) :: rform
      character(LEN=SYS_STRLEN) :: sdiffusionName
      real(DP), dimension(1) :: Dunity = (/1.0_DP/)
      real(DP) :: dalpha
      integer :: i,icomp,idiffusiontype

      ! Retrieve data from parameter list
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'idiffusiontype', idiffusiontype)

      select case(idiffusiontype)
      case (DIFFUSION_ISOTROPIC)

        ! Retrieve name/number of expression describing the diffusion coefficient
        call parlst_getvalue_string(rparlist,&
            ssectionName, 'sdiffusionname', sdiffusionName, isubString=1)
        icomp = fparser_getFunctionNumber(rfparser, sdiffusionName)

        ! Evaluate the constant coefficient from the function parser
        call fparser_evalFunction(rfparser, icomp, Dunity, dalpha)

        ! Assemble the Laplace matrix multiplied by the negative value
        ! of the physical diffusion parameter alpha
        call stdop_assembleLaplaceMatrix(rmatrix, .true., -dalpha)

      case (DIFFUSION_ANISOTROPIC)

        do i = 1, 4
          ! Retrieve name/number of expression describing the diffusion coefficient
          call parlst_getvalue_string(rparlist,&
              ssectionName, 'sdiffusionname', sdiffusionName, isubString=i)
          icomp = fparser_getFunctionNumber(rfparser, sdiffusionName)

          ! Evaluate the constant coefficient from the function parser
          call fparser_evalFunction(rfparser, icomp, Dunity, rform%Dcoefficients(i))
        end do

        ! We have constant coefficients
        rform%ballCoeffConstant = .true.
        rform%BconstantCoeff    = .true.
        rform%Dcoefficients     = -1.0_DP

        ! Initialize the bilinear form
        rform%itermCount = 4
        rform%Idescriptors(1,1) = DER_DERIV2D_X
        rform%Idescriptors(2,1) = DER_DERIV2D_X

        rform%Idescriptors(1,2) = DER_DERIV2D_X
        rform%Idescriptors(2,2) = DER_DERIV2D_Y

        rform%Idescriptors(1,3) = DER_DERIV2D_Y
        rform%Idescriptors(2,3) = DER_DERIV2D_X

        rform%Idescriptors(1,4) = DER_DERIV2D_Y
        rform%Idescriptors(2,4) = DER_DERIV2D_Y

        ! Assemble the anisotropic diffusion matrix
        call bilf_buildMatrixScalar(rform, .true., rmatrix)

      case DEFAULT
        call lsyssc_clearMatrix(rmatrix)
      end select

    end subroutine initDiffusionMatrix2D

    !**************************************************************

    subroutine initDiffusionMatrix3D(rfparser, rmatrix)

      type(t_fparser), intent(in) :: rfparser
      type(t_matrixScalar), intent(inout) :: rmatrix

      ! local variables
      type(t_bilinearform) :: rform
      character(LEN=SYS_STRLEN) :: sdiffusionName
      real(DP), dimension(1) :: Dunity = (/1.0_DP/)
      real(DP) :: dalpha
      integer :: i,icomp,idiffusiontype

      ! Retrieve data from parameter list
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'idiffusiontype', idiffusiontype)

      select case(idiffusiontype)
      case (DIFFUSION_ISOTROPIC)

        ! Retrieve name/number of expression describing the diffusion coefficient
        call parlst_getvalue_string(rparlist,&
            ssectionName, 'sdiffusionname', sdiffusionName, isubString=1)
        icomp = fparser_getFunctionNumber(rfparser, sdiffusionName)

        ! Evaluate the constant coefficient from the function parser
        call fparser_evalFunction(rfparser, icomp, Dunity, dalpha)

        ! Assemble the Laplace matrix multiplied by the negative value
        ! of the physical diffusion parameter alpha
        call stdop_assembleLaplaceMatrix(rmatrix, .true., -dalpha)

      case (DIFFUSION_ANISOTROPIC)

        do i = 1, 9
          ! Retrieve name/number of expression describing the diffusion coefficient
          call parlst_getvalue_string(rparlist,&
              ssectionName, 'sdiffusionname', sdiffusionName, isubString=i)
          icomp = fparser_getFunctionNumber(rfparser, sdiffusionName)

          ! Evaluate the constant coefficient from the function parser
          call fparser_evalFunction(rfparser, icomp, Dunity, rform%Dcoefficients(i))
        end do

        ! We have constant coefficients
        rform%ballCoeffConstant = .true.
        rform%BconstantCoeff    = .true.
        rform%Dcoefficients     = -1.0_DP

        ! Initialize the bilinear form
        rform%itermCount = 9
        rform%Idescriptors(1,1) = DER_DERIV3D_X
        rform%Idescriptors(2,1) = DER_DERIV3D_X

        rform%Idescriptors(1,2) = DER_DERIV3D_X
        rform%Idescriptors(2,2) = DER_DERIV3D_Y

        rform%Idescriptors(1,3) = DER_DERIV3D_X
        rform%Idescriptors(2,3) = DER_DERIV3D_Z

        rform%Idescriptors(1,4) = DER_DERIV3D_Y
        rform%Idescriptors(2,4) = DER_DERIV3D_X

        rform%Idescriptors(1,5) = DER_DERIV3D_Y
        rform%Idescriptors(2,5) = DER_DERIV3D_Y

        rform%Idescriptors(1,6) = DER_DERIV3D_Y
        rform%Idescriptors(2,6) = DER_DERIV3D_Z

        rform%Idescriptors(1,7) = DER_DERIV3D_Z
        rform%Idescriptors(2,7) = DER_DERIV3D_X

        rform%Idescriptors(1,8) = DER_DERIV3D_Z
        rform%Idescriptors(2,8) = DER_DERIV3D_Y

        rform%Idescriptors(1,9) = DER_DERIV3D_Z
        rform%Idescriptors(2,9) = DER_DERIV3D_Z

        ! Assemble the anisotropic diffusion matrix
        call bilf_buildMatrixScalar(rform, .true., rmatrix)

      case DEFAULT
        call lsyssc_clearMatrix(rmatrix)
      end select

    end subroutine initDiffusionMatrix3D

  end subroutine transp_initProblemLevel

  !*****************************************************************************

!<subroutine>

  subroutine transp_initAllProblemLevels(rparlist, ssectionName,&
      rproblem, rcollection)

!<description>
    ! This subroutine initializes the all problem levels attached to
    ! the global problem structure. It generates the discretisation,
    ! the template matrix and the coefficient matrices as duplicates
    ! of the template matrix.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName
!</input>

!<inputoutput>
    ! problem structure
    type(t_problem), intent(inout) :: rproblem

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</intputoutput>
!</subroutine>

    ! pointer to the problem level
    type(t_problemLevel), pointer :: p_rproblemLevel


    ! loop over all problem levels
    p_rproblemLevel => rproblem%p_rproblemLevelMax
    do while(associated(p_rproblemLevel))

      ! Initialize individual problem level
      call transp_initProblemLevel(rparlist, ssectionName,&
          p_rproblemLevel, rcollection)

      ! Switch to next coarser level
      p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
    end do

  end subroutine transp_initAllProblemLevels

  !*****************************************************************************

!<subroutine>

  subroutine transp_initSolution(rparlist, ssectionName,&
    rproblemLevel, dtime, rvector, rcollection)

!<description>
    ! This subroutine initializes the solution vector
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName

    ! problem level
    type(t_problemLevel), intent(in), target :: rproblemLevel

    ! time for solution evaluation
    real(DP), intent(in) :: dtime
!</input>

!<inputoutput>
    ! solution vector
    type(t_vectorBlock), intent(inout) :: rvector

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</intputoutput>
!</subroutine>

    ! local variables
    type(t_pgm) :: rpgm
    type(t_afcstab) :: rafcstab
    type(t_linearForm) :: rform
    type(t_vectorBlock) :: rvectorHigh, rvectorAux
    type(t_matrixScalar), target :: rlumpedMassMatrix, rconsistentMassMatrix
    type(t_matrixScalar), pointer :: p_rlumpedMassMatrix, p_rConsistentMassMatrix
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
    type(t_fparser), pointer :: p_rfparser
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Ddata
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: depsAbsSolution, depsRelSolution, dnorm0, dnorm
    character(LEN=SYS_STRLEN) :: ssolutionname
    integer :: isolutiontype
    integer :: ieq, neq, ndim, icomp, iter
    integer :: lumpedMassMatrix, consistentMassMatrix, systemMatrix
    integer :: nmaxIterationsSolution
    
    ! Get global configuration from parameter list
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'isolutiontype', isolutiontype)

    ! How should the solution be initialised?
    select case(isolutionType)
    case (SOLUTION_ZERO)

      !-------------------------------------------------------------------------
      ! Initialize solution by zeros
      !-------------------------------------------------------------------------

      call lsysbl_clearVector(rvector)
      
      
    case (SOLUTION_GRAYMAP)
      
      !-------------------------------------------------------------------------
      ! Initialize the nodal values by the data of a graymap image
      !-------------------------------------------------------------------------

      ! Get global configuration from parameter list
      call parlst_getvalue_string(rparlist,&
          ssectionName, 'ssolutionname', ssolutionName)

      ! Initialize solution from portable graymap image
      call ppsol_readPGM(0, ssolutionName, rpgm)

      ! Set pointers
      call storage_getbase_double2D(rproblemLevel%rtriangulation&
          %h_DvertexCoords, p_DvertexCoords)
      call storage_getbase_double(rvector%h_Ddata, p_Ddata)

      ! Initialize the solution by the image data
      call ppsol_initArrayPGM_Dble(rpgm, p_DvertexCoords, p_Ddata)

      ! Release portable graymap image
      call ppsol_releasePGM(rpgm)


    case (SOLUTION_ANALYTIC_POINTVALUE)

      !-------------------------------------------------------------------------
      ! Initialize the nodal values by the data of an analytical expression
      !-------------------------------------------------------------------------
      
      ! Get global configuration from parameter list
      call parlst_getvalue_string(rparlist,&
          ssectionName, 'ssolutionname', ssolutionName)

      ! Get function parser from collection structure
      p_rfparser => collct_getvalue_pars(rcollection, 'rfparser')

      ! Set pointer to vertex coordinates
      call storage_getbase_double2D(&
          rproblemLevel%rtriangulation%h_DvertexCoords, p_DvertexCoords)

      ! Get number of spatial dimensions
      ndim = rproblemLevel%rtriangulation%ndim

      ! Initialize variable values
      Dvalue = 0.0_DP

      ! Set pointer to data array
      call lsyssc_getbase_double(rvector%RvectorBlock(1), p_Ddata)

      ! Get number of equations of scalar subvector
      neq = rvector%RvectorBlock(1)%NEQ

      ! Get the number of the component used for evaluating the initial solution
      icomp = fparser_getFunctionNumber(p_rfparser, ssolutionname)

      ! Loop over all equations of scalar subvector
      do ieq = 1, neq

        ! Set coordinates and evalution time
        Dvalue(1:ndim)   = p_DvertexCoords(:,ieq)
        Dvalue(NDIM3D+1) = dtime

        ! Evaluate the function parser
        call fparser_evalFunction(p_rfparser, icomp, Dvalue, p_Ddata(ieq))
      end do

      
    case (SOLUTION_ANALYTIC_L2_CONSISTENT,&
          SOLUTION_ANALYTIC_L2_LUMPED)

      !-------------------------------------------------------------------------
      ! Initialize the FE-function by the L2-projection of the analytical data
      !-------------------------------------------------------------------------

      ! Get global configuration from parameter list
      call parlst_getvalue_string(rparlist,&
          ssectionName, 'ssolutionname', ssolutionName)

      ! Get function parser from collection structure
      p_rfparser => collct_getvalue_pars(rcollection, 'rfparser')

      ! Retrieve the lumped and consistent mass matrices from the
      ! problem level structure or recompute them on-the-fly.
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'consistentMassMatrix', consistentMassMatrix)
      
      if (consistentMassMatrix .gt. 0) then
        p_rconsistentMassMatrix => rproblemLevel%Rmatrix(consistentMassMatrix)
      else
        call bilf_createMatrixStructure(&
            rvector%p_rblockDiscr%RspatialDiscr(1),&
            LSYSSC_MATRIX9, rconsistentMassMatrix)
        call stdop_assembleSimpleMatrix(&
            rconsistentMassMatrix, DER_FUNC, DER_FUNC, 1.0_DP, .true.)
        p_rconsistentMassMatrix => rconsistentMassMatrix
      end if
      
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'lumpedmassmatrix', lumpedMassMatrix)
      
      if (lumpedMassMatrix .gt. 0) then
        p_rlumpedMassMatrix => rproblemLevel%Rmatrix(lumpedMassMatrix)
      else
        call lsyssc_duplicateMatrix(p_rconsistentMassMatrix,&
            rlumpedMassMatrix, LSYSSC_DUP_TEMPLATE, LSYSSC_DUP_TEMPLATE)
        call lsyssc_lumpMatrixScalar(rlumpedMassMatrix, LSYSSC_LUMP_DIAG)
        p_rlumpedMassMatrix => rlumpedMassMatrix
      end if

      ! Set up the linear form
      rform%itermCount = 1
      rform%Idescriptors(1) = DER_FUNC
      
      ! Attach the simulation time and the name of the 
      ! function parser to the collection structure
      rcollection%DquickAccess(1) = dtime
      rcollection%SquickAccess(1) = "rfparser"

      ! Set pointer to spatial discretisation
      p_rspatialDiscr => rvector%RvectorBlock(1)%p_rspatialDiscr

      ! Get the number of the component used for evaluating the initial solution
      rcollection%IquickAccess(1) = fparser_getFunctionNumber(p_rfparser, ssolutionname)
      
      ! Assemble the linear form for the scalar subvector
      call linf_buildVectorScalar2(rform, .true.,&
          rvector%RvectorBlock(1), transp_coeffVectorAnalytic, rcollection)

      ! Store norm of load vector (if required)
      dnorm0 = lsyssc_vectorNorm(&
          rvector%RvectorBlock(1), LINALG_NORML2)

      ! Compute the lumped L2-projection
      call lsyssc_invertedDiagMatVec(p_rlumpedMassMatrix,&
          rvector%RvectorBlock(1), 1.0_DP, rvector%RvectorBlock(1))

      !-------------------------------------------------------------------------
      ! Restore contribution of the consistent mass matrix of the L2-projection
      !-------------------------------------------------------------------------
      if (isolutionType .eq. SOLUTION_ANALYTIC_L2_CONSISTENT) then

        ! Get configuration from parameter list
        call parlst_getvalue_double(rparlist,&
            ssectionName, 'depsAbsSolution', depsAbsSolution, 1e-6_DP)
        call parlst_getvalue_double(rparlist,&
            ssectionName, 'depsRelSolution', depsRelSolution, 1e-4_DP)
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'nmaxIterationsSolution', nmaxIterationsSolution, 100)
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'systemMatrix', systemMatrix)

        ! Compute auxiliary vectors for high-order solution and increment
        call lsysbl_duplicateVector(rvector, rvectorHigh,&
            LSYSSC_DUP_TEMPLATE, LSYSSC_DUP_COPY)
        call lsysbl_duplicateVector(rvector, rvectorAux,&
            LSYSSC_DUP_TEMPLATE, LSYSSC_DUP_EMPTY)
        
        ! Compute the consistent L2-projection by Richardson iteration
        richardson: do iter = 1, nmaxIterationsSolution
          ! Compute the increment for each scalar subvector
          call lsyssc_scalarMatVec(p_rconsistentMassMatrix,&
              rvectorHigh%RvectorBlock(1),&
              rvectorAux%RvectorBlock(1), 1.0_DP, 0.0_DP)
          call lsyssc_invertedDiagMatVec(p_rlumpedMassMatrix,&
              rvectorAux%RvectorBlock(1), 1.0_DP,&
              rvectorAux%RvectorBlock(1))
          call lsyssc_vectorLinearComb(rvector%RvectorBlock(1),&
              rvectorAux%RvectorBlock(1), 1.0_DP, -1.0_DP)
          
          ! Update the scalar subvector of thesolution
          call lsyssc_vectorLinearComb(rvectorAux%RvectorBlock(1),&
              rvectorHigh%RvectorBlock(1), 1.0_DP, 1.0_DP)
          
          ! Check for convergence
          dnorm = lsyssc_vectorNorm(&
              rvectorAux%RvectorBlock(1), LINALG_NORML2)
          if ((dnorm .le. depsAbsSolution) .or.&
              (dnorm .le. depsRelSolution*dnorm0)) exit richardson
        end do richardson
        
        ! Initialise stabilisation structure by hand
        rafcstab%iSpec= AFCSTAB_UNDEFINED
        rafcstab%bprelimiting = .false.
        rafcstab%ctypeAFCstabilisation = AFCSTAB_FEMFCT_MASS
        call gfsc_initStabilisation(rproblemLevel%Rmatrix(systemMatrix), rafcstab)
        call afcstab_generateVerticesAtEdge(rproblemLevel%Rmatrix(systemMatrix), rafcstab)

        ! Compute the raw antidiffusive mass fluxes
        call gfsc_buildFluxFCT(rafcstab, rvectorHigh, rvectorHigh,&
            0.0_DP, 0.0_DP, 1.0_DP, .true., p_rconsistentMassMatrix)

        ! Attach section name to collection structure
        rcollection%SquickAccess(1) = ssectionName
        
        ! Apply flux correction to solution profile
        call gfsc_buildConvectionVectorFCT(p_rlumpedMassMatrix,&
            rafcstab, rvector, 1._DP, .false.,&
            AFCSTAB_FCTALGO_STANDARD+AFCSTAB_FCTALGO_SCALEBYMASS, rvector)

        ! Release stabilisation structure
        call afcstab_releaseStabilisation(rafcstab)

        ! Release auxiliary vectors
        call lsysbl_releaseVector(rvectorHigh)
        call lsysbl_releaseVector(rvectorAux)
      end if

      ! Release temporal matrices (if any)
      call lsyssc_releaseMatrix(rconsistentMassMatrix)
      call lsyssc_releaseMatrix(rlumpedMassMatrix)

      
    case DEFAULT
      call output_line('Invalid type of solution profile!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'transp_initSolution')
      call sys_halt()
    end select

  end subroutine transp_initSolution

  !*****************************************************************************

!<subroutine>

  subroutine transp_initRHS(rparlist, ssectionName, rproblemLevel,&
      dtime, rvector, rcollection)

!<description>
    ! This subroutine initializes the right-hand side vector for the
    ! primal problem based on the data from the parameter list
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName

    ! problem level structure
    type(t_problemLevel), intent(in), target :: rproblemLevel

    ! time for right-hand side evaluation
    real(DP), intent(in) :: dtime
!</input>

!<intputoutput>
    ! right-hand side vector
    type(t_vectorBlock), intent(inout) :: rvector

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_linearForm) :: rform
    character(LEN=SYS_STRLEN) :: srhsname
    integer :: icomp, irhstype


    ! Get global configuration from parameter list
    call parlst_getvalue_int(rparlist, ssectionName,&
                             'irhstype', irhstype)

    ! How should the right-hand side be initialised?
    select case(irhstype)
    case (RHS_ZERO)
      ! Initialize right-hand side by zeros
      call lsysbl_clearVector(rvector)


    case (RHS_ANALYTIC)
      ! Get global configuration from parameter list
      call parlst_getvalue_string(rparlist, ssectionName,&
                                  'srhsname', srhsname)

      ! Get function parser from collection structure
      p_rfparser => collct_getvalue_pars(rcollection, 'rfparser')

      ! Get the number of the component used for evaluating the right-hand side
      icomp = fparser_getFunctionNumber(p_rfparser, srhsname)

      ! Prepare quick access arrays of the collection
      rcollection%DquickAccess(1) = dtime
      rcollection%IquickAccess(1) = icomp
      rcollection%SquickAccess(1) = 'rfparser'

      ! Set up the corresponding linear form
      rform%itermCount      = 1
      rform%Idescriptors(1) = DER_FUNC

      ! Build the discretised right-hand side vector
      call linf_buildVectorScalar2(rform, .true., rvector%RvectorBlock(1),&
                                   transp_coeffVectorAnalytic, rcollection)

    case DEFAULT
      call output_line('Invalid type of target functional!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'transp_initRHS')
      call sys_halt()
    end select

  end subroutine transp_initRHS

  !*****************************************************************************

!<subroutine>

  subroutine transp_initTargetFunc(rparlist, ssectionName,&
      rproblemLevel, dtime, rvector, rcollection)

!<description>
    ! This subroutine calculates the target functional which serves as
    ! right-hand side vector for the dual problem in the framework of
    ! goal-oriented error estimation.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName

    ! problem level structure
    type(t_problemLevel), intent(in), target :: rproblemLevel

    ! time for target function evaluation
    real(DP), intent(in) :: dtime
!</input>

!<intputoutput>
    ! target function vector
    type(t_vectorBlock), intent(inout) :: rvector

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</intputoutput>
!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_linearForm) :: rform
    character(LEN=SYS_STRLEN) :: stargetfuncname
    integer :: itargetfunctype, icomp

    ! Get global configuration from parameter list
    call parlst_getvalue_int(rparlist, ssectionName,&
        'itargetfunctype', itargetfunctype)

    ! How should the target functional be initialised?
    select case(itargetfunctype)
    case (TFUNC_ZERO,&
          TFUNC_SURFINTG)

      ! Initialize target functional by zeros. The contribution of
      ! the surfae integral comes in be weakly impose boundary conditions
      call lsysbl_clearVector(rvector)


    case (TFUNC_VOLINTG,&
          TFUNC_MIXINTG)

      ! Get function parser from collection structure
      p_rfparser => collct_getvalue_pars(rcollection, 'rfparser')

      ! Get global configuration from parameter list
      call parlst_getvalue_string(rparlist, ssectionName,&
          'stargetfuncname', stargetfuncname)

      ! Get the number of the component used for evaluating the target functional
      icomp = fparser_getFunctionNumber(p_rfparser, stargetfuncname)

      ! Prepare quick access arrays of the collection
      rcollection%DquickAccess(1) = dtime
      rcollection%IquickAccess(1) = icomp
      rcollection%SquickAccess(1) = 'rfparser'

      ! Set up the corresponding linear form
      rform%itermCount      = 1
      rform%Idescriptors(1) = DER_FUNC

      ! Build the discretised target functional. The contribution of
      ! the surfae integral comes in be weakly impose boundary conditions
      call linf_buildVectorScalar2(rform, .true., rvector&
          %RvectorBlock(1), transp_coeffVectorAnalytic, rcollection)


    case DEFAULT
      call output_line('Invalid type of target functional!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'transp_initTargetFunc')
      call sys_halt()
    end select

  end subroutine transp_initTargetFunc

  !*****************************************************************************

!<subroutine>

  subroutine transp_outputSolution(rparlist, ssectionName,&
    rproblemLevel, rsolutionPrimal, rsolutionDual, dtime)

!<description>
    ! This subroutine exports the solution vector to file in UCD format
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName

    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! OPTIONAL: solution vector for primal problem
    type(t_vectorBlock), intent(in), optional :: rsolutionPrimal

    ! OPTIONAL: solution vector for dual problem
    type(t_vectorBlock), intent(in), optional :: rsolutionDual

    ! OPTIONAL: simulation time
    real(DP), intent(in), optional :: dtime
!</input>
!</subroutine>

    ! section names
    character(LEN=SYS_STRLEN) :: soutputName
    character(LEN=SYS_STRLEN) :: sucdsolution

    ! persistent variable
    integer, save :: ifilenumber = 1

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata
    type(t_ucdExport) :: rexport
    integer :: iformatUCD


    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist, ssectionName,&
                                'output', soutputName)
    call parlst_getvalue_string(rparlist, trim(soutputName),&
                                'sucdsolution', sucdsolution)
    call parlst_getvalue_int(rparlist, trim(soutputName),&
                             'iformatucd', iformatUCD)

    ! Initialize the UCD exporter
    call flagship_initUCDexport(rproblemLevel, sucdsolution,&
                                iformatUCD, rexport, ifilenumber)

    ! Increase filenumber by one
    ifilenumber = ifilenumber+1

    ! Set simulation time
    if (present(dtime)) call ucd_setSimulationTime(rexport, dtime)

    ! Add primal solution vector
    if (present(rsolutionPrimal)) then
      call lsysbl_getbase_double(rsolutionPrimal, p_Ddata)
      call ucd_addVariableVertexBased (rexport, 'u',&
          UCD_VAR_STANDARD, p_Ddata)
    end if

    ! Add dual solution vector
    if (present(rsolutionDual)) then
      call lsysbl_getbase_double(rsolutionDual, p_Ddata)
      call ucd_addVariableVertexBased (rexport, 'z',&
          UCD_VAR_STANDARD, p_Ddata)
    end if

    ! Write UCD file
    call ucd_write  (rexport)
    call ucd_release(rexport)

  end subroutine transp_outputSolution

  !*****************************************************************************

!<subroutine>

  subroutine transp_outputStatistics(rtimerTotal, rcollection)

!<description>
    ! This subroutine output application statistics
!</description>

!<input>
    ! timer for total time measurement
    type(t_timer), intent(in) :: rtimerTotal
!</input>

!<inputoutput>
    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variable
    type(t_timer), pointer :: p_rtimerSolution
    type(t_timer), pointer :: p_rtimerAdaptation
    type(t_timer), pointer :: p_rtimerErrorEstimation
    type(t_timer), pointer :: p_rtimerTriangulation
    type(t_timer), pointer :: p_rtimerAssemblyCoeff
    type(t_timer), pointer :: p_rtimerAssemblyMatrix
    type(t_timer), pointer :: p_rtimerAssemblyVector
    type(t_timer), pointer :: p_rtimerPrePostprocess
    real(DP) :: dtotalTime, dfraction


    ! Get timer objects from collection
    p_rtimerSolution => collct_getvalue_timer(rcollection, 'rtimerSolution')
    p_rtimerAdaptation => collct_getvalue_timer(rcollection, 'rtimerAdaptation')
    p_rtimerErrorEstimation => collct_getvalue_timer(rcollection, 'rtimerErrorEstimation')
    p_rtimerTriangulation => collct_getvalue_timer(rcollection, 'rtimerTriangulation')
    p_rtimerAssemblyCoeff => collct_getvalue_timer(rcollection, 'rtimerAssemblyCoeff')
    p_rtimerAssemblyMatrix => collct_getvalue_timer(rcollection, 'rtimerAssemblyMatrix')
    p_rtimerAssemblyVector => collct_getvalue_timer(rcollection, 'rtimerAssemblyVector')
    p_rtimerPrePostprocess => collct_getvalue_timer(rcollection, 'rtimerPrePostprocess')

    ! Output statistics
    call output_lbrk()
    call output_line('Time measurement:')
    call output_line('-----------------')

    call stat_subTimers(p_rtimerAssemblyMatrix, p_rtimerSolution)
    call stat_subTimers(p_rtimerAssemblyVector, p_rtimerSolution)

    dtotalTime = max(rtimerTotal%delapsedCPU, rtimerTotal%delapsedReal)
    dfraction  = 100.0_DP/dtotalTime

    call output_line('Time for computing solution:   '//&
                     trim(adjustl(sys_sdE(p_rtimerSolution%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerSolution%delapsedCPU, 5)))//' %')
    call output_line('Time for mesh adaptivity:      '//&
                     trim(adjustl(sys_sdE(p_rtimerAdaptation%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerAdaptation%delapsedCPU, 5)))//' %')
    call output_line('Time for error estimation:     '//&
                     trim(adjustl(sys_sdE(p_rtimerErrorEstimation%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerErrorEstimation%delapsedCPU, 5)))//' %')
    call output_line('Time for triangulation:        '//&
                     trim(adjustl(sys_sdE(p_rtimerTriangulation%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerTriangulation%delapsedCPU, 5)))//' %')
    call output_line('Time for coefficient assembly: '//&
                     trim(adjustl(sys_sdE(p_rtimerAssemblyCoeff%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerAssemblyCoeff%delapsedCPU, 5)))//' %')
    call output_line('Time for matrix assembly:      '//&
                     trim(adjustl(sys_sdE(p_rtimerAssemblyMatrix%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerAssemblyMatrix%delapsedCPU, 5)))//' %')
    call output_line('Time for vector assembly:      '//&
                     trim(adjustl(sys_sdE(p_rtimerAssemblyVector%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerAssemblyVector%delapsedCPU, 5)))//' %')
    call output_line('Time for pre-/post-processing: '//&
                     trim(adjustl(sys_sdE(p_rtimerPrePostprocess%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerPrePostprocess%delapsedCPU, 5)))//' %')
    call output_lbrk()
    call output_line('Time for total simulation:     '//&
                     trim(adjustl(sys_sdE(dtotalTime, 5))))
    call output_lbrk()

  end subroutine transp_outputStatistics

  !*****************************************************************************

!<subroutine>

  subroutine transp_estimateTargetFuncError(rparlist, ssectionName,&
      rproblemLevel, rtimestep, rsolver, rsolutionPrimal,&
      rsolutionDual, rcollection, rtargetError, dtargetError, rrhs)

!<description>
    ! This subroutine estimates the error in the quantity of interest
!</description>

!<input>
    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName

    ! time-stepping algorithm
    type(t_timestep), intent(inout) :: rtimestep

    ! primal solution vector
    type(t_vectorBlock), intent(in), target :: rsolutionPrimal

    ! dual solution vector
    type(t_vectorBlock), intent(in) :: rsolutionDual

    ! OPTIONAL: right-hand side vector
    type(t_vectorBlock), intent(in), optional :: rrhs
!</input>

!<inputoutput>
    ! parameter list
    type(t_parlist), intent(inout) :: rparlist

    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(inout) :: rsolver

    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!<output>
    ! element-wise error distribution
    type(t_vectorScalar), intent(out) :: rtargetError

    ! global error in target qunatity
    real(DP), intent(out) :: dtargetError
!</output>
!</subroutine>

    ! section names
    character(LEN=SYS_STRLEN) :: serrorestimatorName
    character(LEN=SYS_STRLEN) :: sexacttargetfuncName
    character(LEN=SYS_STRLEN) :: sexactsolutionName
    character(LEN=SYS_STRLEN) :: stargetfuncName
    character(LEN=SYS_STRLEN) :: svelocityName

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_vectorBlock) :: rvector1, rvector2
    type(t_matrixScalar) :: rmatrix
    real(DP), dimension(:), pointer :: p_DsolutionDual, p_Dresidual
    real(DP), dimension(:), pointer :: p_DlumpedMassMatrix, p_DtargetError
    integer, dimension(:,:), pointer :: p_IverticesAtElement, p_IneighboursAtElement
    logical, dimension(:), pointer :: p_BisactiveElement
    real(DP) :: dexactTargetError, dexactTargetFunc, dprotectLayerTolerance
    real(DP) :: daux, dtargetFunc, dStep, theta
    integer :: i, icomp, convectionAFC, diffusionAFC, NEQ
    integer :: cconvectionStabilisation, cdiffusionStabilisation
    integer :: lumpedMassMatrix, templateMatrix, velocityfield
    integer :: itargetfunctype, iexactsolutiontype, imasstype
    integer :: iprotectLayer, nprotectLayers, igridindicator
    integer :: h_BisactiveElement


    !---------------------------------------------------------------------------
    ! Goal-oriented error estimation part 1: Galerkin orthogonality error
    !---------------------------------------------------------------------------

    ! Create vector for Galerkin residual
    call lsysbl_createVectorBlock(rsolutionPrimal, rvector1, .true., ST_DOUBLE)
    call lsysbl_createVectorBlock(rsolutionPrimal, rvector2, .true., ST_DOUBLE)

    ! Ok, this is a little bit tricky. We need to compute the residual
    ! vector for the steady-state Galerkin scheme for the zeroth
    ! iteration. To this end, we switch off all types of
    ! stabilisation, mass matrices, time stepping parameter, etc., and
    ! force the velocity field, and hence, the preconditioner to be
    ! updated. Then the steady-state residual vector is evaluated.

    ! Get parameters from parameter list
    call parlst_getvalue_int(rparlist, ssectionName, 'convectionAFC', convectionAFC)
    call parlst_getvalue_int(rparlist, ssectionName, 'diffusionAFC', diffusionAFC)
    call parlst_getvalue_int(rparlist, ssectionName, 'imasstype', imasstype)

    ! Set mass type to 'no mass matrix'
    call parlst_setvalue(rparlist, ssectionName, 'imasstype', '0')

    ! Set time-stepping parameters
    dStep = rtimestep%dStep; rtimestep%dStep = 1.0_DP
    theta = rtimestep%theta; rtimestep%theta = 1.0_DP

    ! Set stabilisation to standard Galerkin
    cconvectionStabilisation =&
        rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation
    rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation = AFCSTAB_GALERKIN

    cdiffusionStabilisation =&
        rproblemLevel%Rafcstab(diffusionAFC)%ctypeAFCstabilisation
    rproblemLevel%Rafcstab(diffusionAFC)%ctypeAFCstabilisation = AFCSTAB_GALERKIN

    ! Set update notification for velocity field/preconditioner
    rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec, PROBLEV_MSPEC_UPDATE)

    ! Calculate the standard Galerkin preconditioner
    ! (required for rhs and residual calculation)
    call transp_calcPrecondThetaScheme(rproblemLevel, rtimestep,&
        rsolver, rsolutionPrimal, rcollection)

    ! Calculate the standard Galerkin right-hand side vector
    call transp_calcRhsThetaScheme(rproblemLevel, rtimestep, rsolver,&
        rsolutionPrimal, rvector1, rcollection, rrhs)

    ! Calculate the standard Galerkin residual
    call transp_calcResidualThetaScheme(rproblemLevel, rtimestep,&
        rsolver, rsolutionPrimal, rsolutionPrimal, rvector1,&
        rvector2, 0, rcollection)

    ! Ok, now we have to switch on all types of stabilisation again
    rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation =&
        cconvectionStabilisation
    rproblemLevel%Rafcstab(diffusionAFC)%ctypeAFCstabilisation =&
        cdiffusionStabilisation

    ! ... and we reset the mass type
    call parlst_setvalue(rparlist, ssectionName, 'imasstype',&
        trim(sys_si(imasstype, 3)))

    ! ... and the time-stepping structure
    rtimestep%dStep = dStep
    rtimestep%theta = theta

    ! Again, set update notification for velocity field/preconditioner
    rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec, PROBLEV_MSPEC_UPDATE)


    ! We need the lumped mass matrix for scaling
    call parlst_getvalue_int(rparlist, ssectionName,&
        'lumpedMassMatrix', lumpedMassMatrix)

    if (lumpedMassMatrix > 0) then
      ! Set pointer to the lumped mass matrix
      call lsyssc_getbase_double(rproblemLevel&
          %Rmatrix(lumpedMassMatrix), p_DlumpedMassMatrix)
      NEQ = rproblemLevel%Rmatrix(lumpedMassMatrix)%NEQ

    else
      ! Compute the lumped mass matrix explicitly
      call parlst_getvalue_int(rparlist, ssectionName,&
          'templatematrix', templateMatrix)
      call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(templateMatrix),&
          rmatrix, LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      call stdop_assembleSimpleMatrix(rmatrix, DER_FUNC, DER_FUNC)
      call lsyssc_lumpMatrixScalar(rmatrix, LSYSSC_LUMP_DIAG)
      call lsyssc_getbase_double(rmatrix, p_DlumpedMassMatrix)
      NEQ = rmatrix%NEQ

    end if   ! lumpedMassMatrix > 0

    ! Set pointers
    call lsysbl_getbase_double(rvector2, p_Dresidual)
    call lsysbl_getbase_double(rsolutionDual, p_DsolutionDual)

    ! Now we compute the global error and its nodal contributions
    dtargetError = 0.0_DP
    do i = 1, NEQ
      dtargetError   = dtargetError + p_Dresidual(i)*p_DsolutionDual(i)
      p_Dresidual(i) = abs(p_Dresidual(i)*p_DsolutionDual(i))/p_DlumpedMassMatrix(i)
    end do
    dtargetError = abs(dtargetError)

    ! Compute the element-wise error.  We create a scalar vector and
    ! compute the local L1-norms of nodal error vector which yields
    ! the local values of the a posteriori error estimator.
    call lsyssc_createVector(rtargetError, rproblemLevel%rtriangulation%NEL, .false.)
    call pperr_scalar(rvector2%RvectorBlock(1), PPERR_L1ERROR, daux,&
        relementError=rtargetError)

    ! Release temporal vectors
    call lsysbl_releaseVector(rvector1)
    call lsysbl_releaseVector(rvector2)

    ! Release temporal mass matrix if required
    if (lumpedMassMatrix .le. 0) call lsyssc_releaseMatrix(rmatrix)


    ! Check if an exact solution is available or if there is an exact
    ! expression for the target functional. If neither is available,
    ! then we are unable to compute the effectivity index.
    call parlst_getvalue_int(rparlist, ssectionName,&
        'iexactsolutiontype', iexactsolutiontype)
    call parlst_getvalue_string(rparlist, ssectionName,&
        'sexacttargetfuncname', sexacttargetfuncname, '')

    if ((iexactsolutiontype .ne. SOLUTION_ANALYTIC_POINTVALUE) .and.&
        (trim(sexacttargetfuncname) .eq. '')) then

      call output_lbrk()
      call output_line('Error Analysis')
      call output_line('--------------')
      call output_line('estimated error in quantity of interest: '//trim(sys_sdEP(dtargetError,15,6)))
      call output_line('exact error in quantity of interest:     '//'n.a.')
      call output_line('effectivity index:                       '//'n.a.')
      call output_line('relative effectivity index:              '//'n.a.')
      call output_lbrk()

    else

      !-------------------------------------------------------------------------
      ! Compute the effectivity index
      !-------------------------------------------------------------------------

      ! Get global configuration from parameter list
      call parlst_getvalue_int(rparlist, ssectionName,&
                               'itargetfunctype', itargetfuncType)
      call parlst_getvalue_string(rparlist, ssectionName,&
                                  'sexactsolutionname', sexactsolutionname, '')

      select case(itargetfunctype)
      case (TFUNC_ZERO)
        call output_line('Zero target functional is not implemented yet!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'transp_estimateTargetFuncError')
        call sys_halt()


      case(TFUNC_VOLINTG)
        ! Get function parser from collection structure
        p_rfparser => collct_getvalue_pars(rcollection, 'rfparser')

        ! Get global configuration from parameter list
        call parlst_getvalue_string(rparlist, ssectionName,&
                                  'stargetfuncname', stargetfuncName)

        ! Prepare quick access arrays of the collection
        rcollection%DquickAccess(1) = rtimestep%dTime
        rcollection%IquickAccess(1) = fparser_getFunctionNumber(p_rfparser, sexactsolutionName)
        rcollection%IquickAccess(2) = fparser_getFunctionNumber(p_rfparser, stargetfuncName)
        rcollection%SquickAccess(1) = 'rfparser'

        if (trim(sexacttargetfuncname) .ne. '') then
          ! Get the number of the component used for evaluating the exact target functional
          icomp = fparser_getFunctionNumber(p_rfparser, sexacttargetfuncName)

          ! Evaluate exact value of the quantity of interest
          call fparser_evalFunction(p_rfparser, icomp,&
                                   (/rtimestep%dTime/), dexactTargetFunc)

          ! Compute the approximate value of the quantity of interest
          call pperr_scalar(PPERR_MEANERROR, dtargetFunc,&
                            rsolutionPrimal%RvectorBlock(1),&
                            rcollection=rcollection,&
                            ffunctionWeight=transp_weightFuncAnalytic)

          ! Compute exact error in target functional. Note that the
          ! routine pperr_scalar computes the value $J(0-u_h)$ so that we
          ! have to change the sign to minus
          dexactTargetError = dexactTargetFunc+dtargetFunc

        else

          ! Compute the exact error of the quantity of interest
          call pperr_scalar(PPERR_MEANERROR, dexactTargetError,&
              rsolutionPrimal%RvectorBlock(1), transp_refFuncAnalytic,&
              rcollection, ffunctionWeight=transp_weightFuncAnalytic)

          ! Compute the exact value of the quantity of interest
          call pperr_scalar(PPERR_MEANERROR, dexactTargetFunc,&
              ffunctionReference=transp_refFuncAnalytic, rcollection&
              =rcollection, rdiscretisation=rsolutionPrimal&
              %RvectorBlock(1)%p_rspatialdiscr, ffunctionWeight&
              =transp_weightFuncAnalytic)
        end if


      case (TFUNC_SURFINTG)
        ! Get function parser from collection structure
        p_rfparser => collct_getvalue_pars(rcollection, 'rfparser')

        ! Get global configuration from parameter list
        call parlst_getvalue_int(rparlist, ssectionName,&
            'velocityfield', velocityfield)
        call parlst_getvalue_string(rparlist, ssectionName,&
            'stargetfuncname', stargetfuncName)

        ! Prepare quick access arrays of the collection
        rcollection%DquickAccess(1) = rtimestep%dTime
        rcollection%IquickAccess(1) = fparser_getFunctionNumber(p_rfparser, sexactsolutionName)
        rcollection%IquickAccess(2) = fparser_getFunctionNumber(p_rfparser, stargetfuncName)
        rcollection%SquickAccess(1) = 'rfparser'

        do i = 1, rproblemLevel%rtriangulation%ndim
          call parlst_getvalue_string(rparlist, ssectionName,&
              'svelocityname', svelocityname, isubString=i)
          rcollection%IquickAccess(i+2) = fparser_getFunctionNumber(p_rfparser, svelocityname)
        end do

        ! Attach vectors to collection structure
        rcollection%p_rvectorQuickAccess1 => rsolutionPrimal
        rcollection%p_rvectorQuickAccess2 => rproblemLevel%RvectorBlock(velocityfield)

        if (trim(sexacttargetfuncname) .ne. '') then
          ! Get the number of the component used for evaluating the exact target functional
          icomp = fparser_getFunctionNumber(p_rfparser, sexacttargetfuncName)

          ! Evaluate exact value of the quantity of interest
          call fparser_evalFunction(p_rfparser, icomp,&
              (/rtimestep%dTime/), dexactTargetFunc)

          ! Prepare quick access arrays of the collection
          rcollection%IquickAccess(1) = fparser_getFunctionNumber(p_rfparser, '@null')

          ! Compute the approximate value of the quantity of interest
          call pperr_scalarBoundary2D(0, CUB_G3_1D, dtargetFunc,&
              ffunctionReference=transp_errorBdrInt2D_sim, rcollection&
              =rcollection, rdiscretisation=rsolutionPrimal&
              %RvectorBlock(1)%p_rspatialdiscr, ffunctionWeight&
              =transp_weightFuncBdrInt2D_sim)

          ! Compute exact error in target functional. Note that the
          ! routine pperr_scalar computes the value $J(0-u_h)$ so that we
          ! have to change the sign to minus
          dexactTargetError = dexactTargetFunc+dtargetFunc

        else

          ! Compute the exact error of the quantity of interest
          call pperr_scalarBoundary2D(0, CUB_G3_1D, dexactTargetError,&
              ffunctionReference=transp_errorBdrInt2D_sim, rcollection&
              =rcollection, rdiscretisation=rsolutionPrimal&
              %RvectorBlock(1)%p_rspatialdiscr, ffunctionWeight&
              =transp_weightFuncBdrInt2D_sim)

          ! Compute the exact value of the quantity of interest
          call pperr_scalarBoundary2D(0, CUB_G3_1D, dexactTargetFunc,&
              ffunctionReference=transp_refFuncBdrInt2D_sim, rcollection&
              =rcollection, rdiscretisation=rsolutionPrimal&
              %RvectorBlock(1)%p_rspatialdiscr, ffunctionWeight&
              =transp_weightFuncBdrInt2D_sim)
        end if


      case (TFUNC_MIXINTG)
        ! Get function parser from collection structure
        p_rfparser => collct_getvalue_pars(rcollection, 'rfparser')

        ! Get global configuration from parameter list
        call parlst_getvalue_int(rparlist, ssectionName,&
            'velocityfield', velocityfield)

        ! Get the name of the function used for evaluating the
        ! volume integral part of the target functional
        if (parlst_querysubstrings(rparlist, ssectionName, 'stargetfuncname') .eq. 0) then
          call parlst_getvalue_string(rparlist, ssectionName,&
              'stargetfuncname', stargetfuncname)
        else
          call parlst_getvalue_string(rparlist, ssectionName,&
              'stargetfuncname', stargetfuncname, isubstring=1)
        end if

        ! Prepare quick access arrays of the collection
        rcollection%DquickAccess(1) = rtimestep%dTime
        rcollection%IquickAccess(1) = fparser_getFunctionNumber(p_rfparser, sexactsolutionName)
        rcollection%IquickAccess(2) = fparser_getFunctionNumber(p_rfparser, stargetfuncName)
        rcollection%SquickAccess(1) = 'rfparser'

        do i = 1, rproblemLevel%rtriangulation%ndim
          call parlst_getvalue_string(rparlist, ssectionName,&
              'svelocityname', svelocityname, isubString=i)
          rcollection%IquickAccess(i+2) = fparser_getFunctionNumber(p_rfparser, svelocityname)
        end do

        ! Attach vectors to collection structure
        rcollection%p_rvectorQuickAccess1 => rsolutionPrimal
        rcollection%p_rvectorQuickAccess2 => rproblemLevel%RvectorBlock(velocityfield)

        if (trim(sexacttargetfuncname) .ne. '') then
          ! Get the number of the component used for evaluating the exact target functional
          icomp = fparser_getFunctionNumber(p_rfparser, sexacttargetfuncName)

          ! Evaluate exact value of the quantity of interest
          call fparser_evalFunction(p_rfparser, icomp,&
                                   (/rtimestep%dTime/), dexactTargetFunc)

          ! Prepare quick access arrays of the collection
          rcollection%IquickAccess(1) = fparser_getFunctionNumber(p_rfparser, '@null')

          ! Compute the approximate value of the quantity of interest
          call pperr_scalar(PPERR_MEANERROR, dtargetFunc,&
              rsolutionPrimal%RvectorBlock(1), rcollection&
              =rcollection, ffunctionWeight=transp_weightFuncAnalytic)

          ! Get the name of the function used for evaluating the
          ! surface integral part of the target functional
          if (parlst_querysubstrings(rparlist, ssectionName, 'stargetfuncname') .eq. 0) then
            call parlst_getvalue_string(rparlist, ssectionName,&
                'stargetfuncname', stargetfuncname)
          else
            call parlst_getvalue_string(rparlist, ssectionName,&
                'stargetfuncname', stargetfuncname, isubstring=2)
          end if

          ! Prepare quick access arrays of the collection
          rcollection%IquickAccess(2) = fparser_getFunctionNumber(p_rfparser, stargetfuncName)

          ! Compute the approximate value of the quantity of interest
          call pperr_scalarBoundary2D(0, CUB_G3_1D, daux,&
              ffunctionReference=transp_errorBdrInt2D_sim, rcollection&
              =rcollection, rdiscretisation=rsolutionPrimal&
              %RvectorBlock(1)%p_rspatialdiscr, ffunctionWeight&
              =transp_weightFuncBdrInt2D_sim)

          ! Add boundary contribution
          dtargetFunc = dtargetFunc + daux

          ! Compute exact error in target functional. Note that the
          ! routine pperr_scalar computes the value $J(0-u_h)$ so that we
          ! have to change the sign to minus
          dexactTargetError = dexactTargetFunc+dtargetFunc

        else

          ! Compute the exact error of the quantity of interest
          call pperr_scalar(PPERR_MEANERROR, dexactTargetError,&
              rsolutionPrimal%RvectorBlock(1), transp_refFuncAnalytic&
              , rcollection, ffunctionWeight=transp_weightFuncAnalytic)

          ! Compute the exact value of the quantity of interest.
          call pperr_scalar(PPERR_MEANERROR, dexactTargetFunc,&
              ffunctionReference=transp_refFuncAnalytic, rcollection&
              =rcollection, rdiscretisation=rsolutionPrimal&
              %RvectorBlock(1)%p_rspatialdiscr, ffunctionWeight&
              =transp_weightFuncAnalytic)

          ! Get the name of the function used for evaluating the
          ! surface integral part of the target functional
          if (parlst_querysubstrings(rparlist, ssectionName, 'stargetfuncname') .eq. 0) then
            call parlst_getvalue_string(rparlist, ssectionName,&
                'stargetfuncname', stargetfuncname)
          else
            call parlst_getvalue_string(rparlist, ssectionName,&
                'stargetfuncname', stargetfuncname, isubstring=2)
          end if

          ! Prepare quick access arrays of the collection
          rcollection%DquickAccess(1) = rtimestep%dTime
          rcollection%IquickAccess(2) = fparser_getFunctionNumber(p_rfparser, stargetfuncName)

          ! Compute the exact error of the quantity of interest at the boundary
          call pperr_scalarBoundary2D(0, CUB_G3_1D, daux,&
              ffunctionReference=transp_errorBdrInt2D_sim, rcollection&
              =rcollection, rdiscretisation=rsolutionPrimal&
              %RvectorBlock(1)%p_rspatialdiscr, ffunctionWeight&
              =transp_weightFuncBdrInt2D_sim)

          ! Add boundary contribution
          dexactTargetError = dexactTargetError + daux

          ! Compute the exact value of the quantity of interest at the boundary
          call pperr_scalarBoundary2D(0, CUB_G3_1D, daux,&
              ffunctionReference=transp_refFuncBdrInt2D_sim, rcollection&
              =rcollection, rdiscretisation=rsolutionPrimal&
              %RvectorBlock(1)%p_rspatialdiscr, ffunctionWeight&
              =transp_weightFuncBdrInt2D_sim)

          ! Add boundary contribution
          dexactTargetFunc = dexactTargetFunc + daux
        end if


      case default
        call output_line('Invalid type of target functional!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'transp_estimateTargetFuncError')
        call sys_halt()
      end select

      call output_lbrk()
      call output_line('Error Analysis')
      call output_line('--------------')
      call output_line('estimated error in quantity of interest: '//trim(sys_sdEP(dtargetError,15,6)))
      call output_line('exact error in quantity of interest:     '//trim(sys_sdEP(abs(dexactTargetError),15,6)))
      call output_line('effectivity index:                       '//trim(sys_sdEP(dtargetError/abs(dexactTargetError),15,6)))
      call output_line('relative effectivity index:              '//trim(sys_sdEP(abs( (dtargetError-abs(dexactTargetError)) /&
                                                                                        dexactTargetFunc ),15,6)))
      call output_lbrk()

    end if


    !---------------------------------------------------------------------------
    ! Apply the adaptation strategy
    !---------------------------------------------------------------------------

    call parlst_getvalue_string(rparlist, ssectionName, 'errorestimator', serrorestimatorName)
    call parlst_getvalue_int(rparlist, trim(serrorestimatorName), 'igridindicator', igridindicator)

    ! What type of grid indicator are we?
    select case(igridIndicator)

    case (ERREST_ASIS)
      ! That is simple, do nothing.

    case DEFAULT
      call output_line('Invalid type of grid indicator!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'transp_estimateTargetFuncError')
      call sys_halt()
    end select


    !---------------------------------------------------------------------------
    ! Calculate protection layers
    !---------------------------------------------------------------------------

    call parlst_getvalue_string(rparlist, ssectionName, 'errorestimator', serrorestimatorName)
    call parlst_getvalue_int(rparlist, trim(serrorestimatorName), 'nprotectLayers', nprotectLayers, 0)
    call parlst_getvalue_double(rparlist, trim(serrorestimatorName),&
                                'dprotectLayerTolerance', dprotectLayerTolerance, 0.0_DP)

    if (nprotectLayers > 0) then

      ! Create auxiliary memory
      h_BisactiveElement = ST_NOHANDLE
      call storage_new('transp_estimateTargetFuncError',' BisactiveElement',&
                       rproblemLevel%rtriangulation%NEL, ST_LOGICAL,&
                       h_BisactiveElement, ST_NEWBLOCK_NOINIT)
      call storage_getbase_logical(h_BisactiveElement, p_BisactiveElement)

      ! Set pointers
      call storage_getbase_int2D(&
          rproblemLevel%rtriangulation%h_IneighboursAtElement, p_IneighboursAtElement)
      call storage_getbase_int2D(&
          rproblemLevel%rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
      call lsyssc_getbase_double(rtargetError, p_DtargetError)

      ! Compute protection layers
      do iprotectLayer = 1, nprotectLayers

        ! Reset activation flag
        p_BisActiveElement = .false.

        ! Compute a single-width protection layer
        call doProtectionLayerUniform(p_IverticesAtElement,&
            p_IneighboursAtElement, rproblemLevel%rtriangulation%NEL,&
            dprotectLayerTolerance, p_DtargetError,&
            p_BisActiveElement)
      end do

      ! Release memory
      call storage_free(h_BisactiveElement)

    end if

  contains

    !**************************************************************
    ! Compute one uniformly distributed protection layer

    subroutine doProtectionLayerUniform(IverticesAtElement,&
        IneighboursAtElement, NEL, dthreshold, Ddata,&
        BisactiveElement)

      integer, dimension(:,:), intent(in) :: IverticesAtElement
      integer, dimension(:,:), intent(in) :: IneighboursAtElement
      real(DP), intent(in) :: dthreshold
      integer, intent(in) :: NEL

      real(DP), dimension(:), intent(inout) :: Ddata
      logical, dimension(:), intent(inout) :: BisactiveElement


      ! local variables
      integer :: iel,jel,ive

      ! Loop over all elements in triangulation
      do iel = 1, NEL

        ! Do nothing if element belongs to active layer
        if (BisactiveElement(iel)) cycle

        ! Do nothing if element indicator does not exceed threshold
        if (Ddata(iel) .le. dthreshold) cycle

        ! Loop over neighbouring elements
        do ive = 1, tria_getNVE(IverticesAtElement, iel)

          ! Get number of neighbouring element
          jel = IneighboursAtElement(ive, iel)

          ! Do nothing at the boundary
          if (jel .eq. 0) cycle

          ! Check if element belongs to active layer
          if (BisactiveElement(jel)) then
            ! If yes, then just update the element indicator
            Ddata(jel) = max(Ddata(jel), Ddata(iel))
          else
            ! Otherwise, we have to check if the neighbouring element
            ! exceeds the prescribed threshold level. If this is the case
            ! it will be processed later or has already been processed
            if (Ddata(jel) .lt. dthreshold) then
              Ddata(jel) = max(Ddata(jel), Ddata(iel))
              BisactiveElement(jel) = .true.
            end if
          end if
        end do
      end do
    end subroutine doProtectionLayerUniform

  end subroutine transp_estimateTargetFuncError

  !*****************************************************************************

!<subroutine>

  subroutine transp_estimateRecoveryError(rparlist, ssectionName,&
      rproblemLevel, rsolution, dtime, rerror, derror, rcollection)

!<description>
    ! This subroutine estimates the error of the discrete solution by
    ! using recovery procedures such as the superconvergent patch
    ! recovery technique or L2-projection. If an exact solution value
    ! is avaialable, it computes the effectivity index of the error
    ! estimator. Moreover it applies the specified strategy to convert
    ! the error estimator into a grid indicator for mesh adaptation.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName

    ! solution vector
    type(t_vectorBlock), intent(in) :: rsolution

    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! simulation time
    real(DP), intent(in) :: dtime
!</input>

!<inputoutput>
    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!<output>
    ! element-wise error distribution
    type(t_vectorScalar), intent(out) :: rerror

    ! global error
    real(DP), intent(out) :: derror
!</output>
!</subroutine>

    ! section names
    character(LEN=SYS_STRLEN) :: serrorestimatorName
    character(LEN=SYS_STRLEN) :: sexactsolutionName

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_vectorScalar) :: rvectorScalar
    real(DP), dimension(:), pointer :: p_Ddata, p_Derror
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:,:), pointer :: p_IneighboursAtElement
    logical, dimension(:), pointer :: p_BisactiveElement
    real(DP) :: dnoiseFilter, dabsFilter, dsolution, dvalue,&
                dexacterror, dprotectLayerTolerance
    integer :: i, ierrorEstimator, igridindicator, iexactsolutiontype, icomp
    integer :: iprotectLayer, nprotectLayers
    integer :: h_BisactiveElement


    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist, ssectionName, 'errorestimator', serrorestimatorName)
    call parlst_getvalue_string(rparlist, ssectionName, 'sexactsolutionname', sexactsolutionName, '')
    call parlst_getvalue_int(rparlist, ssectionName, 'iexactsolutiontype', iexactsolutiontype, 0)
    call parlst_getvalue_int(rparlist, trim(serrorestimatorName), 'ierrorestimator', ierrorestimator)
    call parlst_getvalue_int(rparlist, trim(serrorestimatorName), 'igridindicator', igridindicator)
    call parlst_getvalue_int(rparlist, trim(serrorestimatorName), 'nprotectLayers', nprotectLayers, 0)
    call parlst_getvalue_double(rparlist, trim(serrorestimatorName),&
                                'dprotectLayerTolerance', dprotectLayerTolerance, 0.0_DP)


    !---------------------------------------------------------------------------
    ! Perform recovery-based error estimation
    !---------------------------------------------------------------------------

    ! What type of error estimator are we?
    select case(ierrorEstimator)

    case (ERREST_L2PROJECTION)
      call lsyssc_createVector(rerror, rproblemLevel%rtriangulation%NEL, .false.)
      call ppgrd_calcGradientError(rsolution%RvectorBlock(1), derror,&
                                   PPGRD_INTERPOL, 0, rerror)

    case (ERREST_SPR_VERTEX)
      call lsyssc_createVector(rerror, rproblemLevel%rtriangulation%NEL, .false.)
      call ppgrd_calcGradientError(rsolution%RvectorBlock(1), derror,&
                                   PPGRD_ZZTECHNIQUE, PPGRD_NODEPATCH, rerror)

    case (ERREST_SPR_ELEMENT)
      call lsyssc_createVector(rerror, rproblemLevel%rtriangulation%NEL, .false.)
      call ppgrd_calcGradientError(rsolution%RvectorBlock(1), derror,&
                                   PPGRD_ZZTECHNIQUE, PPGRD_ELEMPATCH, rerror)

    case (ERREST_SPR_FACE)
      call lsyssc_createVector(rerror, rproblemLevel%rtriangulation%NEL, .false.)
      call ppgrd_calcGradientError(rsolution%RvectorBlock(1), derror,&
                                   PPGRD_ZZTECHNIQUE, PPGRD_FACEPATCH, rerror)

    case (ERREST_LIMAVR)
      call lsyssc_createVector(rerror, rproblemLevel%rtriangulation%NEL, .false.)
      call ppgrd_calcGradientError(rsolution%RvectorBlock(1), derror,&
                                   PPGRD_LATECHNIQUE, 0, rerror)

    case (ERREST_SECONDDIFF)
      call parlst_getvalue_double(rparlist,&
          trim(serrorestimatorName), 'dnoiseFilter', dnoiseFilter)
      call parlst_getvalue_double(rparlist,&
          trim(serrorestimatorName), 'dabsFilter', dabsFilter)
      call ppind_secondDifference(rsolution%RvectorBlock(1),&
          dnoiseFilter, dabsFilter, rerror)

      derror = 1.0

    case DEFAULT
      call output_line('Invalid type of error estimator!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'transp_estimateRecoveryError')
      call sys_halt()
    end select


    !---------------------------------------------------------------------------
    ! Compute the effectivity index
    !---------------------------------------------------------------------------

    select case(iexactsolutiontype)
    case (SOLUTION_ANALYTIC_POINTVALUE)

      ! Get function parser from collection
      p_rfparser => collct_getvalue_pars(rcollection, 'rfparser')

      ! Get the number of the component used for evaluating the target functional
      icomp = fparser_getFunctionNumber(p_rfparser, sexactsolutionName)

      ! Prepare quick access arrays of the collection
      rcollection%DquickAccess(1) = dtime
      rcollection%IquickAccess(1) = icomp
      rcollection%SquickAccess(1) = 'rfparser'

      ! Calculate the H1-error of the reference solution
      call pperr_scalar(PPERR_H1ERROR, dexacterror, rsolution%RvectorBlock(1),&
                        transp_refFuncAnalytic, rcollection)

      call output_lbrk()
      call output_line('Error Analysis')
      call output_line('--------------')
      call output_line('estimated H1-error: '//trim(sys_sdEP(derror,15,6)))
      call output_line('exact H1-error:     '//trim(sys_sdEP(dexacterror,15,6)))
      call output_line('effectivity index:  '//trim(sys_sdEP(derror/dexacterror,15,6)))
      call output_lbrk()

    case DEFAULT
      call output_lbrk()
      call output_line('Error Analysis')
      call output_line('--------------')
      call output_line('estimated H1-error: '//trim(sys_sdEP(derror,15,6)))
      call output_lbrk()

    end select


    !---------------------------------------------------------------------------
    ! Apply the adaptation strategy
    !---------------------------------------------------------------------------

    ! What type of grid indicator are we?
    select case(igridIndicator)

    case (ERREST_ASIS)
      ! That is simple, do nothing.


    case (ERREST_EQUIDIST)
      ! Equidistribute the relative percentage error
      call output_lbrk()
      call output_line('Equidistribution strategy')
      call output_line('-------------------------')

      ! We need the global norm of the scalar error variable
      call pperr_scalar(PPERR_H1ERROR, dsolution,&
          rsolution %RvectorBlock(1))

      call output_line('Total percentage error:       '//trim(sys_sdEP(derror/dsolution,15,6)))

      ! Compute global permissible element error
      dvalue = sqrt(dsolution**2 + derror**2)/sqrt(real(rerror%NEQ, DP))

      call output_line('Permissible percentage error: '//trim(sys_sdEP(dvalue,15,6))//' x TOL')
      call output_lbrk()

      ! Scale element error by global permissible error
      call lsyssc_scaleVector(rerror, 1.0_DP/dvalue)


    case (-ERREST_EQUIDIST)
      ! Equidistribute the relative percentage error
      call output_lbrk()
      call output_line('Equidistribution strategy')
      call output_line('-------------------------')

      ! We need the global norm of the scalar error variable
      call lsyssc_createVector(rvectorScalar, rerror%NEQ, .true.)
      call pperr_scalar(PPERR_H1ERROR, dsolution, rsolution&
          %RvectorBlock(1), relementError=rvectorScalar)

      call output_line('Total percentage error: '//trim(sys_sdEP(derror/dsolution,15,6)))
      call output_lbrk()

      ! Compute local permissible element error for each elements
      call lsyssc_getbase_double(rerror, p_Derror)
      call lsyssc_getbase_double(rvectorScalar, p_Ddata)
      dvalue = dsolution/sqrt(real(rerror%NEQ, DP))

      do i = 1, size(p_Derror)
        p_Derror(i) = p_Derror(i)/(0.5*(sqrt(p_Ddata(i)) + dvalue))
      end do


    case (ERREST_LOGEQUIDIST)
      ! Equidistribute the logarithmic error
      call output_lbrk()
      call output_line('Logarithmic equidistribution strategy')
      call output_line('-------------------------------------')

      ! Set pointer
      call lsyssc_getbase_double(rerror, p_Ddata)

      ! Determine largest error value
      dvalue = -SYS_MAXREAL
      do i = 1, size(p_Ddata)
        dvalue = max(dvalue, p_Ddata(i))
      end do

      call output_line('Maximum error: '//trim(sys_sdEP(derror/dvalue,15,6)))

      ! Normalize error by largest value
      call lsyssc_scaleVector(rerror, 1.0_DP/dvalue)

      ! Initialize mean value
      dvalue = 0.0_DP

      ! Loop over all contributions
      do i = 1, size(p_Ddata)
        p_Ddata(i) = log(max(exp(-20.0_DP), p_Ddata(i)))
        dvalue   = dvalue + p_Ddata(i)
      end do

      ! Calculate mean
      dvalue = dvalue/real(size(p_Ddata), DP)

      call output_line('Mean value:    '//trim(sys_sdEP(derror/dvalue,15,6)))
      call output_lbrk()

      ! Subtract mean value from grid indicator
      do i = 1, size(p_Ddata)
        p_Ddata(i) = p_Ddata(i)-dvalue
      end do


    case (ERREST_AUTORMS)
      ! Automatic treshold for RMS
      call output_lbrk()
      call output_line('Automatic treshold for RMS')
      call output_line('--------------------------')

      ! Set pointer
      call lsyssc_getbase_double(rerror, p_Ddata)

      ! Initialize mean value
      dvalue = 0.0_DP

      ! Loop over all  contributions
      do i = 1, size(p_Ddata)
        dvalue = dvalue + p_Ddata(i)**2
      end do

      ! Calculate root mean value
      dvalue = sqrt(dvalue/real(size(p_Ddata), DP))

      call output_line('RMS value: '//trim(sys_sdEP(derror/dvalue,15,6)))
      call output_lbrk()

      ! Normalize grid indicator by RMS
      call lsyssc_scaleVector(rerror, 1.0_DP/dvalue)


    case DEFAULT
      call output_line('Invalid type of grid indicator!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'transp_estimateRecoveryError')
      call sys_halt()
    end select


    !---------------------------------------------------------------------------
    ! Calculate protection layers
    !---------------------------------------------------------------------------

    if (nprotectLayers > 0) then

      ! Create auxiliary memory
      h_BisactiveElement = ST_NOHANDLE
      call storage_new('transp_estimateRecoveryError',' BisactiveEleme' // &
          'nt', rproblemLevel%rtriangulation%NEL, ST_LOGICAL,&
          h_BisactiveElement, ST_NEWBLOCK_NOINIT)
      call storage_getbase_logical(h_BisactiveElement, p_BisactiveElement)

      ! Set pointers
      call storage_getbase_int2D(rproblemLevel%rtriangulation&
          %h_IneighboursAtElement, p_IneighboursAtElement)
      call storage_getbase_int2D(rproblemLevel%rtriangulation&
          %h_IverticesAtElement, p_IverticesAtElement)
      call lsyssc_getbase_double(rerror, p_Ddata)

      ! Compute protection layers
      do iprotectLayer = 1, nprotectLayers

        ! Reset activation flag
        p_BisActiveElement = .false.

        ! Compute a single-width protection layer
        call doProtectionLayerUniform(p_IverticesAtElement,&
            p_IneighboursAtElement, rproblemLevel%rtriangulation%NEL,&
            dprotectLayerTolerance, p_Ddata, p_BisActiveElement)
      end do

      ! Release memory
      call storage_free(h_BisactiveElement)

    end if

  contains

    ! Here, the real working routines follow.

    !**************************************************************
    ! Compute one uniformly distributed protection layer

    subroutine doProtectionLayerUniform(IverticesAtElement,&
        IneighboursAtElement, NEL, dthreshold, Ddata,&
        BisactiveElement)

      integer, dimension(:,:), intent(in) :: IverticesAtElement
      integer, dimension(:,:), intent(in) :: IneighboursAtElement
      real(DP), intent(in) :: dthreshold
      integer, intent(in) :: NEL

      real(DP), dimension(:), intent(inout) :: Ddata
      logical, dimension(:), intent(inout) :: BisactiveElement


      ! local variables
      integer :: iel,jel,ive

      ! Loop over all elements in triangulation
      do iel = 1, NEL

        ! Do nothing if element belongs to active layer
        if (BisactiveElement(iel)) cycle

        ! Do nothing if element indicator does not exceed threshold
        if (Ddata(iel) .le. dthreshold) cycle

        ! Loop over neighbouring elements
        do ive = 1, tria_getNVE(IverticesAtElement, iel)

          ! Get number of neighbouring element
          jel = IneighboursAtElement(ive, iel)

          ! Do nothing at the boundary
          if (jel .eq. 0) cycle

          ! Check if element belongs to active layer
          if (BisactiveElement(jel)) then
            ! If yes, then just update the element indicator
            Ddata(jel) = max(Ddata(jel), Ddata(iel))
          else
            ! Otherwise, we have to check if the neighbouring element
            ! exceeds the prescribed threshold level. If this is the case
            ! it will be processed later or has already been processed
            if (Ddata(jel) .lt. dthreshold) then
              Ddata(jel) = max(Ddata(jel), Ddata(iel))
              BisactiveElement(jel) = .true.
            end if
          end if
        end do
      end do
    end subroutine doProtectionLayerUniform

  end subroutine transp_estimateRecoveryError

  !*****************************************************************************

!<subroutine>

  subroutine transp_adaptTriangulation(rhadapt, rtriangulationSrc,&
      rindicator, rcollection, rtriangulationDest)

!<description>
    ! This subroutine performs h-adaptation for the given triangulation
!</description>

!<inputoutput>
    ! adaptation structure
    type(t_hadapt), intent(inout) :: rhadapt

    ! source triangulation structure
    type(t_triangulation), intent(inout), target :: rtriangulationSrc

    ! element-wise indicator
    type(t_vectorScalar), intent(inout) :: rindicator

    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!<output>
    ! OPTIONAL: destination triangulation structure
    ! If it is not given, the source triangulation is updated
    type(t_triangulation), intent(out), optional, target :: rtriangulationDest
!</output>
!</subroutine>


    ! local variables
    type(t_triangulation), pointer :: p_rtriangulation


    ! Check if adaptation structure has been prepared
    if (rhadapt%iSpec .eq. HADAPT_HAS_PARAMETERS) then

      ! Initialize adaptation structure from triangulation
      call hadapt_initFromTriangulation(rhadapt, rtriangulationSrc)

    else

      ! Refresh adaptation structure
      call hadapt_refreshAdaptation(rhadapt, rtriangulationSrc)

    end if

    ! How many spatial dimensions do we have?
    select case(rtriangulationSrc%ndim)
    case (NDIM1D)
      call transp_hadaptCallback1D(&
          HADAPT_OPR_INITCALLBACK, rcollection)
      call hadapt_performAdaptation(rhadapt, rindicator,&
          rcollection, transp_hadaptCallback1D)

    case (NDIM2D)
      call transp_hadaptCallback2D(&
          HADAPT_OPR_INITCALLBACK, rcollection)
      call hadapt_performAdaptation(rhadapt, rindicator,&
          rcollection, transp_hadaptCallback2D)

    case (NDIM3D)
      call transp_hadaptCallback3D(&
          HADAPT_OPR_INITCALLBACK, rcollection)
      call hadapt_performAdaptation(rhadapt, rindicator,&
          rcollection, transp_hadaptCallback3D)
    end select

    ! Do we have a destination triangulation or should the source
    ! triangulation structure be updated?
    if (present(rtriangulationDest)) then
      p_rtriangulation => rtriangulationDest
    else
      p_rtriangulation => rtriangulationSrc
    end if

    ! Generate raw mesh from adaptation structure
    call hadapt_generateRawMesh(rhadapt, p_rtriangulation)

  end subroutine transp_adaptTriangulation

  !*****************************************************************************

!<subroutine>

    subroutine transp_solveTransientPrimal(rparlist, ssectionName,&
        rbdrCond, rproblem, rtimestep, rsolver, rsolution,&
        rcollection)

!<description>
      ! This subroutine solves the transient primal flow problem
      !
      !  $$ \frac{\partial u}{\partial t}+\nabla\cdot{\bf f}(u)+r(u)=b $$
      !
      ! for the scalar quantity $u$ in the domain $\Omega$. Here,
      ! ${\bf f}(u)$ denotes the flux term and r(u) and b are the
      ! reacitve term and the right-hand side vector, respectively.
!</description>

!<input>
      ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName

    ! boundary condition structure
    type(t_boundaryCondition), intent(in) :: rbdrCond
!</input>

!<inputoutput>
    ! parameter list
    type(t_parlist), intent(inout) :: rparlist

    ! problem structure
    type(t_problem), intent(inout) :: rproblem

    ! time-stepping structure
    type(t_timestep), intent(inout) :: rtimestep

    ! solver struchture
    type(t_solver), intent(inout), target :: rsolver

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!<output>
    ! primal solution vector
    type(t_vectorBlock), intent(out), target :: rsolution
!</subroutine>

    ! Pointer to the multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! Pointer to the discretisation structure
    type(t_blockDiscretisation), pointer :: p_rdiscretisation

    ! Vector for the right-hand side
    type(t_vectorBlock) :: rrhs

    ! Vector for the element-wise error distribution
    type(t_vectorScalar) :: relementError

    ! Structure for h-adaptation
    type(t_hadapt) :: rhadapt

    ! Structure for the sparsity pattern
    type(t_graph) :: rgraph

    ! Timer structures
    type(t_timer), pointer :: p_rtimerPrePostprocess
    type(t_timer), pointer :: p_rtimerSolution
    type(t_timer), pointer :: p_rtimerErrorEstimation
    type(t_timer), pointer :: p_rtimerAdaptation
    type(t_timer), pointer :: p_rtimerTriangulation
    type(t_timer), pointer :: p_rtimerAssemblyCoeff

    ! section names
    character(LEN=SYS_STRLEN) :: sadaptivityName
    character(LEN=SYS_STRLEN) :: soutputName

    ! local variables
    real(dp) :: derror, dstepUCD, dtimeUCD, dstepAdapt, dtimeAdapt
    integer :: templateMatrix, systemMatrix, discretisation
    integer :: nlmin, ipreadapt, npreadapt, irhstype, ivelocitytype
    integer, external :: signal_SIGINT


    ! Get timer structures
    p_rtimerPrePostprocess => collct_getvalue_timer(rcollection, 'rtimerPrePostprocess')
    p_rtimerSolution => collct_getvalue_timer(rcollection, 'rtimerSolution')
    p_rtimerErrorEstimation => collct_getvalue_timer(rcollection, 'rtimerErrorEstimation')
    p_rtimerAdaptation => collct_getvalue_timer(rcollection, 'rtimerAdaptation')
    p_rtimerTriangulation => collct_getvalue_timer(rcollection, 'rtimerTriangulation')
    p_rtimerAssemblyCoeff => collct_getvalue_timer(rcollection, 'rtimerAssemblyCoeff')

    ! Start time measurement for pre-processing
    call stat_startTimer(p_rtimerPrePostprocess, STAT_TIMERSHORT)

    ! Get global parameters
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'irhstype', irhstype)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'ivelocitytype', ivelocitytype)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'discretisation', discretisation)

    ! Set pointer to maximum problem level and discretisation
    p_rproblemLevel   => rproblem%p_rproblemLevelMax
    p_rdiscretisation => p_rproblemLevel%Rdiscretisation(discretisation)

    ! Create the solution vector
    call lsysbl_createVectorBlock(p_rdiscretisation,&
        rsolution, .false., ST_DOUBLE)

    ! Initialize the solution vector and impose boundary conditions
    call transp_initSolution(rparlist, ssectionName, p_rproblemLevel,&
        rtimestep%dinitialTime, rsolution, rcollection)
    call bdrf_filterVectorExplicit(rbdrCond, rsolution,&
        rtimestep%dinitialTime)

    ! Initialize timer for intermediate UCD exporter
    dtimeUCD = rtimestep%dinitialTime
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'output', soutputName)
    call parlst_getvalue_double(rparlist,&
        trim(soutputName), 'dstepUCD', dstepUCD, 0.0_DP)


    !---------------------------------------------------------------------------
    ! Initialize the h-adaptation structure and perform pre-adaptation
    !---------------------------------------------------------------------------

    call parlst_getvalue_string(rparlist,&
        ssectionName, 'adaptivity', sadaptivityName, '')
    if (trim(adjustl(sadaptivityName)) .ne. '') then

      call parlst_getvalue_double(rparlist,&
          trim(sadaptivityName), 'dstepAdapt', dstepAdapt)
      call parlst_getvalue_double(rparlist,&
          trim(sadaptivityName), 'dtimeAdapt', dtimeAdapt)
      call parlst_getvalue_int(rparlist,&
          trim(sadaptivityName), 'npreadapt', npreadapt)


      if ((dstepAdapt > 0.0_DP) .or. (npreadapt > 0)) then

        ! Initialize adaptation structure from parameter list
        call hadapt_initFromParameterlist(rhadapt, rparlist, sadaptivityName)

        ! Generate a dynamic graph for the sparsity pattern and attach
        ! it to the collection structure which is used to pass this
        ! type to the callback function for h-adaptation
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'templateMatrix', templateMatrix)
        call grph_createGraphFromMatrix(p_rproblemLevel&
            %Rmatrix(templateMatrix), rgraph)
        call collct_setvalue_graph(rcollection, 'sparsitypattern',&
            rgraph, .true.)

        ! Perform pre-adaptation?
        if (npreadapt > 0) then

          ! Perform number of pre-adaptation steps
          do ipreadapt = 1, npreadapt

            ! Compute the error estimator using recovery techniques
            call transp_estimateRecoveryError(rparlist, ssectionname,&
                p_rproblemLevel, rsolution, rtimestep%dinitialTime,&
                relementError, derror, rcollection)

            ! Set the name of the template matrix
            rcollection%SquickAccess(1) = 'sparsitypattern'

            ! Attach the primal solution vector to the collection structure
            rcollection%p_rvectorQuickAccess1 => rsolution

            ! Perform h-adaptation and update the triangulation structure
            call transp_adaptTriangulation(rhadapt, p_rproblemLevel&
                %rtriangulation, relementError, rcollection)

            ! Release element-wise error distribution
            call lsyssc_releaseVector(relementError)

            ! Generate standard mesh from raw mesh
            call tria_initStandardMeshFromRaw(&
                p_rproblemLevel%rtriangulation, rproblem%rboundary)

            ! Update the template matrix according to the sparsity pattern
            call grph_generateMatrix(rgraph, p_rproblemLevel&
                %Rmatrix(templateMatrix))

            ! Resize the solution vector accordingly
            call lsysbl_resizeVectorBlock(rsolution, &
                p_rproblemLevel%Rmatrix(templateMatrix)%NEQ, .false.)

            ! Re-generate the initial solution vector and impose
            !  boundary conditions explicitly
            call transp_initSolution(rparlist, ssectionname,&
                p_rproblemLevel, rtimestep%dinitialTime, rsolution,&
                rcollection)
            call bdrf_filterVectorExplicit(rbdrCond, rsolution,&
                rtimestep%dinitialTime)

            ! Re-initialize all constant coefficient matrices
            call transp_initProblemLevel(rparlist, ssectionName,&
                p_rproblemLevel, rcollection)
          end do

          ! Prepare internal data arrays of the solver structure
          call parlst_getvalue_int(rparlist,&
              ssectionName, 'systemMatrix', systemMatrix)
          call flagship_updateSolverMatrix(p_rproblemLevel, rsolver,&
              systemMatrix, SYSTEM_INTERLEAVEFORMAT, UPDMAT_ALL)
          call solver_updateStructure(rsolver)

        end if   ! npreadapt > 0

      end if   ! dstepAdapt > 0

    else

      dstepAdapt = 0.0_DP

    end if

    ! Initialize right-hand side vector
    if (irhstype > 0) then
      call lsysbl_createVectorBlock(rsolution, rrhs)
      call transp_initRHS(rparlist, ssectionName, p_rproblemLevel,&
          rtimestep%dinitialTime, rrhs, rcollection)
    end if

    ! Calculate the initial velocity field
    nlmin = solver_getMinimumMultigridlevel(rsolver)
    call transp_calcVelocityField(rparlist, ssectionName,&
        p_rproblemLevel, rtimestep%dinitialTime, rcollection, nlmin)

    ! Attach the boundary condition
    call solver_setBoundaryCondition(rsolver, rbdrCond, .true.)

    ! Set primal problem mode
    call parlst_addvalue(rparlist, ssectionName, 'mode', 'primal')

    ! Stop time measurement for pre-processing
    call stat_stopTimer(p_rtimerPrePostprocess)


    !---------------------------------------------------------------------------
    ! Infinite time stepping loop
    !---------------------------------------------------------------------------

    timeloop: do

      ! Check for user interaction
      if (signal_SIGINT(-1) > 0 )&
          call transp_outputSolution(rparlist, ssectionName,&
          p_rproblemLevel, rsolution, dtime=rtimestep%dTime)

      !-------------------------------------------------------------------------
      ! Advance solution in time
      !-------------------------------------------------------------------------

      ! Start time measurement for solution procedure
      call stat_startTimer(p_rtimerSolution, STAT_TIMERSHORT)

      ! Prepare quick access arrays of the collection
      rcollection%SquickAccess(1) = ssectionName

      ! What time-stepping scheme should be used?
      select case(rtimestep%ctimestepType)

      case (TSTEP_RK_SCHEME)

        if (irhstype > 0) then
          ! Explicit Runge-Kutta scheme with non-zero right-hand side vector
          call tstep_performRKStep(p_rproblemLevel, rtimestep,&
              rsolver, rsolution, transp_nlsolverCallback,&
              rcollection, rrhs)
        else
          ! Explicit Runge-Kutta scheme without right-hand side vector
          call tstep_performRKStep(p_rproblemLevel, rtimestep,&
              rsolver, rsolution, transp_nlsolverCallback,&
              rcollection)
        end if

      case (TSTEP_THETA_SCHEME)

        if (irhstype > 0) then
          ! Two-level theta-scheme with non-zero right-hand side vector
          call tstep_performThetaStep(p_rproblemLevel, rtimestep,&
              rsolver, rsolution, transp_nlsolverCallback,&
              rcollection, rrhs)
        else
          ! Two-level theta-scheme without right-hand side vector
          call tstep_performThetaStep(p_rproblemLevel, rtimestep,&
              rsolver, rsolution, transp_nlsolverCallback,&
              rcollection)
        end if

      case DEFAULT
        call output_line('Unsupported time-stepping algorithm!',&
            OU_CLASS_ERROR,OU_MODE_STD,'transp_solveTransientPrimal')
        call sys_halt()
      end select

      ! Perform linearised FEM-FCT post-processing
      call transp_calcLinearisedFCT(rbdrCond, p_rproblemLevel,&
          rtimestep, rsolver, rsolution, rcollection)

      ! Stop time measurement for solution procedure
      call stat_stopTimer(p_rtimerSolution)

      ! Reached final time, then exit the infinite time loop?
      if (rtimestep%dTime .ge. rtimestep%dfinalTime) exit timeloop


      !-------------------------------------------------------------------------
      ! Post-process intermediate solution
      !-------------------------------------------------------------------------

      if ((dstepUCD .gt. 0.0_DP) .and. (rtimestep%dTime .ge. dtimeUCD)) then

        ! Set time for next intermediate solution export
        dtimeUCD = dtimeUCD + dstepUCD

        ! Start time measurement for post-processing
        call stat_startTimer(p_rtimerPrepostProcess, STAT_TIMERSHORT)

        ! Export the intermediate solution
        call transp_outputSolution(rparlist, ssectionName,&
            p_rproblemLevel, rsolution, dtime=rtimestep%dTime)

        ! Stop time measurement for post-processing
        call stat_stopTimer(p_rtimerPrepostProcess)

      end if


      !-------------------------------------------------------------------------
      ! Perform mesh adaptation
      !-------------------------------------------------------------------------

      if ((dstepAdapt .gt. 0.0_DP) .and. (rtimestep%dTime .ge. dtimeAdapt)) then

        ! Set time for next adaptation step
        dtimeAdapt = dtimeAdapt + dstepAdapt

        !-----------------------------------------------------------------------
        ! Perform recovery-based error estimation
        !-----------------------------------------------------------------------

        ! Start time measurement for error estimation
        call stat_startTimer(p_rtimerErrorEstimation, STAT_TIMERSHORT)

        ! Compute the error estimator using recovery techniques
        call transp_estimateRecoveryError(rparlist, ssectionname,&
            p_rproblemLevel, rsolution, rtimestep%dTime,&
            relementError, derror, rcollection)

        ! Stop time measurement for error estimation
        call stat_stopTimer(p_rtimerErrorEstimation)


        !-------------------------------------------------------------------------
        ! Perform h-adaptation
        !-------------------------------------------------------------------------

        ! Start time measurement for mesh adaptation
        call stat_startTimer(p_rtimerAdaptation, STAT_TIMERSHORT)

        ! Set the name of the template matrix
        rcollection%SquickAccess(1) = 'sparsitypattern'

        ! Attach the primal solution vector to the collection structure
        rcollection%p_rvectorQuickAccess1 => rsolution

        ! Perform h-adaptation and update the triangulation structure
        call transp_adaptTriangulation(rhadapt, p_rproblemLevel&
            %rtriangulation, relementError, rcollection)

        ! Release element-wise error distribution
        call lsyssc_releaseVector(relementError)

        ! Update the template matrix according to the sparsity pattern
        call grph_generateMatrix(rgraph, p_rproblemLevel&
            %Rmatrix(templateMatrix))

        ! Resize the solution vector accordingly
        call lsysbl_resizeVectorBlock(rsolution, &
            p_rproblemLevel%Rmatrix(templateMatrix)%NEQ, .false.)

        ! Stop time measurement for mesh adaptation
        call stat_stopTimer(p_rtimerAdaptation)


        !-------------------------------------------------------------------------
        ! Re-generate the discretisation and coefficient matrices
        !-------------------------------------------------------------------------

        ! Start time measurement for generation of the triangulation
        call stat_startTimer(p_rtimerTriangulation, STAT_TIMERSHORT)

        ! Generate standard mesh from raw mesh
        call tria_initStandardMeshFromRaw(&
            p_rproblemLevel%rtriangulation, rproblem%rboundary)

        ! Stop time measurement for generation of the triangulation
        call stat_stopTimer(p_rtimerTriangulation)


        ! Start time measurement for generation of constant coefficient matrices
        call stat_startTimer(p_rtimerAssemblyCoeff, STAT_TIMERSHORT)

        ! Re-initialize all constant coefficient matrices
        call transp_initProblemLevel(rparlist, ssectionName,&
            p_rproblemLevel, rcollection)

        ! Prepare internal data arrays of the solver structure
        call parlst_getvalue_int(rparlist, ssectionName, 'systemmatrix', systemMatrix)
        call flagship_updateSolverMatrix(p_rproblemLevel, rsolver,&
            systemMatrix, SYSTEM_INTERLEAVEFORMAT, UPDMAT_ALL)
        call solver_updateStructure(rsolver)

        ! Re-calculate the velocity field
        nlmin = solver_getMinimumMultigridlevel(rsolver)
        call transp_calcVelocityField(rparlist, ssectionName,&
            p_rproblemLevel, rtimestep%dTime, rcollection, nlmin)

        ! Re-initialize the right-hand side vector
        if (irhstype > 0) then
          call lsysbl_resizeVectorBlock(rrhs,&
              p_rproblemLevel%Rmatrix(templateMatrix)%NEQ, .false.)
          call transp_initRHS(rparlist, ssectionName, p_rproblemLevel&
              , rtimestep%dinitialTime, rrhs, rcollection)
        end if

        ! Stop time measurement for generation of constant coefficient matrices
        call stat_stopTimer(p_rtimerAssemblyCoeff)

      elseif(abs(ivelocitytype) .eq. VELOCITY_TIMEDEP) then

        ! Re-calculate the time-dependent velocity field
        nlmin = solver_getMinimumMultigridlevel(rsolver)
        call transp_calcVelocityField(rparlist, ssectionName,&
            p_rproblemLevel, rtimestep%dTime, rcollection, nlmin)

      end if

    end do timeloop


    ! Release adaptation structure
    if (trim(adjustl(sadaptivityName)) .ne. '') then
      if ((dstepAdapt > 0.0_DP) .or. (npreadapt > 0)) then
        call hadapt_releaseAdaptation(rhadapt)
        call grph_releaseGraph(rgraph)
      end if
    end if

    ! Release right-hand side
    if (irhstype > 0) then
      call lsysbl_releaseVector(rrhs)
    end if

  end subroutine transp_solveTransientPrimal

  !*****************************************************************************

!<subroutine>

  subroutine transp_solveTransientPrimalDual(rparlist, ssectionName,&
      rbdrCondPrimal, rbdrCondDual, rproblem, rtimestep, rsolver,&
      rsolutionPrimal, rsolutionDual, rcollection)

!<description>
    ! This subroutine solves the transient primal flow problem
    !
    ! $$\frac{\partial u}{\partial \tau}+\nabla\cdot{\bf f}(u)=s(u)$$
    !
    ! for a scalar quantity $u$ in the domain $\Omega$ both for the
    ! primal and the dual formulation.
!</description>

!<input>
    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName

    ! boundary condition structure for the primal problem
    type(t_boundaryCondition), intent(in) :: rbdrCondPrimal

    ! boundary condition structure for the dual problem
    type(t_boundaryCondition), intent(in) :: rbdrCondDual
!</input>

!<inputoutput>
    ! parameter list
    type(t_parlist), intent(inout) :: rparlist

    ! problem structure
    type(t_problem), intent(inout) :: rproblem

    ! time-stepping structure
    type(t_timestep), intent(inout) :: rtimestep

    ! solver struchture
    type(t_solver), intent(inout), target :: rsolver

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!<output>
    ! primal solution vector
    type(t_vectorBlock), intent(out), target :: rsolutionPrimal

    ! dual solution vector
    type(t_vectorBlock), intent(out), target :: rsolutionDual
!</output>
!</subroutine>

    print *, "Solution of the primal and dual problem for"
    print *, "time-dependent flows is not implemented yet"
    stop

  end subroutine transp_solveTransientPrimalDual

  !*****************************************************************************

!<subroutine>

  subroutine transp_solvePseudoTransientPrimal(rparlist, ssectionName,&
      rbdrCond, rproblem, rtimestep, rsolver, rsolution, rcollection)

!<description>
    ! This subroutine solves the pseudo-transient primal flow problem
    !
    ! $$\frac{\partial u}{\partial \tau}+\nabla\cdot{\bf f}(u)=s(u)$$
    !
    ! for a scalar quantity $u$ in the domain $\Omega$.
!</description>

!<input>
    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName

    ! boundary condition structure
    type(t_boundaryCondition), intent(in) :: rbdrCond
!</input>

!<inputoutput>
    ! parameter list
    type(t_parlist), intent(inout) :: rparlist

    ! problem structure
    type(t_problem), intent(inout) :: rproblem

    ! time-stepping structure
    type(t_timestep), intent(inout) :: rtimestep

    ! solver struchture
    type(t_solver), intent(inout), target :: rsolver

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!<output>
    ! primal solution vector
    type(t_vectorBlock), intent(out), target :: rsolution
!</output>
!</subroutine>

    ! Pointer to the multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! Pointer to the discretisation structure
    type(t_blockDiscretisation), pointer :: p_rdiscretisation

    ! Vector for the right-hand side
    type(t_vectorBlock), target :: rrhs

    ! Vector for the element-wise error distribution
    type(t_vectorScalar) :: relementError

    ! Structure for h-adaptation
    type(t_hadapt) :: rhadapt

    ! Structure for the sparsity pattern
    type(t_graph) :: rgraph

    ! Timer structures
    type(t_timer), pointer :: p_rtimerPrePostprocess
    type(t_timer), pointer :: p_rtimerSolution
    type(t_timer), pointer :: p_rtimerErrorEstimation
    type(t_timer), pointer :: p_rtimerAdaptation
    type(t_timer), pointer :: p_rtimerTriangulation
    type(t_timer), pointer :: p_rtimerAssemblyCoeff

    ! section names
    character(LEN=SYS_STRLEN) :: sadaptivityName

    ! local variables
    real(DP) :: derror
    integer :: templateMatrix, systemMatrix, discretisation
    integer :: nlmin, iadapt, nadapt, irhstype, ivelocitytype


    ! Get timer structures
    p_rtimerPrePostprocess => collct_getvalue_timer(rcollection, 'rtimerPrePostprocess')
    p_rtimerSolution => collct_getvalue_timer(rcollection, 'rtimerSolution')
    p_rtimerErrorEstimation => collct_getvalue_timer(rcollection, 'rtimerErrorEstimation')
    p_rtimerAdaptation => collct_getvalue_timer(rcollection, 'rtimerAdaptation')
    p_rtimerTriangulation => collct_getvalue_timer(rcollection, 'rtimerTriangulation')
    p_rtimerAssemblyCoeff => collct_getvalue_timer(rcollection, 'rtimerAssemblyCoeff')

    ! Start time measurement for pre-processing
    call stat_startTimer(p_rtimerPrePostprocess, STAT_TIMERSHORT)

    ! Get global parameters
    call parlst_getvalue_int(rparlist, ssectionName, 'irhstype', irhstype)
    call parlst_getvalue_int(rparlist, ssectionName, 'ivelocitytype', ivelocitytype)
    call parlst_getvalue_int(rparlist, ssectionName, 'discretisation', discretisation)

    ! Set pointer to maximum problem level
    p_rproblemLevel => rproblem%p_rproblemLevelMax
    p_rdiscretisation => p_rproblemLevel%Rdiscretisation(discretisation)

    ! Initialize the solution vector and impose boundary conditions explicitly
    call lsysbl_createVectorBlock(p_rdiscretisation, rsolution,&
        .false., ST_DOUBLE)
    call transp_initSolution(rparlist, ssectionName, p_rproblemLevel,&
        0.0_DP, rsolution, rcollection)
    call bdrf_filterVectorExplicit(rbdrCond, rsolution, 0.0_DP)

    !---------------------------------------------------------------------------
    ! Initialize the h-adaptation structure
    !---------------------------------------------------------------------------

    call parlst_getvalue_string(rparlist, ssectionName, 'adaptivity', sadaptivityName, '')
    if (trim(adjustl(sadaptivityName)) .ne. '') then

      call parlst_getvalue_int(rparlist, trim(sadaptivityName), 'nadapt', nadapt)

      if (nadapt > 0) then

        ! Initialize adaptation structure from parameter list
        call hadapt_initFromParameterlist(rhadapt, rparlist, sadaptivityName)

        ! Generate a dynamic graph for the sparsity pattern and attach
        ! it to the collection structure which is used to pass this type
        ! to the callback function for h-adaptation
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'templateMatrix', templateMatrix)
        call grph_createGraphFromMatrix(p_rproblemLevel&
            %Rmatrix(templateMatrix), rgraph)
        call collct_setvalue_graph(rcollection, 'sparsitypattern',&
            rgraph, .true.)

        ! Attach the primal solution vector to the collection structure
        rcollection%p_rvectorQuickAccess1 => rsolution

      end if   ! nadapt > 0

    else

      ! No h-adaptation
      nadapt = 0

    end if

    ! Stop time measurement for pre-processing
    call stat_stopTimer(p_rtimerPrePostprocess)


    ! Adaptation loop
    adaptloop: do iadapt = 0, nadapt

      !-------------------------------------------------------------------------
      ! Compute the steady-state solution
      !-------------------------------------------------------------------------

      ! Start time measurement for solution procedure
      call stat_startTimer(p_rtimerSolution, STAT_TIMERSHORT)

      ! Calculate the velocity field
      nlmin = solver_getMinimumMultigridlevel(rsolver)
      call transp_calcVelocityField(rparlist, ssectionName,&
          p_rproblemLevel, 0.0_DP, rcollection, nlmin)

      ! Attach the boundary condition
      call solver_setBoundaryCondition(rsolver, rbdrCond, .true.)

      ! Set primal problem mode
      call parlst_addvalue(rparlist, ssectionName, 'mode', 'primal')

      ! Reset the time-stepping algorithm
      call tstep_resetTimestep(rtimestep, .false.)
      call solver_resetSolver(rsolver, .false.)

      ! Check if right-hand side vector exists
      if (irhstype > 0) then
        call lsysbl_createVectorBlock(rsolution, rrhs)
        call transp_initRHS(rparlist, ssectionName, p_rproblemLevel,&
            0.0_DP, rrhs, rcollection)

        ! Prepare quick access arrays of the collection
        rcollection%SquickAccess(1) = ssectionName

        ! Solve the primal problem with non-zero right-hand side
        call tstep_performPseudoStepping(p_rproblemLevel, rtimestep,&
            rsolver, rsolution, transp_nlsolverCallback, rcollection, rrhs)

        ! Release right-hand side vector
        call lsysbl_releaseVector(rrhs)

      else

        ! Prepare quick access arrays of the collection
        rcollection%SquickAccess(1) = ssectionName

        ! Solve the primal problem without right-hand side
        call tstep_performPseudoStepping(p_rproblemLevel, rtimestep,&
            rsolver, rsolution, transp_nlsolverCallback, rcollection)
      end if

      ! Stop time measurement for solution procedure
      call stat_stopTimer(p_rtimerSolution)


      if (iadapt .eq. nadapt) exit adaptloop

      !-------------------------------------------------------------------------
      ! Perform recovery-based error estimation
      !-------------------------------------------------------------------------

      ! Start time measurement for error estimation
      call stat_startTimer(p_rtimerErrorEstimation, STAT_TIMERSHORT)

      ! Compute the error estimator using recovery techniques
      call transp_estimateRecoveryError(rparlist, ssectionname,&
          p_rproblemLevel, rsolution, 0.0_DP, relementError, derror,&
          rcollection)

      ! Stop time measurement for error estimation
      call stat_stopTimer(p_rtimerErrorEstimation)


      !-------------------------------------------------------------------------
      ! Perform h-adaptation
      !-------------------------------------------------------------------------

      ! Start time measurement for mesh adaptation
      call stat_startTimer(p_rtimerAdaptation, STAT_TIMERSHORT)

      ! Set the name of the template matrix
      rcollection%SquickAccess(1) = 'sparsitypattern'

      ! Attach the primal solution vector to the collection structure
      rcollection%p_rvectorQuickAccess1 => rsolution

      ! Perform h-adaptation and update the triangulation structure
      call transp_adaptTriangulation(rhadapt, p_rproblemLevel&
          %rtriangulation, relementError, rcollection)

      ! Update the template matrix according to the sparsity pattern
      call grph_generateMatrix(rgraph, p_rproblemLevel&
          %Rmatrix(templateMatrix))

      ! Resize the solution vector accordingly
      call lsysbl_resizeVectorBlock(rsolution, p_rproblemLevel&
          %Rmatrix(templateMatrix)%NEQ, .false.)

      ! Release element-wise error distribution
      call lsyssc_releaseVector(relementError)

      ! Stop time measurement for mesh adaptation
      call stat_stopTimer(p_rtimerAdaptation)


      !-------------------------------------------------------------------------
      ! Re-generate the discretisation and coefficient matrices
      !-------------------------------------------------------------------------

      ! Start time measurement for generation of the triangulation
      call stat_startTimer(p_rtimerTriangulation, STAT_TIMERSHORT)

      ! Generate standard mesh from raw mesh
      call tria_initStandardMeshFromRaw(&
          p_rproblemLevel%rtriangulation, rproblem%rboundary)

      ! Stop time measurement for generation of the triangulation
      call stat_stopTimer(p_rtimerTriangulation)


      ! Start time measurement for generation of constant coefficient matrices
      call stat_startTimer(p_rtimerAssemblyCoeff, STAT_TIMERSHORT)

      ! Re-initialize all constant coefficient matrices
      call transp_initProblemLevel(rparlist, ssectionName,&
          p_rproblemLevel, rcollection)

      ! Prepare internal data arrays of the solver structure
      call parlst_getvalue_int(rparlist, ssectionName, 'systemmatrix', systemMatrix)
      call flagship_updateSolverMatrix(p_rproblemLevel, rsolver,&
          systemMatrix, SYSTEM_INTERLEAVEFORMAT, UPDMAT_ALL)
      call solver_updateStructure(rsolver)

      ! Stop time measurement for generation of constant coefficient matrices
      call stat_stopTimer(p_rtimerAssemblyCoeff)

    end do adaptloop


    ! Release adaptation structure
    if (nadapt > 0) then
      call hadapt_releaseAdaptation(rhadapt)
      call grph_releaseGraph(rgraph)
    end if

  end subroutine transp_solvePseudoTransientPrimal

  !*****************************************************************************

!<subroutine>

  subroutine transp_solvePseudoTransientPrimalDual(rparlist, ssectionName,&
      rbdrCondPrimal, rbdrCondDual, rproblem, rtimestep, rsolver,&
      rsolutionPrimal, rsolutionDual, rcollection)

!<description>
    ! This subroutine solves the pseudo-transient primal flow problem
    !
    ! $$\frac{\partial u}{\partial \tau}+\nabla\cdot{\bf f}(u)=s(u)$$
    !
    ! for a scalar quantity $u$ in the domain $\Omega$ both for the
    ! primal and the dual formulation.
!</description>

!<input>
    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName

    ! boundary condition structure for the primal problem
    type(t_boundaryCondition), intent(in) :: rbdrCondPrimal

    ! boundary condition structure for the dual problem
    type(t_boundaryCondition), intent(in) :: rbdrCondDual
!</input>

!<inputoutput>
    ! parameter list
    type(t_parlist), intent(inout) :: rparlist

    ! problem structure
    type(t_problem), intent(inout) :: rproblem

    ! time-stepping structure
    type(t_timestep), intent(inout) :: rtimestep

    ! solver struchture
    type(t_solver), intent(inout), target :: rsolver

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!<output>
    ! primal solution vector
    type(t_vectorBlock), intent(out), target :: rsolutionPrimal

    ! dual solution vector
    type(t_vectorBlock), intent(out), target :: rsolutionDual
!</output>
!</subroutine>

    ! Pointer to the multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! Pointer to the discretisation structure
    type(t_blockDiscretisation), pointer :: p_rdiscretisation

    ! Vector for the linear target functional and the right-hand side
    type(t_vectorBlock), target :: rtargetFunc, rrhs

    ! Vector for the element-wise error distribution
    type(t_vectorScalar) :: relementError

    ! Structure for h-adaptation
    type(t_hadapt) :: rhadapt

    ! Structure for the sparsity pattern
    type(t_graph) :: rgraph

    ! Timer structures
    type(t_timer), pointer :: p_rtimerPrePostprocess
    type(t_timer), pointer :: p_rtimerSolution
    type(t_timer), pointer :: p_rtimerErrorEstimation
    type(t_timer), pointer :: p_rtimerAdaptation
    type(t_timer), pointer :: p_rtimerTriangulation
    type(t_timer), pointer :: p_rtimerAssemblyCoeff

    ! section names
    character(LEN=SYS_STRLEN) :: sadaptivityName

    ! local variables
    real(dp) :: derror
    integer :: templateMatrix, systemMatrix, discretisation
    integer :: nlmin, iadapt, nadapt, irhstype, ivelocitytype


    ! Get timer structures
    p_rtimerPrePostprocess => collct_getvalue_timer(rcollection, 'rtimerPrePostprocess')
    p_rtimerSolution => collct_getvalue_timer(rcollection, 'rtimerSolution')
    p_rtimerErrorEstimation => collct_getvalue_timer(rcollection, 'rtimerErrorEstimation')
    p_rtimerAdaptation => collct_getvalue_timer(rcollection, 'rtimerAdaptation')
    p_rtimerTriangulation => collct_getvalue_timer(rcollection, 'rtimerTriangulation')
    p_rtimerAssemblyCoeff => collct_getvalue_timer(rcollection, 'rtimerAssemblyCoeff')

    ! Start time measurement for pre-processing
    call stat_startTimer(p_rtimerPrePostprocess, STAT_TIMERSHORT)

    ! Get global parameters
    call parlst_getvalue_int(rparlist, ssectionName, 'irhstype', irhstype)
    call parlst_getvalue_int(rparlist, ssectionName, 'ivelocitytype', ivelocitytype)
    call parlst_getvalue_int(rparlist, ssectionName, 'discretisation', discretisation)

    ! Set pointer to maximum problem level
    p_rproblemLevel   => rproblem%p_rproblemLevelMax
    p_rdiscretisation => p_rproblemLevel%Rdiscretisation(discretisation)

    ! Initialize the solution vector and impose boundary conditions explicitly
    call lsysbl_createVectorBlock(p_rdiscretisation, rsolutionPrimal,&
        .false., ST_DOUBLE)
    call transp_initSolution(rparlist, ssectionName, p_rproblemLevel,&
        0.0_DP, rsolutionPrimal, rcollection)
    call bdrf_filterVectorExplicit(rbdrCondPrimal, rsolutionPrimal, 0.0_DP)

    !---------------------------------------------------------------------------
    ! Initialize the h-adaptation structure
    !---------------------------------------------------------------------------

    call parlst_getvalue_string(rparlist, ssectionName, 'adaptivity', sadaptivityName, '')
    if (trim(adjustl(sadaptivityName)) .ne. '') then

      call parlst_getvalue_int(rparlist, trim(sadaptivityName), 'nadapt', nadapt)

      if (nadapt > 0) then

        ! Initialize adaptation structure from parameter list
        call hadapt_initFromParameterlist(rhadapt, rparlist, sadaptivityName)

        ! Generate a dynamic graph for the sparsity pattern and attach
        ! it to the collection structure which is used to pass this type
        ! to the callback function for h-adaptation
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'templateMatrix', templateMatrix)
        call grph_createGraphFromMatrix(p_rproblemLevel&
            %Rmatrix(templateMatrix), rgraph)
        call collct_setvalue_graph(rcollection, 'sparsitypattern',&
            rgraph, .true.)

        ! Attach the primal solution vector to the collection structure
        rcollection%p_rvectorQuickAccess1 => rsolutionPrimal

      end if   ! nadapt > 0

    else

      ! No h-adaptation
      nadapt = 0

    end if

    ! Stop time measurement for pre-processing
    call stat_stopTimer(p_rtimerPrePostprocess)

    adaptloop: do iadapt = 0, nadapt

      !-------------------------------------------------------------------------
      ! Compute steady-state solution for the primal problem
      !-------------------------------------------------------------------------

      ! Start time measurement for solution procedure
      call stat_startTimer(p_rtimerSolution, STAT_TIMERSHORT)

      ! Calculate the velocity field
      nlmin = solver_getMinimumMultigridlevel(rsolver)
      call transp_calcVelocityField(rparlist, ssectionName,&
          p_rproblemLevel, rtimestep%dTime, rcollection, nlmin)

      ! Attach the boundary condition
      call solver_setBoundaryCondition(rsolver, rbdrCondPrimal, .true.)

      ! Set primal problem mode
      call parlst_addvalue(rparlist, ssectionName, 'mode', 'primal')

      ! Reset the time-stepping algorithm
      call tstep_resetTimestep(rtimestep, .false.)
      call solver_resetSolver(rsolver, .false.)

      ! Check if right-hand side vector exists
      if (irhstype > 0) then
        call lsysbl_createVectorBlock(rsolutionPrimal, rrhs)
        call transp_initRHS(rparlist, ssectionName, p_rproblemLevel,&
            0.0_DP, rrhs, rcollection)

        ! Prepare quick access arrays of the collection
        rcollection%SquickAccess(1) = ssectionName

        ! Solve the primal problem with non-zero right-hand side
        call tstep_performPseudoStepping(p_rproblemLevel, rtimestep,&
            rsolver, rsolutionPrimal, transp_nlsolverCallback,&
            rcollection, rrhs)

        ! Release right-hand side vector
        call lsysbl_releaseVector(rrhs)

      else

        ! Prepare quick access arrays of the collection
        rcollection%SquickAccess(1) = ssectionName

        ! Solve the primal problem without right-hand side
        call tstep_performPseudoStepping(p_rproblemLevel, rtimestep,&
            rsolver, rsolutionPrimal, transp_nlsolverCallback,&
            rcollection)
      end if

      ! Stop time measurement for solution procedure
      call stat_stopTimer(p_rtimerSolution)


      !-------------------------------------------------------------------------
      ! Compute the right-hand side for the dual problem
      !-------------------------------------------------------------------------

      ! Start time measurement for error estimation
      call stat_startTimer(p_rtimerErrorEstimation, STAT_TIMERSHORT)

      ! Initialize target functional
      call lsysbl_createVectorBlock(rsolutionPrimal, rtargetFunc)
      call transp_initTargetFunc(rparlist, ssectionName,&
          p_rproblemLevel, 0.0_DP, rtargetFunc, rcollection)

      ! Stop time measurement for error estimation
      call stat_stopTimer(p_rtimerErrorEstimation)


      !-------------------------------------------------------------------------
      ! Compute steady-state solution for the dual problem
      !-------------------------------------------------------------------------

      ! Start time measurement for solution procedure
      call stat_startTimer(p_rtimerSolution, STAT_TIMERSHORT)

      ! Attach the boundary condition
      call solver_setBoundaryCondition(rsolver, rbdrCondDual, .true.)

      ! Set dual problem mode
      call parlst_addvalue(rparlist, ssectionName, 'mode', 'dual')

      ! Create dual solution vector initialised by zeros
      call lsysbl_releaseVector(rsolutionDual)
      call lsysbl_createVectorBlock(rsolutionPrimal, rsolutionDual, .true.)

      ! Reset the time-stepping and solver algorithms
      call tstep_resetTimestep(rtimestep, .false.)
      call solver_resetSolver(rsolver, .false.)

      ! Prepare quick access arrays of the collection
      rcollection%SquickAccess(1) = ssectionName

      ! Solve the dual problem
      call tstep_performPseudoStepping(p_rproblemLevel, rtimestep,&
          rsolver, rsolutionDual, transp_nlsolverCallback,&
          rcollection, rtargetFunc)

      ! Release discretised target functional
      call lsysbl_releaseVector(rtargetFunc)

      ! Stop time measurement for solution procedure
      call stat_stopTimer(p_rtimerSolution)


      if (iadapt .eq. nadapt) exit adaptloop

      !-------------------------------------------------------------------------
      ! Perform goal-oriented error estimation
      !-------------------------------------------------------------------------

      ! Start time measurement for error estimation
      call stat_startTimer(p_rtimerErrorEstimation, STAT_TIMERSHORT)

      ! Calculate the velocity field
      nlmin = solver_getMinimumMultigridlevel(rsolver)
      call transp_calcVelocityField(rparlist, ssectionName,&
          p_rproblemLevel, rtimestep%dTime, rcollection, nlmin)

      ! Attach the boundary condition
      call solver_setBoundaryCondition(rsolver, rbdrCondPrimal, .true.)

      ! Set primal problem mode
      call parlst_addvalue(rparlist, ssectionName, 'mode', 'primal')

      ! Check if right-hand side vector exists
      if (irhstype > 0) then

        ! Initialize right-hand side vector
        call lsysbl_createVectorBlock(rsolutionPrimal, rrhs)
        call transp_initRHS(rparlist, ssectionName, p_rproblemLevel,&
            0.0_DP, rrhs, rcollection)

        ! Prepare quick access arrays of the collection
        rcollection%SquickAccess(1) = ssectionName

        ! Compute the error in the quantity of interest
        call transp_estimateTargetFuncError(rparlist, ssectionName,&
            p_rproblemLevel, rtimestep, rsolver, rsolutionPrimal,&
            rsolutionDual, rcollection, relementError, derror, rrhs)

        ! Release right-hand side vector
        call lsysbl_releaseVector(rrhs)

      else

        ! Prepare quick access arrays of the collection
        rcollection%SquickAccess(1) = ssectionName

        ! Compute the error in the quantity of interest
        call transp_estimateTargetFuncError(rparlist, ssectionName,&
            p_rproblemLevel, rtimestep, rsolver, rsolutionPrimal,&
            rsolutionDual, rcollection, relementError, derror)

      end if

      ! Stop time measurement for error estimation
      call stat_stopTimer(p_rtimerErrorEstimation)


      !-------------------------------------------------------------------------
      ! Perform h-adaptation
      !-------------------------------------------------------------------------

      ! Start time measurement for mesh adaptation
      call stat_startTimer(p_rtimerAdaptation, STAT_TIMERSHORT)

      ! Set the names of the template matrix and the solution vector
      rcollection%SquickAccess(1) = 'sparsitypattern'

      ! Attach the primal solution vector to the collection structure
      rcollection%p_rvectorQuickAccess1 => rsolutionPrimal

      ! Perform h-adaptation and update the triangulation structure
      call transp_adaptTriangulation(rhadapt, p_rproblemLevel&
          %rtriangulation, relementError, rcollection)

      ! Update the template matrix according to the sparsity pattern
      call grph_generateMatrix(rgraph, p_rproblemLevel&
          %Rmatrix(templateMatrix))

      ! Resize the solution vector accordingly
      call lsysbl_resizeVectorBlock(rsolutionPrimal, p_rproblemLevel&
          %Rmatrix(templateMatrix)%NEQ, .false.)

      ! Release element-wise error distribution
      call lsyssc_releaseVector(relementError)

      ! Stop time measurement for mesh adaptation
      call stat_stopTimer(p_rtimerAdaptation)


      !-------------------------------------------------------------------------
      ! Re-generate the discretisation and coefficient matrices
      !-------------------------------------------------------------------------

      ! Start time measurement for generation of the triangulation
      call stat_startTimer(p_rtimerTriangulation, STAT_TIMERSHORT)

      ! Generate standard mesh from raw mesh
      call tria_initStandardMeshFromRaw(&
          p_rproblemLevel%rtriangulation, rproblem%rboundary)

      ! Stop time measurement for generation of the triangulation
      call stat_stopTimer(p_rtimerTriangulation)

      !-------------------------------------------------------------------

      ! Start time measurement for generation of constant coefficient matrices
      call stat_startTimer(p_rtimerAssemblyCoeff, STAT_TIMERSHORT)

      ! Re-initialize all constant coefficient matrices
      call transp_initProblemLevel(rparlist, ssectionName, p_rproblemLevel, rcollection)

      ! Prepare internal data arrays of the solver structure
      call parlst_getvalue_int(rparlist, ssectionName, 'systemMatrix', systemMatrix)
      call flagship_updateSolverMatrix(p_rproblemLevel, rsolver,&
          systemMatrix, SYSTEM_INTERLEAVEFORMAT, UPDMAT_ALL)
      call solver_updateStructure(rsolver)

      ! Stop time measurement for generation of constant coefficient matrices
      call stat_stopTimer(p_rtimerAssemblyCoeff)

    end do adaptloop


    ! Release adaptation structure
    if (nadapt > 0) then
      call hadapt_releaseAdaptation(rhadapt)
      call grph_releaseGraph(rgraph)
    end if

  end subroutine transp_solvePseudoTransientPrimalDual

  !*****************************************************************************

!<subroutine>

  subroutine transp_solveSteadyStatePrimal(rparlist, ssectionName,&
      rbdrCond, rproblem, rtimestep, rsolver, rsolution, rcollection)

!<description>
    ! This subroutine solves the steady-state primal flow problem
    !
    ! $$\nabla\cdot{\bf f}(u)=s(u)$$
    !
    ! for a scalar quantity $u$ in the domain $\Omega$.
!</description>

!<input>
    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName

    ! boundary condition structure
    type(t_boundaryCondition), intent(in) :: rbdrCond
!</input>

!<inputoutput>
    ! parameter list
    type(t_parlist), intent(inout) :: rparlist

    ! problem structure
    type(t_problem), intent(inout) :: rproblem

    ! time-stepping structure
    type(t_timestep), intent(inout) :: rtimestep

    ! solver struchture
    type(t_solver), intent(inout), target :: rsolver

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!<output>
    ! primal solution vector
    type(t_vectorBlock), intent(out), target :: rsolution
!</output>
!</subroutine>

    ! Pointer to the multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! Pointer to the discretisation structure
    type(t_blockDiscretisation), pointer :: p_rdiscretisation

    ! Vector for the right-hand side
    type(t_vectorBlock), target :: rrhs

    ! Vector for the element-wise error distribution
    type(t_vectorScalar) :: relementError

    ! Structure for h-adaptation
    type(t_hadapt) :: rhadapt

    ! Structure for the sparsity pattern
    type(t_graph) :: rgraph

    ! Timer structures
    type(t_timer), pointer :: p_rtimerPrePostprocess
    type(t_timer), pointer :: p_rtimerSolution
    type(t_timer), pointer :: p_rtimerErrorEstimation
    type(t_timer), pointer :: p_rtimerAdaptation
    type(t_timer), pointer :: p_rtimerTriangulation
    type(t_timer), pointer :: p_rtimerAssemblyCoeff

    ! section names
    character(LEN=SYS_STRLEN) :: sadaptivityName

    ! local variables
    real(dp) :: derror
    integer :: templateMatrix, systemMatrix, discretisation
    integer :: nlmin, iadapt, nadapt, irhstype, ivelocitytype


    ! Get timer structures
    p_rtimerPrePostprocess => collct_getvalue_timer(rcollection, 'rtimerPrePostprocess')
    p_rtimerSolution => collct_getvalue_timer(rcollection, 'rtimerSolution')
    p_rtimerErrorEstimation => collct_getvalue_timer(rcollection, 'rtimerErrorEstimation')
    p_rtimerAdaptation => collct_getvalue_timer(rcollection, 'rtimerAdaptation')
    p_rtimerTriangulation => collct_getvalue_timer(rcollection, 'rtimerTriangulation')
    p_rtimerAssemblyCoeff => collct_getvalue_timer(rcollection, 'rtimerAssemblyCoeff')

    ! Start time measurement for pre-processing
    call stat_startTimer(p_rtimerPrePostprocess, STAT_TIMERSHORT)

    ! Adjust time stepping scheme
    rtimestep%ctimestepType = TSTEP_THETA_SCHEME
    rtimestep%dinitialTime  = 0.0_DP
    rtimestep%dinitialStep  = 1.0_DP
    rtimestep%dfinalTime    = 1.0_DP
    rtimestep%theta         = 1.0_DP

    ! Get global parameters
    call parlst_getvalue_int(rparlist, ssectionName, 'irhstype', irhstype)
    call parlst_getvalue_int(rparlist, ssectionName, 'ivelocitytype', ivelocitytype)
    call parlst_getvalue_int(rparlist, ssectionName, 'discretisation', discretisation)

    ! Set pointer to maximum problem level
    p_rproblemLevel   => rproblem%p_rproblemLevelMax
    p_rdiscretisation => p_rproblemLevel%Rdiscretisation(discretisation)

    ! Initialize the solution vector and impose boundary conditions explicitly
    call lsysbl_createVectorBlock(p_rdiscretisation, rsolution, .false., ST_DOUBLE)
    call transp_initSolution(rparlist, ssectionName, p_rproblemLevel,&
        0.0_DP, rsolution, rcollection)
    call bdrf_filterVectorExplicit(rbdrCond, rsolution, 0.0_DP)

    !---------------------------------------------------------------------------
    ! Initialize the h-adaptation structure
    !---------------------------------------------------------------------------

    call parlst_getvalue_string(rparlist, ssectionName, 'adaptivity', sadaptivityName, '')
    if (trim(adjustl(sadaptivityName)) .ne. '') then

      call parlst_getvalue_int(rparlist, trim(sadaptivityName), 'nadapt', nadapt)

      if (nadapt > 0) then

        ! Initialize adaptation structure from parameter list
        call hadapt_initFromParameterlist(rhadapt, rparlist, sadaptivityName)

        ! Generate a dynamic graph for the sparsity pattern and attach
        ! it to the collection structure which is used to pass this type
        ! to the callback function for h-adaptation
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'templateMatrix', templateMatrix)
        call grph_createGraphFromMatrix(p_rproblemLevel&
            %Rmatrix(templateMatrix), rgraph)
        call collct_setvalue_graph(rcollection, 'sparsitypattern',&
            rgraph, .true.)

        ! Attach the primal solution vector to the collection structure
        rcollection%p_rvectorQuickAccess1 => rsolution

      end if   ! nadapt > 0

    else

      ! No h-adaptation
      nadapt = 0

    end if

    ! Stop time measurement for pre-processing
    call stat_stopTimer(p_rtimerPrePostprocess)


    adaptloop: do iadapt = 0, nadapt

      !-------------------------------------------------------------------------
      ! Compute steady-state solution for the primal problem
      !-------------------------------------------------------------------------

      ! Start time measurement for solution procedure
      call stat_startTimer(p_rtimerSolution, STAT_TIMERSHORT)

      ! Calculate the velocity field
      nlmin = solver_getMinimumMultigridlevel(rsolver)
      call transp_calcVelocityField(rparlist, ssectionName, p_rproblemLevel,&
                                    rtimestep%dTime, rcollection, nlmin)

      ! Attach the boundary condition
      call solver_setBoundaryCondition(rsolver, rbdrCond, .true.)

      ! Set primal problem mode
      call parlst_addvalue(rparlist, ssectionName, 'mode', 'primal')

      ! Reset the time-stepping algorithm
      call tstep_resetTimestep(rtimestep, .false.)
      call solver_resetSolver(rsolver, .false.)

      ! Check if right-hand side vector exists
      if (irhstype > 0) then
        call lsysbl_createVectorblock(rsolution, rrhs)
        call transp_initRHS(rparlist, ssectionName, p_rproblemLevel,&
            0.0_DP, rrhs, rcollection)

        ! Prepare quick access arrays of the collection
        rcollection%SquickAccess(1) = ssectionName

        ! Solve the primal problem with non-zero right-hand side
        call tstep_performThetaStep(p_rproblemLevel, rtimestep,&
            rsolver, rsolution, transp_nlsolverCallback, rcollection,&
            rrhs)

        ! Release right-hand side vector
        call lsysbl_releaseVector(rrhs)

      else

        ! Prepare quick access arrays of the collection
        rcollection%SquickAccess(1) = ssectionName

        ! Solve the primal problem without right-hand side
        call tstep_performThetaStep(p_rproblemLevel, rtimestep,&
            rsolver, rsolution, transp_nlsolverCallback, rcollection)
      end if

      ! Stop time measurement for solution procedure
      call stat_stopTimer(p_rtimerSolution)


      if (iadapt .eq. nadapt) exit adaptloop

      !-------------------------------------------------------------------------
      ! Perform recovery-based error estimation
      !-------------------------------------------------------------------------

      ! Start time measurement for error estimation
      call stat_startTimer(p_rtimerErrorEstimation, STAT_TIMERSHORT)

      ! Compute the error estimator using recovery techniques
      call transp_estimateRecoveryError(rparlist, ssectionname,&
          p_rproblemLevel, rsolution, 0.0_DP, relementError, derror,&
          rcollection)

      ! Stop time measurement for error estimation
      call stat_stopTimer(p_rtimerErrorEstimation)


      !-------------------------------------------------------------------------
      ! Perform h-adaptation
      !-------------------------------------------------------------------------

      ! Start time measurement for mesh adaptation
      call stat_startTimer(p_rtimerAdaptation, STAT_TIMERSHORT)

      ! Set the names of the template matrix and the solution vector
      rcollection%SquickAccess(1) = 'sparsitypattern'

      ! Attach the primal solution vector to the collection structure
      rcollection%p_rvectorQuickAccess1 => rsolution

      ! Perform h-adaptation and update the triangulation structure
      call transp_adaptTriangulation(rhadapt, p_rproblemLevel&
          %rtriangulation, relementError, rcollection)

      ! Update the template matrix according to the sparsity pattern
      call grph_generateMatrix(rgraph, p_rproblemLevel&
          %Rmatrix(templateMatrix))

      ! Resize the solution vector accordingly
      call lsysbl_resizeVectorBlock(rsolution, p_rproblemLevel&
          %Rmatrix(templateMatrix)%NEQ, .false.)

      ! Release element-wise error distribution
      call lsyssc_releaseVector(relementError)

      ! Stop time measurement for mesh adaptation
      call stat_stopTimer(p_rtimerAdaptation)


      !-------------------------------------------------------------------------
      ! Re-generate the discretisation and coefficient matrices
      !-------------------------------------------------------------------------

      ! Start time measurement for generation of the triangulation
      call stat_startTimer(p_rtimerTriangulation, STAT_TIMERSHORT)

      ! Generate standard mesh from raw mesh
      call tria_initStandardMeshFromRaw(&
          p_rproblemLevel%rtriangulation, rproblem%rboundary)

      ! Stop time measurement for generation of the triangulation
      call stat_stopTimer(p_rtimerTriangulation)


      ! Start time measurement for generation of constant coefficient matrices
      call stat_startTimer(p_rtimerAssemblyCoeff, STAT_TIMERSHORT)

      ! Re-initialize all constant coefficient matrices
      call transp_initProblemLevel(rparlist, ssectionName, p_rproblemLevel, rcollection)

      ! Prepare internal data arrays of the solver structure
      call parlst_getvalue_int(rparlist, ssectionName, 'systemmatrix', systemMatrix)
      call flagship_updateSolverMatrix(p_rproblemLevel, rsolver,&
          systemMatrix, SYSTEM_INTERLEAVEFORMAT, UPDMAT_ALL)
      call solver_updateStructure(rsolver)

      ! Stop time measurement for generation of constant coefficient matrices
      call stat_stopTimer(p_rtimerAssemblyCoeff)

    end do adaptloop


    ! Release adaptation structure
    if (nadapt > 0) then
      call hadapt_releaseAdaptation(rhadapt)
      call grph_releaseGraph(rgraph)
    end if

  end subroutine transp_solveSteadyStatePrimal

  !*****************************************************************************

!<subroutine>

  subroutine transp_solveSteadyStatePrimalDual(rparlist, ssectionName,&
      rbdrCondPrimal, rbdrCondDual, rproblem, rtimestep, rsolver,&
      rsolutionPrimal, rsolutionDual, rcollection)

!<description>
    ! This subroutine solves the steady-state primal flow problem
    !
    ! $$\nabla\cdot{\bf f}(u)=s(u)$$
    !
    ! for a scalar quantity $u$ in the domain $\Omega$ both for the
    ! primal and the dual formulation.
!</description>

!<input>
    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName

    ! boundary condition structure for the primal problem
    type(t_boundaryCondition), intent(in) :: rbdrCondPrimal

    ! boundary condition structure for the dual problem
    type(t_boundaryCondition), intent(in) :: rbdrCondDual
!</input>

!<inputoutput>
    ! parameter list
    type(t_parlist), intent(inout) :: rparlist

    ! problem structure
    type(t_problem), intent(inout) :: rproblem

    ! time-stepping structure
    type(t_timestep), intent(inout) :: rtimestep

    ! solver struchture
    type(t_solver), intent(inout), target :: rsolver

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!<output>
    ! primal solution vector
    type(t_vectorBlock), intent(out), target :: rsolutionPrimal

    ! dual solution vector
    type(t_vectorBlock), intent(out), target :: rsolutionDual
!</output>
!</subroutine>

    ! Pointer to the multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! Pointer to the discretisation structure
    type(t_blockDiscretisation), pointer :: p_rdiscretisation

    ! Vector for the linear target functional and the right-hand side
    type(t_vectorBlock), target :: rtargetFunc, rrhs

    ! Vector for the element-wise error distribution
    type(t_vectorScalar) :: relementError

    ! Structure for h-adaptation
    type(t_hadapt) :: rhadapt

    ! Structure for the sparsity pattern
    type(t_graph) :: rgraph

    ! Timer structures
    type(t_timer), pointer :: p_rtimerPrePostprocess
    type(t_timer), pointer :: p_rtimerSolution
    type(t_timer), pointer :: p_rtimerErrorEstimation
    type(t_timer), pointer :: p_rtimerAdaptation
    type(t_timer), pointer :: p_rtimerTriangulation
    type(t_timer), pointer :: p_rtimerAssemblyCoeff

    ! section names
    character(LEN=SYS_STRLEN) :: sadaptivityName

    ! local variables
    real(dp) :: derror
    integer :: templateMatrix, systemMatrix, discretisation
    integer :: nlmin, iadapt, nadapt, irhstype, ivelocitytype


    ! Get timer structures
    p_rtimerPrePostprocess => collct_getvalue_timer(rcollection, 'rtimerPrePostprocess')
    p_rtimerSolution => collct_getvalue_timer(rcollection, 'rtimerSolution')
    p_rtimerErrorEstimation => collct_getvalue_timer(rcollection, 'rtimerErrorEstimation')
    p_rtimerAdaptation => collct_getvalue_timer(rcollection, 'rtimerAdaptation')
    p_rtimerTriangulation => collct_getvalue_timer(rcollection, 'rtimerTriangulation')
    p_rtimerAssemblyCoeff => collct_getvalue_timer(rcollection, 'rtimerAssemblyCoeff')

    ! Start time measurement for pre-processing
    call stat_startTimer(p_rtimerPrePostprocess, STAT_TIMERSHORT)

    ! Adjust time stepping scheme
    rtimestep%ctimestepType = TSTEP_THETA_SCHEME
    rtimestep%dinitialTime  = 0.0_DP
    rtimestep%dinitialStep  = 1.0_DP
    rtimestep%dfinalTime    = 1.0_DP
    rtimestep%theta         = 1.0_DP

    ! Get global parameters
    call parlst_getvalue_int(rparlist, ssectionName, 'irhstype', irhstype)
    call parlst_getvalue_int(rparlist, ssectionName, 'ivelocitytype', ivelocitytype)
    call parlst_getvalue_int(rparlist, ssectionName, 'discretisation', discretisation)

    ! Set pointer to maximum problem level
    p_rproblemLevel   => rproblem%p_rproblemLevelMax
    p_rdiscretisation => p_rproblemLevel%Rdiscretisation(discretisation)

    ! Initialize the solution vector and impose boundary conditions explicitly
    call lsysbl_createVectorBlock(p_rdiscretisation, rsolutionPrimal,&
        .false., ST_DOUBLE)
    call transp_initSolution(rparlist, ssectionName, p_rproblemLevel,&
        0.0_DP, rsolutionPrimal, rcollection)
    call bdrf_filterVectorExplicit(rbdrCondPrimal, rsolutionPrimal, 0.0_DP)

    !---------------------------------------------------------------------------
    ! Initialize the h-adaptation structure
    !---------------------------------------------------------------------------

    call parlst_getvalue_string(rparlist, ssectionName, 'adaptivity', sadaptivityName, '')
    if (trim(adjustl(sadaptivityName)) .ne. '') then

      call parlst_getvalue_int(rparlist, trim(sadaptivityName), 'nadapt', nadapt)

      if (nadapt > 0) then

        ! Initialize adaptation structure from parameter list
        call hadapt_initFromParameterlist(rhadapt, rparlist, sadaptivityName)

        ! Generate a dynamic graph for the sparsity pattern and attach
        ! it to the collection structure which is used to pass this type
        ! to the callback function for h-adaptation
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'templateMatrix', templateMatrix)
        call grph_createGraphFromMatrix(p_rproblemLevel&
            %Rmatrix(templateMatrix), rgraph)
        call collct_setvalue_graph(rcollection, 'sparsitypattern',&
            rgraph, .true.)

        ! Attach the primal solution vector to the collection structure
        rcollection%p_rvectorQuickAccess1 => rsolutionPrimal

      end if   ! nadapt > 0

    else

      ! No h-adaptation
      nadapt = 0

    end if

    ! Stop time measurement for pre-processing
    call stat_stopTimer(p_rtimerPrePostprocess)


    adaptloop: do iadapt = 0, nadapt

      !-------------------------------------------------------------------------
      ! Compute steady-state solution for the primal problem
      !-------------------------------------------------------------------------

      ! Start time measurement for solution procedure
      call stat_startTimer(p_rtimerSolution, STAT_TIMERSHORT)

      ! Calculate the velocity field
      nlmin = solver_getMinimumMultigridlevel(rsolver)
      call transp_calcVelocityField(rparlist, ssectionName,&
          p_rproblemLevel, rtimestep%dTime, rcollection, nlmin)

      ! Attach the boundary condition
      call solver_setBoundaryCondition(rsolver, rbdrCondPrimal, .true.)

      ! Set primal problem mode
      call parlst_addvalue(rparlist, ssectionName, 'mode', 'primal')

      ! Reset the time-stepping algorithm
      call tstep_resetTimestep(rtimestep, .false.)
      call solver_resetSolver(rsolver, .false.)

      ! Check if right-hand side vector exists
      if (irhstype > 0) then
        call lsysbl_createVectorBlock(rsolutionPrimal, rrhs)
        call transp_initRHS(rparlist, ssectionName, p_rproblemLevel,&
            0.0_DP, rrhs, rcollection)

        ! Prepare quick access arrays of the collection
        rcollection%SquickAccess(1) = ssectionName

        ! Solve the primal problem with non-zero right-hand side
        call tstep_performThetaStep(p_rproblemLevel, rtimestep,&
            rsolver, rsolutionPrimal, transp_nlsolverCallback,&
            rcollection, rrhs)

        ! Release right-hand side vector
        call lsysbl_releaseVector(rrhs)

      else

        ! Prepare quick access arrays of the collection
        rcollection%SquickAccess(1) = ssectionName

        ! Solve the primal problem without right-hand side
        call tstep_performThetaStep(p_rproblemLevel, rtimestep,&
            rsolver, rsolutionPrimal, transp_nlsolverCallback,&
            rcollection)
      end if

      ! Stop time measurement for solution procedure
      call stat_stopTimer(p_rtimerSolution)


      !-------------------------------------------------------------------------
      ! Compute the right-hand side for the dual problem
      !-------------------------------------------------------------------------

      ! Start time measurement for error estimation
      call stat_startTimer(p_rtimerErrorEstimation, STAT_TIMERSHORT)

      ! Initialize target functional
      call lsysbl_createVectorBlock(rsolutionPrimal, rtargetFunc)
      call transp_initTargetFunc(rparlist, ssectionName,&
          p_rproblemLevel, 0.0_DP, rtargetFunc, rcollection)

      ! Stop time measurement for error estimation
      call stat_stopTimer(p_rtimerErrorEstimation)


      !-------------------------------------------------------------------------
      ! Compute steady-state solution for the dual problem
      !-------------------------------------------------------------------------

      ! Start time measurement for solution procedure
      call stat_startTimer(p_rtimerSolution, STAT_TIMERSHORT)

      ! Attach the boundary condition
      call solver_setBoundaryCondition(rsolver, rbdrCondDual, .true.)

      ! Set dual problem mode
      call parlst_addvalue(rparlist, ssectionName, 'mode', 'dual')

      ! Create dual solution vector initialised by zeros
      call lsysbl_releaseVector(rsolutionDual)
      call lsysbl_createVectorBlock(rsolutionPrimal, rsolutionDual, .true.)

      ! Reset the time-stepping and solver algorithms
      call tstep_resetTimestep(rtimestep, .false.)
      call solver_resetSolver(rsolver, .false.)

      ! Prepare quick access arrays of the collection
      rcollection%SquickAccess(1) = ssectionName

      ! Solve the dual problem
      call tstep_performThetaStep(p_rproblemLevel, rtimestep, rsolver,&
          rsolutionDual, transp_nlsolverCallback, rcollection,&
          rtargetFunc)

      ! Release discretised target functional
      call lsysbl_releaseVector(rtargetFunc)

      ! Stop time measurement for solution procedure
      call stat_stopTimer(p_rtimerSolution)


      if (iadapt .eq. nadapt) exit adaptloop

      !-------------------------------------------------------------------------
      ! Perform goal-oriented error estimation
      !-------------------------------------------------------------------------

      ! Start time measurement for error estimation
      call stat_startTimer(p_rtimerErrorEstimation, STAT_TIMERSHORT)

      ! Calculate the velocity field
      nlmin = solver_getMinimumMultigridlevel(rsolver)
      call transp_calcVelocityField(rparlist, ssectionName,&
          p_rproblemLevel, rtimestep%dTime, rcollection, nlmin)

      ! Attach the boundary condition
      call solver_setBoundaryCondition(rsolver, rbdrCondPrimal, .true.)

      ! Set primal problem mode
      call parlst_addvalue(rparlist, ssectionName, 'mode', 'primal')

      ! Check if right-hand side vector exists
      if (irhstype > 0) then

        ! Initialize right-hand side vector
        call lsysbl_createVectorBlock(rsolutionPrimal, rrhs)
        call transp_initRHS(rparlist, ssectionName, p_rproblemLevel,&
            0.0_DP, rrhs, rcollection)

        ! Prepare quick access arrays of the collection
        rcollection%SquickAccess(1) = ssectionName

        ! Compute the error in the quantity of interest
        call transp_estimateTargetFuncError(rparlist, ssectionName,&
            p_rproblemLevel, rtimestep, rsolver, rsolutionPrimal,&
            rsolutionDual, rcollection, relementError, derror, rrhs)

        ! Release right-hand side vector
        call lsysbl_releaseVector(rrhs)

      else

        ! Prepare quick access arrays of the collection
        rcollection%SquickAccess(1) = ssectionName

        ! Compute the error in the quantity of interest
        call transp_estimateTargetFuncError(rparlist, ssectionName,&
            p_rproblemLevel, rtimestep, rsolver, rsolutionPrimal,&
            rsolutionDual, rcollection, relementError, derror)

      end if

      ! Stop time measurement for error estimation
      call stat_stopTimer(p_rtimerErrorEstimation)


      !-------------------------------------------------------------------------
      ! Perform h-adaptation
      !-------------------------------------------------------------------------

      ! Start time measurement for mesh adaptation
      call stat_startTimer(p_rtimerAdaptation, STAT_TIMERSHORT)

      ! Set the names of the template matrix and the solution vector
      rcollection%SquickAccess(1) = 'sparsitypattern'

      ! Attach the primal solution vector to the collection structure
      rcollection%p_rvectorQuickAccess1 => rsolutionPrimal

      ! Perform h-adaptation and update the triangulation structure
      call transp_adaptTriangulation(rhadapt, p_rproblemLevel&
          %rtriangulation, relementError, rcollection)

      ! Update the template matrix according to the sparsity pattern
      call grph_generateMatrix(rgraph, p_rproblemLevel&
          %Rmatrix(templateMatrix))

      ! Resize the solution vector accordingly
      call lsysbl_resizeVectorBlock(rsolutionPrimal, p_rproblemLevel&
          %Rmatrix(templateMatrix)%NEQ, .false.)

      ! Release element-wise error distribution
      call lsyssc_releaseVector(relementError)

      ! Stop time measurement for mesh adaptation
      call stat_stopTimer(p_rtimerAdaptation)


      !-------------------------------------------------------------------------
      ! Re-generate the discretisation and coefficient matrices
      !-------------------------------------------------------------------------

      ! Start time measurement for generation of the triangulation
      call stat_startTimer(p_rtimerTriangulation, STAT_TIMERSHORT)

      ! Generate standard mesh from raw mesh
      call tria_initStandardMeshFromRaw(&
          p_rproblemLevel%rtriangulation, rproblem%rboundary)

      ! Stop time measurement for generation of the triangulation
      call stat_stopTimer(p_rtimerTriangulation)

      !-------------------------------------------------------------------

      ! Start time measurement for generation of constant coefficient matrices
      call stat_startTimer(p_rtimerAssemblyCoeff, STAT_TIMERSHORT)

      ! Re-initialize all constant coefficient matrices
      call transp_initProblemLevel(rparlist, ssectionName, p_rproblemLevel, rcollection)

      ! Prepare internal data arrays of the solver structure
      call parlst_getvalue_int(rparlist, ssectionName, 'systemMatrix', systemMatrix)
      call flagship_updateSolverMatrix(p_rproblemLevel, rsolver,&
          systemMatrix, SYSTEM_INTERLEAVEFORMAT, UPDMAT_ALL)
      call solver_updateStructure(rsolver)

      ! Stop time measurement for generation of constant coefficient matrices
      call stat_stopTimer(p_rtimerAssemblyCoeff)

    end do adaptloop


    ! Release adaptation structure
    if (nadapt > 0) then
      call hadapt_releaseAdaptation(rhadapt)
      call grph_releaseGraph(rgraph)
    end if

  end subroutine transp_solveSteadyStatePrimalDual

  !*****************************************************************************
  ! AUXILIARY ROUTINES
  !*****************************************************************************

!<subroutine>

  subroutine transp_parseCmdlArguments(rparlist)

!<description>
    ! This subroutine parses the commandline arguments and modifies the
    ! parameter values in the global parameter list.
!</description>

!<inputoutput>
    ! parameter list
    type(t_parlist), intent(inout) :: rparlist
!</inputoutput>
!</subroutine>

    ! local variables
    character(LEN=SYS_STRLEN) :: cbuffer
    integer :: iarg, narg

    iarg = 1; narg = command_argument_count()

    cmdarg: do
      ! Retrieve next command line argument
      call get_command_argument(iarg,cbuffer)

      if ((trim(adjustl(cbuffer)) .eq. '-I') .or.&
          (trim(adjustl(cbuffer)) .eq. '--inviscid')) then
        
        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        call parlst_setvalue(rparlist, '', 'inviscid', trim(adjustl(cbuffer)))
        
      elseif ((trim(adjustl(cbuffer)) .eq. '-T') .or.&
              (trim(adjustl(cbuffer)) .eq. '--timestep')) then 
        
        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        call parlst_setvalue(rparlist, '', 'timestep', trim(adjustl(cbuffer)))

      elseif ((trim(adjustl(cbuffer)) .eq. '-S') .or.&
              (trim(adjustl(cbuffer)) .eq. '--solver')) then 

        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        call parlst_setvalue(rparlist, '', 'solver', trim(adjustl(cbuffer)))

      elseif ((trim(adjustl(cbuffer)) .eq. '-O') .or.&
              (trim(adjustl(cbuffer)) .eq. '--output')) then 

        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        call parlst_setvalue(rparlist, '', 'output', trim(adjustl(cbuffer)))

      elseif ((trim(adjustl(cbuffer)) .eq. '-E') .or.&
              (trim(adjustl(cbuffer)) .eq. '--errorestimator')) then 

        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        call parlst_setvalue(rparlist, '', 'errorestimator', trim(adjustl(cbuffer)))

      elseif ((trim(adjustl(cbuffer)) .eq. '-A') .or.&
              (trim(adjustl(cbuffer)) .eq. '--adaptivity')) then 

        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        call parlst_setvalue(rparlist, '', 'adaptivity', trim(adjustl(cbuffer)))

      else
        iarg = iarg+1
        if (iarg .ge. narg) exit cmdarg
      end if
    end do cmdarg

  end subroutine transp_parseCmdlArguments

  !*****************************************************************************

!<subroutine>

  subroutine transp_adjustParameterlist(rparlist, ssectionName)

!<description>
    ! This subroutine adjusts the content of the parameter list
    ! depending on internal data, i.e., the dimension
!</description>

!<inputoutput>
    ! parameter list
    type(t_parlist), intent(inout) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: ndimension
    integer :: imasstype
    integer :: imassantidiffusiontype
    integer :: idiffusiontype
    integer :: ivelocitytype


    ! Check if mass matrix needs to be built
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'imasstype', imasstype)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'imassantidiffusiontype', imassantidiffusiontype)

    if ((imasstype .eq. MASS_ZERO) .and. &
        (imassantidiffusiontype .eq. MASS_ZERO)) then
      call parlst_setvalue(rparlist,&
          ssectionName, 'ConsistentMassMatrix', '0')
      call parlst_setvalue(rparlist,&
          ssectionName, 'LumpedMassMatrix', '0')
    end if


    ! Check if diffusion matrix needs to be built
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'idiffusiontype', idiffusiontype)

    if (idiffusiontype .eq. DIFFUSION_ZERO) then
      call parlst_setvalue(rparlist,&
          ssectionName, 'CoeffMatrix_S', '0')
    end if


    ! Check if convection matrix needs to be build
    call parlst_getvalue_int(rparlist, ssectionName, 'ivelocitytype', ivelocitytype)

    if (ivelocitytype .eq. VELOCITY_ZERO) then
      call parlst_setvalue(rparlist,&
          ssectionName, 'CoeffMatrix_CX', '0')
      call parlst_setvalue(rparlist,&
          ssectionName, 'CoeffMatrix_CY', '0')
      call parlst_setvalue(rparlist,&
          ssectionName, 'CoeffMatrix_CZ', '0')
    else

      call parlst_getvalue_int(rparlist,&
          ssectionName, 'ndimension', ndimension)

      select case(ndimension)
      case (NDIM1D)
        call parlst_setvalue(rparlist,&
            ssectionName, 'CoeffMatrix_CY', '0')
        call parlst_setvalue(rparlist,&
            ssectionName, 'CoeffMatrix_CZ', '0')
      case (NDIM2D)
        call parlst_setvalue(rparlist,&
            ssectionName, 'CoeffMatrix_CZ', '0')
      case (NDIM3D)
        ! We actually need all three matrices
      case default
        call parlst_setvalue(rparlist,&
            ssectionName, 'CoeffMatrix_CX', '0')
        call parlst_setvalue(rparlist,&
            ssectionName, 'CoeffMatrix_CY', '0')
        call parlst_setvalue(rparlist,&
            ssectionName, 'CoeffMatrix_CZ', '0')
      end select

    end if

  end subroutine transp_adjustParameterlist

end module transport_application
