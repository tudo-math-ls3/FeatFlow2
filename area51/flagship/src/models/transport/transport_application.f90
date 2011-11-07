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
!# 3.) transp_initProblemDescriptor
!#     -> Initializes the abstract problem descriptor based on the
!#        parameter settings given by the parameter list.
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
!# 15.) transp_solveTransientPrimDual
!#      -> Solves the primal and the dual formulation of the time-dependent
!#         convection-diffusion-reaction equation using pseudo time-stepping.
!#
!# 16.) transp_solvePseudoTransPrimal
!#      -> Solves the primal formulation of the steady
!#         convection-diffusion-reaction equation using pseudo time-stepping.
!#
!# 17.) transp_solvePseudoTransPrimDual
!#      -> Solves the primal and the dual formulation of the steady
!#         convection-diffusion-reaction equation using pseudo time-stepping.
!#
!# 18.) transp_solveSteadyStatePrimal
!#      -> Solves the primal formulation of the steady
!#         convection-diffusion-reaction equation directly
!#
!# 19.) transp_solveSteadyStatePrimDual
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

  use afcstabbase
  use afcstabscalar
  use bilinearformevaluation
  use boundary
  use boundarycondaux
  use boundaryfilter
  use collection
  use cubature
  use derivatives
  use dofmapping
  use elementbase
  use element
  use flagship_basic
  use fparser
  use fsystem
  use genoutput
  use graph
  use groupfembase
  use groupfemscalar
  use hadaptaux
  use hadaptivity
  use linearformevaluation
  use lineariser
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
  use transport_basic
  use transport_callback
  use transport_callback1d
  use transport_callback2d
  use transport_callback3d
  use triangulation
  use ucd

  implicit none

  private
  public :: transp_app
  public :: transp_adaptTriangulation
  public :: transp_estimateRecoveryError
  public :: transp_estimateTargetFuncError
  public :: transp_initAllProblemLevels
  public :: transp_initProblemDescriptor
  public :: transp_initProblemLevel
  public :: transp_initRHS
  public :: transp_initSolution
  public :: transp_initSolvers
  public :: transp_initTargetFunc
  public :: transp_outputSolution
  public :: transp_outputStatistics
  public :: transp_solveTransientPrimal
  public :: transp_solveTransientPrimDual
  public :: transp_solvePseudoTransPrimal
  public :: transp_solvePseudoTransPrimDual
  public :: transp_solveSteadyStatePrimal
  public :: transp_solveSteadyStatePrimDual

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

    ! Abstract problem descriptor
    type(t_problemDescriptor) :: rproblemDescriptor

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

    ! Initialize global collection structure
    call collct_init(rcollection)

    ! Create a separate section for the scalar transport model
    call collct_addsection(rcollection, ssectionName)
    
    ! Define section name of this application
    call collct_setvalue_string(rcollection, 'ssectionName',&
        ssectionName, .true.)

    ! Attach the parameter list and the timers to the collection
    call collct_setvalue_parlst(rcollection,&
        'rparlist', rparlist, .true.,&
        ssectionName=ssectionName)
    call collct_setvalue_timer(rcollection,&
        'rtimerSolution', rtimerSolution, .true.,&
        ssectionName=ssectionName)
    call collct_setvalue_timer(rcollection,&
        'rtimerAdaptation', rtimerAdaptation, .true.,&
        ssectionName=ssectionName)
    call collct_setvalue_timer(rcollection,&
        'rtimerErrorEstimation', rtimerErrorEstimation, .true.,&
        ssectionName=ssectionName)
    call collct_setvalue_timer(rcollection,&
        'rtimerTriangulation', rtimerTriangulation, .true.,&
        ssectionName=ssectionName)
    call collct_setvalue_timer(rcollection,&
        'rtimerAssemblyCoeff', rtimerAssemblyCoeff, .true.,&
        ssectionName=ssectionName)
    call collct_setvalue_timer(rcollection,&
        'rtimerAssemblyMatrix', rtimerAssemblyMatrix, .true.,&
        ssectionName=ssectionName)
    call collct_setvalue_timer(rcollection,&
        'rtimerAssemblyVector', rtimerAssemblyVector, .true.,&
        ssectionName=ssectionName)
    call collct_setvalue_timer(rcollection,&
        'rtimerPrePostprocess', rtimerPrePostprocess, .true.,&
        ssectionName=ssectionName)

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
    call collct_setvalue_pars(rcollection, 'rfparser', rfparser, .true.,&
        ssectionName=ssectionName)

    ! Initialize the solver structures
    call transp_initSolvers(rparlist, ssectionName, rtimestep, rsolver)

    ! Initialize the boundary conditions for the primal problem
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'ndimension', ndimension)
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'sprimalbdrcondname', sbdrcondName)
    call bdrc_readBoundaryCondition(rbdrCondPrimal,&
        sindatfileName, '['//trim(sbdrcondName)//']',&
        ndimension, transp_parseBoundaryCondition)

    ! Initialize the boundary conditions for the dual problem
    call parlst_getvalue_string(rparlist, ssectionName,&
        'sdualbdrcondname', sbdrcondName, '')
    if (sbdrcondName .ne. '') then
      call bdrc_readBoundaryCondition(rbdrCondDual,&
          sindatfileName, '['//trim(sbdrcondName)//']',&
          ndimension, transp_parseBoundaryCondition)
    end if
    
    ! Initialize the abstract problem structure
    call transp_initProblemDescriptor(rparlist, ssectionName,&
        solver_getMinimumMultigridlevel(rsolver),&
        solver_getMaximumMultigridlevel(rsolver),&
        rproblemDescriptor)
    call problem_initProblem(rproblemDescriptor, rproblem)

    ! Initialize all individual problem levels with primal and dual
    ! boundary conditions (if available)
    if (sbdrcondName .ne. '') then
      call transp_initAllProblemLevels(rparlist,&
          ssectionName, rproblem, rcollection,&
          rbdrCondPrimal, rbdrCondDual)
    else
      call transp_initAllProblemLevels(rparlist,&
          ssectionName, rproblem, rcollection,&
          rbdrCondPrimal)
    end if
    
    ! Prepare internal data arrays of the solver structure
    call parlst_getvalue_int(rparlist, ssectionName,&
        'systemMatrix', systemMatrix)
    call flagship_updateSolverMatrix(rproblem%p_rproblemLevelMax,&
        rsolver, systemMatrix, SYSTEM_INTERLEAVEFORMAT, UPDMAT_ALL)
    call solver_updateStructure(rsolver)

    ! Stop time measurement for pre-processing
    call stat_stopTimer(rtimerPrePostprocess)

    !---------------------------------------------------------------------------
    ! Solution algorithm
    !---------------------------------------------------------------------------

    if (rtimestep%dfinalTime .ge. rtimestep%dinitialTime) then

      ! Get global configuration from parameter list
      call parlst_getvalue_string(rparlist,&
          ssectionName, 'algorithm', algorithm)
      
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

        call transp_solveTransientPrimDual(rparlist, ssectionName,&
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
        call transp_solvePseudoTransPrimal(rparlist, ssectionName,&
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
        call transp_solvePseudoTransPrimDual(rparlist, ssectionName,&
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
        call transp_solveSteadyStatePrimDual(rparlist, ssectionName,&
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
    call transp_outputStatistics(rtimerTotal, ssectionName, rcollection)

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
  
  subroutine transp_initProblemDescriptor(rparlist, ssectionName,&
      nlmin, nlmax, rproblemDescriptor)

!<description>
    ! This subroutine initializes the abstract problem descriptor
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
    ! problem descriptor
    type(t_problemDescriptor), intent(out) :: rproblemDescriptor
!</output>
!</subroutine>

    ! local variables
    integer :: discretisation,velocityfield,dofCoords
    integer :: massAFC,templateGFEM
    integer :: convectionAFC,convectionGFEM
    integer :: diffusionAFC,diffusionGFEM
    integer :: primalBdrGFEM,dualBdrGFEM
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
    integer :: iconvToTria
    
    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'trifile', rproblemDescriptor%trifile)
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'prmfile', rproblemDescriptor%prmfile, '')
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'ndimension', rproblemDescriptor%ndimension)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'iconvtotria', iconvToTria, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'discretisation', discretisation, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'dofCoords', dofCoords, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'velocityfield', velocityfield, 0)

    ! Get global positions of matrices
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'templateMatrix', templateMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'systemMatrix', systemMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'jacobianMatrix', jacobianMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'transportMatrix', transportMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'consistentMassMatrix', consistentMassMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'lumpedMassMatrix', lumpedMassMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffMatrix_CX', coeffMatrix_CX, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffMatrix_CY', coeffMatrix_CY, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffMatrix_CZ', coeffMatrix_CZ, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffMatrix_S', coeffMatrix_S, 0)

    ! Default is no stabilization
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'massAFC', massAFC, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'convectionAFC', convectionAFC, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'diffusionAFC', diffusionAFC, 0)

    ! By default the same identifier is used for the group finite
    ! element formulation and the stabilization structure
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'templateGFEM', templateGFEM, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'convectionGFEM', convectionGFEM, convectionAFC)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'diffusionGFEM', diffusionGFEM, diffusionAFC)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'primalbdrGFEM', primalBdrGFEM, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'dualbdrGFEM', dualBdrGFEM, 0)

    ! Consistency check
    if (templateGFEM .eq. 0) then
      convectionGFEM = 0
      diffusionGFEM  = 0
      primalBdrGFEM  = 0
      dualBdrGFEM    = 0
    end if

    ! Set additional problem descriptor
    rproblemDescriptor%ndiscretisation = max(0, discretisation)
    rproblemDescriptor%nafcstab        = max(0, massAFC,&
                                                convectionAFC,&
                                                diffusionAFC)
    rproblemDescriptor%ngroupfemBlock  = max(0, templateGFEM,&
                                                convectionGFEM,&
                                                diffusionGFEM,&
                                                primalBdrGFEM,&
                                                dualBdrGFEM)
    rproblemDescriptor%nlmin           = nlmin
    rproblemDescriptor%nlmax           = nlmax
    rproblemDescriptor%nmatrixScalar   = max(0, templateMatrix,&
                                                systemMatrix,&
                                                jacobianMatrix,&
                                                transportMatrix,&
                                                consistentMassMatrix,&
                                                lumpedMassMatrix,&
                                                coeffMatrix_CX,&
                                                coeffMatrix_CY,&
                                                coeffMatrix_CZ,&
                                                coeffMatrix_S)
    rproblemDescriptor%nmatrixBlock    = 0
    rproblemDescriptor%nvectorScalar   = 0
    rproblemDescriptor%nvectorBlock    = max(0, velocityfield, dofCoords)

    ! Check if quadrilaterals should be converted to triangles
    if (iconvToTria .ne. 0) then
      rproblemDescriptor%iproblemSpec = rproblemDescriptor%iproblemSpec &
                                      + PROBDESC_MSPEC_CONVTRIANGLES
    end if

  end subroutine transp_initProblemDescriptor

  !*****************************************************************************

!<subroutine>

  subroutine transp_initProblemLevel(rparlist, ssectionName,&
      rproblemLevel, rcollection, rbdrCondPrimal, rbdrCondDual)

!<description>
    ! This subroutine initielizes the individual problem level. It
    ! generates the discretisation, the template matrix and the
    ! coefficient matrices as duplicates of the template matrix.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: boundary condition for primal problem
    type(t_boundaryCondition), intent(in), optional :: rbdrCondPrimal

    ! OPTIONAL: boundary condition for dual problem
    type(t_boundaryCondition), intent(in), optional :: rbdrCondDual
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
    integer :: massAFC,templateGFEM
    integer :: convectionAFC,convectionGFEM
    integer :: diffusionAFC,diffusionGFEM
    integer :: discretisation
    integer :: dofCoords
    integer :: ijacobianFormat
    integer :: imatrixFormat
    integer :: primalbdrGFEM,dualbdrGFEM

    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    type(t_triangulation) , pointer :: p_rtriangulation
    type(t_groupFEMSet), pointer :: p_rgroupFEMSet
    type(t_boundary) , pointer :: p_rboundary
    type(t_fparser), pointer :: p_rfparser
    type(t_boundaryRegion) :: rboundaryRegion
    type(t_matrixScalar) :: rmatrixSX,rmatrixSY,rmatrixSZ
    integer, dimension(:), allocatable :: Celement
    integer, dimension(:), pointer :: p_IbdrCondCpIdx
    integer :: nsumcubRefBilForm,nsumcubRefLinForm,nsumcubRefEval
    integer :: ccubType,i,ibdc,isegment,j,nmatrices,nsubstrings,neq
    character(len=SYS_STRLEN) :: selemName,smass,sconvection,sdiffusion

    ! Retrieve application specific parameters from the parameter list
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'templatematrix', templateMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'systemmatrix', systemMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'jacobianmatrix', jacobianMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'transportmatrix', transportMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'consistentmassmatrix', consistentMassMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'lumpedmassmatrix', lumpedMassMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_s', coeffMatrix_S, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_cx', coeffMatrix_CX, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_cy', coeffMatrix_CY, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_cz', coeffMatrix_CZ, 0)

    ! Default is no stabilization
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'massAFC', massAFC, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'convectionAFC', convectionAFC, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'diffusionAFC', diffusionAFC, 0)

    ! By default the same identifier is used for the group finite
    ! element formulation and the stabilization structure
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'templateGFEM', templateGFEM, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'convectionGFEM', convectionGFEM, convectionAFC)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'diffusionGFEM', diffusionGFEM, diffusionAFC)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'primalBdrGFEM', primalBdrGFEM, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'dualBdrGFEM', dualBdrGFEM, 0)

    ! Default no summed cubature
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'discretisation', discretisation, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'dofCoords', dofCoords, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'nsumcubRefBilForm', nsumcubRefBilForm, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'nsumcubRefLinForm', nsumcubRefLinForm, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'nsumcubRefEval', nsumcubRefEval, 0)

    ! Default is empty section, i.e. no configuration
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'mass', smass, '')
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'diffusion', sdiffusion, '')
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'convection', sconvection, '')

    ! Set pointers to triangulation and boundary structure
    p_rtriangulation => rproblemLevel%rtriangulation
    p_rboundary      => rproblemLevel%p_rproblem%rboundary

    !---------------------------------------------------------------------------
    ! Create discretisation structure
    !---------------------------------------------------------------------------
    if (discretisation > 0) then

      ! Initialize the discretisation structure
      p_rdiscretisation => rproblemLevel%Rdiscretisation(discretisation)
      if (p_rdiscretisation%ndimension .eq. 0) then
        call spdiscr_initBlockDiscr(p_rdiscretisation, 1, p_rtriangulation)
      end if

      ! Allocate temporal memory
      nsubstrings = max(1, parlst_querysubstrings(rparlist,&
                           ssectionName, 'celement'))
      allocate(Celement(nsubstrings))

      ! Get IDs of element types
      do i = 1, nsubstrings
        call parlst_getvalue_string(rparlist,&
            ssectionName, 'celement', selemName, isubstring=i)
        Celement(i) = elem_igetID(selemName)
      end do
      
      ! Get spatial dimension
      select case(p_rdiscretisation%ndimension)
      case (NDIM1D)
        call spdiscr_initDiscr_simple(&
            p_rdiscretisation%RspatialDiscr(1),&
            Celement(1), SPDISC_CUB_AUTOMATIC,&
            p_rtriangulation, p_rboundary)

      case (NDIM2D)
        if (size(Celement) .eq. 1) then
          call spdiscr_initDiscr_simple(&
              p_rdiscretisation%RspatialDiscr(1),&
              Celement(1), SPDISC_CUB_AUTOMATIC,&
              p_rtriangulation, p_rboundary)
        else
          call spdiscr_initDiscr_triquad(&
              p_rdiscretisation%RspatialDiscr(1), Celement(1), Celement(2),&
              SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC,&
              p_rtriangulation, p_rboundary)
        end if

      case (NDIM3D)
        call spdiscr_initDiscr_simple(&
            p_rdiscretisation%RspatialDiscr(1),&
            Celement(1), SPDISC_CUB_AUTOMATIC,&
            p_rtriangulation, p_rboundary)
      
      case default
        call output_line('Invalid number of spatial dimensions',&
            OU_CLASS_ERROR,OU_MODE_STD,'transp_initProblemLevel')
        call sys_halt()
      end select
      
      ! Deallocate temporal memory
      deallocate(Celement)

      !-------------------------------------------------------------------------
      ! Configure evaluation of (bi-)linear formes
      !-------------------------------------------------------------------------

      if (parlst_queryvalue(rparlist, ssectionName, 'ccubTypeBilForm') .ne. 0) then
        ! Check if special cubature formula for evaluating integral
        ! terms of the bilinear form are requested by the user
        nsubstrings = max(1, parlst_querysubstrings(rparlist,&
                             ssectionName, 'ccubTypeBilForm'))
        if (nsubstrings .ne. 0) then
          if (nsubstrings .eq. p_rdiscretisation%RspatialDiscr(1)%inumFESpaces) then
            do i = 1, nsubstrings
              call parlst_getvalue_int(rparlist,&
                  ssectionName, 'ccubTypeBilForm', ccubType, 0, i)
              if (ccubType .ne. 0)&
                  p_rdiscretisation%RspatialDiscr(1)%RelementDistr(i)%ccubTypeBilForm = ccubType
            end do
          else
            call output_line('Number of substrings does not match number of FE spaces!',&
                OU_CLASS_ERROR,OU_MODE_STD,'transp_initProblemLevel')
            call sys_halt()
          end if
        end if
      end if
      
      if (parlst_queryvalue(rparlist, ssectionName, 'ccubTypeLinForm') .ne. 0) then
        ! Check if special cubature formula for evaluating integral
        ! terms of the linear form are requested by the user
        nsubstrings = max(1, parlst_querysubstrings(rparlist,&
                             ssectionName, 'ccubTypeLinForm'))
        if (nsubstrings .ne. 0) then
          if (nsubstrings .eq. p_rdiscretisation%RspatialDiscr(1)%inumFESpaces) then
            do i = 1, nsubstrings
              call parlst_getvalue_int(rparlist,&
                  ssectionName, 'ccubTypeLinForm', ccubType, 0, i)
              if (ccubType .ne. 0)&
                  p_rdiscretisation%RspatialDiscr(1)%RelementDistr(i)%ccubTypeLinForm = ccubType
            end do
          else
            call output_line('Number of substrings does not match number of FE spaces!',&
                OU_CLASS_ERROR,OU_MODE_STD,'transp_initProblemLevel')
            call sys_halt()
          end if
        end if
      end if

      if (parlst_queryvalue(rparlist, ssectionName, 'ccubTypeEval') .ne. 0) then
        ! Check if special cubature formula for evaluating integral
        ! terms of allother terms are requested by the user
        nsubstrings = max(1, parlst_querysubstrings(rparlist,&
                             ssectionName, 'ccubTypeEval'))
        if (nsubstrings .ne. 0) then
          if (nsubstrings .eq. p_rdiscretisation%RspatialDiscr(1)%inumFESpaces) then
            do i = 1, nsubstrings
              call parlst_getvalue_int(rparlist,&
                  ssectionName, 'ccubTypeEval', ccubType, 0, i)
              if (ccubType .ne. 0)&
                  p_rdiscretisation%RspatialDiscr(1)%RelementDistr(i)%ccubTypeEval = ccubType
            end do
          else
            call output_line('Number of substrings does not match number of FE spaces!',&
                OU_CLASS_ERROR,OU_MODE_STD,'transp_initProblemLevel')
            call sys_halt()
          end if
        end if
      end if

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
      
      !-------------------------------------------------------------------------
      ! Calculate coordinates of the global DOF`s
      if (dofCoords > 0) then
        ! Check if block vector has been initialized
        if (rproblemLevel%RvectorBlock(dofCoords)%nblocks .eq. 0) then
          call lsysbl_createVectorBlock(p_rdiscretisation,&
              p_rdiscretisation%ndimension,&
              rproblemLevel%RvectorBlock(dofCoords), .false.)
        else
          neq = dof_igetNDofGlobBlock(p_rdiscretisation)
          if (rproblemLevel%RvectorBlock(dofCoords)%NEQ .ne. neq) then
            call lsysbl_resizeVectorBlock(p_rdiscretisation,&
                rproblemLevel%RvectorBlock(dofCoords), .false.)
          end if
        end if
        ! Calculate coordinates
        call lin_calcDofCoordsBlock(p_rdiscretisation,&
            rproblemLevel%RvectorBlock(dofCoords))
      end if

      !-------------------------------------------------------------------------
      ! Create finite element matrices
      !-------------------------------------------------------------------------
      
      ! If the template matrix has no structure data then generate the
      ! finite element matrix sparsity structure based on the spatial
      ! descretisation and store it as the template matrix. Otherwise we
      ! assume that the template matrix has been generated externally.
      if (.not.lsyssc_hasMatrixStructure(rproblemLevel%Rmatrix(templateMatrix))) then
        call parlst_getvalue_int(rparlist, ssectionName, 'imatrixFormat', imatrixFormat)
        call bilf_createMatrixStructure(p_rdiscretisation%RspatialDiscr(1),&
            imatrixFormat, rproblemLevel%Rmatrix(templateMatrix))
      end if
      
      !-------------------------------------------------------------------------
      ! Create system matrix as duplicate of the template matrix
      if (systemMatrix > 0)&
          call initMatrixStructure(rproblemLevel%Rmatrix(templateMatrix),&
                                   rproblemLevel%Rmatrix(systemMatrix))

      !-------------------------------------------------------------------------
      ! Create transport matrix as duplicate of the template matrix
      if (transportMatrix > 0)&
          call initMatrixStructure(rproblemLevel%Rmatrix(templateMatrix),&
                                   rproblemLevel%Rmatrix(transportMatrix))

      !-------------------------------------------------------------------------
      ! Create Jacobian matrix. This is a little bit tricky. If the
      ! Jacobian matrix has the same sparsity pattern as the template
      ! matrix, we can just create the Jacobian matrix as a duplicate of
      ! the template matrix. If the Jacobian matrix has an extended
      ! sparsity pattern we must create it by using the template matrix
      if (jacobianMatrix > 0) then
        
        ! What format do we have for the Jacobian matrix?
        call parlst_getvalue_int(rparlist, ssectionName,&
            'ijacobianFormat', ijacobianFormat)
        
        if (lsyssc_hasMatrixStructure(rproblemLevel%Rmatrix(jacobianMatrix))) then
          if (ijacobianFormat .eq. 0) then
            call lsyssc_resizeMatrix(&
                rproblemLevel%Rmatrix(jacobianMatrix),&
                rproblemLevel%Rmatrix(templateMatrix),&
                .false., .false., bforce=.true.)
          else
            call afcstab_genExtSparsity(&
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
            call afcstab_genExtSparsity(&
                rproblemLevel%Rmatrix(templateMatrix),&
                rproblemLevel%Rmatrix(jacobianMatrix))
          end if
        end if
      end if

      !-------------------------------------------------------------------------
      ! Create consistent mass matrix as duplicate of the template matrix
      if (consistentMassMatrix > 0) then
        call initMatrixStructure(rproblemLevel%Rmatrix(templateMatrix),&
                                 rproblemLevel%Rmatrix(consistentMassMatrix))
        call stdop_assembleSimpleMatrix(&
            rproblemLevel%Rmatrix(consistentMassMatrix), DER_FUNC, DER_FUNC)
        
        ! Create lumped mass matrix
        if (lumpedMassMatrix > 0) then
          call lsyssc_duplicateMatrix(&
              rproblemLevel%Rmatrix(consistentMassMatrix),&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)
          call lsyssc_lumpMatrixScalar(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), LSYSSC_LUMP_DIAG)
        end if
      elseif (lumpedMassMatrix > 0) then
        ! Create lumped mass matrix
        call lsyssc_duplicateMatrix(&
            rproblemLevel%Rmatrix(templateMatrix),&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
        call stdop_assembleSimpleMatrix(&
            rproblemLevel%Rmatrix(lumpedMassMatrix), DER_FUNC, DER_FUNC)
        call lsyssc_lumpMatrixScalar(&
            rproblemLevel%Rmatrix(lumpedMassMatrix), LSYSSC_LUMP_DIAG)
      end if
      
      !-------------------------------------------------------------------------
      ! Create diffusion matrix as duplicate of the template matrix
      if (coeffMatrix_S > 0) then
        call initMatrixStructure(rproblemLevel%Rmatrix(templateMatrix),&
                                 rproblemLevel%Rmatrix(coeffMatrix_S))
        
        ! Get function parser from collection
        p_rfparser => collct_getvalue_pars(rcollection,&
            'rfparser', ssectionName=ssectionName)
        
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
        case default
          call lsyssc_releaseMatrix(rproblemLevel%Rmatrix(coeffMatrix_S))
        end select
      end if
      
      !-------------------------------------------------------------------------
      ! Create coefficient matrix (phi, dphi/dx) as duplicate of the
      ! template matrix
      if (coeffMatrix_CX > 0) then
        call initMatrixStructure(rproblemLevel%Rmatrix(templateMatrix),&
                                 rproblemLevel%Rmatrix(coeffMatrix_CX))
        call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(coeffMatrix_CX),&
                                       DER_DERIV3D_X, DER_FUNC)
      end if

      !-------------------------------------------------------------------------
      ! Create coefficient matrix (phi, dphi/dy) as duplicate of the
      ! template matrix
      if (coeffMatrix_CY > 0) then
        call initMatrixStructure(rproblemLevel%Rmatrix(templateMatrix),&
                                 rproblemLevel%Rmatrix(coeffMatrix_CY))
        call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(coeffMatrix_CY),&
                                        DER_DERIV3D_Y, DER_FUNC)
      end if
      
      !-------------------------------------------------------------------------
      ! Create coefficient matrix (phi, dphi/dz) as duplicate of the
      ! template matrix
      if (coeffMatrix_CZ > 0) then
        call initMatrixStructure(rproblemLevel%Rmatrix(templateMatrix),&
                                 rproblemLevel%Rmatrix(coeffMatrix_CZ))
        call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(coeffMatrix_CZ),&
                                        DER_DERIV3D_Z, DER_FUNC)
      end if
      
      !-------------------------------------------------------------------------
      ! Create group finite element structures and AFC-stabilisations
      !-------------------------------------------------------------------------
      
      ! Initialize/resize template group finite elment structure and
      ! generate the edge structure derived from the template matrix
      if (templateGFEM > 0) then
        ! Check if structure has been initialized
        if (rproblemLevel%RgroupFEMBlock(templateGFEM)%nblocks .eq. 0)&
            call gfem_initGroupFEMBlock(rproblemLevel%RgroupFEMBlock(templateGFEM), 1)
        
        ! Set pointer to first group finite element set of this block
        p_rgroupFEMSet =>&
            rproblemLevel%RgroupFEMBlock(templateGFEM)%RgroupFEMBlock(1)
        
        if (p_rgroupFEMSet%isetSpec .eq. GFEM_UNDEFINED) then
          ! Initialize first group finite element set for edge-based assembly
          call gfem_initGroupFEMSet(p_rgroupFEMSet,&
              rproblemLevel%Rmatrix(templateMatrix), 0, 0, 0, GFEM_EDGEBASED)
        else
          ! Resize first group finite element set
          call gfem_resizeGroupFEMSet(p_rgroupFEMSet,&
              rproblemLevel%Rmatrix(templateMatrix))
        end if
        
        ! Generate diagonal and edge structure derived from template matrix
        call gfem_genDiagList(rproblemLevel%Rmatrix(templateMatrix),&
            p_rgroupFEMSet)
        call gfem_genEdgeList(rproblemLevel%Rmatrix(templateMatrix),&
            p_rgroupFEMSet)
      else
        convectionGFEM = 0
        diffusionGFEM  = 0
        primalBdrGFEM  = 0
        dualBdrGFEM    = 0
      end if
      
      !-------------------------------------------------------------------------
      ! Initialize/resize group finite element structure as duplicate of
      ! the template group finite element structure and fill it with the
      ! precomputed matrix coefficients for the convective term
      if (convectionGFEM > 0) then
        ! Check if structure has been initialized
        if (rproblemLevel%RgroupFEMBlock(convectionGFEM)%nblocks .eq. 0)&
            call gfem_initGroupFEMBlock(rproblemLevel%RgroupFEMBlock(convectionGFEM), 1)
        
        ! Set pointer to first group finite element set of this block
        p_rgroupFEMSet =>&
            rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1)
        
        if (p_rgroupFEMSet%isetSpec .eq. GFEM_UNDEFINED) then
          ! Initialize first group finite element set for edge-based
          ! assembly as aduplicate of the template structure
          call gfem_duplicateGroupFEMSet(&
              rproblemLevel%RgroupFEMBlock(templateGFEM)%RgroupFEMBlock(1),&
              p_rgroupFEMSet, GFEM_DUP_STRUCTURE, .false.)
          
          ! Compute number of matrices to be copied
          nmatrices = 0
          if (coeffMatrix_CX > 0) nmatrices = nmatrices+1
          if (coeffMatrix_CY > 0) nmatrices = nmatrices+1
          if (coeffMatrix_CZ > 0) nmatrices = nmatrices+1
          
          ! Allocate memory for matrix entries
          call gfem_allocCoeffs(p_rgroupFEMSet, nmatrices, 0, nmatrices)
        else
          ! Resize first group finite element set
          call gfem_resizeGroupFEMSet(p_rgroupFEMSet,&
              rproblemLevel%Rmatrix(templateMatrix))
        end if
        
        ! Duplicate edge-based structure from template
        call gfem_duplicateGroupFEMSet(&
            rproblemLevel%RgroupFEMBlock(templateGFEM)%RgroupFEMBlock(1),&
            p_rgroupFEMSet, GFEM_DUP_DIAGLIST+GFEM_DUP_EDGELIST, .true.)
        
        ! Copy constant coefficient matrices to group finite element set
        nmatrices = 0
        if (coeffMatrix_CX > 0) then
          nmatrices = nmatrices+1
          call gfem_initCoeffsFromMatrix(p_rgroupFEMSet,&
              rproblemLevel%Rmatrix(coeffMatrix_CX), nmatrices)
        end if
        if (coeffMatrix_CY > 0) then
          nmatrices = nmatrices+1
          call gfem_initCoeffsFromMatrix(p_rgroupFEMSet,&
              rproblemLevel%Rmatrix(coeffMatrix_CY), nmatrices)
        end if
        if (coeffMatrix_CZ > 0) then
          nmatrices = nmatrices+1
          call gfem_initCoeffsFromMatrix(p_rgroupFEMSet,&
              rproblemLevel%Rmatrix(coeffMatrix_CZ), nmatrices)
        end if
      end if
      
      !-------------------------------------------------------------------------
      ! Initialise/Resize stabilisation structure for the convective
      ! term by duplicating parts of the template group finite element set
      if (convectionAFC > 0) then
        if (rproblemLevel%Rafcstab(convectionAFC)%istabilisationSpec&
            .eq. AFCSTAB_UNDEFINED) then
          call afcstab_initFromParameterlist(rparlist, sconvection,&
              rproblemLevel%Rafcstab(convectionAFC))
          call afcsc_initStabilisation(&
              rproblemLevel%RgroupFEMBlock(templateGFEM)%RgroupFEMBlock(1),&
              rproblemLevel%Rafcstab(convectionAFC), p_rdiscretisation)
        else
          call afcstab_resizeStabilisation(rproblemLevel%Rafcstab(convectionAFC),&
              rproblemLevel%RgroupFEMBlock(templateGFEM)%RgroupFEMBlock(1))
        end if
      end if
      
      !-------------------------------------------------------------------------
      ! Initialize/resize group finite element structure as duplicate of
      ! the template group finite element structure and fill it with the
      ! precomputed matrix coefficients for the diffusice term
      if (diffusionGFEM > 0) then
        ! Check if structure has been initialized
        if (rproblemLevel%RgroupFEMBlock(diffusionGFEM)%nblocks .eq. 0)&
            call gfem_initGroupFEMBlock(rproblemLevel%RgroupFEMBlock(diffusionGFEM), 1)
        
        ! Set pointer to first group finite element set of this block
        p_rgroupFEMSet =>&
            rproblemLevel%RgroupFEMBlock(diffusionGFEM)%RgroupFEMBlock(1)
        
        if (p_rgroupFEMSet%isetSpec .eq. GFEM_UNDEFINED) then
          ! Initialize first group finite element set for edge-based
          ! assembly as aduplicate of the template structure
          call gfem_duplicateGroupFEMSet(&
              rproblemLevel%RgroupFEMBlock(templateGFEM)%RgroupFEMBlock(1),&
              p_rgroupFEMSet, GFEM_DUP_STRUCTURE, .false.)
          
          ! Compute number of matrices to by copied
          nmatrices = 0
          if (coeffMatrix_S > 0) nmatrices = nmatrices+1
          
          ! Allocate memory for matrix entries
          call gfem_allocCoeffs(p_rgroupFEMSet, nmatrices, 0, nmatrices)
        else
          ! Resize first group finite element set
          call gfem_resizeGroupFEMSet(p_rgroupFEMSet,&
              rproblemLevel%Rmatrix(templateMatrix))
        end if
        
        ! Duplicate edge-based structure from template
        call gfem_duplicateGroupFEMSet(&
            rproblemLevel%RgroupFEMBlock(templateGFEM)%RgroupFEMBlock(1),&
            p_rgroupFEMSet, GFEM_DUP_DIAGLIST+GFEM_DUP_EDGELIST, .true.)
        
        ! Copy constant coefficient matrices to group finite element set
        nmatrices = 0
        if (coeffMatrix_S > 0) then
          nmatrices = nmatrices+1
          call gfem_initCoeffsFromMatrix(p_rgroupFEMSet,&
              rproblemLevel%Rmatrix(coeffMatrix_S), nmatrices)
        end if
      end if
      
      !-------------------------------------------------------------------------
      ! Initialise/Resize stabilisation structure for the diffusive
      ! term by duplicating parts of the template group finite element set
      if (diffusionAFC > 0) then
        if (rproblemLevel%Rafcstab(diffusionAFC)%istabilisationSpec&
            .eq. AFCSTAB_UNDEFINED) then
          call afcstab_initFromParameterlist(rparlist, sdiffusion,&
              rproblemLevel%Rafcstab(diffusionAFC))
          call afcsc_initStabilisation(&
              rproblemLevel%RgroupFEMBlock(templateGFEM)%RgroupFEMBlock(1),&
              rproblemLevel%Rafcstab(diffusionAFC), p_rdiscretisation)
        else
          call afcstab_resizeStabilisation(rproblemLevel%Rafcstab(diffusionAFC),&
              rproblemLevel%RgroupFEMBlock(templateGFEM)%RgroupFEMBlock(1))
        end if
      end if
      
      !-------------------------------------------------------------------------
      ! Initialize/Resize stabilisation structure for the mass matrix by
      ! duplicating parts of the template group finite element set
      if (massAFC > 0) then
        if (rproblemLevel%Rafcstab(massAFC)%istabilisationSpec&
            .eq. AFCSTAB_UNDEFINED) then
          call afcstab_initFromParameterlist(rparlist, smass,&
              rproblemLevel%Rafcstab(massAFC))
          call afcsc_initStabilisation(&
              rproblemLevel%RgroupFEMBlock(templateGFEM)%RgroupFEMBlock(1),&
              rproblemLevel%Rafcstab(massAFC), p_rdiscretisation)
        else
          call afcstab_resizeStabilisation(rproblemLevel%Rafcstab(massAFC),&
              rproblemLevel%RgroupFEMBlock(templateGFEM)%RgroupFEMBlock(1))
        end if
      end if
      
      !-------------------------------------------------------------------------
      ! Create group finite element structures on the boundary
      !-------------------------------------------------------------------------
      
      if ((primalBdrGFEM > 0) .and. present(rbdrCondPrimal) .or.&
          (dualBdrGFEM   > 0) .and. present(rbdrCondDual) ) then
        
        ! Compute X-component of mass matrices on the boundary
        if (coeffMatrix_CX > 0) then
          call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(coeffMatrix_CX),&
              rmatrixSX, LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
          call lsyssc_createMatrixSymmPart(rproblemLevel%Rmatrix(coeffMatrix_CX),&
              1.0_DP, rmatrixSX)
        end if
        
        ! Compute Y-component of mass matrices on the boundary
        if (coeffMatrix_CY > 0) then
          call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(coeffMatrix_CY),&
              rmatrixSY, LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
          call lsyssc_createMatrixSymmPart(rproblemLevel%Rmatrix(coeffMatrix_CY),&
              1.0_DP, rmatrixSY)
        end if
        
        ! Compute Z-component of mass matrices on the boundary
        if (coeffMatrix_CZ > 0) then
          call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(coeffMatrix_CZ),&
              rmatrixSZ, LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
          call lsyssc_createMatrixSymmPart(rproblemLevel%Rmatrix(coeffMatrix_CZ),&
              1.0_DP, rmatrixSZ)
        end if
      end if
      
      !-------------------------------------------------------------------------
      ! Primal boundary condition
      !-------------------------------------------------------------------------
      
      if ((primalBdrGFEM > 0) .and. present(rbdrCondPrimal)) then
        
        ! Check if structure has been initialized
        if (rproblemLevel%RgroupFEMBlock(primalBdrGFEM)%nblocks .eq. 0)&
            call gfem_initGroupFEMBlock(rproblemLevel%RgroupFEMBlock(primalBdrGFEM),&
            bdrc_getNumberOfRegions(rbdrCondPrimal))
        
        ! Compute number of matrices to be copied
        nmatrices = 0
        if (coeffMatrix_CX > 0) nmatrices = nmatrices+1
        if (coeffMatrix_CY > 0) nmatrices = nmatrices+1
        if (coeffMatrix_CZ > 0) nmatrices = nmatrices+1
        
        ! Set pointer
        call storage_getbase_int(rbdrCondPrimal%h_IbdrCondCpIdx, p_IbdrCondCpIdx)
        
        ! Loop over all boundary components
        do ibdc = 1, rbdrCondPrimal%iboundarycount
      
          ! Loop over all boundary segments
          do isegment = p_IbdrCondCpIdx(ibdc), p_IbdrCondCpIdx(ibdc+1)-1
            
            ! Create boundary region
            call bdrc_createRegion(rbdrCondPrimal, ibdc,&
                isegment-p_IbdrCondCpIdx(ibdc)+1, rboundaryRegion)
            
            ! Set pointer to group finite element set
            p_rgroupFEMSet =>&
                rproblemLevel%RgroupFEMBlock(primalBdrGFEM)%RgroupFEMBlock(isegment)
            
            ! Initialize group finite element
            call initGroupFEMSetBoundary(rboundaryRegion,&
                rproblemLevel%Rmatrix(templateMatrix), nmatrices, p_rgroupFEMSet)
            
            ! Copy constant coefficient matrices to group finite element set
            if (coeffMatrix_CX > 0)&
                call gfem_initCoeffsFromMatrix(p_rgroupFEMSet, rmatrixSX, 1)
            if (coeffMatrix_CY > 0)&
                call gfem_initCoeffsFromMatrix(p_rgroupFEMSet, rmatrixSY, 2)
            if (coeffMatrix_CZ > 0)&
                call gfem_initCoeffsFromMatrix(p_rgroupFEMSet, rmatrixSZ, 3)
          end do
        end do
      end if
      
      !-------------------------------------------------------------------------
      ! Dual boundary condition
      !-------------------------------------------------------------------------
      
      if ((dualBdrGFEM > 0) .and. present(rbdrCondDual)) then
        
        ! Check if structure has been initialized
        if (rproblemLevel%RgroupFEMBlock(dualBdrGFEM)%nblocks .eq. 0)&
            call gfem_initGroupFEMBlock(rproblemLevel%RgroupFEMBlock(dualBdrGFEM),&
            bdrc_getNumberOfRegions(rbdrCondDual))
        
        ! Compute number of matrices to be copied
        nmatrices = 0
        if (coeffMatrix_CX > 0) nmatrices = nmatrices+1
        if (coeffMatrix_CY > 0) nmatrices = nmatrices+1
        if (coeffMatrix_CZ > 0) nmatrices = nmatrices+1
        
        ! Set pointer
        call storage_getbase_int(rbdrCondDual%h_IbdrCondCpIdx, p_IbdrCondCpIdx)
        
        ! Loop over all boundary components
        do ibdc = 1, rbdrCondDual%iboundarycount
          
          ! Loop over all boundary segments
          do isegment = p_IbdrCondCpIdx(ibdc), p_IbdrCondCpIdx(ibdc+1)-1
            
            ! Create boundary region
            call bdrc_createRegion(rbdrCondDual, ibdc,&
                isegment-p_IbdrCondCpIdx(ibdc)+1, rboundaryRegion)
            
            ! Set pointer to group finite element set
            p_rgroupFEMSet =>&
                rproblemLevel%RgroupFEMBlock(dualBdrGFEM)%RgroupFEMBlock(isegment)
            
            ! Initialize group finite element
            call initGroupFEMSetBoundary(rboundaryRegion,&
                rproblemLevel%Rmatrix(templateMatrix), nmatrices, p_rgroupFEMSet)
            
            ! Copy constant coefficient matrices to group finite element set
            if (coeffMatrix_CX > 0)&
                call gfem_initCoeffsFromMatrix(p_rgroupFEMSet, rmatrixSX, 1)
            if (coeffMatrix_CY > 0)&
                call gfem_initCoeffsFromMatrix(p_rgroupFEMSet, rmatrixSY, 2)
            if (coeffMatrix_CZ > 0)&
                call gfem_initCoeffsFromMatrix(p_rgroupFEMSet, rmatrixSZ, 3)
          end do
        end do
      end if
      
      ! The auxiliary matrices are not needed any more
      if (coeffMatrix_CX > 0)&
          call lsyssc_releaseMatrix(rmatrixSX)
      if (coeffMatrix_CY > 0)&
          call lsyssc_releaseMatrix(rmatrixSY)
      if (coeffMatrix_CZ > 0)&
          call lsyssc_releaseMatrix(rmatrixSZ)

    end if   ! discretisation > 0

    !---------------------------------------------------------------------------
    ! Set update notifiers for the discrete transport operator and the
    ! preconditioned in the problem level structure
    !---------------------------------------------------------------------------
    rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                     TRANSP_TROPER_UPDATE)
    rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                     TRANSP_PRECOND_UPDATE)

  contains

    !**************************************************************
    ! Initialize the matrix structure by duplicating the template matrix

    subroutine initMatrixStructure(rmatrixTemplate, rmatrix)

      type(t_matrixScalar), intent(in) :: rmatrixTemplate
      type(t_matrixScalar), intent(inout) :: rmatrix

      if (lsyssc_isMatrixStructureShared(rmatrix, rmatrixTemplate)) then
        call lsyssc_resizeMatrix(rmatrix, rmatrixTemplate,&
            .false., .false., bforce=.true.)
      else
        call lsyssc_duplicateMatrix(rmatrixTemplate, rmatrix,&
            LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      end if

    end subroutine initMatrixStructure  
    
    !**************************************************************
    ! Initialize the diffusion matrix in 1D

    subroutine initDiffusionMatrix1D(rfparser, rmatrix)

      type(t_fparser), intent(in) :: rfparser
      type(t_matrixScalar), intent(inout) :: rmatrix

      ! local variables
      character(LEN=SYS_STRLEN) :: sdiffusionName
      real(DP), dimension(1) :: Dunity = (/1.0_DP/)
      real(DP) :: dalpha
      integer :: idiffusiontype

      ! Retrieve data from parameter list
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'idiffusiontype', idiffusiontype)

      select case(idiffusiontype)
      case (DIFFUSION_ISOTROPIC,&
            DIFFUSION_ANISOTROPIC)

        ! Retrieve name/number of expression describing the diffusion coefficient
        call parlst_getvalue_string(rparlist, ssectionName, 'sdiffusionname',&
                                    sdiffusionName, isubString=1)

        ! Evaluate the constant coefficient from the function parser
        call fparser_evalFunction(rfparser, sdiffusionName, Dunity, dalpha)

        ! Assemble the Laplace matrix multiplied by the negative value
        ! of the physical diffusion parameter alpha
        call stdop_assembleLaplaceMatrix(rmatrix, .true., -dalpha)

      case default
        call lsyssc_clearMatrix(rmatrix)
      end select

    end subroutine initDiffusionMatrix1D

    !**************************************************************
    ! Initialize the diffusion matrix in 2D

    subroutine initDiffusionMatrix2D(rfparser, rmatrix)

      type(t_fparser), intent(in) :: rfparser
      type(t_matrixScalar), intent(inout) :: rmatrix

      ! local variables
      type(t_bilinearform) :: rform
      character(LEN=SYS_STRLEN) :: sdiffusionName
      real(DP), dimension(1) :: Dunity = (/1.0_DP/)
      real(DP) :: dalpha
      integer :: i,idiffusiontype

      ! Retrieve data from parameter list
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'idiffusiontype', idiffusiontype)

      select case(idiffusiontype)
      case (DIFFUSION_ISOTROPIC)

        ! Retrieve name/number of expression describing the diffusion coefficient
        call parlst_getvalue_string(rparlist,&
            ssectionName, 'sdiffusionname', sdiffusionName, isubString=1)

        ! Evaluate the constant coefficient from the function parser
        call fparser_evalFunction(rfparser, sdiffusionName, Dunity, dalpha)

        ! Assemble the Laplace matrix multiplied by the negative value
        ! of the physical diffusion parameter alpha
        call stdop_assembleLaplaceMatrix(rmatrix, .true., -dalpha)

      case (DIFFUSION_ANISOTROPIC)

        do i = 1, 4
          ! Retrieve name/number of expression describing the diffusion coefficient
          call parlst_getvalue_string(rparlist,&
              ssectionName, 'sdiffusionname', sdiffusionName, isubString=i)

          ! Evaluate the constant coefficient from the function parser
          call fparser_evalFunction(rfparser, sdiffusionName, Dunity,&
              rform%Dcoefficients(i))
        end do

        ! We have constant coefficients
        rform%ballCoeffConstant   = .true.
        rform%BconstantCoeff(1:4) = .true.
        rform%Dcoefficients(1:4)  = -rform%Dcoefficients(1:4)

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

      case default
        call lsyssc_clearMatrix(rmatrix)
      end select

    end subroutine initDiffusionMatrix2D

    !**************************************************************
    ! Initialize the diffusion matrix in 3D

    subroutine initDiffusionMatrix3D(rfparser, rmatrix)

      type(t_fparser), intent(in) :: rfparser
      type(t_matrixScalar), intent(inout) :: rmatrix

      ! local variables
      type(t_bilinearform) :: rform
      character(LEN=SYS_STRLEN) :: sdiffusionName
      real(DP), dimension(1) :: Dunity = (/1.0_DP/)
      real(DP) :: dalpha
      integer :: i,idiffusiontype

      ! Retrieve data from parameter list
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'idiffusiontype', idiffusiontype)

      select case(idiffusiontype)
      case (DIFFUSION_ISOTROPIC)

        ! Retrieve name/number of expression describing the diffusion coefficient
        call parlst_getvalue_string(rparlist,&
            ssectionName, 'sdiffusionname', sdiffusionName, isubString=1)

        ! Evaluate the constant coefficient from the function parser
        call fparser_evalFunction(rfparser, sdiffusionName, Dunity, dalpha)

        ! Assemble the Laplace matrix multiplied by the negative value
        ! of the physical diffusion parameter alpha
        call stdop_assembleLaplaceMatrix(rmatrix, .true., -dalpha)

      case (DIFFUSION_ANISOTROPIC)

        do i = 1, 9
          ! Retrieve name/number of expression describing the diffusion coefficient
          call parlst_getvalue_string(rparlist,&
              ssectionName, 'sdiffusionname', sdiffusionName, isubString=i)

          ! Evaluate the constant coefficient from the function parser
          call fparser_evalFunction(rfparser, sdiffusionName, Dunity,&
              rform%Dcoefficients(i))
        end do

        ! We have constant coefficients
        rform%ballCoeffConstant   = .true.
        rform%BconstantCoeff(1:9) = .true.
        rform%Dcoefficients(1:9)  = -rform%Dcoefficients(1:9)

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

      case default
        call lsyssc_clearMatrix(rmatrix)
      end select

    end subroutine initDiffusionMatrix3D

    !**************************************************************
    ! Initialize the group finite element set for evaluating the
    ! bilinear and linear forms on the boundary node-by-node. 

    subroutine initGroupFEMSetBoundary(rregion, rmatrix, nmatrices, rgroupFEMSet)

      ! input parameters
      type(t_boundaryRegion), intent(in) :: rregion
      type(t_matrixScalar), intent(in) :: rmatrix
      integer, intent(in) :: nmatrices

      ! input/output parameters
      type(t_groupFEMSet), intent(inout) :: rgroupFEMSet
      
      if (rgroupFEMSet%isetSpec .eq. GFEM_UNDEFINED) then
        ! Initialize finite element set for node-based assembly
        call gfem_initGroupFEMSetBoundary(rgroupFEMSet, rmatrix,&
            0, 0, 0, GFEM_NODEBASED, rregionTest=rregion,&
            brestrictToBoundary=.true.)
        
        ! Allocate memory for matrix entries
        call gfem_allocCoeffs(rgroupFEMSet, 0, nmatrices,0)
      else
        ! Resize first group finite element set
        call gfem_resizeGroupFEMSetBoundary(rgroupFEMSet, rmatrix,&
            rregionTest=rregion, rregionTrial=rregion)
      end if

      ! Generate node structure derived from template matrix
      call gfem_genNodeList(rmatrix, rgroupFEMSet)

    end subroutine initGroupFEMSetBoundary
    
  end subroutine transp_initProblemLevel

  !*****************************************************************************

!<subroutine>

  subroutine transp_initAllProblemLevels(rparlist, ssectionName,&
      rproblem, rcollection, rbdrCondPrimal, rbdrCondDual)

!<description>
    ! This subroutine initializes the all problem levels attached to
    ! the global problem structure. It generates the discretisation,
    ! the template matrix and the coefficient matrices as duplicates
    ! of the template matrix.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: boundary condition for primal problem
    type(t_boundaryCondition), intent(in), optional :: rbdrCondPrimal

    ! OPTIONAL: boundary condition for dual problem
    type(t_boundaryCondition), intent(in), optional :: rbdrCondDual
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
          p_rproblemLevel, rcollection, rbdrCondPrimal, rbdrCondDual)

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

    ! section name in parameter list and collection structure
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
    type(t_collection), intent(inout), target :: rcollection
!</intputoutput>
!</subroutine>

    ! local variables
    type(t_afcstab) :: rafcstab
    type(t_collection) :: rcollectionTmp
    type(t_fparser), pointer :: p_rfparser
    type(t_linearForm) :: rform
    type(t_matrixScalar), pointer :: p_rlumpedMassMatrix,p_rConsistentMassMatrix
    type(t_matrixScalar), target :: rlumpedMassMatrix,rconsistentMassMatrix
    type(t_pgm) :: rpgm
    type(t_vectorBlock) :: rvectorHigh,rvectorAux
    real(DP), dimension(:), pointer :: p_Ddata,p_DdofCoords
    real(DP) :: depsAbsSolution,depsRelSolution,dnorm0,dnorm
    character(LEN=SYS_STRLEN) :: ssolutionname
    integer :: iter,isolutiontype,nmaxIterationsSolution
    integer :: lumpedMassMatrix,consistentMassMatrix,systemMatrix,dofCoords

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
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'dofCoords', dofCoords, 0)
      
      if (dofCoords > 0) then
        ! Set pointers
        call lsyssc_getbase_double(rvector%RvectorBlock(1), p_Ddata)
        call lsyssc_getbase_double(&
            rproblemLevel%RvectorBlock(dofCoords)%RvectorBlock(1), p_DdofCoords)
              
        ! Initialize solution from portable graymap image
        call ppsol_readPGM(0, ssolutionName, rpgm)
        
        ! Initialize the solution by the image data
        call ppsol_initArrayPGMDble(rpgm, p_DdofCoords, p_Ddata)
        
        ! Release portable graymap image
        call ppsol_releasePGM(rpgm)
      else
        call output_line('Coordinates of DOFs not available!',&
            OU_CLASS_ERROR,OU_MODE_STD,'transp_initSolution')
        call sys_halt()
      end if


    case (SOLUTION_ANALYTIC_POINTVALUE)

      !-------------------------------------------------------------------------
      ! Initialize the nodal values by the data of an analytical expression
      !-------------------------------------------------------------------------
      
      ! Get global configuration from parameter list
      call parlst_getvalue_string(rparlist,&
          ssectionName, 'ssolutionname', ssolutionName)
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'dofCoords', dofCoords, 0)
      
      ! Get function parser from collection structure
      p_rfparser => collct_getvalue_pars(rcollection,&
          'rfparser', ssectionName=ssectionName)
      
      if (dofCoords > 0) then
        ! Set pointers
        call lsyssc_getbase_double(rvector%RvectorBlock(1), p_Ddata)
        call lsyssc_getbase_double(&
            rproblemLevel%RvectorBlock(dofCoords)%RvectorBlock(1), p_DdofCoords)

        ! Evaluate solution values in the positions of the degrees of freedom
        call fparser_evalFuncBlockByName2(p_rfparser, ssolutionname,&
            rproblemLevel%RvectorBlock(dofCoords)%RvectorBlock(1)%NVAR,&
            rproblemLevel%RvectorBlock(dofCoords)%RvectorBlock(1)%NEQ,&
            p_DdofCoords, rvector%RvectorBlock(1)%NEQ, p_Ddata, (/dtime/))
      else
        call output_line('Coordinates of DOFs not available!',&
            OU_CLASS_ERROR,OU_MODE_STD,'transp_initSolution')
        call sys_halt()
      end if

      
    case (SOLUTION_ANALYTIC_L2_CONSISTENT,&
          SOLUTION_ANALYTIC_L2_LUMPED)

      !-------------------------------------------------------------------------
      ! Initialize the FE-function by the L2-projection of the analytical data
      !-------------------------------------------------------------------------

      ! Get global configuration from parameter list
      call parlst_getvalue_string(rparlist,&
          ssectionName, 'ssolutionname', ssolutionName)

      ! Get function parser from collection structure
      p_rfparser => collct_getvalue_pars(rcollection,&
          'rfparser', ssectionName=ssectionName)

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
      
      ! Initialize temporal collection structure
      call collct_init(rcollectionTmp)
      
      ! Prepare quick access arrays of the temporal collection structure
      rcollectionTmp%SquickAccess(1) = ''
      rcollectionTmp%SquickAccess(2) = 'rfparser'
      rcollectionTmp%DquickAccess(1) = dtime
      rcollectionTmp%IquickAccess(1) = fparser_getFunctionNumber(p_rfparser, ssolutionname)
      
      ! Attach user-defined collection structure to temporal collection
      ! structure (may be required by the callback function)
      rcollectionTmp%p_rnextCollection => rcollection

      ! Attach function parser from boundary conditions to collection
      ! structure and specify its name in quick access string array
      call collct_setvalue_pars(rcollectionTmp, 'rfparser', p_rfparser, .true.)

      ! Set up the linear form
      rform%itermCount = 1
      rform%Idescriptors(1) = DER_FUNC
      
      ! Assemble the linear form for the scalar subvector
      call linf_buildVectorScalar2(rform, .true.,&
          rvector%RvectorBlock(1), transp_coeffVectorAnalytic, rcollectionTmp)

      ! Release temporal collection structure
      call collct_done(rcollectionTmp)

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
        rafcstab%istabilisationSpec = AFCSTAB_UNDEFINED
        rafcstab%cprelimitingType   = AFCSTAB_PRELIMITING_NONE
        rafcstab%cafcstabType = AFCSTAB_LINFCT_MASS
        call afcsc_initStabilisation(rproblemLevel%Rmatrix(systemMatrix), rafcstab)

        ! Compute the raw antidiffusive mass fluxes
        call afcsc_buildFluxFCT(rafcstab, rvectorHigh,&
            0.0_DP, 0.0_DP, 1.0_DP, .true., .true.,&
            AFCSTAB_FCTFLUX_EXPLICIT,&
            rmatrix=p_rconsistentMassMatrix, rxTimeDeriv=rvectorHigh)

        ! Apply flux correction to solution profile
        call afcsc_buildVectorFCT(rafcstab,&
            p_rlumpedMassMatrix, rvector, 1.0_DP, .false.,&
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


    case default
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

    ! section name in parameter list and collection structure
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
    type(t_collection), intent(inout), target :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_linearForm) :: rform
    type(t_collection) :: rcollectionTmp
    character(LEN=SYS_STRLEN) :: srhsname
    integer :: irhstype


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
      p_rfparser => collct_getvalue_pars(rcollection,&
          'rfparser', ssectionName=ssectionName)

      ! Initialize temporal collection structure
      call collct_init(rcollectionTmp)

      ! Prepare quick access arrays of the temporal collection structure
      rcollectionTmp%SquickAccess(1) = ''
      rcollectionTmp%SquickAccess(2) = 'rfparser'
      rcollectionTmp%DquickAccess(1) = dtime
      rcollectionTmp%IquickAccess(1) = fparser_getFunctionNumber(p_rfparser, srhsname)
      
      ! Attach user-defined collection structure to temporal collection
      ! structure (may be required by the callback function)
      rcollectionTmp%p_rnextCollection => rcollection

      ! Attach function parser from boundary conditions to collection
      ! structure and specify its name in quick access string array
      call collct_setvalue_pars(rcollectionTmp, 'rfparser', p_rfparser, .true.)
      
      ! Set up the corresponding linear form
      rform%itermCount      = 1
      rform%Idescriptors(1) = DER_FUNC

      ! Build the discretised right-hand side vector
      call linf_buildVectorScalar2(rform, .true., rvector%RvectorBlock(1),&
                                   transp_coeffVectorAnalytic, rcollectionTmp)

      ! Release temporal collection structure
      call collct_done(rcollectionTmp)


    case default
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

    ! section name in parameter list and collection structure
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
    type(t_collection), intent(inout), target :: rcollection
!</intputoutput>
!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_linearForm) :: rform
    type(t_collection) :: rcollectionTmp
    character(LEN=SYS_STRLEN) :: stargetfuncname
    integer :: itargetfunctype

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
      ! Get global configuration from parameter list
      call parlst_getvalue_string(rparlist, ssectionName,&
          'stargetfuncname', stargetfuncname)
      
      ! Get function parser from collection structure
      p_rfparser => collct_getvalue_pars(rcollection,&
          'rfparser', ssectionName=ssectionName)
      
      ! Initialize temporal collection structure
      call collct_init(rcollectionTmp)
      
      ! Prepare quick access arrays of the temporal collection structure
      rcollectionTmp%SquickAccess(1) = ''
      rcollectionTmp%SquickAccess(2) = 'rfparser'
      rcollectionTmp%DquickAccess(1) = dtime
      rcollectionTmp%IquickAccess(1) = fparser_getFunctionNumber(p_rfparser, stargetfuncname)
      
      ! Attach user-defined collection structure to temporal collection
      ! structure (may be required by the callback function)
      rcollectionTmp%p_rnextCollection => rcollection

      ! Attach function parser from boundary conditions to collection
      ! structure and specify its name in quick access string array
      call collct_setvalue_pars(rcollectionTmp, 'rfparser', p_rfparser, .true.)


      ! Set up the corresponding linear form
      rform%itermCount      = 1
      rform%Idescriptors(1) = DER_FUNC

      ! Build the discretised target functional. The contribution of
      ! the surfae integral comes in be weakly impose boundary conditions
      call linf_buildVectorScalar2(rform, .true., rvector%RvectorBlock(1),&
          transp_coeffVectorAnalytic, rcollectionTmp)

      ! Release temporal collection structure
      call collct_done(rcollectionTmp)


    case default
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
    type(t_ucdExport) :: rexport
    type(t_triangulation) :: rtriangulationPrimal,rtriangulationDual
    type(t_blockDiscretisation) :: rdiscretisationPrimal
    type(t_blockDiscretisation) :: rdiscretisationDual
    type(t_vectorBlock) :: rvectorPrimal,rvectorDual
    real(DP), dimension(:), pointer :: p_DdataPrimal, p_DdataDual
    integer :: iformatUCD,ilineariseUCD,nrefineUCD
    logical :: bexportMeshOnly,bdiscontinuous

    ! Initialisation
    bexportMeshOnly = .true.
    if (present(rsolutionPrimal) .or.&
        present(rsolutionDual)) bexportMeshOnly=.false.

    nullify(p_DdataPrimal, p_DdataDual)

    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist, ssectionName,&
                                'output', soutputName)
    call parlst_getvalue_string(rparlist, trim(soutputName),&
                                'sucdsolution', sucdsolution)
    call parlst_getvalue_int(rparlist, trim(soutputName),&
                             'iformatucd', iformatUCD)
    call parlst_getvalue_int(rparlist, trim(soutputName),&
                             'ilineariseucd', ilineariseUCD, UCDEXPORT_STD)
    call parlst_getvalue_int(rparlist, trim(soutputName),&
                             'nrefineucd', nrefineUCD, 0)
    
    ! Initialize the UCD exporter
    select case(ilineariseUCD)
    case (UCDEXPORT_STD)
      call flagship_initUCDexport(rproblemLevel,&
          sucdsolution, iformatUCD, rexport, ifilenumber)
      
      ! Set pointers to solution(s)
      if (present(rsolutionPrimal))&
          call lsysbl_getbase_double(rsolutionPrimal, p_DdataPrimal)
      if (present(rsolutionDual))&
          call lsysbl_getbase_double(rsolutionDual, p_DdataDual)
      
    case (UCDEXPORT_P1CONTINUOUS,&
          UCDEXPORT_P1DISCONTINUOUS)
      bdiscontinuous = (ilineariseUCD .eq. UCDEXPORT_P1DISCONTINUOUS)
      
      if (present(rsolutionPrimal)) then
        call lin_lineariseVectorGlobal(rsolutionPrimal, rdiscretisationPrimal,&
            rtriangulationPrimal, rvectorPrimal, nrefineUCD, 0, bdiscontinuous)
        call lsysbl_getbase_double(rvectorPrimal, p_DdataPrimal)
        
        if (present(rsolutionDual)) then
          call lin_lineariseVectorGlobal(rsolutionDual, rdiscretisationDual,&
              rtriangulationDual, rvectorDual, nrefineUCD, 0, bdiscontinuous)
          call lsysbl_getbase_double(rvectorDual, p_DdataDual)
        end if
        
        ! We assume that both primal and dual solutions are based on
        ! the same triangulation, thus both vectors can be exported
        ! using the same triangulation.       
        call flagship_initUCDexport(rproblemLevel, sucdsolution,&
            iformatUCD, rexport, ifilenumber, rtriangulationPrimal)

      elseif (present(rsolutionDual)) then
        call lin_lineariseVectorGlobal(rsolutionDual, rdiscretisationDual,&
            rtriangulationDual, rvectorDual, nrefineUCD, 0, bdiscontinuous)
        call lsysbl_getbase_double(rvectorDual, p_DdataDual)
        
        call flagship_initUCDexport(rproblemLevel, sucdsolution,&
            iformatUCD, rexport, ifilenumber, rtriangulationDual)
      end if

    case default
      call output_line('Unsupported type of solution output!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_outputSolution')
      call sys_halt()
    end select

    ! Increase filenumber by one
    ifilenumber = ifilenumber+1

    ! Set simulation time
    if (present(dtime)) call ucd_setSimulationTime(rexport, dtime)

    ! Add primal/dual solution vectors
    if (associated(p_DdataPrimal))&
        call ucd_addVariableVertexBased (rexport, 'u',&
        UCD_VAR_STANDARD, p_DdataPrimal)
    if (associated(p_DdataDual))&
        call ucd_addVariableVertexBased (rexport, 'z',&
        UCD_VAR_STANDARD, p_DdataDual)

    ! Write UCD file
    call ucd_write  (rexport)
    call ucd_release(rexport)

    ! Release temporal memory
    call lsysbl_releaseVector(rvectorPrimal)
    call lsysbl_releaseVector(rvectorDual)
    call spdiscr_releaseBlockDiscr(rdiscretisationPrimal)
    call spdiscr_releaseBlockDiscr(rdiscretisationDual)
    call tria_done(rtriangulationPrimal)
    call tria_done(rtriangulationDual)

  end subroutine transp_outputSolution

  !*****************************************************************************

!<subroutine>

  subroutine transp_outputStatistics(rtimerTotal, ssectionName, rcollection)

!<description>
    ! This subroutine output application statistics
!</description>

!<input>
    ! timer for total time measurement
    type(t_timer), intent(in) :: rtimerTotal

    ! section name in parameter collection structure
    character(LEN=*), intent(in) :: ssectionName
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
    p_rtimerSolution => collct_getvalue_timer(rcollection,&
        'rtimerSolution', ssectionName=ssectionName)
    p_rtimerAdaptation => collct_getvalue_timer(rcollection,&
        'rtimerAdaptation', ssectionName=ssectionName)
    p_rtimerErrorEstimation => collct_getvalue_timer(rcollection,&
        'rtimerErrorEstimation', ssectionName=ssectionName)
    p_rtimerTriangulation => collct_getvalue_timer(rcollection,&
        'rtimerTriangulation', ssectionName=ssectionName)
    p_rtimerAssemblyCoeff => collct_getvalue_timer(rcollection,&
        'rtimerAssemblyCoeff', ssectionName=ssectionName)
    p_rtimerAssemblyMatrix => collct_getvalue_timer(rcollection,&
        'rtimerAssemblyMatrix', ssectionName=ssectionName)
    p_rtimerAssemblyVector => collct_getvalue_timer(rcollection,&
        'rtimerAssemblyVector', ssectionName=ssectionName)
    p_rtimerPrePostprocess => collct_getvalue_timer(rcollection,&
        'rtimerPrePostprocess', ssectionName=ssectionName)

    ! Output statistics
    call output_lbrk()
    call output_line('Time measurement: '//trim(adjustl(ssectionName)))
    call output_line('-----------------')

    call stat_subTimers(p_rtimerAssemblyMatrix, p_rtimerSolution)
    call stat_subTimers(p_rtimerAssemblyVector, p_rtimerSolution)

    dtotalTime = max(rtimerTotal%delapsedCPU, rtimerTotal%delapsedReal)
    dfraction  = 100.0_DP/dtotalTime

    call output_line('Time for computing solution   : '//&
                     trim(adjustl(sys_sdE(p_rtimerSolution%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerSolution%delapsedCPU, 5)))//' %')
    call output_line('Time for mesh adaptivity      : '//&
                     trim(adjustl(sys_sdE(p_rtimerAdaptation%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerAdaptation%delapsedCPU, 5)))//' %')
    call output_line('Time for error estimation     : '//&
                     trim(adjustl(sys_sdE(p_rtimerErrorEstimation%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerErrorEstimation%delapsedCPU, 5)))//' %')
    call output_line('Time for triangulation        : '//&
                     trim(adjustl(sys_sdE(p_rtimerTriangulation%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerTriangulation%delapsedCPU, 5)))//' %')
    call output_line('Time for coefficient assembly : '//&
                     trim(adjustl(sys_sdE(p_rtimerAssemblyCoeff%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerAssemblyCoeff%delapsedCPU, 5)))//' %')
    call output_line('Time for matrix assembly      : '//&
                     trim(adjustl(sys_sdE(p_rtimerAssemblyMatrix%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerAssemblyMatrix%delapsedCPU, 5)))//' %')
    call output_line('Time for vector assembly      : '//&
                     trim(adjustl(sys_sdE(p_rtimerAssemblyVector%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerAssemblyVector%delapsedCPU, 5)))//' %')
    call output_line('Time for pre-/post-processing : '//&
                     trim(adjustl(sys_sdE(p_rtimerPrePostprocess%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerPrePostprocess%delapsedCPU, 5)))//' %')
    call output_lbrk()
    call output_line('Time for total simulation     : '//&
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
    ! section name in parameter list and collection structure
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

    ! collection structure
    type(t_collection), intent(inout), target :: rcollection
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
    type(t_collection) :: rcollectionTmp
    real(DP), dimension(:), pointer :: p_DsolutionDual, p_Dresidual
    real(DP), dimension(:), pointer :: p_DlumpedMassMatrix, p_DtargetError
    integer, dimension(:,:), pointer :: p_IverticesAtElement, p_IneighboursAtElement
    logical, dimension(:), pointer :: p_BisactiveElement
    real(DP) :: dexactTargetError, dexactTargetFunc, dprotectLayerTolerance
    real(DP) :: daux, dtargetFunc, dStep, theta
    integer :: i, convectionAFC, diffusionAFC, NEQ
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
        rproblemLevel%Rafcstab(convectionAFC)%cafcstabType
    rproblemLevel%Rafcstab(convectionAFC)%cafcstabType = AFCSTAB_GALERKIN

    cdiffusionStabilisation =&
        rproblemLevel%Rafcstab(diffusionAFC)%cafcstabType
    rproblemLevel%Rafcstab(diffusionAFC)%cafcstabType = AFCSTAB_GALERKIN

    ! Set update notifiers for the discrete transport operator and the
    ! preconditioner in the problem level structure
    rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                     TRANSP_TROPER_UPDATE)
    rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                     TRANSP_PRECOND_UPDATE)

    ! Calculate the standard Galerkin preconditioner
    ! (required for rhs and residual calculation)
    call transp_calcPrecondThetaScheme(rproblemLevel, rtimestep,&
        rsolver, rsolutionPrimal, ssectionName, rcollection)

    ! Calculate the standard Galerkin right-hand side vector
    call transp_calcRhsThetaScheme(rproblemLevel, rtimestep, rsolver,&
        rsolutionPrimal, rvector1, ssectionName, rcollection, rrhs)

    ! Calculate the standard Galerkin residual
    call transp_calcResidualThetaScheme(rproblemLevel, rtimestep,&
        rsolver, rsolutionPrimal, rsolutionPrimal, rvector1,&
        rvector2, 0, ssectionName, rcollection)

    ! Ok, now we have to switch on all types of stabilisation again
    rproblemLevel%Rafcstab(convectionAFC)%cafcstabType =&
        cconvectionStabilisation
    rproblemLevel%Rafcstab(diffusionAFC)%cafcstabType =&
        cdiffusionStabilisation

    ! ... and we reset the mass type
    call parlst_setvalue(rparlist, ssectionName, 'imasstype',&
        trim(sys_si(imasstype, 3)))

    ! ... and the time-stepping structure
    rtimestep%dStep = dStep
    rtimestep%theta = theta

    ! Again, set update notifiers for the discrete transport operator
    ! and the preconditioner in the problem level structure
    rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                     TRANSP_TROPER_UPDATE)
    rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                     TRANSP_PRECOND_UPDATE)

    ! We need the lumped mass matrix for scaling
    call parlst_getvalue_int(rparlist, ssectionName,&
        'lumpedMassMatrix', lumpedMassMatrix)

    if (lumpedMassMatrix > 0) then
      ! Set pointer to the lumped mass matrix
      call lsyssc_getbase_double(&
          rproblemLevel%Rmatrix(lumpedMassMatrix), p_DlumpedMassMatrix)
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
        ! Get global configuration from parameter list
        call parlst_getvalue_string(rparlist, ssectionName,&
            'stargetfuncname', stargetfuncName)
        
        ! Get function parser from collection structure
        p_rfparser => collct_getvalue_pars(rcollection,&
            'rfparser', ssectionName=ssectionName)

        ! Initialize temporal collection structure
        call collct_init(rcollectionTmp)
        
        ! Prepare quick access arrays of the temporal collection structure
        rcollectionTmp%SquickAccess(1) = ''
        rcollectionTmp%SquickAccess(2) = 'rfparser'
        rcollectionTmp%DquickAccess(1) = rtimestep%dTime
        rcollectionTmp%IquickAccess(1) = fparser_getFunctionNumber(p_rfparser, sexactsolutionName)
        rcollectionTmp%IquickAccess(2) = fparser_getFunctionNumber(p_rfparser, stargetfuncName)

        ! Attach user-defined collection structure to temporal collection
        ! structure (may be required by the callback function)
        rcollectionTmp%p_rnextCollection => rcollection
        
        ! Attach function parser from boundary conditions to collection
        ! structure and specify its name in quick access string array
        call collct_setvalue_pars(rcollectionTmp, 'rfparser', p_rfparser, .true.)
        

        if (trim(sexacttargetfuncname) .ne. '') then
          ! Evaluate exact value of the quantity of interest
          call fparser_evalFunction(p_rfparser, sexacttargetfuncName,&
              (/rtimestep%dTime/), dexactTargetFunc)

          ! Compute the approximate value of the quantity of interest
          call pperr_scalar(PPERR_MEANERROR, dtargetFunc,&
              rsolutionPrimal%RvectorBlock(1), rcollection=rcollectionTmp,&
              ffunctionWeight=transp_weightFuncAnalytic)
          
          ! Compute exact error in target functional. Note that the
          ! routine pperr_scalar computes the value $J(0-u_h)$ so that we
          ! have to change the sign to minus
          dexactTargetError = dexactTargetFunc+dtargetFunc
          
        else
          
          ! Compute the exact error of the quantity of interest
          call pperr_scalar(PPERR_MEANERROR, dexactTargetError,&
              rsolutionPrimal%RvectorBlock(1), transp_refFuncAnalytic,&
              rcollectionTmp, ffunctionWeight=transp_weightFuncAnalytic)
          
          ! Compute the exact value of the quantity of interest
          call pperr_scalar(PPERR_MEANERROR, dexactTargetFunc,&
              ffunctionReference=transp_refFuncAnalytic,&
              rcollection=rcollectionTmp,&
              rdiscretisation=rsolutionPrimal%RvectorBlock(1)%p_rspatialdiscr,&
              ffunctionWeight=transp_weightFuncAnalytic)
        end if

        ! Release temporal collection structure
        call collct_done(rcollectionTmp)


      case (TFUNC_SURFINTG)
        ! Get global configuration from parameter list
        call parlst_getvalue_int(rparlist, ssectionName,&
            'velocityfield', velocityfield)
        call parlst_getvalue_string(rparlist, ssectionName,&
            'stargetfuncname', stargetfuncName)

        ! Get function parser from collection structure
        p_rfparser => collct_getvalue_pars(rcollection,&
            'rfparser', ssectionName=ssectionName)
        
        ! Initialize temporal collection structure
        call collct_init(rcollectionTmp)

        ! Prepare quick access arrays of the temporal collection structure
        rcollectionTmp%SquickAccess(1) = ''
        rcollectionTmp%SquickAccess(2) = 'rfparser'
        rcollectionTmp%DquickAccess(1) = rtimestep%dTime
        rcollectionTmp%IquickAccess(1) = fparser_getFunctionNumber(p_rfparser, sexactsolutionName)
        rcollectionTmp%IquickAccess(2) = fparser_getFunctionNumber(p_rfparser, stargetfuncName)

        ! ... and also the numbers of the exact velocity function
        do i = 1, rproblemLevel%rtriangulation%ndim
          call parlst_getvalue_string(rparlist, ssectionName,&
              'svelocityname', svelocityname, isubString=i)
          rcollectionTmp%IquickAccess(i+2) = fparser_getFunctionNumber(p_rfparser, svelocityname)
        end do

        ! Attach primal solution vector and velocity fied to first and
        ! second quick access vector of the temporal collection structure
        rcollectionTmp%p_rvectorQuickAccess1 => rsolutionPrimal
        rcollectionTmp%p_rvectorQuickAccess2 => rproblemLevel%RvectorBlock(velocityfield)
        
        ! Attach user-defined collection structure to temporal collection
        ! structure (may be required by the callback function)
        rcollectionTmp%p_rnextCollection => rcollection
        
        ! Attach function parser from boundary conditions to collection
        ! structure and specify its name in quick access string array
        call collct_setvalue_pars(rcollectionTmp, 'rfparser', p_rfparser, .true.)
        

        if (trim(sexacttargetfuncname) .ne. '') then
          ! Evaluate exact value of the quantity of interest
          call fparser_evalFunction(p_rfparser, sexacttargetfuncName,&
              (/rtimestep%dTime/), dexactTargetFunc)
          
          ! Prepare quick access arrays of the temporal collection structure
          rcollectionTmp%IquickAccess(1) = fparser_getFunctionNumber(p_rfparser, '@null')
          
          ! Compute the approximate value of the quantity of interest
          call pperr_scalarBoundary2D(0, CUB_G3_1D, dtargetFunc,&
              ffunctionReference=transp_errorBdrInt2D_sim,&
              rcollection=rcollectionTmp,&
              rdiscretisation=rsolutionPrimal%RvectorBlock(1)%p_rspatialdiscr,&
              ffunctionWeight=transp_weightFuncBdrInt2D_sim)

          ! Compute exact error in target functional. Note that the
          ! routine pperr_scalar computes the value $J(0-u_h)$ so that we
          ! have to change the sign to minus
          dexactTargetError = dexactTargetFunc+dtargetFunc

        else

          ! Compute the exact error of the quantity of interest
          call pperr_scalarBoundary2D(0, CUB_G3_1D, dexactTargetError,&
              ffunctionReference=transp_errorBdrInt2D_sim,&
              rcollection=rcollectionTmp,&
              rdiscretisation=rsolutionPrimal%RvectorBlock(1)%p_rspatialdiscr,&
              ffunctionWeight=transp_weightFuncBdrInt2D_sim)

          ! Compute the exact value of the quantity of interest
          call pperr_scalarBoundary2D(0, CUB_G3_1D, dexactTargetFunc,&
              ffunctionReference=transp_refFuncBdrInt2D_sim,&
              rcollection=rcollectionTmp,&
              rdiscretisation=rsolutionPrimal%RvectorBlock(1)%p_rspatialdiscr,&
              ffunctionWeight=transp_weightFuncBdrInt2D_sim)
        end if

        ! Release temporal collection structure
        call collct_done(rcollectionTmp)


      case (TFUNC_MIXINTG)
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

        ! Get function parser from collection structure
        p_rfparser => collct_getvalue_pars(rcollection,&
            'rfparser', ssectionName=ssectionName)
        
        ! Initialize temporal collection structure
        call collct_init(rcollectionTmp)

        ! Prepare quick access arrays of the temporal collection structure
        rcollectionTmp%SquickAccess(1) = ''
        rcollectionTmp%SquickAccess(2) = 'rfparser'
        rcollectionTmp%DquickAccess(1) = rtimestep%dTime
        rcollectionTmp%IquickAccess(1) = fparser_getFunctionNumber(p_rfparser, sexactsolutionName)
        rcollectionTmp%IquickAccess(2) = fparser_getFunctionNumber(p_rfparser, stargetfuncName)
        
        ! ... and also the numbers of the exact velocity function
        do i = 1, rproblemLevel%rtriangulation%ndim
          call parlst_getvalue_string(rparlist, ssectionName,&
              'svelocityname', svelocityname, isubString=i)
          rcollectionTmp%IquickAccess(i+2) = fparser_getFunctionNumber(p_rfparser, svelocityname)
        end do

        ! Attach primal solution vector and velocity fied to first and
        ! second quick access vector of the temporal collection structure
        rcollectionTmp%p_rvectorQuickAccess1 => rsolutionPrimal
        rcollectionTmp%p_rvectorQuickAccess2 => rproblemLevel%RvectorBlock(velocityfield)

        ! Attach user-defined collection structure to temporal collection
        ! structure (may be required by the callback function)
        rcollectionTmp%p_rnextCollection => rcollection
        
        ! Attach function parser from boundary conditions to collection
        ! structure and specify its name in quick access string array
        call collct_setvalue_pars(rcollectionTmp, 'rfparser', p_rfparser, .true.)


        if (trim(sexacttargetfuncname) .ne. '') then
          ! Evaluate exact value of the quantity of interest
          call fparser_evalFunction(p_rfparser, sexacttargetfuncName,&
              (/rtimestep%dTime/), dexactTargetFunc)

          ! Prepare quick access arrays of the collection
          rcollectionTmp%IquickAccess(1) = fparser_getFunctionNumber(p_rfparser, '@null')

          ! Compute the approximate value of the quantity of interest
          call pperr_scalar(PPERR_MEANERROR, dtargetFunc,&
              rsolutionPrimal%RvectorBlock(1), rcollection=rcollectionTmp,&
              ffunctionWeight=transp_weightFuncAnalytic)
          
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
          rcollectionTmp%IquickAccess(2) = fparser_getFunctionNumber(p_rfparser, stargetfuncName)

          ! Compute the approximate value of the quantity of interest
          call pperr_scalarBoundary2D(0, CUB_G3_1D, daux,&
              ffunctionReference=transp_errorBdrInt2D_sim,&
              rcollection=rcollectionTmp,&
              rdiscretisation=rsolutionPrimal%RvectorBlock(1)%p_rspatialdiscr,&
              ffunctionWeight=transp_weightFuncBdrInt2D_sim)

          ! Add boundary contribution
          dtargetFunc = dtargetFunc + daux

          ! Compute exact error in target functional. Note that the
          ! routine pperr_scalar computes the value $J(0-u_h)$ so that we
          ! have to change the sign to minus
          dexactTargetError = dexactTargetFunc+dtargetFunc

        else

          ! Compute the exact error of the quantity of interest
          call pperr_scalar(PPERR_MEANERROR, dexactTargetError,&
              rsolutionPrimal%RvectorBlock(1), transp_refFuncAnalytic,&
              rcollectionTmp, ffunctionWeight=transp_weightFuncAnalytic)

          ! Compute the exact value of the quantity of interest.
          call pperr_scalar(PPERR_MEANERROR, dexactTargetFunc,&
              ffunctionReference=transp_refFuncAnalytic,&
              rcollection=rcollectionTmp,&
              rdiscretisation=rsolutionPrimal%RvectorBlock(1)%p_rspatialdiscr,&
              ffunctionWeight=transp_weightFuncAnalytic)
          
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
          rcollectionTmp%DquickAccess(1) = rtimestep%dTime
          rcollectionTmp%IquickAccess(2) = fparser_getFunctionNumber(p_rfparser, stargetfuncName)

          ! Compute the exact error of the quantity of interest at the boundary
          call pperr_scalarBoundary2D(0, CUB_G3_1D, daux,&
              ffunctionReference=transp_errorBdrInt2D_sim,&
              rcollection=rcollectionTmp,&
              rdiscretisation=rsolutionPrimal%RvectorBlock(1)%p_rspatialdiscr,&
              ffunctionWeight=transp_weightFuncBdrInt2D_sim)

          ! Add boundary contribution
          dexactTargetError = dexactTargetError + daux

          ! Compute the exact value of the quantity of interest at the boundary
          call pperr_scalarBoundary2D(0, CUB_G3_1D, daux,&
              ffunctionReference=transp_refFuncBdrInt2D_sim,&
              rcollection=rcollectionTmp,&
              rdiscretisation=rsolutionPrimal%RvectorBlock(1)%p_rspatialdiscr,&
              ffunctionWeight=transp_weightFuncBdrInt2D_sim)

          ! Add boundary contribution
          dexactTargetFunc = dexactTargetFunc + daux
        end if

        ! Release temporal collection structure
        call collct_done(rcollectionTmp)


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

    case default
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

    ! section name in parameter list and collection structure
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
    type(t_collection), intent(inout), target :: rcollection
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
    type(t_collection) :: rcollectionTmp
    real(DP), dimension(:), pointer :: p_Ddata, p_Derror
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:,:), pointer :: p_IneighboursAtElement
    logical, dimension(:), pointer :: p_BisactiveElement
    real(DP) :: dnoiseFilter, dabsFilter, dsolution, dvalue,&
                dexacterror, dprotectLayerTolerance
    integer :: i, ierrorEstimator, igridindicator, iexactsolutiontype
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

    case default
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
      p_rfparser => collct_getvalue_pars(rcollection,&
          'rfparser', ssectionName=ssectionName)

      ! Initialize temporal collection structure
      call collct_init(rcollectionTmp)

      ! Prepare quick access arrays of the temporal collection structure
      rcollectionTmp%SquickAccess(1) = ''
      rcollectionTmp%SquickAccess(2) = 'rfparser'
      rcollectionTmp%DquickAccess(1) = dtime
      rcollectionTmp%IquickAccess(1) = fparser_getFunctionNumber(p_rfparser, sexactsolutionName)
      
      ! Attach user-defined collection structure to temporal collection
      ! structure (may be required by the callback function)
      rcollectionTmp%p_rnextCollection => rcollection

      ! Attach function parser from boundary conditions to collection
      ! structure and specify its name in quick access string array
      call collct_setvalue_pars(rcollectionTmp, 'rfparser', p_rfparser, .true.)

      ! Calculate the H1-error of the reference solution
      call pperr_scalar(PPERR_H1ERROR, dexacterror, rsolution%RvectorBlock(1),&
                        transp_refFuncAnalytic, rcollectionTmp)

      call output_lbrk()
      call output_line('Error Analysis')
      call output_line('--------------')
      call output_line('estimated H1-error: '//trim(sys_sdEP(derror,15,6)))
      call output_line('exact H1-error:     '//trim(sys_sdEP(dexacterror,15,6)))
      call output_line('effectivity index:  '//trim(sys_sdEP(derror/dexacterror,15,6)))
      call output_lbrk()

      ! Release temporal collection structure
      call collct_done(rcollectionTmp)

    case default
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
          rsolution%RvectorBlock(1))

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
      call pperr_scalar(PPERR_H1ERROR, dsolution,&
          rsolution%RvectorBlock(1), relementError=rvectorScalar)

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
      dvalue = -SYS_MAXREAL_DP
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


    case default
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
      call storage_getbase_int2D(&
          rproblemLevel%rtriangulation%h_IneighboursAtElement, p_IneighboursAtElement)
      call storage_getbase_int2D(&
          rproblemLevel%rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
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

    ! collection structure
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
      ! section name in parameter list and collection structure
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

    ! Matrix for the linearised FCT algorithm
    type(t_matrixScalar) :: rmatrix1

    ! Vectors for the linearised FCT algorithm
    type(t_vectorBlock) :: rvector1, rvector2, rvector3

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
    p_rtimerPrePostprocess => collct_getvalue_timer(rcollection,&
        'rtimerPrePostprocess', ssectionName=ssectionName)
    p_rtimerSolution => collct_getvalue_timer(rcollection,&
        'rtimerSolution', ssectionName=ssectionName)
    p_rtimerErrorEstimation => collct_getvalue_timer(rcollection,&
        'rtimerErrorEstimation', ssectionName=ssectionName)
    p_rtimerAdaptation => collct_getvalue_timer(rcollection,&
        'rtimerAdaptation', ssectionName=ssectionName)
    p_rtimerTriangulation => collct_getvalue_timer(rcollection,&
        'rtimerTriangulation', ssectionName=ssectionName)
    p_rtimerAssemblyCoeff => collct_getvalue_timer(rcollection,&
        'rtimerAssemblyCoeff', ssectionName=ssectionName)

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
        call grph_createGraphFromMatrix(&
            p_rproblemLevel%Rmatrix(templateMatrix), rgraph)
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
            call transp_adaptTriangulation(rhadapt,&
                p_rproblemLevel%rtriangulation, relementError, rcollection)

            ! Release element-wise error distribution
            call lsyssc_releaseVector(relementError)

            ! Generate standard mesh from raw mesh
            call tria_initStandardMeshFromRaw(&
                p_rproblemLevel%rtriangulation, rproblem%rboundary)

            ! Update the template matrix according to the sparsity pattern
            call grph_generateMatrix(rgraph,&
                p_rproblemLevel%Rmatrix(templateMatrix))

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
                p_rproblemLevel, rcollection, rbdrCond)
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

    ! Force initialisation of the discrete transport operator and the
    ! preconditioner. This may be necessary if the velocity is
    ! constant, and hence, no repeated updates would be performed.
    call problem_setSpec(rproblem, TRANSP_TROPER_INIT+&
                                   TRANSP_PRECOND_INIT,'ior')

    ! Initialize right-hand side vector
    if (irhstype > 0) then
      call lsysbl_createVectorBlock(rsolution, rrhs)
      call transp_initRHS(rparlist, ssectionName, p_rproblemLevel,&
          rtimestep%dinitialTime, rrhs, rcollection)
    end if

    ! Calculate the initial velocity field on all levels
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

      case default
        call output_line('Unsupported time-stepping algorithm!',&
            OU_CLASS_ERROR,OU_MODE_STD,'transp_solveTransientPrimal')
        call sys_halt()
      end select

      ! Perform linearised FEM-FCT post-processing
      call transp_calcLinearisedFCT(p_rproblemLevel, rtimestep, rsolver,&
          rsolution, ssectionName, rcollection, rmatrix=rmatrix1,&
          rvector1=rvector1, rvector2=rvector2, rvector3=rvector3)

      ! Perform linearised FEM-LPT post-processing
      call transp_calcLinearisedLPT(p_rproblemLevel, rtimestep, rsolver,&
          rsolution, ssectionName, rcollection, rmatrix=rmatrix1,&
          rvector1=rvector1, rvector2=rvector2, rvector3=rvector3)

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


        !-----------------------------------------------------------------------
        ! Perform h-adaptation
        !-----------------------------------------------------------------------

        ! Start time measurement for mesh adaptation
        call stat_startTimer(p_rtimerAdaptation, STAT_TIMERSHORT)

        ! Set the name of the template matrix
        rcollection%SquickAccess(1) = 'sparsitypattern'

        ! Attach the primal solution vector to the collection structure
        rcollection%p_rvectorQuickAccess1 => rsolution

        ! Perform h-adaptation and update the triangulation structure
        call transp_adaptTriangulation(rhadapt,&
            p_rproblemLevel%rtriangulation, relementError, rcollection)

        ! Release element-wise error distribution
        call lsyssc_releaseVector(relementError)

        ! Update the template matrix according to the sparsity pattern
        call grph_generateMatrix(rgraph,&
            p_rproblemLevel%Rmatrix(templateMatrix))

        ! Resize the solution vector accordingly
        call lsysbl_resizeVectorBlock(rsolution, &
            p_rproblemLevel%Rmatrix(templateMatrix)%NEQ, .false.)

        ! Stop time measurement for mesh adaptation
        call stat_stopTimer(p_rtimerAdaptation)


        !-----------------------------------------------------------------------
        ! Re-generate the discretisation and coefficient matrices
        !-----------------------------------------------------------------------

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
            p_rproblemLevel, rcollection, rbdrCond)

        ! Prepare internal data arrays of the solver structure
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'systemmatrix', systemMatrix)
        call flagship_updateSolverMatrix(p_rproblemLevel, rsolver,&
            systemMatrix, SYSTEM_INTERLEAVEFORMAT, UPDMAT_ALL)
        call solver_updateStructure(rsolver)

        ! Re-calculate the velocity field ...
        if (abs(ivelocitytype) .eq. VELOCITY_TIMEDEP) then
          ! ... on all problem levels if it depends on time
          nlmin = solver_getMinimumMultigridlevel(rsolver)
          call transp_calcVelocityField(rparlist, ssectionName,&
              p_rproblemLevel, rtimestep%dTime, rcollection, nlmin)
        else
          ! ... on the finest mesh only if the velocity is constant
          call transp_calcVelocityField(rparlist, ssectionName,&
              p_rproblemLevel, rtimestep%dTime, rcollection, nlmin)
        end if

        ! Re-initialize the right-hand side vector
        if (irhstype > 0) then
          call lsysbl_resizeVectorBlock(rrhs,&
              p_rproblemLevel%Rmatrix(templateMatrix)%NEQ, .false.)
          call transp_initRHS(rparlist, ssectionName, p_rproblemLevel,&
              rtimestep%dinitialTime, rrhs, rcollection)
        end if

        ! Force initialisation of the discrete transport operator and
        ! the preconditioner. This may be necessary if the velocity is
        ! constant, and hence, no repeated updates would be performed.
        p_rproblemLevel%iproblemSpec = ior(p_rproblemLevel%iproblemSpec,&
            TRANSP_TROPER_INIT+TRANSP_PRECOND_INIT)
        
        ! Stop time measurement for generation of constant coefficient matrices
        call stat_stopTimer(p_rtimerAssemblyCoeff)

      elseif(abs(ivelocitytype) .eq. VELOCITY_TIMEDEP) then

        ! Re-calculate the time-dependent velocity field on all levels
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

    ! Release vectors and matrices for the linearised FCT algorithm
    call lsyssc_releaseMatrix(rmatrix1)
    call lsysbl_releaseVector(rvector1)
    call lsysbl_releaseVector(rvector2)
    call lsysbl_releaseVector(rvector3)

  end subroutine transp_solveTransientPrimal

  !*****************************************************************************

!<subroutine>

  subroutine transp_solveTransientPrimDual(rparlist, ssectionName,&
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
    ! section name in parameter list and collection structure
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

  end subroutine transp_solveTransientPrimDual

  !*****************************************************************************

!<subroutine>

  subroutine transp_solvePseudoTransPrimal(rparlist, ssectionName,&
      rbdrCond, rproblem, rtimestep, rsolver, rsolution, rcollection)

!<description>
    ! This subroutine solves the pseudo-transient primal flow problem
    !
    ! $$\frac{\partial u}{\partial \tau}+\nabla\cdot{\bf f}(u)=s(u)$$
    !
    ! for a scalar quantity $u$ in the domain $\Omega$.
!</description>

!<input>
    ! section name in parameter list and collection structure
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
    p_rtimerPrePostprocess => collct_getvalue_timer(rcollection,&
        'rtimerPrePostprocess', ssectionName=ssectionName)
    p_rtimerSolution => collct_getvalue_timer(rcollection,&
        'rtimerSolution', ssectionName=ssectionName)
    p_rtimerErrorEstimation => collct_getvalue_timer(rcollection,&
        'rtimerErrorEstimation', ssectionName=ssectionName)
    p_rtimerAdaptation => collct_getvalue_timer(rcollection,&
        'rtimerAdaptation', ssectionName=ssectionName)
    p_rtimerTriangulation => collct_getvalue_timer(rcollection,&
        'rtimerTriangulation', ssectionName=ssectionName)
    p_rtimerAssemblyCoeff => collct_getvalue_timer(rcollection,&
        'rtimerAssemblyCoeff', ssectionName=ssectionName)

    ! Start time measurement for pre-processing
    call stat_startTimer(p_rtimerPrePostprocess, STAT_TIMERSHORT)

    ! Get global parameters
    call parlst_getvalue_int(rparlist, ssectionName,&
        'irhstype', irhstype)
    call parlst_getvalue_int(rparlist, ssectionName,&
        'ivelocitytype', ivelocitytype)
    call parlst_getvalue_int(rparlist, ssectionName,&
        'discretisation', discretisation)

    ! Set pointer to maximum problem level
    p_rproblemLevel => rproblem%p_rproblemLevelMax
    p_rdiscretisation => p_rproblemLevel%Rdiscretisation(discretisation)

    ! Create the solution vector
    call lsysbl_createVectorBlock(p_rdiscretisation, rsolution,&
        .false., ST_DOUBLE)

    ! Initialize the solution vector and impose boundary conditions
    call transp_initSolution(rparlist, ssectionName, p_rproblemLevel,&
        rtimestep%dinitialTime, rsolution, rcollection)
    call bdrf_filterVectorExplicit(rbdrCond, rsolution,&
        rtimestep%dinitialTime)

    !---------------------------------------------------------------------------
    ! Initialize the h-adaptation structure
    !---------------------------------------------------------------------------

    call parlst_getvalue_string(rparlist,&
        ssectionName, 'adaptivity', sadaptivityName, '')
    if (trim(adjustl(sadaptivityName)) .ne. '') then

      call parlst_getvalue_int(rparlist,&
          trim(sadaptivityName), 'nadapt', nadapt)

      if (nadapt > 0) then

        ! Initialize adaptation structure from parameter list
        call hadapt_initFromParameterlist(rhadapt, rparlist, sadaptivityName)

        ! Generate a dynamic graph for the sparsity pattern and attach
        ! it to the collection structure which is used to pass this type
        ! to the callback function for h-adaptation
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'templateMatrix', templateMatrix)
        call grph_createGraphFromMatrix(&
            p_rproblemLevel%Rmatrix(templateMatrix), rgraph)
        call collct_setvalue_graph(rcollection, 'sparsitypattern',&
            rgraph, .true.)

        ! Attach the primal solution vector to the collection structure
        rcollection%p_rvectorQuickAccess1 => rsolution

      end if   ! nadapt > 0

    else

      ! No h-adaptation
      nadapt = 0

    end if

    ! Force initialisation of the discrete transport operator and the
    ! preconditioner. This may be necessary if the velocity is
    ! constant, and hence, no repeated updates would be performed.
    call problem_setSpec(rproblem, TRANSP_TROPER_INIT+&
                                   TRANSP_PRECOND_INIT,'ior')

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

        ! Solve the primal problem with non-zero right-hand side
        call tstep_performPseudoStepping(p_rproblemLevel, rtimestep,&
            rsolver, rsolution, transp_nlsolverCallback, rcollection, rrhs)

        ! Release right-hand side vector
        call lsysbl_releaseVector(rrhs)

      else

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
      call transp_adaptTriangulation(rhadapt,&
          p_rproblemLevel%rtriangulation, relementError, rcollection)

      ! Update the template matrix according to the sparsity pattern
      call grph_generateMatrix(rgraph,&
          p_rproblemLevel%Rmatrix(templateMatrix))

      ! Resize the solution vector accordingly
      call lsysbl_resizeVectorBlock(rsolution,&
          p_rproblemLevel%Rmatrix(templateMatrix)%NEQ, .false.)

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
          p_rproblemLevel, rcollection, rbdrCond)

      ! Prepare internal data arrays of the solver structure
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'systemmatrix', systemMatrix)
      call flagship_updateSolverMatrix(p_rproblemLevel, rsolver,&
          systemMatrix, SYSTEM_INTERLEAVEFORMAT, UPDMAT_ALL)
      call solver_updateStructure(rsolver)

      ! Force initialisation of the discrete transport operator and
      ! the preconditioner. This may be necessary if the velocity is
      ! constant, and hence, no repeated updates would be performed.
      p_rproblemLevel%iproblemSpec = ior(p_rproblemLevel%iproblemSpec,&
          TRANSP_TROPER_INIT+TRANSP_PRECOND_INIT)
      
      ! Stop time measurement for generation of constant coefficient matrices
      call stat_stopTimer(p_rtimerAssemblyCoeff)

    end do adaptloop


    ! Release adaptation structure
    if (nadapt > 0) then
      call hadapt_releaseAdaptation(rhadapt)
      call grph_releaseGraph(rgraph)
    end if

  end subroutine transp_solvePseudoTransPrimal

  !*****************************************************************************

!<subroutine>

  subroutine transp_solvePseudoTransPrimDual(rparlist, ssectionName,&
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
    ! section name in parameter list and collection structure
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
    p_rtimerPrePostprocess => collct_getvalue_timer(rcollection,&
        'rtimerPrePostprocess', ssectionName=ssectionName)
    p_rtimerSolution => collct_getvalue_timer(rcollection,&
        'rtimerSolution', ssectionName=ssectionName)
    p_rtimerErrorEstimation => collct_getvalue_timer(rcollection,&
        'rtimerErrorEstimation', ssectionName=ssectionName)
    p_rtimerAdaptation => collct_getvalue_timer(rcollection,&
        'rtimerAdaptation', ssectionName=ssectionName)
    p_rtimerTriangulation => collct_getvalue_timer(rcollection,&
        'rtimerTriangulation', ssectionName=ssectionName)
    p_rtimerAssemblyCoeff => collct_getvalue_timer(rcollection,&
        'rtimerAssemblyCoeff', ssectionName=ssectionName)

    ! Start time measurement for pre-processing
    call stat_startTimer(p_rtimerPrePostprocess, STAT_TIMERSHORT)

    ! Get global parameters
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'irhstype', irhstype)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'ivelocitytype', ivelocitytype)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'discretisation', discretisation)

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

    call parlst_getvalue_string(rparlist,&
        ssectionName, 'adaptivity', sadaptivityName, '')
    if (trim(adjustl(sadaptivityName)) .ne. '') then

      call parlst_getvalue_int(rparlist,&
          trim(sadaptivityName), 'nadapt', nadapt)

      if (nadapt > 0) then

        ! Initialize adaptation structure from parameter list
        call hadapt_initFromParameterlist(rhadapt, rparlist, sadaptivityName)

        ! Generate a dynamic graph for the sparsity pattern and attach
        ! it to the collection structure which is used to pass this type
        ! to the callback function for h-adaptation
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'templateMatrix', templateMatrix)
        call grph_createGraphFromMatrix(&
            p_rproblemLevel%Rmatrix(templateMatrix), rgraph)
        call collct_setvalue_graph(rcollection, 'sparsitypattern',&
            rgraph, .true.)

        ! Attach the primal solution vector to the collection structure
        rcollection%p_rvectorQuickAccess1 => rsolutionPrimal

      end if   ! nadapt > 0

    else

      ! No h-adaptation
      nadapt = 0

    end if

    ! Force initialisation of the discrete transport operator and the
    ! preconditioner. This may be necessary if the velocity is
    ! constant, and hence, no repeated updates would be performed.
    call problem_setSpec(rproblem, TRANSP_TROPER_INIT+&
                                   TRANSP_PRECOND_INIT,'ior')

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

        ! Solve the primal problem with non-zero right-hand side
        call tstep_performPseudoStepping(p_rproblemLevel, rtimestep,&
            rsolver, rsolutionPrimal, transp_nlsolverCallback,&
            rcollection, rrhs)

        ! Release right-hand side vector
        call lsysbl_releaseVector(rrhs)

      else

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

        ! Compute the error in the quantity of interest
        call transp_estimateTargetFuncError(rparlist, ssectionName,&
            p_rproblemLevel, rtimestep, rsolver, rsolutionPrimal,&
            rsolutionDual, rcollection, relementError, derror, rrhs)

        ! Release right-hand side vector
        call lsysbl_releaseVector(rrhs)

      else

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
      call transp_adaptTriangulation(rhadapt,&
          p_rproblemLevel%rtriangulation, relementError, rcollection)

      ! Update the template matrix according to the sparsity pattern
      call grph_generateMatrix(rgraph,&
          p_rproblemLevel%Rmatrix(templateMatrix))

      ! Resize the solution vector accordingly
      call lsysbl_resizeVectorBlock(rsolutionPrimal,&
          p_rproblemLevel%Rmatrix(templateMatrix)%NEQ, .false.)

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
      call transp_initProblemLevel(rparlist, ssectionName,&
          p_rproblemLevel, rcollection, rbdrCondPrimal, rbdrCondDual)

      ! Prepare internal data arrays of the solver structure
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'systemMatrix', systemMatrix)
      call flagship_updateSolverMatrix(p_rproblemLevel, rsolver,&
          systemMatrix, SYSTEM_INTERLEAVEFORMAT, UPDMAT_ALL)
      call solver_updateStructure(rsolver)

      ! Force initialisation of the discrete transport operator and
      ! the preconditioner. This may be necessary if the velocity is
      ! constant, and hence, no repeated updates would be performed.
      p_rproblemLevel%iproblemSpec = ior(p_rproblemLevel%iproblemSpec,&
          TRANSP_TROPER_INIT+TRANSP_PRECOND_INIT)
      
      ! Stop time measurement for generation of constant coefficient matrices
      call stat_stopTimer(p_rtimerAssemblyCoeff)

    end do adaptloop


    ! Release adaptation structure
    if (nadapt > 0) then
      call hadapt_releaseAdaptation(rhadapt)
      call grph_releaseGraph(rgraph)
    end if

  end subroutine transp_solvePseudoTransPrimDual

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
    ! section name in parameter list and collection structure
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
    p_rtimerPrePostprocess => collct_getvalue_timer(rcollection,&
        'rtimerPrePostprocess', ssectionName=ssectionName)
    p_rtimerSolution => collct_getvalue_timer(rcollection,&
        'rtimerSolution', ssectionName=ssectionName)
    p_rtimerErrorEstimation => collct_getvalue_timer(rcollection,&
        'rtimerErrorEstimation', ssectionName=ssectionName)
    p_rtimerAdaptation => collct_getvalue_timer(rcollection,&
        'rtimerAdaptation', ssectionName=ssectionName)
    p_rtimerTriangulation => collct_getvalue_timer(rcollection,&
        'rtimerTriangulation', ssectionName=ssectionName)
    p_rtimerAssemblyCoeff => collct_getvalue_timer(rcollection,&
        'rtimerAssemblyCoeff', ssectionName=ssectionName)

    ! Start time measurement for pre-processing
    call stat_startTimer(p_rtimerPrePostprocess, STAT_TIMERSHORT)

    ! Adjust time stepping scheme
    rtimestep%ctimestepType = TSTEP_THETA_SCHEME
    rtimestep%dinitialTime  = 0.0_DP
    rtimestep%dinitialStep  = 1.0_DP
    rtimestep%dfinalTime    = 1.0_DP
    rtimestep%theta         = 1.0_DP

    ! Get global parameters
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'irhstype', irhstype)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'ivelocitytype', ivelocitytype)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'discretisation', discretisation)

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

    call parlst_getvalue_string(rparlist,&
        ssectionName, 'adaptivity', sadaptivityName, '')
    if (trim(adjustl(sadaptivityName)) .ne. '') then

      call parlst_getvalue_int(rparlist,&
          trim(sadaptivityName), 'nadapt', nadapt)

      if (nadapt > 0) then

        ! Initialize adaptation structure from parameter list
        call hadapt_initFromParameterlist(rhadapt, rparlist, sadaptivityName)

        ! Generate a dynamic graph for the sparsity pattern and attach
        ! it to the collection structure which is used to pass this type
        ! to the callback function for h-adaptation
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'templateMatrix', templateMatrix)
        call grph_createGraphFromMatrix(&
            p_rproblemLevel%Rmatrix(templateMatrix), rgraph)
        call collct_setvalue_graph(rcollection, 'sparsitypattern',&
            rgraph, .true., ssectionName=ssectionName)

        ! Attach the primal solution vector to the collection structure
        rcollection%p_rvectorQuickAccess1 => rsolution

      end if   ! nadapt > 0

    else

      ! No h-adaptation
      nadapt = 0

    end if

    ! Force initialisation of the discrete transport operator and the
    ! preconditioner. This may be necessary if the velocity is
    ! constant, and hence, no repeated updates would be performed.
    call problem_setSpec(rproblem, TRANSP_TROPER_INIT+&
                                   TRANSP_PRECOND_INIT,'ior')

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

        ! Solve the primal problem with non-zero right-hand side
        call tstep_performThetaStep(p_rproblemLevel, rtimestep,&
            rsolver, rsolution, transp_nlsolverCallback, rcollection,&
            rrhs)

        ! Release right-hand side vector
        call lsysbl_releaseVector(rrhs)

      else

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
      call transp_adaptTriangulation(rhadapt,&
          p_rproblemLevel%rtriangulation, relementError, rcollection)

      ! Update the template matrix according to the sparsity pattern
      call grph_generateMatrix(rgraph,&
          p_rproblemLevel%Rmatrix(templateMatrix))

      ! Resize the solution vector accordingly
      call lsysbl_resizeVectorBlock(rsolution,&
          p_rproblemLevel%Rmatrix(templateMatrix)%NEQ, .false.)

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
          p_rproblemLevel, rcollection, rbdrCond)

      ! Prepare internal data arrays of the solver structure
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'systemmatrix', systemMatrix)
      call flagship_updateSolverMatrix(p_rproblemLevel, rsolver,&
          systemMatrix, SYSTEM_INTERLEAVEFORMAT, UPDMAT_ALL)
      call solver_updateStructure(rsolver)

      ! Force initialisation of the discrete transport operator and
      ! the preconditioner. This may be necessary if the velocity is
      ! constant, and hence, no repeated updates would be performed.
      p_rproblemLevel%iproblemSpec = ior(p_rproblemLevel%iproblemSpec,&
          TRANSP_TROPER_INIT+TRANSP_PRECOND_INIT)

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

  subroutine transp_solveSteadyStatePrimDual(rparlist, ssectionName,&
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
    ! section name in parameter list and collection structure
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
    p_rtimerPrePostprocess => collct_getvalue_timer(rcollection,&
        'rtimerPrePostprocess', ssectionName=ssectionName)
    p_rtimerSolution => collct_getvalue_timer(rcollection,&
        'rtimerSolution', ssectionName=ssectionName)
    p_rtimerErrorEstimation => collct_getvalue_timer(rcollection,&
        'rtimerErrorEstimation', ssectionName=ssectionName)
    p_rtimerAdaptation => collct_getvalue_timer(rcollection,&
        'rtimerAdaptation', ssectionName=ssectionName)
    p_rtimerTriangulation => collct_getvalue_timer(rcollection,&
        'rtimerTriangulation', ssectionName=ssectionName)
    p_rtimerAssemblyCoeff => collct_getvalue_timer(rcollection,&
        'rtimerAssemblyCoeff', ssectionName=ssectionName)

    ! Start time measurement for pre-processing
    call stat_startTimer(p_rtimerPrePostprocess, STAT_TIMERSHORT)

    ! Adjust time stepping scheme
    rtimestep%ctimestepType = TSTEP_THETA_SCHEME
    rtimestep%dinitialTime  = 0.0_DP
    rtimestep%dinitialStep  = 1.0_DP
    rtimestep%dfinalTime    = 1.0_DP
    rtimestep%theta         = 1.0_DP

    ! Get global parameters
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'irhstype', irhstype)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'ivelocitytype', ivelocitytype)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'discretisation', discretisation)

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

    call parlst_getvalue_string(rparlist,&
        ssectionName, 'adaptivity', sadaptivityName, '')
    if (trim(adjustl(sadaptivityName)) .ne. '') then

      call parlst_getvalue_int(rparlist,&
          trim(sadaptivityName), 'nadapt', nadapt)

      if (nadapt > 0) then

        ! Initialize adaptation structure from parameter list
        call hadapt_initFromParameterlist(rhadapt, rparlist, sadaptivityName)

        ! Generate a dynamic graph for the sparsity pattern and attach
        ! it to the collection structure which is used to pass this type
        ! to the callback function for h-adaptation
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'templateMatrix', templateMatrix)
        call grph_createGraphFromMatrix(&
            p_rproblemLevel%Rmatrix(templateMatrix), rgraph)
        call collct_setvalue_graph(rcollection, 'sparsitypattern',&
            rgraph, .true.)

        ! Attach the primal solution vector to the collection structure
        rcollection%p_rvectorQuickAccess1 => rsolutionPrimal

      end if   ! nadapt > 0

    else

      ! No h-adaptation
      nadapt = 0

    end if

    ! Force initialisation of the discrete transport operator and the
    ! preconditioner. This may be necessary if the velocity is
    ! constant, and hence, no repeated updates would be performed.
    call problem_setSpec(rproblem, TRANSP_TROPER_INIT+&
                                   TRANSP_PRECOND_INIT,'ior')

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

        ! Solve the primal problem with non-zero right-hand side
        call tstep_performThetaStep(p_rproblemLevel, rtimestep,&
            rsolver, rsolutionPrimal, transp_nlsolverCallback,&
            rcollection, rrhs)

        ! Release right-hand side vector
        call lsysbl_releaseVector(rrhs)

      else

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

      ! Force update of the discrete transport operator and the
      ! preconditioner. This may be necessary if the velocity is
      ! constant, and hence, no repeated updates would be performed.
      call problem_setSpec(rproblem, TRANSP_TROPER_UPDATE+&
                                     TRANSP_PRECOND_UPDATE,'ior')

      ! Reset the time-stepping and solver algorithms
      call tstep_resetTimestep(rtimestep, .false.)
      call solver_resetSolver(rsolver, .false.)

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

        ! Compute the error in the quantity of interest
        call transp_estimateTargetFuncError(rparlist, ssectionName,&
            p_rproblemLevel, rtimestep, rsolver, rsolutionPrimal,&
            rsolutionDual, rcollection, relementError, derror, rrhs)

        ! Release right-hand side vector
        call lsysbl_releaseVector(rrhs)

      else

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
      call transp_adaptTriangulation(rhadapt,&
          p_rproblemLevel%rtriangulation, relementError, rcollection)

      ! Update the template matrix according to the sparsity pattern
      call grph_generateMatrix(rgraph, p_rproblemLevel%Rmatrix(templateMatrix))

      ! Resize the solution vector accordingly
      call lsysbl_resizeVectorBlock(rsolutionPrimal,&
          p_rproblemLevel%Rmatrix(templateMatrix)%NEQ, .false.)

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
      call transp_initProblemLevel(rparlist, ssectionName,&
          p_rproblemLevel, rcollection, rbdrCondPrimal, rbdrCondDual)

      ! Prepare internal data arrays of the solver structure
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'systemMatrix', systemMatrix)
      call flagship_updateSolverMatrix(p_rproblemLevel, rsolver,&
          systemMatrix, SYSTEM_INTERLEAVEFORMAT, UPDMAT_ALL)
      call solver_updateStructure(rsolver)

      ! Force initialisation of the discrete transport operator and
      ! the preconditioner. This may be necessary if the velocity is
      ! constant, and hence, no repeated updates would be performed.
      p_rproblemLevel%iproblemSpec = ior(p_rproblemLevel%iproblemSpec,&
          TRANSP_TROPER_INIT+TRANSP_PRECOND_INIT)
      
      ! Stop time measurement for generation of constant coefficient matrices
      call stat_stopTimer(p_rtimerAssemblyCoeff)

    end do adaptloop


    ! Release adaptation structure
    if (nadapt > 0) then
      call hadapt_releaseAdaptation(rhadapt)
      call grph_releaseGraph(rgraph)
    end if

  end subroutine transp_solveSteadyStatePrimDual

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

end module transport_application
