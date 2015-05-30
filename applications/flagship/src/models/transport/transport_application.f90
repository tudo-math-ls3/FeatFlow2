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
!# 2.) transp_solveTransientPrimal
!#     -> Solves the primal formulation of the time-dependent
!#        convection-diffusion-reaction equation.
!#
!# 3.) transp_solveTransientPrimDual
!#     -> Solves the primal and the dual formulation of the time-dependent
!#        convection-diffusion-reaction equation using pseudo time-stepping.
!#
!# 4.) transp_solvePseudoTransPrimal
!#     -> Solves the primal formulation of the steady
!#        convection-diffusion-reaction equation using pseudo time-stepping.
!#
!# 5.) transp_solvePseudoTransPrimDual
!#     -> Solves the primal and the dual formulation of the steady
!#        convection-diffusion-reaction equation using pseudo time-stepping.
!#
!# 6.) transp_solveSteadyStatePrimal
!#     -> Solves the primal formulation of the steady
!#        convection-diffusion-reaction equation directly
!#
!# 7.) transp_solveSteadyStatePrimDual
!#     -> Solves the primal and the dual formulation of the steady
!#        convection-diffusion-reaction equation directly
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

#include "flagship.h"

!$ use omp_lib
  use boundarycondaux
  use boundaryfilter
  use collection
  use flagship_basic
  use flagship_signals
  use fparser
  use fsystem
  use genoutput
  use graph
  use hadaptaux
  use hadaptivity
  use linearalgebra
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use problem
  use solverbase
  use spatialdiscretisation
  use statistics
  use storage
  use timestep
  use timestepbase
  use triangulation
  use ucd

  ! Modules from transport model
  use transport_basic
  use transport_callback
  use transport_errorestimation
  use transport_meshadaptation
  use transport_postprocessing
  use transport_preprocessing

  implicit none

  private

  public :: transp_app
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

    ! Timer for file IO
    type(t_timer) :: rtimerFileIO

    ! Abstract problem descriptor
    type(t_problemDescriptor) :: rproblemDescriptor

    ! Parameter file and section names
    character(LEN=SYS_STRLEN) :: sindatfileName
    character(LEN=SYS_STRLEN) :: sbdrcondName
    character(LEN=SYS_STRLEN) :: algorithm
    character(LEN=SYS_STRLEN) :: benchmark

    ! local variables
    real(DP) :: derrorL1,derrorL2,derrorH1,derrorDispersion
    integer :: systemMatrix, ndimension, version


    ! Check if configuration file is valid
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'version', version, 0)
    if (version .lt. TRANSP_VERSION) then
      call output_lbrk()
      call output_separator(OU_SEP_DOLLAR,OU_CLASS_ERROR,OU_MODE_STD)
      call output_line('This configuration file is deprecated!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_app')
      call output_line('Please wait for an updated version of Featflow2 or contact',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_app')
      call output_line('matthias.moeller@math.tu-dortmund.de',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_app')
      call output_line('if running this benchmark is urgent!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_app')
      call output_separator(OU_SEP_DOLLAR,OU_CLASS_ERROR,OU_MODE_STD)
      call output_lbrk()
      call sys_halt()
    end if

    ! Start total time measurement
    call stat_startTimer(rtimerTotal)

    !---------------------------------------------------------------------------
    ! Pre-processing
    !---------------------------------------------------------------------------

    ! Start time measurement
    call stat_startTimer(rtimerPrePostprocess, STAT_TIMERSHORT)

    ! Initialise global collection structure
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
    call collct_setvalue_timer(rcollection,&
        'rtimerFileIO', rtimerFileIO, .true.,&
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

    ! Initialise the solver structures
    call transp_initSolvers(rparlist, ssectionName, rtimestep, rsolver)

    ! Initialise the boundary conditions for the primal problem
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'ndimension', ndimension)
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'sprimalbdrcondname', sbdrcondName)
    call bdrc_readBoundaryCondition(rbdrCondPrimal,&
        sindatfileName, sbdrcondName,&
        ndimension, transp_parseBoundaryCondition)

    ! Initialise the boundary conditions for the dual problem
    call parlst_getvalue_string(rparlist, ssectionName,&
        'sdualbdrcondname', sbdrcondName, '')
    if (sbdrcondName .ne. '') then
      call bdrc_readBoundaryCondition(rbdrCondDual,&
          sindatfileName, sbdrcondName,&
          ndimension, transp_parseBoundaryCondition)
    end if
    
    ! Initialise the abstract problem structure
    call transp_initProblemDescriptor(rparlist, ssectionName,&
        solver_getMinimumMultigridlevel(rsolver),&
        solver_getMaximumMultigridlevel(rsolver),&
        rproblemDescriptor)
    call problem_initProblem(rproblem, rproblemDescriptor)

    ! Initialise all individual problem levels with primal and dual
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
      call parlst_getvalue_string(rparlist,&
          ssectionName, 'benchmark', benchmark, '')
      
      ! What solution algorithm should be applied?
      if (trim(algorithm) .eq. 'transient_primal') then
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Solve the primal formulation for
        ! the time-dependent problem
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call transp_solveTransientPrimal(rparlist, ssectionName,&
            rbdrCondPrimal, rproblem, rtimestep, rsolver,&
            rsolutionPrimal, rcollection)

        ! Estimate the error to the exact solution
        call transp_errestExact(rparlist, ssectionName,&
            rproblem%p_rproblemLevelMax, rsolutionPrimal,&
            rtimestep%dTime, derrorL1=derrorL1, derrorL2=derrorL2,&
            derrorH1=derrorH1, rcollection=rcollection)

        ! Start time measurement for post-processing
        call stat_startTimer(rtimerFileIO, STAT_TIMERSHORT)

        ! Output solution to file
        call transp_outputSolution(rparlist, ssectionName,&
            rproblem%p_rproblemLevelMax,&
            rsolutionPrimal=rsolutionPrimal,&
            dtime=rtimestep%dTime)

        ! Stop time measurement for post-processing
        call stat_stopTimer(rtimerFileIO)

      elseif (trim(algorithm) .eq. 'transient_primaldual') then
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Solve the primal and dual formulation for
        ! the time-dependent problem
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call transp_solveTransientPrimDual(rparlist, ssectionName,&
            rbdrCondPrimal, rbdrCondDual, rproblem, rtimestep,&
            rsolver, rsolutionPrimal, rsolutionDual, rcollection)

        ! Estimate the error to the exact solution
        call transp_errestExact(rparlist, ssectionName,&
            rproblem%p_rproblemLevelMax, rsolutionPrimal,&
            rtimestep%dTime, derrorL1=derrorL1, derrorL2=derrorL2,&
            derrorH1=derrorH1, rcollection=rcollection)

        ! Start time measurement for post-processing
        call stat_startTimer(rtimerFileIO, STAT_TIMERSHORT)

        ! Output solution to file
        call transp_outputSolution(rparlist, ssectionName,&
            rproblem%p_rproblemLevelMax,&
            rsolutionPrimal, rsolutionDual, rtimestep%dTime)

        ! Stop time measurement for post-processing
        call stat_stopTimer(rtimerFileIO)

      elseif (trim(algorithm) .eq. 'pseudotransient_primal') then
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Solve the primal formulation for
        ! the pseudo time-dependent problem
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call transp_solvePseudoTransPrimal(rparlist, ssectionName,&
            rbdrCondPrimal, rproblem, rtimestep, rsolver,&
            rsolutionPrimal, rcollection)

        ! Estimate the error to the exact solution
        call transp_errestExact(rparlist, ssectionName,&
            rproblem%p_rproblemLevelMax, rsolutionPrimal,&
            rtimestep%dTime, derrorL1=derrorL1, derrorL2=derrorL2,&
            derrorH1=derrorH1, rcollection=rcollection)

        ! Start time measurement for post-processing
        call stat_startTimer(rtimerFileIO, STAT_TIMERSHORT)

        ! Output solution to file
        call transp_outputSolution(rparlist, ssectionName,&
            rproblem%p_rproblemLevelMax,&
            rsolutionPrimal=rsolutionPrimal,&
            dtime=rtimestep%dTime)

        ! Stop time measurement for post-processing
        call stat_stopTimer(rtimerFileIO)

      elseif (trim(algorithm) .eq. 'pseudotransient_primaldual') then
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Solve the primal and dual formulation for
        ! the pseudo time-dependent problem
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call transp_solvePseudoTransPrimDual(rparlist, ssectionName,&
            rbdrCondPrimal, rbdrCondDual, rproblem, rtimestep,&
            rsolver, rsolutionPrimal, rsolutionDual, rcollection)

        ! Estimate the error to the exact solution
        call transp_errestExact(rparlist, ssectionName,&
            rproblem%p_rproblemLevelMax, rsolutionPrimal,&
            rtimestep%dTime, derrorL1=derrorL1, derrorL2=derrorL2,&
            derrorH1=derrorH1, rcollection=rcollection)

        ! Start time measurement for post-processing
        call stat_startTimer(rtimerFileIO, STAT_TIMERSHORT)

        ! Output solution to file
        call transp_outputSolution(rparlist, ssectionName,&
            rproblem%p_rproblemLevelMax,&
            rsolutionPrimal, rsolutionDual, rtimestep%dTime)

        ! Stop time measurement for post-processing
        call stat_stopTimer(rtimerFileIO)

      elseif (trim(algorithm) .eq. 'stationary_primal') then
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Solve the primal formulation for
        ! the stationary problem
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call transp_solveSteadyStatePrimal(rparlist, ssectionName,&
            rbdrCondPrimal, rproblem, rtimestep, rsolver,&
            rsolutionPrimal, rcollection)
        
        ! Estimate the error to the exact solution
        call transp_errestExact(rparlist, ssectionName,&
            rproblem%p_rproblemLevelMax, rsolutionPrimal,&
            rtimestep%dTime, derrorL1=derrorL1, derrorL2=derrorL2,&
            derrorH1=derrorH1, rcollection=rcollection)

!!$        ! Estimate the error to the exact solution
!!$        call transp_errestExactBoundary2D(rparlist, ssectionName,&
!!$            rproblem%p_rproblemLevelMax, rsolutionPrimal,&
!!$            rtimestep%dTime, derrorL1=derrorL1, derrorL2=derrorL2,&
!!$            derrorH1=derrorH1, rcollection=rcollection)

        ! Start time measurement for post-processing
        call stat_startTimer(rtimerFileIO, STAT_TIMERSHORT)

        ! Output solution to file
        call transp_outputSolution(rparlist, ssectionName,&
            rproblem%p_rproblemLevelMax,&
            rsolutionPrimal=rsolutionPrimal,&
            dtime=rtimestep%dTime)

        ! Stop time measurement for post-processing
        call stat_stopTimer(rtimerFileIO)

      elseif (trim(algorithm) .eq. 'stationary_primaldual') then
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Solve the primal and dual formulation for
        ! the stationary problem
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call transp_solveSteadyStatePrimDual(rparlist, ssectionName,&
            rbdrCondPrimal, rbdrCondDual, rproblem, rtimestep,&
            rsolver, rsolutionPrimal, rsolutionDual, rcollection)

        ! Estimate the error to the exact solution
        call transp_errestExact(rparlist, ssectionName,&
            rproblem%p_rproblemLevelMax, rsolutionPrimal,&
            rtimestep%dTime, derrorL1=derrorL1, derrorL2=derrorL2,&
            derrorH1=derrorH1, rcollection=rcollection)

        ! Start time measurement for post-processing
        call stat_startTimer(rtimerFileIO, STAT_TIMERSHORT)

        ! Output solution to file
        call transp_outputSolution(rparlist, ssectionName,&
            rproblem%p_rproblemLevelMax,&
            rsolutionPrimal, rsolutionDual, rtimestep%dTime)

        ! Stop time measurement for post-processing
        call stat_stopTimer(rtimerFileIO)

      else
        call output_line(trim(algorithm)//' is not a valid solution algorithm!',&
            OU_CLASS_ERROR,OU_MODE_STD,'transp_app')
        call sys_halt()
      end if
      
      ! Perform benchmark specific checks
      if (trim(adjustl(benchmark)) .eq. 'Gaussian-Hill-2D') then
        call transp_errestDispersionGHill(rparlist, ssectionName,&
            rproblem%p_rproblemLevelMax, rsolutionPrimal,&
            rtimestep%dTime, derrorDispersion, rcollection)
      end if

    else

      ! Start time measurement for post-processing
      call stat_startTimer(rtimerFileIO, STAT_TIMERSHORT)
        
      ! Just output the computational mesh and exit
      call transp_outputSolution(rparlist, ssectionName,&
          rproblem%p_rproblemLevelMax)

      ! Stop time measurement for post-processing
      call stat_stopTimer(rtimerFileIO)

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

    ! solver structure
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
    type(t_timer), pointer :: p_rtimerAssemblyVector
    type(t_timer), pointer :: p_rtimerFileIO

    ! section names
    character(LEN=SYS_STRLEN) :: sadaptivityName
    character(LEN=SYS_STRLEN) :: soutputName
    character(LEN=SYS_STRLEN) :: sucdimport
    
    ! local variables
    real(DP) :: dnormL1, dnormL2
    type(t_ucdExport) :: rimport
    real(dp) :: derror, dstepUCD, dtimeUCD, dstepAdapt, dtimeAdapt
    integer :: templateMatrix, systemMatrix, discretisation
    integer :: nlmin, ipreadapt, npreadapt, irhstype, ivelocitytype

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
    p_rtimerAssemblyVector => collct_getvalue_timer(rcollection,&
        'rtimerAssemblyVector', ssectionName=ssectionName)
    p_rtimerFileIO => collct_getvalue_timer(rcollection,&
        'rtimerFileIO', ssectionName=ssectionName)

    ! Start time measurement for pre-processing
    call stat_startTimer(p_rtimerPrePostprocess, STAT_TIMERSHORT)

    ! Get global parameters
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'irhstype', irhstype)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'ivelocitytype', ivelocitytype)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'discretisation', discretisation)

    ! Set pointer to maximum problem level and block discretisation
    p_rproblemLevel   => rproblem%p_rproblemLevelMax
    p_rdiscretisation => p_rproblemLevel%RblockDiscretisation(discretisation)

    ! Create the solution vector
    call lsysbl_createVectorBlock(p_rdiscretisation,&
        rsolution, .false., ST_DOUBLE)

    ! Initialise the solution vector and impose boundary conditions
    call transp_initSolution(rparlist, ssectionName, p_rproblemLevel,&
        rtimestep%dinitialTime, rsolution, rcollection)
    call bdrf_filterVectorExplicit(rbdrCond, rsolution,&
        rtimestep%dinitialTime)

    ! Initialise timer for intermediate UCD exporter
    dtimeUCD = rtimestep%dinitialTime
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'output', soutputName)
    call parlst_getvalue_double(rparlist,&
        trim(soutputName), 'dstepUCD', dstepUCD, 0.0_DP)


    !---------------------------------------------------------------------------
    ! Initialise the h-adaptation structure and perform pre-adaptation
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

        ! Initialise adaptation structure from parameter list
        call hadapt_initFromParameterlist(rhadapt, rparlist, sadaptivityName)

        ! Generate a dynamic graph for the sparsity pattern and attach
        ! it to the collection structure which is used to pass this
        ! type to the callback function for h-adaptation
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'templateMatrix', templateMatrix)
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'systemMatrix', systemMatrix)
        call lsyssc_createGraphFromMatrix(&
            p_rproblemLevel%RmatrixScalar(templateMatrix), rgraph)
        call collct_setvalue_graph(rcollection, 'sparsitypattern',&
            rgraph, .true., ssectionName=ssectionName)

        ! Perform pre-adaptation?
        if (npreadapt > 0) then

          ! Perform number of pre-adaptation steps
          do ipreadapt = 1, npreadapt

            ! Compute the error estimator using recovery techniques
            call transp_errestRecovery(rparlist, ssectionname,&
                p_rproblemLevel, rsolution, rtimestep%dinitialTime,&
                relementError, derror, rcollection)

            ! Set the name of the template matrix
            rcollection%SquickAccess(1) = 'sparsitypattern'

            ! Attach the primal solution vector to the collection structure
            rcollection%p_rvectorQuickAccess1 => rsolution

            ! Perform h-adaptation and update the triangulation structure
            call transp_adaptTriangulation(&
                rhadapt, p_rproblemLevel%rtriangulation,&
                relementError, rcollection)

            ! Release element-wise error distribution
            call lsyssc_releaseVector(relementError)

            ! Generate standard mesh from raw mesh
            call tria_initStandardMeshFromRaw(&
                p_rproblemLevel%rtriangulation, rproblem%rboundary)

            ! Update the template matrix according to the sparsity pattern
            call lsyssc_createMatrixFromGraph(rgraph,&
                p_rproblemLevel%RmatrixScalar(templateMatrix))
            
            ! Re-initialise all constant coefficient matrices
            call transp_initProblemLevel(rparlist, ssectionName,&
                p_rproblemLevel, rcollection, rbdrCond)
            
            ! Resize the solution vector accordingly
            call lsysbl_resizeVectorBlock(p_rdiscretisation, &
                rsolution, .false.)
            
            ! Re-generate the initial solution vector and impose
            !  boundary conditions explicitly
            call transp_initSolution(rparlist, ssectionname,&
                p_rproblemLevel, rtimestep%dinitialTime, rsolution,&
                rcollection)
            call bdrf_filterVectorExplicit(rbdrCond, rsolution,&
                rtimestep%dinitialTime)
          end do

          ! Prepare internal data arrays of the solver structure
          call flagship_updateSolverMatrix(p_rproblemLevel, rsolver,&
              systemMatrix, SYSTEM_INTERLEAVEFORMAT, UPDMAT_ALL) 
         call solver_updateStructure(rsolver)

        end if   ! npreadapt > 0

      end if   ! dstepAdapt > 0

    else

      dstepAdapt = 0.0_DP

    end if

    ! Get name of import file (if any)
    call parlst_getvalue_string(rparlist,&
        trim(soutputName), 'sucdimport', sucdimport, '')
    
    ! Do we have to read in a pre-computed solution?
    if (trim(sucdimport) .ne. '') then
      call ucd_readGMV(sucdimport, rimport, p_rproblemLevel%rtriangulation)
      call ucd_getSimulationTime(rimport, rtimestep%dinitialTime)
      call ucd_getSimulationTime(rimport, rtimestep%dTime)
      call transp_setVariable(rimport, 'u', rsolution)
      call ucd_release(rimport)
      
      ! Set time for solution output
      dtimeUCD = rtimestep%dinitialTime
    end if

    ! Force initialisation of the discrete transport operator and the
    ! preconditioner. This may be necessary if the velocity is
    ! constant, and hence, no repeated updates would be performed.
    call problem_setSpec(rproblem, TRANSP_TROPER_INIT+&
                                   TRANSP_PRECOND_INIT,'ior')

    ! Initialise right-hand side vector
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
      if (flagship_SIGINT(-1) > 0 ) then

        ! Start time measurement for post-processing
        call stat_startTimer(p_rtimerFileIO, STAT_TIMERSHORT)

        call transp_outputSolution(rparlist, ssectionName,&
            p_rproblemLevel, rsolution, dtime=rtimestep%dTime)

        ! Stop time measurement for post-processing
        call stat_stopTimer(p_rtimerFileIO)
      end if

      !-------------------------------------------------------------------------
      ! Advance solution in time
      !-------------------------------------------------------------------------

      ! Start time measurement for solution procedure
      call stat_startTimer(p_rtimerSolution, STAT_TIMERSHORT)

      ! Perform a single time step
      if (irhstype > 0) then
        call tstep_performTimestep(p_rproblemLevel, rtimestep,&
            rsolver, rsolution, transp_nlsolverCallback,&
            rcollection, 1, rrhs)
      else
        call tstep_performTimestep(p_rproblemLevel, rtimestep,&
            rsolver, rsolution, transp_nlsolverCallback,&
            rcollection, 1)
      end if
      
      ! Stop time measurement for solution procedure
      call stat_stopTimer(p_rtimerSolution)

      ! Start time measurement for vector assembly
      call stat_startTimer(p_rtimerAssemblyVector)

      ! Perform linearised FEM-FCT post-processing
      call transp_calcLinearisedFCT(p_rproblemLevel, rtimestep, rsolver,&
          rsolution, ssectionName, rcollection, rmatrix=rmatrix1,&
          rvector1=rvector1, rvector2=rvector2, rvector3=rvector3)

      ! Perform linearised FEM-LPT post-processing
      call transp_calcLinearisedLPT(p_rproblemLevel, rtimestep, rsolver,&
          rsolution, ssectionName, rcollection, rmatrix=rmatrix1,&
          rvector1=rvector1, rvector2=rvector2, rvector3=rvector3)

      ! Stop time measurement for vector assembly
      call stat_stopTimer(p_rtimerAssemblyVector)

      ! Write norm of solution to benchmark logfile
      if (OU_BENCHLOG .ne. 0) then
        dnormL1 = lsysbl_vectorNorm(rsolution, LINALG_NORML1)
        dnormL2 = lsysbl_vectorNorm(rsolution, LINALG_NORML2)
        call output_line('T = '//trim(sys_sdEL(rtimestep%dTime,5))//&
            '   ||u||_1 = '//trim(sys_sdEL(dnormL1,5))//&
            '   ||u||_2 = '//trim(sys_sdEL(dnormL2,5)),&
            OU_CLASS_MSG, OU_MODE_BENCHLOG)
      end if

      ! Reached final time, then exit the infinite time loop?
      if (rtimestep%dTime .ge. rtimestep%dfinalTime) exit timeloop


      !-------------------------------------------------------------------------
      ! Post-process intermediate solution
      !-------------------------------------------------------------------------

      if ((dstepUCD .gt. 0.0_DP) .and. (rtimestep%dTime .ge. dtimeUCD)) then

        ! Set time for next intermediate solution export
        dtimeUCD = dtimeUCD + dstepUCD

        ! Start time measurement for post-processing
        call stat_startTimer(p_rtimerFileIO, STAT_TIMERSHORT)

        ! Export the intermediate solution
        call transp_outputSolution(rparlist, ssectionName,&
            p_rproblemLevel, rsolution, dtime=rtimestep%dTime)

        ! Stop time measurement for post-processing
        call stat_stopTimer(p_rtimerFileIO)

      end if


      !-------------------------------------------------------------------------
      ! Perform h-adaptation
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
        call transp_errestRecovery(rparlist, ssectionname,&
            p_rproblemLevel, rsolution, rtimestep%dTime,&
            relementError, derror, rcollection)

        ! Stop time measurement for error estimation
        call stat_stopTimer(p_rtimerErrorEstimation)


        !-----------------------------------------------------------------------
        ! Perform h-adaptation
        !-----------------------------------------------------------------------

        ! Start time measurement for h-adaptation
        call stat_startTimer(p_rtimerAdaptation, STAT_TIMERSHORT)

        ! Set the name of the template matrix
        rcollection%SquickAccess(1) = 'sparsitypattern'

        ! Attach the primal solution vector to the collection structure
        rcollection%p_rvectorQuickAccess1 => rsolution

        ! Perform h-adaptation and update the triangulation structure
        call transp_adaptTriangulation(&
	   rhadapt, p_rproblemLevel%rtriangulation,&
	   relementError, rcollection)

        ! Release element-wise error distribution
        call lsyssc_releaseVector(relementError)

        ! Update the template matrix according to the sparsity pattern
        call lsyssc_createMatrixFromGraph(rgraph,&
            p_rproblemLevel%RmatrixScalar(templateMatrix))

        ! Resize the solution vector accordingly
        call lsysbl_resizeVectorBlock(p_rdiscretisation,&
            rsolution, .false.)

        ! Stop time measurement for h-adaptation
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

        ! Re-initialise all constant coefficient matrices
        call transp_initProblemLevel(rparlist, ssectionName,&
            p_rproblemLevel, rcollection, rbdrCond)

        ! Prepare internal data arrays of the solver structure
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

        ! Re-initialise the right-hand side vector
        if (irhstype > 0) then
          call lsysbl_resizeVectorBlock(p_rdiscretisation,&
              rrhs, .false.)
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
    real(DP) :: derror, dnormL1, dnormL2
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

    ! Set pointer to maximum problem level and block discretisation
    p_rproblemLevel => rproblem%p_rproblemLevelMax
    p_rdiscretisation => p_rproblemLevel%RblockDiscretisation(discretisation)

    ! Create the solution vector
    call lsysbl_createVectorBlock(p_rdiscretisation, rsolution,&
        .false., ST_DOUBLE)

    ! Initialise the solution vector and impose boundary conditions
    call transp_initSolution(rparlist, ssectionName, p_rproblemLevel,&
        rtimestep%dinitialTime, rsolution, rcollection)
    call bdrf_filterVectorExplicit(rbdrCond, rsolution,&
        rtimestep%dinitialTime)

    !---------------------------------------------------------------------------
    ! Initialise the h-adaptation structure
    !---------------------------------------------------------------------------

    call parlst_getvalue_string(rparlist,&
        ssectionName, 'adaptivity', sadaptivityName, '')
    if (trim(adjustl(sadaptivityName)) .ne. '') then

      call parlst_getvalue_int(rparlist,&
          trim(sadaptivityName), 'nadapt', nadapt)

      if (nadapt > 0) then

        ! Initialise adaptation structure from parameter list
        call hadapt_initFromParameterlist(rhadapt, rparlist, sadaptivityName)

        ! Generate a dynamic graph for the sparsity pattern and attach
        ! it to the collection structure which is used to pass this type
        ! to the callback function for h-adaptation
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'templateMatrix', templateMatrix)
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'systemMatrix', systemMatrix)
        call lsyssc_createGraphFromMatrix(&
            p_rproblemLevel%RmatrixScalar(templateMatrix), rgraph)
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
        call tstep_performTimestep(p_rproblemLevel, rtimestep,&
            rsolver, rsolution, transp_nlsolverCallback,&
            rcollection, SYS_INFINITY_INT, rrhs)

        ! Release right-hand side vector
        call lsysbl_releaseVector(rrhs)

      else

        ! Solve the primal problem without right-hand side
        call tstep_performTimestep(p_rproblemLevel, rtimestep,&
            rsolver, rsolution, transp_nlsolverCallback,&
            rcollection, SYS_INFINITY_INT)
      end if

      ! Stop time measurement for solution procedure
      call stat_stopTimer(p_rtimerSolution)

      ! Write norm of solution to benchmark logfile
      if (OU_BENCHLOG .ne. 0) then
        dnormL1 = lsysbl_vectorNorm(rsolution, LINALG_NORML1)
        dnormL2 = lsysbl_vectorNorm(rsolution, LINALG_NORML2)
        call output_line('T = '//trim(sys_sdEL(rtimestep%dTime,5))//&
            '   ||u||_1 = '//trim(sys_sdEL(dnormL1,5))//&
            '   ||u||_2 = '//trim(sys_sdEL(dnormL2,5)),&
            OU_CLASS_MSG, OU_MODE_BENCHLOG)
      end if

      if (iadapt .eq. nadapt) exit adaptloop

      !-------------------------------------------------------------------------
      ! Perform recovery-based error estimation
      !-------------------------------------------------------------------------

      ! Start time measurement for error estimation
      call stat_startTimer(p_rtimerErrorEstimation, STAT_TIMERSHORT)

      ! Compute the error estimator using recovery techniques
      call transp_errestRecovery(rparlist, ssectionname,&
          p_rproblemLevel, rsolution, 0.0_DP, relementError, derror,&
          rcollection)

      ! Stop time measurement for error estimation
      call stat_stopTimer(p_rtimerErrorEstimation)


      !-------------------------------------------------------------------------
      ! Perform h-adaptation
      !-------------------------------------------------------------------------

      ! Start time measurement for h-adaptation
      call stat_startTimer(p_rtimerAdaptation, STAT_TIMERSHORT)

      ! Set the name of the template matrix
      rcollection%SquickAccess(1) = 'sparsitypattern'

      ! Attach the primal solution vector to the collection structure
      rcollection%p_rvectorQuickAccess1 => rsolution

      ! Perform h-adaptation and update the triangulation structure
      call transp_adaptTriangulation(&
          rhadapt, p_rproblemLevel%rtriangulation,&
          relementError, rcollection)

      ! Update the template matrix according to the sparsity pattern
      call lsyssc_createMatrixFromGraph(rgraph,&
          p_rproblemLevel%RmatrixScalar(templateMatrix))

      ! Resize the solution vector accordingly
      call lsysbl_resizeVectorBlock(p_rdiscretisation,&
          rsolution, .false.)

      ! Release element-wise error distribution
      call lsyssc_releaseVector(relementError)

      ! Stop time measurement for h-adaptation
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

      ! Re-initialise all constant coefficient matrices
      call transp_initProblemLevel(rparlist, ssectionName,&
          p_rproblemLevel, rcollection, rbdrCond)

      ! Prepare internal data arrays of the solver structure
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
    real(dp) :: derror, dnormL1, dnormL2
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

    ! Set pointer to maximum problem level and block discretisation
    p_rproblemLevel   => rproblem%p_rproblemLevelMax
    p_rdiscretisation => p_rproblemLevel%RblockDiscretisation(discretisation)

    ! Initialise the solution vector and impose boundary conditions explicitly
    call lsysbl_createVectorBlock(p_rdiscretisation, rsolutionPrimal,&
        .false., ST_DOUBLE)
    call transp_initSolution(rparlist, ssectionName, p_rproblemLevel,&
        0.0_DP, rsolutionPrimal, rcollection)
    call bdrf_filterVectorExplicit(rbdrCondPrimal, rsolutionPrimal, 0.0_DP)

    !---------------------------------------------------------------------------
    ! Initialise the h-adaptation structure
    !---------------------------------------------------------------------------

    call parlst_getvalue_string(rparlist,&
        ssectionName, 'adaptivity', sadaptivityName, '')
    if (trim(adjustl(sadaptivityName)) .ne. '') then

      call parlst_getvalue_int(rparlist,&
          trim(sadaptivityName), 'nadapt', nadapt)

      if (nadapt > 0) then

        ! Initialise adaptation structure from parameter list
        call hadapt_initFromParameterlist(rhadapt, rparlist, sadaptivityName)

        ! Generate a dynamic graph for the sparsity pattern and attach
        ! it to the collection structure which is used to pass this type
        ! to the callback function for h-adaptation
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'templateMatrix', templateMatrix)
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'systemMatrix', systemMatrix)
        call lsyssc_createGraphFromMatrix(&
            p_rproblemLevel%RmatrixScalar(templateMatrix), rgraph)
        call collct_setvalue_graph(rcollection, 'sparsitypattern',&
            rgraph, .true., ssectionName=ssectionName)

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
        call tstep_performTimestep(p_rproblemLevel, rtimestep,&
            rsolver, rsolutionPrimal, transp_nlsolverCallback,&
            rcollection, SYS_INFINITY_INT, rrhs)

        ! Release right-hand side vector
        call lsysbl_releaseVector(rrhs)

      else

        ! Solve the primal problem without right-hand side
        call tstep_performTimestep(p_rproblemLevel, rtimestep,&
            rsolver, rsolutionPrimal, transp_nlsolverCallback,&
            rcollection, SYS_INFINITY_INT)
      end if

      ! Stop time measurement for solution procedure
      call stat_stopTimer(p_rtimerSolution)

      ! Write norm of primal solution to benchmark logfile
      if (OU_BENCHLOG .ne. 0) then
        dnormL1 = lsysbl_vectorNorm(rsolutionPrimal, LINALG_NORML1)
        dnormL2 = lsysbl_vectorNorm(rsolutionPrimal, LINALG_NORML2)
        call output_line('T = '//trim(sys_sdEL(rtimestep%dTime,5))//&
            '   ||u||_1 = '//trim(sys_sdEL(dnormL1,5))//&
            '   ||u||_2 = '//trim(sys_sdEL(dnormL2,5)),&
            OU_CLASS_MSG, OU_MODE_BENCHLOG)
      end if


      !-------------------------------------------------------------------------
      ! Compute the right-hand side for the dual problem
      !-------------------------------------------------------------------------

      ! Start time measurement for error estimation
      call stat_startTimer(p_rtimerErrorEstimation, STAT_TIMERSHORT)

      ! Initialise target functional
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
      call tstep_performTimestep(p_rproblemLevel, rtimestep,&
          rsolver, rsolutionDual, transp_nlsolverCallback,&
          rcollection, SYS_INFINITY_INT, rtargetFunc)

      ! Release discretised target functional
      call lsysbl_releaseVector(rtargetFunc)

      ! Stop time measurement for solution procedure
      call stat_stopTimer(p_rtimerSolution)

      ! Write norm of dual solution to benchmark logfile
      if (OU_BENCHLOG .ne. 0) then
        dnormL1 = lsysbl_vectorNorm(rsolutionDual, LINALG_NORML1)
        dnormL2 = lsysbl_vectorNorm(rsolutionDual, LINALG_NORML2)
        call output_line('T = '//trim(sys_sdEL(rtimestep%dTime,5))//&
            '   ||z||_1 = '//trim(sys_sdEL(dnormL1,5))//&
            '   ||z||_2 = '//trim(sys_sdEL(dnormL2,5)),&
            OU_CLASS_MSG, OU_MODE_BENCHLOG)
      end if


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

        ! Initialise right-hand side vector
        call lsysbl_createVectorBlock(rsolutionPrimal, rrhs)
        call transp_initRHS(rparlist, ssectionName, p_rproblemLevel,&
            0.0_DP, rrhs, rcollection)

        ! Compute the error in the quantity of interest
        call transp_errestTargetFunc(rparlist, ssectionName,&
            p_rproblemLevel, rtimestep, rsolver, rsolutionPrimal,&
            rsolutionDual, rcollection, relementError, derror, rrhs)

        ! Release right-hand side vector
        call lsysbl_releaseVector(rrhs)

      else

        ! Compute the error in the quantity of interest
        call transp_errestTargetFunc(rparlist, ssectionName,&
            p_rproblemLevel, rtimestep, rsolver, rsolutionPrimal,&
            rsolutionDual, rcollection, relementError, derror)

      end if

      ! Stop time measurement for error estimation
      call stat_stopTimer(p_rtimerErrorEstimation)


      !-------------------------------------------------------------------------
      ! Perform h-adaptation
      !-------------------------------------------------------------------------

      ! Start time measurement for h-adaptation
      call stat_startTimer(p_rtimerAdaptation, STAT_TIMERSHORT)

      ! Set the name of the template matrix
      rcollection%SquickAccess(1) = 'sparsitypattern'

      ! Attach the primal solution vector to the collection structure
      rcollection%p_rvectorQuickAccess1 => rsolutionPrimal

      ! Perform h-adaptation and update the triangulation structure
      call transp_adaptTriangulation(&
          rhadapt, p_rproblemLevel%rtriangulation,&
          relementError, rcollection)

      ! Update the template matrix according to the sparsity pattern
      call lsyssc_createMatrixFromGraph(rgraph,&
          p_rproblemLevel%RmatrixScalar(templateMatrix))

      ! Resize the solution vector accordingly
      call lsysbl_resizeVectorBlock(p_rdiscretisation,&
          rsolutionPrimal, .false.)

      ! Release element-wise error distribution
      call lsyssc_releaseVector(relementError)

      ! Stop time measurement for h-adaptation
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

      ! Re-initialise all constant coefficient matrices
      call transp_initProblemLevel(rparlist, ssectionName,&
          p_rproblemLevel, rcollection, rbdrCondPrimal, rbdrCondDual)

      ! Prepare internal data arrays of the solver structure
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
    real(dp) :: derror, dnormL1, dnormL2
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

    ! Set pointer to maximum problem level and block discretisation
    p_rproblemLevel   => rproblem%p_rproblemLevelMax
    p_rdiscretisation => p_rproblemLevel%RblockDiscretisation(discretisation)

    ! Initialise the solution vector and impose boundary conditions explicitly
    call lsysbl_createVectorBlock(p_rdiscretisation, rsolution, .false., ST_DOUBLE)
    call transp_initSolution(rparlist, ssectionName, p_rproblemLevel,&
        0.0_DP, rsolution, rcollection)
    call bdrf_filterVectorExplicit(rbdrCond, rsolution, 0.0_DP)

    !---------------------------------------------------------------------------
    ! Initialise the h-adaptation structure
    !---------------------------------------------------------------------------

    call parlst_getvalue_string(rparlist,&
        ssectionName, 'adaptivity', sadaptivityName, '')
    if (trim(adjustl(sadaptivityName)) .ne. '') then

      call parlst_getvalue_int(rparlist,&
          trim(sadaptivityName), 'nadapt', nadapt)

      if (nadapt > 0) then

        ! Initialise adaptation structure from parameter list
        call hadapt_initFromParameterlist(rhadapt, rparlist, sadaptivityName)

        ! Generate a dynamic graph for the sparsity pattern and attach
        ! it to the collection structure which is used to pass this type
        ! to the callback function for h-adaptation
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'templateMatrix', templateMatrix)
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'systemMatrix', systemMatrix)
        call lsyssc_createGraphFromMatrix(&
            p_rproblemLevel%RmatrixScalar(templateMatrix), rgraph)
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
        call tstep_performTimestep(p_rproblemLevel, rtimestep,&
            rsolver, rsolution, transp_nlsolverCallback,&
            rcollection, 1, rrhs)

        ! Release right-hand side vector
        call lsysbl_releaseVector(rrhs)

      else

        ! Solve the primal problem without right-hand side
        call tstep_performTimestep(p_rproblemLevel, rtimestep,&
            rsolver, rsolution, transp_nlsolverCallback,&
            rcollection, 1)
      end if

      ! Stop time measurement for solution procedure
      call stat_stopTimer(p_rtimerSolution)

      ! Write norm of solution to benchmark logfile
      if (OU_BENCHLOG .ne. 0) then
        dnormL1 = lsysbl_vectorNorm(rsolution, LINALG_NORML1)
        dnormL2 = lsysbl_vectorNorm(rsolution, LINALG_NORML2)
        call output_line('T = '//trim(sys_sdEL(rtimestep%dTime,5))//&
            '   ||u||_1 = '//trim(sys_sdEL(dnormL1,5))//&
            '   ||u||_2 = '//trim(sys_sdEL(dnormL2,5)),&
            OU_CLASS_MSG, OU_MODE_BENCHLOG)
      end if


      if (iadapt .eq. nadapt) exit adaptloop

      !-------------------------------------------------------------------------
      ! Perform recovery-based error estimation
      !-------------------------------------------------------------------------

      ! Start time measurement for error estimation
      call stat_startTimer(p_rtimerErrorEstimation, STAT_TIMERSHORT)

      ! Compute the error estimator using recovery techniques
      call transp_errestRecovery(rparlist, ssectionname,&
          p_rproblemLevel, rsolution, 0.0_DP, relementError, derror,&
          rcollection)

      ! Stop time measurement for error estimation
      call stat_stopTimer(p_rtimerErrorEstimation)


      !-------------------------------------------------------------------------
      ! Perform h-adaptation
      !-------------------------------------------------------------------------

      ! Start time measurement for h-adaptation
      call stat_startTimer(p_rtimerAdaptation, STAT_TIMERSHORT)

      ! Set the name of the template matrix
      rcollection%SquickAccess(1) = 'sparsitypattern'

      ! Attach the primal solution vector to the collection structure
      rcollection%p_rvectorQuickAccess1 => rsolution

      ! Perform h-adaptation and update the triangulation structure
      call transp_adaptTriangulation(&
          rhadapt, p_rproblemLevel%rtriangulation,&
          relementError, rcollection)

      ! Update the template matrix according to the sparsity pattern
      call lsyssc_createMatrixFromGraph(rgraph,&
          p_rproblemLevel%RmatrixScalar(templateMatrix))

      ! Resize the solution vector accordingly
      call lsysbl_resizeVectorBlock(p_rdiscretisation,&
          rsolution, .false.)

      ! Release element-wise error distribution
      call lsyssc_releaseVector(relementError)

      ! Stop time measurement for h-adaptation
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

      ! Re-initialise all constant coefficient matrices
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
    real(dp) :: derror, dnormL1, dnormL2
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
    allocate(rtimestep%p_rthetaScheme)
    rtimestep%p_rthetaScheme%theta = 1.0_DP

    ! Get global parameters
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'irhstype', irhstype)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'ivelocitytype', ivelocitytype)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'discretisation', discretisation)

    ! Set pointer to maximum problem level and block discretisation
    p_rproblemLevel   => rproblem%p_rproblemLevelMax
    p_rdiscretisation => p_rproblemLevel%RblockDiscretisation(discretisation)

    ! Initialise the solution vector and impose boundary conditions explicitly
    call lsysbl_createVectorBlock(p_rdiscretisation, rsolutionPrimal,&
        .false., ST_DOUBLE)
    call transp_initSolution(rparlist, ssectionName, p_rproblemLevel,&
        0.0_DP, rsolutionPrimal, rcollection)
    call bdrf_filterVectorExplicit(rbdrCondPrimal, rsolutionPrimal, 0.0_DP)

    !---------------------------------------------------------------------------
    ! Initialise the h-adaptation structure
    !---------------------------------------------------------------------------

    call parlst_getvalue_string(rparlist,&
        ssectionName, 'adaptivity', sadaptivityName, '')
    if (trim(adjustl(sadaptivityName)) .ne. '') then

      call parlst_getvalue_int(rparlist,&
          trim(sadaptivityName), 'nadapt', nadapt)

      if (nadapt > 0) then

        ! Initialise adaptation structure from parameter list
        call hadapt_initFromParameterlist(rhadapt, rparlist, sadaptivityName)

        ! Generate a dynamic graph for the sparsity pattern and attach
        ! it to the collection structure which is used to pass this type
        ! to the callback function for h-adaptation
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'templateMatrix', templateMatrix)
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'systemMatrix', systemMatrix)
        call lsyssc_createGraphFromMatrix(&
            p_rproblemLevel%RmatrixScalar(templateMatrix), rgraph)
        call collct_setvalue_graph(rcollection, 'sparsitypattern',&
            rgraph, .true., ssectionName=ssectionName)

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
        call tstep_performTimestep(p_rproblemLevel, rtimestep,&
            rsolver, rsolutionPrimal, transp_nlsolverCallback,&
            rcollection, 1, rrhs)

        ! Release right-hand side vector
        call lsysbl_releaseVector(rrhs)

      else

        ! Solve the primal problem without right-hand side
        call tstep_performTimestep(p_rproblemLevel, rtimestep,&
            rsolver, rsolutionPrimal, transp_nlsolverCallback,&
            rcollection, 1)
      end if

      ! Stop time measurement for solution procedure
      call stat_stopTimer(p_rtimerSolution)

      ! Write norm of primal solution to benchmark logfile
      if (OU_BENCHLOG .ne. 0) then
        dnormL1 = lsysbl_vectorNorm(rsolutionPrimal, LINALG_NORML1)
        dnormL2 = lsysbl_vectorNorm(rsolutionPrimal, LINALG_NORML2)
        call output_line('T = '//trim(sys_sdEL(rtimestep%dTime,5))//&
            '   ||u||_1 = '//trim(sys_sdEL(dnormL1,5))//&
            '   ||u||_2 = '//trim(sys_sdEL(dnormL2,5)),&
            OU_CLASS_MSG, OU_MODE_BENCHLOG)
      end if

      !-------------------------------------------------------------------------
      ! Compute the right-hand side for the dual problem
      !-------------------------------------------------------------------------

      ! Start time measurement for error estimation
      call stat_startTimer(p_rtimerErrorEstimation, STAT_TIMERSHORT)

      ! Initialise target functional
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
      call tstep_performTimestep(p_rproblemLevel, rtimestep,&
          rsolver, rsolutionDual, transp_nlsolverCallback,&
          rcollection, 1, rtargetFunc)

      ! Release discretised target functional
      call lsysbl_releaseVector(rtargetFunc)

      ! Stop time measurement for solution procedure
      call stat_stopTimer(p_rtimerSolution)

      ! Write norm of dual solution to benchmark logfile
      if (OU_BENCHLOG .ne. 0) then
        dnormL1 = lsysbl_vectorNorm(rsolutionDual, LINALG_NORML1)
        dnormL2 = lsysbl_vectorNorm(rsolutionDual, LINALG_NORML2)
        call output_line('T = '//trim(sys_sdEL(rtimestep%dTime,5))//&
            '   ||u||_1 = '//trim(sys_sdEL(dnormL1,5))//&
            '   ||u||_2 = '//trim(sys_sdEL(dnormL2,5)),&
            OU_CLASS_MSG, OU_MODE_BENCHLOG)
      end if


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

        ! Initialise right-hand side vector
        call lsysbl_createVectorBlock(rsolutionPrimal, rrhs)
        call transp_initRHS(rparlist, ssectionName, p_rproblemLevel,&
            0.0_DP, rrhs, rcollection)

        ! Compute the error in the quantity of interest
        call transp_errestTargetFunc(rparlist, ssectionName,&
            p_rproblemLevel, rtimestep, rsolver, rsolutionPrimal,&
            rsolutionDual, rcollection, relementError, derror, rrhs)

        ! Release right-hand side vector
        call lsysbl_releaseVector(rrhs)

      else

        ! Compute the error in the quantity of interest
        call transp_errestTargetFunc(rparlist, ssectionName,&
            p_rproblemLevel, rtimestep, rsolver, rsolutionPrimal,&
            rsolutionDual, rcollection, relementError, derror)

      end if

      ! Stop time measurement for error estimation
      call stat_stopTimer(p_rtimerErrorEstimation)


      !-------------------------------------------------------------------------
      ! Perform h-adaptation
      !-------------------------------------------------------------------------

      ! Start time measurement for h-adaptation
      call stat_startTimer(p_rtimerAdaptation, STAT_TIMERSHORT)

      ! Set the name of the template matrix
      rcollection%SquickAccess(1) = 'sparsitypattern'

      ! Attach the primal solution vector to the collection structure
      rcollection%p_rvectorQuickAccess1 => rsolutionPrimal

      ! Perform h-adaptation and update the triangulation structure
      call transp_adaptTriangulation(&
          rhadapt, p_rproblemLevel%rtriangulation,&
          relementError, rcollection)

      ! Update the template matrix according to the sparsity pattern
      call lsyssc_createMatrixFromGraph(rgraph, p_rproblemLevel%RmatrixScalar(templateMatrix))

      ! Resize the solution vector accordingly
      call lsysbl_resizeVectorBlock(p_rdiscretisation,&
          rsolutionPrimal, .false.)

      ! Release element-wise error distribution
      call lsyssc_releaseVector(relementError)

      ! Stop time measurement for h-adaptation
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

      ! Re-initialise all constant coefficient matrices
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

end module transport_application
