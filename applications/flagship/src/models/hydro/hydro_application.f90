!##############################################################################
!# ****************************************************************************
!# <name> hydro_application </name>
!# ****************************************************************************
!#
!# <purpose>
!# This application solves the time-dependent compressible Euler equations
!#
!#   $$\frac{\partial U}{\partial t}+\nabla\cdot{\bf F}(U)=0$$
!#
!# for the vector of conservative variables
!#
!#   $$U=[\rho,\rho{\bf v},\rho E]^T$$
!#
!# and the set of inviscid fluxes ${\bf F}=(F^1,F^2,F^3)$
!# for each coordinate direction
!#
!#   $${\bf F}=[\rho{\bf v},
!#              \rho{\bf v}\otimes{\bf v}+p{\mathcal I},
!#              (\rho E+p){\bf v}]^T$$
!#
!# in the one-, two- or three-dimensional domain $\Omega$.
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
!# 1.) hydro_app
!#     -> The main routine of the application called from the main
!#        program. The routine gets all required information from the
!#        parameter list which needs to be initialised and filled in
!#        the main program. It then works black-box, that is, it
!#        determines the solution algorithm to be used and performs
!#        the simulation. The user should only have to modify this
!#        routine if another solution algorithm is implemented.
!#
!# 2.) hydro_solveTransientPrimal
!#     -> Solves the primal formulation of the time-dependent
!#        compressible Euler equations.
!#
!# </purpose>
!##############################################################################

module hydro_application

#include "../../flagship.h"

!$use omp_lib
  use basicgeometry
  use boundarycondaux
  use boundaryfilter
  use collection
  use flagship_basic
  use fparser
  use genoutput
  use graph
  use hadaptaux
  use hadaptivity
  use linearalgebra
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use problem
  use solveraux
  use spatialdiscretisation
  use statistics
  use storage
  use timestep
  use timestepaux
  use triangulation
  use ucd

  ! Modules from hydrodynamic model
  use hydro_basic
  use hydro_callback
  use hydro_callback1d
  use hydro_callback2d
  use hydro_callback3d
  use hydro_errorestimation
  use hydro_meshadaptation
  use hydro_postprocessing
  use hydro_preprocessing

  implicit none

  private

  public :: hydro_app
  public :: hydro_solveTransientPrimal


contains

  !*****************************************************************************

!<subroutine>

  subroutine hydro_app(rparlist, ssectionName)

!<description>
    ! This is the main application for the compressible Euler
    ! equations.  It is a so-called driver routine which can be used
    ! to start a standalone hydrodynamic simulation.
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

    ! Abstract problem descriptor
    type(t_problemDescriptor) :: rproblemDescriptor

    ! Parameter file and section names
    character(LEN=SYS_STRLEN) :: sindatfileName
    character(LEN=SYS_STRLEN) :: sbdrcondName
    character(LEN=SYS_STRLEN) :: algorithm
    character(LEN=SYS_STRLEN) :: benchmark

    ! local variables
    integer :: isystemFormat, systemMatrix, ndimension, version


    ! Check if configuration file is valid
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'version', version, 0)
    if (version .lt. HYDRO_VERSION) then
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

    ! Create a separate section for the hydrodynamic model
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
    
    ! Initialise the solver structures
    call hydro_initSolvers(rparlist, ssectionName, rtimestep, rsolver)

    ! Initialise the boundary conditions for the primal problem
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'ndimension', ndimension)
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'sprimalbdrcondname', sbdrcondName)
    call bdrc_readBoundaryCondition(rbdrCondPrimal,&
        sindatfileName, sbdrcondName,&
        ndimension, hydro_parseBoundaryCondition)

    ! Initialise the boundary conditions for the dual problem
    call parlst_getvalue_string(rparlist, ssectionName,&
        'sdualbdrcondname', sbdrcondName, '')
    if (sbdrcondName .ne. '') then
      call bdrc_readBoundaryCondition(rbdrCondDual,&
          sindatfileName, sbdrcondName,&
          ndimension, hydro_parseBoundaryCondition)
    end if

    ! Initialise the abstract problem structure
    call hydro_initProblemDescriptor(rparlist, ssectionName,&
        solver_getMinimumMultigridlevel(rsolver),&
        solver_getMaximumMultigridlevel(rsolver),&
        rproblemDescriptor)
    call problem_initProblem(rproblemDescriptor, rproblem)

    ! Initialise all individual problem levels with primal and dual
    ! boundary conditions (if available)
    if (sbdrcondName .ne. '') then
      call hydro_initAllProblemLevels(rparlist,&
          ssectionName, rproblem, rcollection,&
          rbdrCondPrimal, rbdrCondDual)
    else
      call hydro_initAllProblemLevels(rparlist,&
          ssectionName, rproblem, rcollection,&
          rbdrCondPrimal)
    end if

    ! Prepare internal data arrays of the solver structure
    call parlst_getvalue_int(rparlist, ssectionName,&
        'systemMatrix', systemMatrix)
    call parlst_getvalue_int(rparlist, ssectionName,&
        'isystemFormat', isystemFormat)
    call flagship_updateSolverMatrix(rproblem%p_rproblemLevelMax,&
        rsolver, systemMatrix, isystemFormat, UPDMAT_ALL)
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
        ! Solve the primal formulation for the time-dependent problem
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call hydro_solveTransientPrimal(rparlist, ssectionName,&
            rbdrCondPrimal, rproblem, rtimestep, rsolver,&
            rsolutionPrimal, rcollection)
        
        ! Output solution to file
        call hydro_outputSolution(rparlist, ssectionName,&
            rproblem%p_rproblemLevelMax, rsolutionPrimal,&
            dtime=rtimestep%dTime)

      else
        call output_line(trim(algorithm)//' is not a valid solution algorithm!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_app')
        call sys_halt()
      end if
      
    else

      ! Just output the computational mesh and exit
      call hydro_outputSolution(rparlist, ssectionName,&
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
    call hydro_outputStatistics(rtimerTotal, ssectionName, rcollection)

    ! Release collection
    call collct_done(rcollection)
 
  end subroutine hydro_app
  
!*****************************************************************************

!<subroutine>

  subroutine hydro_solveTransientPrimal(rparlist, ssectionName,&
      rbdrCond, rproblem, rtimestep, rsolver, rsolution, rcollection)

!<description>
      ! This subroutine solves the transient primal compressible Euler equations
      !
      !  $$ \frac{\partial U}{\partial t}+\nabla\cdot{\bf F}(U)=b $$
      !
      ! for the vector of conservaive variables $U$ in the domain
      ! $\Omega$. Here, ${\bf F}(u)$ denotes the inviscid fluxes b is
      ! the right-hand side vector.
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
    character(LEN=SYS_STRLEN) :: sucdimport
    
    ! local variables
    real(DP) :: dnormL1, dnormL2
    type(t_ucdExport) :: rimport
    real(dp) :: derror, dstepUCD, dtimeUCD, dstepAdapt, dtimeAdapt
    integer :: templateMatrix, systemMatrix, isystemFormat, discretisation
    integer :: ipreadapt, npreadapt, irhstype, ndimension
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
        ssectionName, 'irhstype', irhstype, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'isystemFormat', isystemFormat)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'ndimension', ndimension)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'discretisation', discretisation)

    ! Set pointer to maximum problem level and discretisation
    p_rproblemLevel   => rproblem%p_rproblemLevelMax
    p_rdiscretisation => p_rproblemLevel%Rdiscretisation(discretisation)

    ! Create the solution vector
    if (p_rdiscretisation%ncomponents .eq. hydro_getNVAR(p_rproblemLevel)) then
      call lsysbl_createVectorBlock(p_rdiscretisation,&
          rsolution, .false., ST_DOUBLE)
    else
      call lsysbl_createVectorBlock(p_rdiscretisation,&
          hydro_getNVAR(p_rproblemLevel),&
          rsolution, .false., ST_DOUBLE)
    end if

    ! Initialise the solution vector and impose boundary conditions
    call hydro_initSolution(rparlist, ssectionName, p_rproblemLevel,&
        rtimestep%dinitialTime, rsolution, rcollection)

    select case(ndimension)
    case (NDIM1D)
      call bdrf_filterVectorExplicit(rbdrCond, rsolution,&
          rtimestep%dinitialTime, hydro_calcBoundaryvalues1d)

    case (NDIM2D)
      call bdrf_filterVectorExplicit(rbdrCond, rsolution,&
          rtimestep%dinitialTime, hydro_calcBoundaryvalues2d)

    case (NDIM3D)
      call bdrf_filterVectorExplicit(rbdrCond, rsolution,&
          rtimestep%dinitialTime, hydro_calcBoundaryvalues3d)
    end select

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
            p_rproblemLevel%Rmatrix(templateMatrix), rgraph)
        call collct_setvalue_graph(rcollection, 'sparsitypattern',&
            rgraph, .true., ssectionName=ssectionName)

        ! Perform pre-adaptation?
        if (npreadapt > 0) then

          ! Perform number of pre-adaptation steps
          do ipreadapt = 1, npreadapt

            ! Compute the error estimator using recovery techniques
            call hydro_errestRecovery(rparlist, ssectionname,&
                p_rproblemLevel, rsolution, rtimestep%dinitialTime,&
                relementError, derror, rcollection)

            ! Set the name of the template matrix
            rcollection%SquickAccess(1) = 'sparsitypattern'

            ! Attach the primal solution vector to the collection structure
            rcollection%p_rvectorQuickAccess1 => rsolution

            ! Perform h-adaptation and update the triangulation structure
            call hydro_adaptTriangulation(rparlist, ssectionname,&
                rhadapt, p_rproblemLevel%rtriangulation,&
                relementError, rcollection)

            ! Release element-wise error distribution
            call lsyssc_releaseVector(relementError)

            ! Generate standard mesh from raw mesh
            call tria_initStandardMeshFromRaw(&
                p_rproblemLevel%rtriangulation, rproblem%rboundary)

            ! Update the template matrix according to the sparsity pattern
            call lsyssc_createMatrixFromGraph(rgraph,&
                p_rproblemLevel%Rmatrix(templateMatrix))

            ! Re-initialise all constant coefficient matrices
            call hydro_initProblemLevel(rparlist, ssectionName,&
                p_rproblemLevel, rcollection, rbdrCond)
            
            ! Resize the solution vector accordingly
            call lsysbl_resizeVectorBlock(&
                p_rproblemLevel%RmatrixBlock(systemMatrix), rsolution, .false.)
            
            ! Re-generate the initial solution vector
            call hydro_initSolution(rparlist, ssectionname,&
                p_rproblemLevel, rtimestep%dinitialTime, rsolution,&
                rcollection)

            ! Impose boundary conditions explicitly
            select case(ndimension)
            case (NDIM1D)
              call bdrf_filterVectorExplicit(rbdrCond, rsolution,&

                  rtimestep%dinitialTime, hydro_calcBoundaryvalues1d)
            case (NDIM2D)
              call bdrf_filterVectorExplicit(rbdrCond, rsolution,&
                  rtimestep%dinitialTime, hydro_calcBoundaryvalues2d)

            case (NDIM3D)
              call bdrf_filterVectorExplicit(rbdrCond, rsolution,&
                  rtimestep%dinitialTime, hydro_calcBoundaryvalues3d)
            end select    
          end do

          ! Prepare internal data arrays of the solver structure
          call flagship_updateSolverMatrix(p_rproblemLevel, rsolver,&
              systemMatrix, isystemFormat, UPDMAT_ALL)
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
      call hydro_setVariables(rimport, rsolution)
      call ucd_release(rimport)
      
      ! Set time for solution output
      dtimeUCD = rtimestep%dinitialTime
    end if
    
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
          call hydro_outputSolution(rparlist, ssectionName,&
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
              rsolver, rsolution, hydro_nlsolverCallback,&
              rcollection, rrhs)
        else
          ! Explicit Runge-Kutta scheme without right-hand side vector
          call tstep_performRKStep(p_rproblemLevel, rtimestep,&
              rsolver, rsolution, hydro_nlsolverCallback,&
              rcollection)
        end if

      case (TSTEP_THETA_SCHEME)

        if (irhstype > 0) then
          ! Two-level theta-scheme with non-zero right-hand side vector
          call tstep_performThetaStep(p_rproblemLevel, rtimestep,&
              rsolver, rsolution, hydro_nlsolverCallback,&
              rcollection, rrhs)
        else
          ! Two-level theta-scheme without right-hand side vector
          call tstep_performThetaStep(p_rproblemLevel, rtimestep,&
              rsolver, rsolution, hydro_nlsolverCallback,&
              rcollection)
        end if

      case default
        call output_line('Unsupported time-stepping algorithm!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_solveTransientPrimal')
        call sys_halt()
      end select

      ! Perform linearised FEM-FCT post-processing
      call hydro_calcLinearisedFCT(p_rproblemLevel, rtimestep,&
          rsolver, rsolution, ssectionName, rcollection,&
          rvector1=rvector1, rvector2=rvector2, rvector3=rvector3)

      ! Stop time measurement for solution procedure
      call stat_stopTimer(p_rtimerSolution)
      
      ! Write norm of solution to benchmark logfile
      if (OU_BENCHLOG .ne. 0) then
        dnormL1 = lsysbl_vectorNorm(rsolution, LINALG_NORML1)
        dnormL2 = lsysbl_vectorNorm(rsolution, LINALG_NORML2)
        call output_line('T = '//trim(sys_sdEL(rtimestep%dTime,5))//&
            '   ||U||_1 = '//trim(sys_sdEL(dnormL1,5))//&
            '   ||U||_2 = '//trim(sys_sdEL(dnormL2,5)),&
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
        call stat_startTimer(p_rtimerPrepostProcess, STAT_TIMERSHORT)

        ! Export the intermediate solution
        call hydro_outputSolution(rparlist, ssectionName,&
            p_rproblemLevel, rsolution, dtime=rtimestep%dTime)

        ! Stop time measurement for post-processing
        call stat_stopTimer(p_rtimerPrepostProcess)

      end if


      !-------------------------------------------------------------------------
      ! Perform adaptation
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
        call hydro_errestRecovery(rparlist, ssectionname,&
            p_rproblemLevel, rsolution, rtimestep%dTime,&
            relementError, derror, rcollection)

        ! Stop time measurement for error estimation
        call stat_stopTimer(p_rtimerErrorEstimation)


        !-----------------------------------------------------------------------
        ! Perform h-adaptation
        !-----------------------------------------------------------------------

        ! Start time measurement for h-adaptation
        call stat_startTimer(p_rtimerAdaptation, STAT_TIMERSHORT)

        ! Set the names of the template matrix
        rcollection%SquickAccess(1) = 'sparsitypattern'

        ! Attach the primal solution vector to the collection structure
        rcollection%p_rvectorQuickAccess1 => rsolution

        ! Perform h-adaptation and update the triangulation structure
        call hydro_adaptTriangulation(rparlist, ssectionname,&
            rhadapt, p_rproblemLevel%rtriangulation,&
            relementError, rcollection)

        ! Release element-wise error distribution
        call lsyssc_releaseVector(relementError)

        ! Update the template matrix according to the sparsity pattern
        call lsyssc_createMatrixFromGraph(rgraph,&
            p_rproblemLevel%Rmatrix(templateMatrix))

        ! Resize the solution vector accordingly
        call lsysbl_resizeVectorBlock(&
            p_rproblemLevel%RmatrixBlock(systemMatrix), rsolution, .false.)
        
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
        call hydro_initProblemLevel(rparlist, ssectionName,&
            p_rproblemLevel, rcollection, rbdrCond)

        ! Prepare internal data arrays of the solver structure
        call flagship_updateSolverMatrix(p_rproblemLevel, rsolver,&
            systemMatrix, isystemFormat, UPDMAT_ALL)
        call solver_updateStructure(rsolver)

        ! Stop time measurement for generation of constant
        ! coefficient matrices
        call stat_stopTimer(p_rtimerAssemblyCoeff)

      end if

    end do timeloop


    ! Release adaptation structure
    if (trim(adjustl(sadaptivityName)) .ne. '') then
      if ((dstepAdapt > 0.0_DP) .or. (npreadapt > 0)) then
        call hadapt_releaseAdaptation(rhadapt)
        call grph_releaseGraph(rgraph)
      end if
    end if

    ! Release vectors for the linearised FCT algorithm
    call lsysbl_releaseVector(rvector1)
    call lsysbl_releaseVector(rvector2)
    call lsysbl_releaseVector(rvector3)

  end subroutine hydro_solveTransientPrimal

end module hydro_application
