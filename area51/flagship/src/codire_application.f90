!##############################################################################
!# ****************************************************************************
!# <name> codire_application </name>
!# ****************************************************************************
!#
!# <purpose>
!# This application solves the generic time-dependent conservation law
!#
!#   $$\frac{\partial u}{\partial t}+\nabla\cdot{\bf f}(u)=s(u)$$
!#
!# for a scalar quantity $u$ in the domain $\Omega$ in 1D, 2D and 3D.
!# The (possibly nonlinear) flux function ${\bf f}(u)$ can be
!#
!#   $${\bf f}(u):={\bf v}u-D\nabla u$$
!#
!# whereby $\bf v$ denotes an externally given velocity profile. 
!# The physical diffusion tensor
!#
!#   $$D=\left[\begin{array}{ccc}
!#       d_{11} & d_{12} & d_{13}\\
!#       d_{21} & d_{22} & d_{23}\\
!#       d_{31} & d_{32} & d_{33}
!#       \end{array}\right]$$
!#
!# may by anitotropic and it can also reduce to the scalar coefficient
!#
!#   $$D=dI,\quad d_{ij}=0\,\forall j\ne i$$
!#
!# Moreover, ${\bf f}(u)$ can also be defined such that
!#
!#   $$\nabl\cdot{\bf f}(u):=u_t+\frac{\partial f(u)}{\partial x}$$
!#
!# that is, the above equation stands for a 1d or 2d conservation
!# law which is solved in the space-time domain $\Omega=$.
!#
!# The spatial discretization is perform by means of the algebraic
!# flux correction (AFC) paradigm by Kuzmin, Moeller and Turek. In
!# particular, high-resolution finite element schemes of TVD- and
!# FCT-type are available. For the temporal discretization, the 
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
!# 1.) codire
!#     -> The application's main routine called from the main problem
!#
!# 2.) codire_initApplication
!#     -> Initializes the application descriptor by the parameter
!#        settings given by the parameter list
!#
!# 3.) codire_doneApplication
!#     -> Finalizes the application descriptor
!# 
!# 4.) codire_initCollection
!#     -> Initializes the collection structure based on the parameter
!#        settings given by the application descriptor
!#
!# 5.) codire_initSolvers
!#     -> Initializes the solve structures from the parameter list
!#
!# 6.) codire_initProblem
!#     -> Initializes the global problem structure based on the
!#        parameter settings given by the parameter list
!#
!# 7.) codire_initProblemLevel
!#     -> Initializes the individual problem level based on the
!#        parameter settings given by the application descriptor
!#
!# 8.) codire_initAllProblemLevels
!#     -> Initializes all problem levels attached to the global
!#        problem structure based on the parameter settings
!#        given by the parameter list
!#
!# 9.) codire_initSolution
!#     -> Initializes the solution vector based on the parameter
!#        settings given by the parameter list
!#
!# 10.) codire_initTargetFunc
!#      -> Initializes the target functional for the dual problem
!#
!# 11.) codire_outputSolution
!#      -> Outputs the solution vector to file in UCD format
!#
!# 12.) codire_outputStatistics
!#      -> Outputs the application statitics
!#
!# 13.) codire_estimateTargetFuncError
!#      -> Estimates the error in the quantity of interest
!#
!# 14.) codire_adaptTriangulation
!#      -> Performs h-adaptation for the given triangulation
!#
!# 15.) codire_solveTransientPrimal
!#      -> Solves the primal formulation of the time-dependent 
!#         convection-diffusion-reaction equation.
!#
!# 16.) codire_solvePseudoTransientPrimal
!#      -> Solves the primal formulation of the steady convection-
!#         diffusion-reaction equation using pseudo time-stepping.
!#
!# 17.) codire_solveSteadyStatePrimal
!#      -> Solves the primal formulation of the steady convection-
!#         diffusion-reaction equation directly
!#
!#
!# The following auxiliary routines are available:
!#
!# 1.) codire_parseCmdlArguments
!#     -> Parses the list of commandline arguments and overwrites
!#        parameter values from the parameter files
!# 
!# </purpose>
!##############################################################################

module codire_application

  use afcstabilisation
  use bilinearformevaluation
  use boundaryfilter
  use codire_basic
  use codire_callback
  use collection
  use flagship_basic
  use fparser
  use fsystem
  use genoutput
  use graph
  use hadaptaux1d
  use hadaptaux2d
  use hadaptaux3d
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
  use solver
  use spatialdiscretisation
  use statistics
  use stdoperators
  use storage
  use timestep
  use ucd

  implicit none
  
  private
  public :: codire

contains

  !*****************************************************************************

!<subroutine>

  subroutine codire(rparlist)

!<description>
    ! This is the main application for the convection-diffusion-reaction 
    ! benchmark problem. It is a so-called driver routine which can be 
    ! used to start a standalone convection-diffusion-reaction application.
!</description>

!<inputoutput>
    ! parameter list
    type(t_parlist), intent(INOUT) :: rparlist
!</inputoutput>
!</subroutine>

    !****************************************************************************
    ! Structures required for this application
    
    ! Global collection which is used to pass arguments to callback routines
    type(t_collection) :: rcollection

    ! Application descriptor which holds all internal information
    ! for the convection-diffusion-reaction application
    type(t_codire) :: rappDescriptor

    ! Boundary condition structure for the primal problem
    type(t_boundaryCondition) :: rbdrCondPrimal

    ! Boundary condition structure for the dual problem
    type(t_boundaryCondition) :: rbdrCondDual

    ! Problem structure which holds all internal data (vectors/matrices)
    ! for the convection-diffusion-reaction application
    type(t_problem) :: rproblem
    
    ! Time-stepping structure for the convection-
    ! diffusion-reaction application
    type(t_timestep) :: rtimestep
    
    ! Global solver structure for the convection-
    ! diffusion-reaction application
    type(t_solver) :: rsolver

    ! Solution vector for the primal problem
    type(t_vectorBlock) :: rsolutionPrimal

    ! Solution vector for the dual problem
    type(t_vectorBlock) :: rsolutionDual

    ! Timer for the total solution process
    type(t_timer) :: rtimerTotal

    ! Parameter file and section names
    character(LEN=SYS_STRLEN) :: sindatfileName
    character(LEN=SYS_STRLEN) :: sbdrcondName
    character(LEN=SYS_STRLEN) :: algorithm
    

    ! Start total time measurement
    call stat_clearTimer(rtimerTotal)
    call stat_startTimer(rtimerTotal)
    
    !---------------------------------------------------------------------------
    ! Pre-processing
    !---------------------------------------------------------------------------
    
    ! Overwrite global configuration from command line arguments. After
    ! this subroutine has been called, the parameter list remains unchanged.
    call codire_parseCmdlArguments(rparlist)

    ! Initialize global collection structure
    call collct_init(rcollection)

    ! Initialize the application descriptor
    call codire_initApplication(rparlist, 'codire', rappDescriptor)

    ! Start time measurement for pre-processing
    call stat_startTimer(rappDescriptor%rtimerPrepostProcess, STAT_TIMERSHORT)
    
    ! Initialize the global collection
    call codire_initCollection(rappDescriptor, rparlist, rcollection)

    ! Initialize the solver structures
    call codire_initSolvers(rparlist, 'codire', rtimestep, rsolver)

    ! Initialize the abstract problem structure
    call codire_initProblem(rparlist, 'codire',&
                            solver_getMinimumMultigridlevel(rsolver),&
                            solver_getMaximumMultigridlevel(rsolver),&
                            rproblem, rcollection)

    ! Initialize the individual problem levels
    call codire_initAllProblemLevels(rappDescriptor, rproblem, rcollection)

    ! Initialize the primal solution vector
    call codire_initSolution(rparlist, 'codire', rproblem%p_rproblemLevelMax,&
                             0.0_DP, rsolutionPrimal)
        

    ! Prepare internal data arrays of the solver structure
    call flagship_updateSolverMatrix(rproblem%p_rproblemLevelMax, rsolver,&
                                     1, SYSTEM_INTERLEAVEFORMAT, UPDMAT_ALL)
    call solver_updateStructure(rsolver)

    ! Stop time measurement for pre-processing
    call stat_stopTimer(rappDescriptor%rtimerPrePostprocess)
    

    !---------------------------------------------------------------------------
    ! Solution algorithm
    !---------------------------------------------------------------------------
    
    if (rtimestep%dfinalTime > 0) then
      
      ! Get global configuration from parameter list
      call parlst_getvalue_string(rparlist, 'codire', 'algorithm', algorithm)
      call parlst_getvalue_string(rparlist, 'codire', 'indatfile', sindatfileName)
      
      ! The boundary condition for the primal problem is required for all 
      ! solution strategies so initialize it from the parameter file
      call parlst_getvalue_string(rparlist, 'codire', 'sprimalbdrcondname', sbdrcondName)
      call bdrf_readBoundaryCondition(rbdrCondPrimal, sindatfileName,&
                                      '['//trim(sbdrcondName)//']', rappDescriptor%ndimension)

      ! Impose primal boundary conditions explicitely
      call bdrf_filterVectorExplicit(rbdrCondPrimal,&
          rproblem%p_rproblemLevelMax%rtriangulation, rsolutionPrimal, 0.0_DP)

      
      ! What solution algorithm should be applied?
      select case(trim(algorithm))

      case ('transient_primal')
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Solve the primal formulation for the time-dependent problem
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call codire_solveTransientPrimal(rappDescriptor, rparlist, 'codire',&
                                         rbdrCondPrimal, rproblem, rtimestep,&
                                         rsolver, rsolutionPrimal, rcollection)
        call codire_outputSolution(rparlist, 'codire', rproblem%p_rproblemLevelMax,&
                                   rsolutionPrimal, dtime=rtimestep%dTime)
        
      case ('transient_primaldual')
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Solve the primal and dual formulation for the time-dependent problem
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


      case ('pseudotransient_primal')
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Solve the primal formulation for the pseudo time-dependent problem
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call codire_solvePseudoTransientPrimal(rappDescriptor, rparlist, 'codire',&
                                               rbdrCondPrimal, rproblem, rtimestep,&
                                               rsolver, rsolutionPrimal, rcollection)
        call codire_outputSolution(rparlist, 'codire', rproblem%p_rproblemLevelMax,&
                                   rsolutionPrimal, dtime=rtimestep%dTime)


      case ('pseudotransient_primaldual')
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Solve the primal and dual formulation for the pseudo time-dependent problem
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        
      case ('stationary_primal')
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Solve the primal formulation for the stationary problem
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call codire_solveSteadyStatePrimal(rappDescriptor, rparlist, 'codire',&
                                           rbdrCondPrimal, rproblem, rtimestep,&
                                           rsolver, rsolutionPrimal, rcollection)
        call codire_outputSolution(rparlist, 'codire', rproblem%p_rproblemLevelMax,&
                                   rsolutionPrimal, dtime=rtimestep%dTime)


      case ('stationary_primaldual')
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Solve the primal and dual formulation for the stationary problem
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call parlst_getvalue_string(rparlist, 'codire', 'sdualbdrcondname', sbdrcondName)
        call bdrf_readBoundaryCondition(rbdrCondDual, sindatfileName,&
                                      '['//trim(sbdrcondName)//']', rappDescriptor%ndimension)

        call codire_solveSteadyStatePrimalDual(rappDescriptor, rparlist, 'codire',&
                                               rbdrCondPrimal, rbdrCondDual,&
                                               rproblem, rtimestep, rsolver,&
                                               rsolutionPrimal, rsolutionDual, rcollection)
        call codire_outputSolution(rparlist, 'codire', rproblem%p_rproblemLevelMax,&
                                   rsolutionPrimal, rsolutionDual, rtimestep%dTime)


      case DEFAULT
        call output_line(trim(algorithm)//' is not a valid solution algorithm!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'codire_start')
        call sys_halt()
      end select
    end if

    
    !---------------------------------------------------------------------------
    ! Post-processing
    !---------------------------------------------------------------------------

    ! Start time measurement for pre-processing
    call stat_startTimer(rappDescriptor%rtimerPrepostProcess, STAT_TIMERSHORT)
    
    ! Release solvers
    call solver_releaseTimestep(rtimestep)
    call solver_releaseSolver(rsolver)
    
    ! Release problem structure
    call problem_releaseProblem(rproblem)

    ! Release application descriptor
    call codire_doneApplication(rappDescriptor)
    
    ! Release boundary conditions
    call bdrf_release(rbdrCondPrimal)
    call bdrf_release(rbdrCondDual)

    ! Release vectors
    call lsysbl_releaseVector(rsolutionPrimal)
    call lsysbl_releaseVector(rsolutionDual)

    ! Release collection
    call collct_done(rcollection)

    ! Release function parser
    call fparser_done()

    ! Stop time measurement for pre-processing
    call stat_stopTimer(rappDescriptor%rtimerPrePostprocess)

    ! Stop time measurement for total time measurement
    call stat_stopTimer(rtimerTotal)

    ! Output statistics
    call codire_outputStatistics(rappDescriptor, rtimerTotal)
    
  end subroutine codire
  
  !*****************************************************************************

!<subroutine>

  subroutine codire_initApplication(rparlist, ssectionName, rappDescriptor)

!<description>
    ! This subroutine initializes the application descriptor
    ! by the parameter settings given by the parameter list
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(IN) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(IN) :: ssectionName
!</input>

!<output>
    ! application descriptor
    type(t_codire), intent(OUT) :: rappDescriptor
!</output>
!</subroutine>

    ! section names
    character(LEN=SYS_STRLEN) :: sindatfileName
    character(LEN=SYS_STRLEN) :: svelocityName
    character(LEN=SYS_STRLEN) :: sdiffusionName
    character(LEN=SYS_STRLEN) :: sreactionName
    character(LEN=SYS_STRLEN) :: srhsName
    character(LEN=SYS_STRLEN) :: stargetfuncName

    ! symbolic variable names
    character(LEN=*), dimension(4), parameter ::&
                      cvariables = (/ (/'x'/), (/'y'/), (/'z'/), (/'t'/) /)


    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist, ssectionName, 'indatfile', sindatfileName)

    ! Initialize function parser
    call fparser_init()
    call fparser_parseFileForKeyword(sindatfileName, 'defconst', FPAR_CONSTANT)
    call fparser_parseFileForKeyword(sindatfileName, 'defexpr',  FPAR_EXPRESSION)

    ! Get application specifig parameters from the parameterlist
    call parlst_getvalue_int(rparlist, ssectionName, 'ndimension', rappDescriptor%ndimension)
    call parlst_getvalue_int(rparlist, ssectionName, 'imasstype', rappDescriptor%imasstype)
    call parlst_getvalue_int(rparlist, ssectionName, 'imassantidiffusion', rappDescriptor%imassantidiffusion)
    call parlst_getvalue_int(rparlist, ssectionName, 'ivelocitytype', rappDescriptor%ivelocitytype)
    call parlst_getvalue_int(rparlist, ssectionName, 'idiffusiontype', rappDescriptor%idiffusiontype)
    call parlst_getvalue_int(rparlist, ssectionName, 'ireactiontype', rappDescriptor%ireactiontype)
    call parlst_getvalue_int(rparlist, ssectionName, 'irhstype', rappDescriptor%irhstype)
    call parlst_getvalue_int(rparlist, ssectionName, 'itargetfunctype', rappDescriptor%itargetfunctype)
    call parlst_getvalue_int(rparlist, ssectionName, 'ieltype', rappDescriptor%ieltype)
    call parlst_getvalue_int(rparlist, ssectionName, 'imatrixformat', rappDescriptor%imatrixformat)
    call parlst_getvalue_int(rparlist, ssectionName, 'ijacobianformat', rappDescriptor%ijacobianformat)

    
    ! Initialize the function parser for the velocity field if required
    if (rappDescriptor%ivelocitytype .ne. VELOCITY_ZERO) then
      call parlst_getvalue_string(rparlist, ssectionName, 'svelocityname', svelocityName)
      call flagship_readParserFromFile(sindatfileName, '['//trim(svelocityName)//']',&
                                       cvariables, rappDescriptor%rfparserVelocityField)
    end if

    ! Initialize the function parser for the diffusion tensor if required
    if (rappDescriptor%idiffusiontype .ne. DIFFUSION_ZERO) then
      call parlst_getvalue_string(rparlist, ssectionName, 'sdiffusionname', sdiffusionName)
      call flagship_readParserFromFile(sindatfileName, '['//trim(sdiffusionName)//']',&
                                       cvariables, rappDescriptor%rfparserDiffusionTensor)
    end if

    ! Initialize the function parser for the reactive term if required
    if (rappDescriptor%ireactiontype .ne. REACTION_ZERO) then
      call parlst_getvalue_string(rparlist, ssectionName, 'sreactionname', sreactionName)
      call flagship_readParserFromFile(sindatfileName, '['//trim(sreactionName)//']',&
                                       cvariables, rappDescriptor%rfparserReaction)
    end if

    ! Initialize the function parser for the right-hand side if required
    if (rappDescriptor%irhstype .ne. RHS_ZERO) then
      call parlst_getvalue_string(rparlist, ssectionName, 'srhsname', srhsName)
      call flagship_readParserFromFile(sindatfileName, '['//trim(srhsName)//']',&
                                       cvariables, rappDescriptor%rfparserRHS)
    end if

    ! Initialize the function parser for the target functional if required
    if (rappDescriptor%itargetfunctype .ne. TFUNC_ZERO) then
      call parlst_getvalue_string(rparlist, ssectionName, 'stargetfuncname', stargetfuncName)
      call flagship_readParserFromFile(sindatfileName, '['//trim(stargetfuncName)//']',&
                                       cvariables, rappDescriptor%rfparserTargetFunc)
    end if

  end subroutine codire_initApplication

  !*****************************************************************************

!<subroutine>

  subroutine codire_doneApplication(rappDescriptor)

!<description>
    ! This subroutine finalizes the application descriptor
!</description>

!<inputoutput>
    ! application descriptor
    type(t_codire), intent(INOUT) :: rappDescriptor
!</inputoutput>
!</subroutine>

    if (rappDescriptor%ivelocitytype .ne. VELOCITY_ZERO) then
      call fparser_release(rappDescriptor%rfparserVelocityField)
    end if

    if (rappDescriptor%idiffusiontype .ne. DIFFUSION_ZERO) then
      call fparser_release(rappDescriptor%rfparserDiffusionTensor)
    end if

    if (rappDescriptor%ireactiontype .ne. REACTION_ZERO) then
      call fparser_release(rappDescriptor%rfparserReaction)
    end if

    if (rappDescriptor%ireactiontype .ne. RHS_ZERO) then
      call fparser_release(rappDescriptor%rfparserRHS)
    end if

    if (rappDescriptor%itargetfunctype .ne. TFUNC_ZERO) then
      call fparser_release(rappDescriptor%rfparserTargetFunc)
    end if

  end subroutine codire_doneApplication

  !*****************************************************************************

!<subroutine>

  subroutine codire_initCollection(rappDescriptor, rparlist, rcollection)

!<description>
    ! This subroutine initializes the collection based on
    ! the parameter settings of the application descriptor
!</description>

!<input>
    ! application descriptor
    type(t_codire), intent(IN) :: rappDescriptor

    ! parameter list
    type(t_parlist), intent(IN) :: rparlist
!</input>

!<inputoutput>
    ! collection structure
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: i


    ! Add parameter settings from application descriptor
    call collct_setvalue_int(rcollection, 'ndimension',&
                             rappDescriptor%ndimension, .true.)
    call collct_setvalue_int(rcollection, 'imasstype',&
                             rappDescriptor%imasstype, .true.)
    call collct_setvalue_int(rcollection, 'imassantidiffusion',&
                             rappDescriptor%imassantidiffusion, .true.)
    call collct_setvalue_int(rcollection, 'ivelocitytype',&
                             rappDescriptor%ivelocitytype, .true.)
    call collct_setvalue_int(rcollection, 'idiffusiontype',&
                             rappDescriptor%idiffusiontype, .true.)
    call collct_setvalue_int(rcollection, 'ireactiontype',&
                             rappDescriptor%ireactiontype, .true.)
    call collct_setvalue_int(rcollection, 'irhstype',&
                             rappDescriptor%irhstype, .true.)
    call collct_setvalue_int(rcollection, 'ieltype',&
                             rappDescriptor%ieltype, .true.)
    call collct_setvalue_int(rcollection, 'imatrixformat',&
                             rappDescriptor%imatrixFormat, .true.)
    call collct_setvalue_int(rcollection, 'ijacobianformat',&
                             rappDescriptor%ijacobianFormat, .true.)
    

    ! Add timer structures
    call collct_setvalue_timer(rcollection, 'timerSolution',&
                               rappDescriptor%rtimerSolution, .true.)
    call collct_setvalue_timer(rcollection, 'timerAdaptation',&
                               rappDescriptor%rtimerAdaptation, .true.)
    call collct_setvalue_timer(rcollection, 'timerErrorEstimation',&
                               rappDescriptor%rtimerErrorEstimation, .true.)
    call collct_setvalue_timer(rcollection, 'timerTriangulation',&
                               rappDescriptor%rtimerTriangulation, .true.)
    call collct_setvalue_timer(rcollection, 'timerAssemblyCoeff',&
                               rappDescriptor%rtimerAssemblyCoeff, .true.)
    call collct_setvalue_timer(rcollection, 'timerAssemblyMatrix',&
                               rappDescriptor%rtimerAssemblyMatrix, .true.)
    call collct_setvalue_timer(rcollection, 'timerAssemblyVector',&
                               rappDescriptor%rtimerAssemblyVector, .true.)
    call collct_setvalue_timer(rcollection, 'timerPrePostprocess',&
                               rappDescriptor%rtimerPrePostprocess, .true.)

    ! Add internal information about the position of 
    ! constant coefficient matrices and auxiliary vectors
    call parlst_getvalue_int(rparlist, '', 'VelocityField', i, 1)
    call collct_setvalue_int(rcollection,  'VelocityField', i, .true.)
    
    call parlst_getvalue_int(rparlist, '', 'ConvectionAFC', i, 1)
    call collct_setvalue_int(rcollection,  'ConvectionAFC', i, .true.)

    call parlst_getvalue_int(rparlist, '', 'DiffusionAFC', i, 2)
    call collct_setvalue_int(rcollection,  'DiffusionAFC', i, .true.)

    call parlst_getvalue_int(rparlist, '', 'TemplateMatrix', i, 1)
    call collct_setvalue_int(rcollection,  'TemplateMatrix', i, .true.)

    call parlst_getvalue_int(rparlist, '', 'SystemMatrix', i, 2)
    call collct_setvalue_int(rcollection,  'SystemMatrix', i, .true.)
    
    call parlst_getvalue_int(rparlist, '', 'JacobianMatrix', i, 3)
    call collct_setvalue_int(rcollection,  'JacobianMatrix', i, .true.)

    call parlst_getvalue_int(rparlist, '', 'TransportMatrix', i, 4)
    call collct_setvalue_int(rcollection,  'TransportMatrix', i, .true.)

    if ((rappDescriptor%imasstype .ne. MASS_ZERO) .or.&
        (rappDescriptor%imassantidiffusion .ne. MASS_ZERO)) then
      
      call parlst_getvalue_int(rparlist, '', 'ConsistentMassMatrix', i, 5)
      call collct_setvalue_int(rcollection,  'ConsistentMassMatrix', i, .true.)

      call parlst_getvalue_int(rparlist, '', 'LumpedMassMatrix', i, 6)
      call collct_setvalue_int(rcollection,  'LumpedMassMatrix', i, .true.)
      
    else
      
      call collct_setvalue_int(rcollection, 'ConsistentMassMatrix', 0, .true.)
      call collct_setvalue_int(rcollection, 'LumpedMassMatrix',     0, .true.)

    end if

    if (rappDescriptor%idiffusiontype .ne. 0)then
      call parlst_getvalue_int(rparlist, '', 'CoeffMatrix_S', i, 7)
      call collct_setvalue_int(rcollection,  'CoeffMatrix_S', i, .true.)
    else
      call collct_setvalue_int(rcollection,  'CoeffMatrix_S', 0, .true.)
    end if

    if (rappDescriptor%ivelocitytype .ne. 0)then
      select case(rappDescriptor%ndimension)
      case (NDIM1D)
        call parlst_getvalue_int(rparlist, '', 'CoeffMatrix_CX', i, 8)
        call collct_setvalue_int(rcollection,  'CoeffMatrix_CX', i, .true.)
        call collct_setvalue_int(rcollection,  'CoeffMatrix_CY', 0, .true.)
        call collct_setvalue_int(rcollection,  'CoeffMatrix_CZ', 0, .true.)
      case (NDIM2D)
        call parlst_getvalue_int(rparlist, '', 'CoeffMatrix_CX', i, 8)
        call collct_setvalue_int(rcollection,  'CoeffMatrix_CX', i, .true.)
        call parlst_getvalue_int(rparlist, '', 'CoeffMatrix_CY', i, 9)
        call collct_setvalue_int(rcollection,  'CoeffMatrix_CY', i, .true.)
        call collct_setvalue_int(rcollection,  'CoeffMatrix_CZ', 0, .true.)
      case (NDIM3D)
        call parlst_getvalue_int(rparlist, '', 'CoeffMatrix_CX', i, 8)
        call collct_setvalue_int(rcollection,  'CoeffMatrix_CX', i, .true.)
        call parlst_getvalue_int(rparlist, '', 'CoeffMatrix_CY', i, 9)
        call collct_setvalue_int(rcollection,  'CoeffMatrix_CY', i, .true.)
        call parlst_getvalue_int(rparlist, '', 'CoeffMatrix_CZ', i, 10)
        call collct_setvalue_int(rcollection,  'CoeffMatrix_CZ', i, .true.)
      case DEFAULT
        call collct_setvalue_int(rcollection, 'CoeffMatrix_CX',  0, .true.)
        call collct_setvalue_int(rcollection, 'CoeffMatrix_CY',  0, .true.)
        call collct_setvalue_int(rcollection, 'CoeffMatrix_CZ',  0, .true.)
      end select
    else
      call collct_setvalue_int(rcollection, 'CoeffMatrix_CX',  0, .true.)
      call collct_setvalue_int(rcollection, 'CoeffMatrix_CY',  0, .true.)
      call collct_setvalue_int(rcollection, 'CoeffMatrix_CZ',  0, .true.)
    end if

  end subroutine codire_initCollection

  !*****************************************************************************

!<subroutine>

  subroutine codire_initSolvers(rparlist, ssectionName, rtimestep, rsolver)

!<description>
    ! This subroutine initializes the time-stepping structure and
    ! the top-level solver structure from the parameter list
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(IN) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(IN) :: ssectionName
!</input>

!<output>
    ! time-stepping structure
    type(t_timestep), intent(OUT) :: rtimestep

    ! solver struchture
    type(t_solver), intent(OUT) :: rsolver
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
    call solver_createTimestep(rparlist, stimestepName, rtimestep)
    
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
    
  end subroutine codire_initSolvers

  !*****************************************************************************

!<subroutine>

  subroutine codire_initProblem(rparlist, ssectionName, nlmin, nlmax,&
                                rproblem, rcollection)

!<description>
    ! This subroutine initializes the abstract problem structure 
    ! based on the parameters settings given by the parameter list
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(IN) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(IN) :: ssectionName

    ! minimum/maximum problem level
    integer, intent(IN) :: nlmin, nlmax
!</input>

!<inputoutput>
    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</intputoutput>

!<output>
    ! problem structure
    type(t_problem), intent(OUT) :: rproblem
!</output>
!</subroutine>

    ! section names
    character(LEN=SYS_STRLEN) :: sconvectionName
    character(LEN=SYS_STRLEN) :: sdiffusionName

    ! abstract problem descriptor
    type(t_problemDescriptor) :: rproblemDescriptor

    ! pointer to the problem level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! local variables
    integer :: convectionAFC
    integer :: diffusionAFC
    integer :: iconvToTria
    

    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist, ssectionName, 'diffusion', sdiffusionName)
    call parlst_getvalue_string(rparlist, ssectionName, 'convection', sconvectionName)
    

    call parlst_getvalue_int(rparlist, ssectionName, 'ndimension', rproblemDescriptor%ndimension)
    call parlst_getvalue_string(rparlist, ssectionName, 'trifile', rproblemDescriptor%trifile)
    call parlst_getvalue_string(rparlist, ssectionName, 'prmfile', rproblemDescriptor%prmfile)
    
    ! Set additional problem descriptor
    rproblemDescriptor%nafcstab      = 2   ! for convective and diffusive stabilization
    rproblemDescriptor%nlmin         = nlmin
    rproblemDescriptor%nlmax         = nlmax
    rproblemDescriptor%nmatrixScalar = rproblemDescriptor%ndimension + 7
    rproblemDescriptor%nmatrixBlock  = 0
    rproblemDescriptor%nvectorScalar = 0
    rproblemDescriptor%nvectorBlock  = 1   ! for velocity field

    ! Check if quadrilaterals should be converted to triangles
    call parlst_getvalue_int(rparlist, ssectionName, 'iconvtotria', iconvToTria, 0)
    if (iconvToTria .eq. 1)&
        rproblemDescriptor%iproblemSpec = rproblemDescriptor%iproblemSpec &
                                        + PROBDESC_MSPEC_CONVTRIANGLES   

    ! Initialize problem structure
    call problem_initProblem(rproblemDescriptor, rproblem)

    
    ! Initialize the stabilisation structure
    convectionAFC = collct_getvalue_int(rcollection, 'convectionAFC')
    diffusionAFC  = collct_getvalue_int(rcollection, 'diffusionAFC')
    
    ! loop over all problem levels
    p_rproblemLevel => rproblem%p_rproblemLevelMax
    do while(associated(p_rproblemLevel))

      if (convectionAFC > 0) then
        call afcstab_initFromParameterlist(rparlist, sconvectionName,&
                                           p_rproblemLevel%Rafcstab(convectionAFC))
      end if

      if (diffusionAFC > 0) then
        call afcstab_initFromParameterlist(rparlist, sdiffusionName,&
                                           p_rproblemLevel%Rafcstab(diffusionAFC))
      end if

      ! Switch to next coarser level
      p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
    end do

  end subroutine codire_initProblem

  !*****************************************************************************

!<subroutine>

  subroutine codire_initProblemLevel(rappDescriptor, rproblemLevel, rcollection)

!<description>
    ! This subroutine initielizes the individual problem level. It
    ! generates the discretization, the template matrix and the
    ! coefficient matrices as duplicates of the template matrix.
!</description>

!<input>
    ! application descriptor
    type(t_codire), intent(IN) :: rappDescriptor
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! collection structure
    type(t_collection), intent(INOUT) :: rcollection
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
    integer :: ijacobianFormat
    

    ! retrieve application specific parameters from the collection
    templateMatrix       = collct_getvalue_int(rcollection, 'templatematrix')
    systemMatrix         = collct_getvalue_int(rcollection, 'systemmatrix')
    jacobianMatrix       = collct_getvalue_int(rcollection, 'jacobianmatrix')
    transportMatrix      = collct_getvalue_int(rcollection, 'transportmatrix')
    consistentMassMatrix = collct_getvalue_int(rcollection, 'consistentmassmatrix')
    lumpedMassMatrix     = collct_getvalue_int(rcollection, 'lumpedmassmatrix')
    coeffMatrix_S        = collct_getvalue_int(rcollection, 'coeffmatrix_s')
    coeffMatrix_CX       = collct_getvalue_int(rcollection, 'coeffmatrix_cx')
    coeffMatrix_CY       = collct_getvalue_int(rcollection, 'coeffmatrix_cy')
    coeffMatrix_CZ       = collct_getvalue_int(rcollection, 'coeffmatrix_cz')
    convectionAFC        = collct_getvalue_int(rcollection, 'convectionAFC')
    diffusionAFC         = collct_getvalue_int(rcollection, 'diffusionAFC')

    
    ! Initialize the discretization structure
    if (rproblemLevel%rdiscretisation%ndimension .eq. 0) then
      call spdiscr_initBlockDiscr(rproblemLevel%rdiscretisation, 1,&
                                  rproblemLevel%rtriangulation)
    end if
    
    ! Get spatial dimension
    select case(rproblemLevel%rdiscretisation%ndimension)
    case (NDIM1D)
      select case(rappDescriptor%ieltype)
      case (-1,1,11)
        ! P1=Q1 finite elements
        call spdiscr_initDiscr_simple(rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                      EL_E001_1D, SPDISC_CUB_AUTOMATIC,&
                                      rproblemLevel%rtriangulation,&
                                      rproblemLevel%p_rproblem%rboundary)
        ! P2=Q2 finite elements
        call spdiscr_initDiscr_simple(rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                      EL_E002_1D, SPDISC_CUB_AUTOMATIC,&
                                      rproblemLevel%rtriangulation,&
                                      rproblemLevel%p_rproblem%rboundary)
      case DEFAULT
        call output_line('Unsupproted element type!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'codire_initProblemLevel')
        call sys_halt()
      end select
      
    case (NDIM2D)
      select case(rappDescriptor%ieltype)
      case (1)
        ! P1 finite elements
        call spdiscr_initDiscr_simple(rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                      EL_E001, SPDISC_CUB_AUTOMATIC,&
                                      rproblemLevel%rtriangulation,&
                                      rproblemLevel%p_rproblem%rboundary)
      case (2)
        ! P2 finite elements
        call spdiscr_initDiscr_simple(rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                      EL_E002, SPDISC_CUB_AUTOMATIC,&
                                      rproblemLevel%rtriangulation,&
                                      rproblemLevel%p_rproblem%rboundary)
      case (11)
        ! Q1 finite elements
        call spdiscr_initDiscr_simple(rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                      EL_E011, SPDISC_CUB_AUTOMATIC,&
                                      rproblemLevel%rtriangulation,&
                                      rproblemLevel%p_rproblem%rboundary)
      case (12)
        ! Q2 finite elements
        call spdiscr_initDiscr_simple(rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                      EL_E013, SPDISC_CUB_AUTOMATIC,&
                                      rproblemLevel%rtriangulation,&
                                      rproblemLevel%p_rproblem%rboundary)
      case (-1)
        ! mixed P1/Q1 finite elements
        call spdiscr_initDiscr_triquad(rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                       EL_E001, EL_E011, SPDISC_CUB_AUTOMATIC,&
                                       SPDISC_CUB_AUTOMATIC,&
                                       rproblemLevel%rtriangulation,&
                                       rproblemLevel%p_rproblem%rboundary)
      case (-2)
        ! mixed P2/Q2 finite elements
        call spdiscr_initDiscr_triquad(rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                       EL_E002, EL_E013, SPDISC_CUB_AUTOMATIC,&
                                       SPDISC_CUB_AUTOMATIC,&
                                       rproblemLevel%rtriangulation,&
                                       rproblemLevel%p_rproblem%rboundary)
      case DEFAULT
        call output_line('Unsupproted element type!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'codire_initProblemLevel')
        call sys_halt()
      end select
      
    case (NDIM3D)
      select case(rappDescriptor%ieltype)
      case (1)
        ! P1 finite elements
        call spdiscr_initDiscr_simple(rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                      EL_E001_3D, SPDISC_CUB_AUTOMATIC,&
                                      rproblemLevel%rtriangulation,&
                                      rproblemLevel%p_rproblem%rboundary)
      case (11)
        ! Q1 finite elements
        call spdiscr_initDiscr_simple(rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                      EL_E010_3D, SPDISC_CUB_AUTOMATIC,&
                                      rproblemLevel%rtriangulation,&
                                      rproblemLevel%p_rproblem%rboundary)
      case DEFAULT
        call output_line('Unsupproted element type!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'codire_initProblemLevel')
        call sys_halt()
      end select
        
    case DEFAULT
      call output_line('Invalid number of spatial dimensions',&
                       OU_CLASS_ERROR,OU_MODE_STD,'codire_initProblemLevel')
      call sys_halt()
    end select


    ! If the template matrix has no structure data then generate the
    ! finite element matrix sparsity structure based on the spatial
    ! descretization and store it as the template matrix. Otherwise we
    ! assume that the template matrix has been generated externally.
    if (.not.lsyssc_hasMatrixStructure(rproblemLevel%Rmatrix(templateMatrix))) then
      call bilf_createMatrixStructure(rproblemLevel%rdiscretisation%RspatialDiscr(1),&
                                      rappDescriptor%imatrixFormat,&
                                      rproblemLevel%Rmatrix(templateMatrix))
    end if

    ! Create system matrix as duplicate of the template matrix
    if (systemMatrix > 0) then
      if (lsyssc_isMatrixStructureShared(rproblemLevel%Rmatrix(systemMatrix),&
                                         rproblemLevel%Rmatrix(templateMatrix))) then
        call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(systemMatrix),&
                                 rproblemLevel%Rmatrix(templateMatrix),&
                                 .false., .false., .true.)
      else
        call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(templateMatrix),&
                                    rproblemLevel%Rmatrix(systemMatrix),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      end if
    end if

    ! Create transport matrix as duplicate of the template matrix
    if (transportMatrix > 0) then
      if (lsyssc_isMatrixStructureShared(rproblemLevel%Rmatrix(transportMatrix),&
                                         rproblemLevel%Rmatrix(templateMatrix))) then
        call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(transportMatrix),&
                                 rproblemLevel%Rmatrix(templateMatrix),&
                                 .false., .false., .true.)
      else
        call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(templateMatrix),&
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
      ijacobianFormat = collct_getvalue_int(rcollection, 'ijacobianFormat')

      if (lsyssc_hasMatrixStructure(rproblemLevel%Rmatrix(jacobianMatrix))) then
        if (ijacobianFormat .eq. 0) then
          call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(jacobianMatrix),&
                                   rproblemLevel%Rmatrix(templateMatrix),&
                                   .false., .false., .true.)
        else
          call afcstab_generateExtSparsity(rproblemLevel%Rmatrix(templateMatrix),&
                                           rproblemLevel%Rmatrix(jacobianMatrix))
        end if
      else
        if (ijacobianFormat .eq. 0) then
          call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(templateMatrix),&
                                      rproblemLevel%Rmatrix(jacobianMatrix),&
                                      LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
        else
          call afcstab_generateExtSparsity(rproblemLevel%Rmatrix(templateMatrix),&
                                           rproblemLevel%Rmatrix(jacobianMatrix))
        end if
      end if     
    end if
    
    ! Create consistent (and lumped) mass matrix as duplicate of the template matrix
    if (consistentMassMatrix > 0) then
      if (lsyssc_isMatrixStructureShared(rproblemLevel%Rmatrix(consistentMassMatrix),&
                                         rproblemLevel%Rmatrix(templateMatrix))) then
        call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(consistentMassMatrix),&
                                 rproblemLevel%Rmatrix(templateMatrix),&
                                 .false., .false., .true.)
      else
        call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(templateMatrix),&
                                    rproblemLevel%Rmatrix(consistentMassMatrix),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      end if
      call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(consistentMassMatrix),&
                                      DER_FUNC, DER_FUNC) 
      if (lumpedMassMatrix > 0) then
        call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(consistentMassMatrix),&
                                    rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)
        call lsyssc_lumpMatrixScalar(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                     LSYSSC_LUMP_DIAG)
      end if
    elseif (lumpedMassMatrix > 0) then
      call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(templateMatrix),&
                                  rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                      DER_FUNC, DER_FUNC) 
      call lsyssc_lumpMatrixScalar(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                   LSYSSC_LUMP_DIAG)
    end if
    
    ! Create diffusion matrix as duplicate of the template matrix
    if (coeffMatrix_S > 0) then
      if (lsyssc_isMatrixStructureShared(rproblemLevel%Rmatrix(coeffMatrix_S),&
                                         rproblemLevel%Rmatrix(templateMatrix))) then
        call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(coeffMatrix_S),&
                                 rproblemLevel%Rmatrix(templateMatrix),&
                                 .false., .false., .true.)
      else
        call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(templateMatrix),&
                                    rproblemLevel%Rmatrix(coeffMatrix_S),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      end if

      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        call initDiffusionMatrix1D(rappDescriptor%rfparserDiffusionTensor,&
                                   rappDescriptor%idiffusiontype,&
                                   rproblemLevel%Rmatrix(coeffMatrix_S))
      case (NDIM2D)
        call initDiffusionMatrix2D(rappDescriptor%rfparserDiffusionTensor,&
                                   rappDescriptor%idiffusiontype,&
                                   rproblemLevel%Rmatrix(coeffMatrix_S))
      case (NDIM3D)
        call initDiffusionMatrix3D(rappDescriptor%rfparserDiffusionTensor,&
                                   rappDescriptor%idiffusiontype,&
                                   rproblemLevel%Rmatrix(coeffMatrix_S))
      case DEFAULT
        call lsyssc_releaseMatrix(rproblemLevel%Rmatrix(coeffMatrix_S))
      end select
    end if
    
    ! Create coefficient matrix (phi, dphi/dx) duplicate of the template matrix
    if (coeffMatrix_CX > 0) then
      if (lsyssc_isMatrixStructureShared(rproblemLevel%Rmatrix(coeffMatrix_CX),&
                                         rproblemLevel%Rmatrix(templateMatrix))) then
        call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(coeffMatrix_CX),&
                                 rproblemLevel%Rmatrix(templateMatrix),&
                                 .false., .false., .true.)
      else
        call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(templateMatrix),&
                                    rproblemLevel%Rmatrix(coeffMatrix_CX),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      end if
      call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(coeffMatrix_CX),&
                                      DER_DERIV3D_X, DER_FUNC)
    end if
    
    ! Create coefficient matrix (phi, dphi/dy) duplicate of the template matrix
    if (coeffMatrix_CY > 0) then
      if (lsyssc_isMatrixStructureShared(rproblemLevel%Rmatrix(coeffMatrix_CY),&
                                         rproblemLevel%Rmatrix(templateMatrix))) then
        call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(coeffMatrix_CY),&
                                 rproblemLevel%Rmatrix(templateMatrix),&
                                 .false., .false., .true.)
      else
        call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(templateMatrix),&
                                    rproblemLevel%Rmatrix(coeffMatrix_CY),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      end if
      call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(coeffMatrix_CY),&
                                      DER_DERIV3D_Y, DER_FUNC)
    end if
    
    ! Create coefficient matrix (phi, dphi/dz) duplicate of the template matrix
    if (coeffMatrix_CZ > 0) then
      if (lsyssc_isMatrixStructureShared(rproblemLevel%Rmatrix(coeffMatrix_CZ),&
                                         rproblemLevel%Rmatrix(templateMatrix))) then
        call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(coeffMatrix_CZ),&
                                 rproblemLevel%Rmatrix(templateMatrix),&
                                 .false., .false., .true.)
      else
        call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(templateMatrix),&
                                    rproblemLevel%Rmatrix(coeffMatrix_CZ),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      end if
      call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(coeffMatrix_CZ),&
                                      DER_DERIV3D_Z, DER_FUNC)
    end if

    ! Resize stabilization structures if necessary and remove the
    ! indicator for the subdiagonal edge structure. If they are
    ! needed, then they are re-generated on-the-fly.
    if (convectionAFC > 0) then
      if (rproblemLevel%Rafcstab(convectionAFC)%iSpec .ne. AFCSTAB_UNDEFINED) then
        call afcstab_resizeStabilisation(rproblemLevel%Rafcstab(convectionAFC),&
                                         rproblemLevel%Rmatrix(templateMatrix))
        rproblemLevel%Rafcstab(convectionAFC)%iSpec =&
            iand(rproblemLevel%Rafcstab(convectionAFC)%iSpec,&
            not(AFCSTAB_SUBDIAGONALEDGES))
      end if
    end if

    if (diffusionAFC > 0) then
      if (rproblemLevel%Rafcstab(diffusionAFC)%iSpec .ne. AFCSTAB_UNDEFINED) then
        call afcstab_resizeStabilisation(rproblemLevel%Rafcstab(diffusionAFC),&
                                         rproblemLevel%Rmatrix(templateMatrix))
        rproblemLevel%Rafcstab(diffusionAFC)%iSpec =&
            iand(rproblemLevel%Rafcstab(diffusionAFC)%iSpec,&
            not(AFCSTAB_SUBDIAGONALEDGES))
      end if
    end if

    ! Set update notification in problem level structure
    rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                     PROBLEV_MSPEC_UPDATE)
    
  contains

    ! Here, the initialization routines for the diffusion matrix follow.
    
    !**************************************************************

    subroutine initDiffusionMatrix1D(rfparser, idiffusiontype, rmatrix)
      type(t_fparser), intent(IN) :: rfparser
      integer, intent(IN) :: idiffusiontype
      type(t_matrixScalar), intent(INOUT) :: rmatrix
      
      ! local variables
      real(DP), dimension(1) :: Dunity = (/1.0_DP/)
      real(DP) :: dalpha
      
      select case(idiffusiontype)
      case (DIFFUSION_ISOTROPIC,&
            DIFFUSION_ANISOTROPIC)
        ! Evaluate the constant coefficient from the function parser
        call fparser_evalFunction(rfparser, 1, Dunity, dalpha)
        
        ! Assemble the Laplace matrix multiplied by the negative value
        ! of the parameter alpha
        call stdop_assembleLaplaceMatrix(rmatrix, .true., -dalpha)
        
      case DEFAULT
        call lsyssc_clearMatrix(rmatrix)
      end select
      
    end subroutine initDiffusionMatrix1D

    !**************************************************************
    
    subroutine initDiffusionMatrix2D(rfparser, idiffusiontype, rmatrix)
      type(t_fparser), intent(IN) :: rfparser
      integer, intent(IN) :: idiffusiontype
      type(t_matrixScalar), intent(INOUT) :: rmatrix

      ! local variables
      type(t_bilinearform) :: rform
      real(DP), dimension(1) :: Dunity = (/1.0_DP/)
      real(DP) :: dalpha
      
      select case(idiffusiontype)
      case (DIFFUSION_ISOTROPIC)
        ! Evaluate the constant coefficient from the function parser
        call fparser_evalFunction(rfparser, 1, Dunity, dalpha)
        
        ! Assemble the Laplace matrix multiplied by the negative value
        ! of the parameter alpha
        call stdop_assembleLaplaceMatrix(rmatrix, .true., -dalpha)
        
      case (DIFFUSION_ANISOTROPIC)
        ! Evaluate the constant coefficient from the function parser
        call fparser_evalFunction(rfparser, 1, Dunity, rform%Dcoefficients(1))
        call fparser_evalFunction(rfparser, 2, Dunity, rform%Dcoefficients(2))
        call fparser_evalFunction(rfparser, 3, Dunity, rform%Dcoefficients(3))
        call fparser_evalFunction(rfparser, 4, Dunity, rform%Dcoefficients(4))
        
        ! We have constant coefficients
        rform%ballCoeffConstant = .true.
        rform%BconstantCoeff    = .true.
        rform%Dcoefficients     = -rform%Dcoefficients
        
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
    
    subroutine initDiffusionMatrix3D(rfparser, idiffusiontype, rmatrix)
      type(t_fparser), intent(IN) :: rfparser
      integer, intent(IN) :: idiffusiontype
      type(t_matrixScalar), intent(INOUT) :: rmatrix
      
      ! local variables
      type(t_bilinearform) :: rform
      real(DP), dimension(1) :: Dunity = (/1.0_DP/)
      real(DP) :: dalpha
      
      select case(idiffusiontype)
      case (DIFFUSION_ISOTROPIC)
        ! Evaluate the constant coefficient from the function parser
        call fparser_evalFunction(rfparser, 1, Dunity, dalpha)
        
        ! Assemble the Laplace matrix multiplied by the negative value
        ! of the parameter alpha
        call stdop_assembleLaplaceMatrix(rmatrix, .true., -dalpha)
        
      case (DIFFUSION_ANISOTROPIC)
        ! Evaluate the constant coefficient from the function parser
        call fparser_evalFunction(rfparser, 1, Dunity, rform%Dcoefficients(1))
        call fparser_evalFunction(rfparser, 2, Dunity, rform%Dcoefficients(2))
        call fparser_evalFunction(rfparser, 3, Dunity, rform%Dcoefficients(3))
        call fparser_evalFunction(rfparser, 4, Dunity, rform%Dcoefficients(4))
        call fparser_evalFunction(rfparser, 5, Dunity, rform%Dcoefficients(5))
        call fparser_evalFunction(rfparser, 6, Dunity, rform%Dcoefficients(6))
        call fparser_evalFunction(rfparser, 7, Dunity, rform%Dcoefficients(7))
        call fparser_evalFunction(rfparser, 8, Dunity, rform%Dcoefficients(8))
        call fparser_evalFunction(rfparser, 9, Dunity, rform%Dcoefficients(9))
        
        ! We have constant coefficients
        rform%ballCoeffConstant = .true.
        rform%BconstantCoeff    = .true.
        rform%Dcoefficients     = -rform%Dcoefficients
        
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
    
  end subroutine codire_initProblemLevel

  !*****************************************************************************

!<subroutine>

  subroutine codire_initAllProblemLevels(rappDescriptor, rproblem, rcollection)

!<description>
    ! This subroutine initializes the all problem levels attached to
    ! the global problem structure. It generates the discretization,
    ! the template matrix and the coefficient matrices as duplicates
    ! of the template matrix.
!</description>

!<input>
    ! application descriptor
    type(t_codire), intent(IN) :: rappDescriptor
!</input>

!<inputoutput>
    ! problem structure
    type(t_problem), intent(INOUT) :: rproblem

    ! collection structure
    type(t_collection), intent(INOUT) :: rcollection
!</intputoutput>
!</subroutine>

    ! pointer to the problem level
    type(t_problemLevel), pointer :: p_rproblemLevel

    
    ! loop over all problem levels
    p_rproblemLevel => rproblem%p_rproblemLevelMax
    do while(associated(p_rproblemLevel))
      
      ! Initialize individual problem level
      call codire_initProblemLevel(rappDescriptor, p_rproblemLevel, rcollection)
      
      ! Switch to next coarser level
      p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
    end do

  end subroutine codire_initAllProblemLevels

  !*****************************************************************************

!<subroutine>

  subroutine codire_initSolution(rparlist, ssectionName, rproblemLevel,&
                                 dtime, rvector)

!<description>
    ! This subroutine initializes the solution vector
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(IN) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(IN) :: ssectionName

    ! problem level
    type(t_problemLevel), intent(IN), target :: rproblemLevel

    ! time for solution evaluation
    real(DP), intent(IN) :: dtime
!</input>

!<output>
    ! solution vector
    type(t_vectorBlock), intent(OUT) :: rvector
!</output>
!</subroutine>

    ! section names
    character(LEN=SYS_STRLEN) :: sindatfileName
    character(LEN=SYS_STRLEN) :: ssolutionName

    ! local variables
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    type(t_fparser) :: rfparser
    type(t_pgm) :: rpgm
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Ddata
    real(DP), dimension(NDIM3D+1) :: Dvalue
    integer :: isolutiontype,iblock,ieq,neq,ndim

    ! symbolic variable names
    character(LEN=*), dimension(4), parameter ::&
                      cvariables = (/ (/'x'/), (/'y'/), (/'z'/), (/'t'/) /)
    
    
    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist, ssectionName, 'indatfile', sindatfileName)
    call parlst_getvalue_string(rparlist, ssectionName, 'ssolutionname', ssolutionName)
    call parlst_getvalue_int(rparlist, ssectionName, 'isolutiontype', isolutiontype)
    
    ! Create new solution vector based on the spatial discretisation
    p_rdiscretisation => rproblemLevel%rdiscretisation
    call lsysbl_createVectorBlock(p_rdiscretisation, rvector, .false., ST_DOUBLE)

    ! How should the solution be initialized?
    select case(isolutionType)
    case (0)
      ! Initialize solution by zeros
      call lsysbl_clearVector(rvector)

    case (1)
      ! Initialize solution from analytic profile
      call flagship_readParserFromFile(sindatfileName, '['//trim(ssolutionName)//']',&
                                       cvariables, rfparser)
      
      ! Set pointers
      call storage_getbase_double2D(&
          rproblemLevel%rtriangulation%h_DvertexCoords, p_DvertexCoords)
      
      ndim = rproblemLevel%rtriangulation%ndim
      
      ! Initialize variable values
      Dvalue = 0.0_DP
      
      ! Loop over all blocks of the solution vector
      do iblock = 1, rvector%nblocks
        
        ! Get number of equations of scalar subvector
        neq = rvector%RvectorBlock(iblock)%NEQ
        call lsyssc_getbase_double(rvector%RvectorBlock(iblock), p_Ddata)
        
        ! Loop over all equations of scalar subvector
        do ieq = 1, neq
          Dvalue(1:ndim)   = p_DvertexCoords(:,ieq)
          Dvalue(NDIM3D+1) = dtime
          call fparser_evalFunction(rfparser, iblock, Dvalue, p_Ddata(ieq))
        end do
      end do
      
      call fparser_release(rfparser)
      

    case (2)
      ! Initialize solution from portable graymap image
      call ppsol_readPGM(0, ssolutionName, rpgm)
      
      ! Set pointers
      call storage_getbase_double2D(&
          rproblemLevel%rtriangulation%h_DvertexCoords, p_DvertexCoords)
      call storage_getbase_double(rvector%h_Ddata, p_Ddata)
    
      ! Initialize the solution by the image data
      call ppsol_initArrayPGM_Dble(rpgm, p_DvertexCoords, p_Ddata)
      
      ! Release portable graymap image
      call ppsol_releasePGM(rpgm)
      
      
    case DEFAULT
      call output_line('Invalid type of solution profile!',&
                       OU_CLASS_ERROR, OU_MODE_STD, 'codire_initSolution')
      call sys_halt()
    end select

  end subroutine codire_initSolution

  !*****************************************************************************

!<subroutine>

  subroutine codire_initRHS(rappDescriptor, rproblemLevel,&
                            dtime, rvector)

!<description>
    ! This subroutine initializes the right-hand side vector.
!</description>

!<input>
    ! application descriptor
    type(t_codire), intent(IN) :: rappDescriptor

    ! problem level structure
    type(t_problemLevel), intent(IN), target :: rproblemLevel

    ! time for right-hand side evaluation
    real(DP), intent(IN) :: dtime
!</input>

!<output>
    ! right-hand side vector
    type(t_vectorBlock), intent(OUT) :: rvector
!</output>
!</subroutine>

    ! local variables
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    type(t_collection) :: rcollection
    type(t_linearForm) :: rform


    ! Create new block vector based on the spatial discretisation
    p_rdiscretisation => rproblemLevel%rdiscretisation
    call lsysbl_createVectorBlock(p_rdiscretisation, rvector, .false., ST_DOUBLE)
    
    ! How should the target functional be initialized?
    select case(rappDescriptor%irhstype)
    case (RHS_ZERO)
      call lsysbl_clearVector(rvector)


    case (RHS_ANALYTIC)
      ! Create a collection and attach the function parser
      call collct_init(rcollection)
      call collct_setvalue_pars(rcollection, 'rhsParser',&
                                rappDescriptor%rfparserRHS, .true.)
      rcollection%SquickAccess(1) = 'rhsParser'
      rcollection%DquickAccess(1) = dtime
      
      ! Set up the corresponding linear form
      rform%itermCount      = 1
      rform%Idescriptors(1) = DER_FUNC

      ! Build the discretized target functional
      call linf_buildVectorScalar(rvector%RvectorBlock(1)%p_rspatialDiscr,&
                                  rform, .true., rvector%RvectorBlock(1),&
                                  codire_coeffRHS, rcollection)

      ! Release collection structure
      call collct_done(rcollection)

      
    case DEFAULT
      call output_line('Invalid type of target functional!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'codire_initRHS')
      call sys_halt()
    end select

  end subroutine codire_initRHS

  !*****************************************************************************

!<subroutine>

  subroutine codire_initTargetFunc(rappDescriptor, rproblemLevel,&
                                   dtime, rvector)

!<description>
    ! This subroutine initializes the target functional which serves as
    ! right-hand side vector for the dual problem in the framework of
    ! goal-oriented error estimation.
!</description>

!<input>
    ! application descriptor
    type(t_codire), intent(IN) :: rappDescriptor

    ! problem level structure
    type(t_problemLevel), intent(IN), target :: rproblemLevel

    ! time for target function evaluation
    real(DP), intent(IN) :: dtime
!</input>

!<output>
    ! target function vector
    type(t_vectorBlock), intent(OUT) :: rvector
!</output>
!</subroutine>

    ! local variables
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    type(t_collection) :: rcollection
    type(t_linearForm) :: rform


    ! Create new block vector based on the spatial discretisation
    p_rdiscretisation => rproblemLevel%rdiscretisation
    call lsysbl_createVectorBlock(p_rdiscretisation, rvector, .false., ST_DOUBLE)
    
    ! How should the target functional be initialized?
    select case(rappDescriptor%itargetfunctype)
    case (TFUNC_ZERO)
      call lsysbl_clearVector(rvector)


    case (TFUNC_VOLINTG)
      ! Create a collection and attach the function parser
      call collct_init(rcollection)
      call collct_setvalue_pars(rcollection, 'targetFuncParser',&
                                rappDescriptor%rfparserTargetFunc, .true.)
      rcollection%SquickAccess(1) = 'targetFuncParser'
      rcollection%DquickAccess(1) = dtime
      
      ! Set up the corresponding linear form
      rform%itermCount      = 1
      rform%Idescriptors(1) = DER_FUNC

      ! Build the discretized target functional
      call linf_buildVectorScalar(rvector%RvectorBlock(1)%p_rspatialDiscr,&
                                  rform, .true., rvector%RvectorBlock(1),&
                                  codire_coeffTargetFuncVolInt, rcollection)

      ! Release collection structure
      call collct_done(rcollection)

      
    case DEFAULT
      call output_line('Invalid type of target functional!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'codire_initTargetFunc')
      call sys_halt()
    end select
    
  end subroutine codire_initTargetFunc

  !*****************************************************************************

!<subroutine>

  subroutine codire_outputSolution(rparlist, ssectionName, rproblemLevel,&
                                   rsolutionPrimal, rsolutionDual, dtime)

!<description>
    ! This subroutine exports the solution vector to file in UCD format
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(IN) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(IN) :: ssectionName

    ! problem level structure
    type(t_problemLevel), intent(IN) :: rproblemLevel

    ! solution vector for primal problem
    type(t_vectorBlock), intent(IN) :: rsolutionPrimal

    ! OPTIONAL: solution vector for dual problem
    type(t_vectorBlock), intent(IN), optional :: rsolutionDual

    ! OPTIONAL: simulation time
    real(DP), intent(IN), optional :: dtime
!</input>
!</subroutine>

    ! section names
    character(LEN=SYS_STRLEN) :: soutputName
    character(LEN=SYS_STRLEN) :: ucdsolution

    ! persistent variable
    integer, save :: ifilenumber = 1
    
    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata
    type(t_ucdExport) :: rexport
    integer :: iformatUCD


    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist, ssectionName, 'output', soutputName)
    call parlst_getvalue_string(rparlist, trim(soutputName), 'ucdsolution', ucdsolution)
    call parlst_getvalue_int(rparlist, trim(soutputName), 'iformatucd', iformatUCD)

    ! Initialize the UCD exporter
    call flagship_initUCDexport(rproblemLevel, ucdsolution,&
                                iformatUCD, rexport, ifilenumber)

    ! Increase filenumber by one
    ifilenumber = ifilenumber+1

    ! Set simulation time
    if (present(dtime)) call ucd_setSimulationTime(rexport, dtime)

    ! Add primal solution vector
    call lsysbl_getbase_double(rsolutionPrimal, p_Ddata)
    call ucd_addVariableVertexBased (rexport, 'u', UCD_VAR_STANDARD, p_Ddata)

    ! Add dual solution vector
    if (present(rsolutionDual)) then
      call lsysbl_getbase_double(rsolutionDual, p_Ddata)
      call ucd_addVariableVertexBased (rexport, 'z', UCD_VAR_STANDARD, p_Ddata)
    end if

    ! Write UCD file
    call ucd_write  (rexport)
    call ucd_release(rexport)

  end subroutine codire_outputSolution
  
  !*****************************************************************************

!<subroutine>

  subroutine codire_outputStatistics(rappDescriptor, rtimerTotal)

!<description>
    ! This subroutine output application statistics
!</description>

!<input>
    ! application descriptor
    type(t_codire), intent(INOUT) :: rappDescriptor

    ! timer for total time measurement
    type(t_timer), intent(IN) :: rtimerTotal
!</input>
!</subroutine>

    ! local variable
    type(t_timer) :: rtimerSolution
    real(DP) :: dtotalTime, dfraction

    call output_lbrk(nlbrk=5)
    call output_separator (OU_SEP_STAR)
    write(*,FMT='(24X,A)') '*** TIME MEASUREMENT ***'
    call output_separator (OU_SEP_STAR)
    
    rtimerSolution = rappDescriptor%rtimerSolution
    call stat_subTimers(rappDescriptor%rtimerAssemblyMatrix, rtimerSolution)
    call stat_subTimers(rappDescriptor%rtimerAssemblyVector, rtimerSolution)

    dtotalTime = max(rtimerTotal%delapsedCPU, rtimerTotal%delapsedReal)
    dfraction  = 100.0_DP/dtotalTime

    call output_line('  Time for computing solution   : '//&
                     trim(sys_sdL(rtimerSolution%delapsedCPU, 2))//'  '//&
                     trim(sys_sdL(dfraction*rtimerSolution%delapsedCPU, 2))//' %')
    call output_line('  Time for mesh adaptivity      : '//&
                     trim(sys_sdL(rappDescriptor%rtimerAdaptation%delapsedCPU, 2))//'  '//&
                     trim(sys_sdL(dfraction*rappDescriptor%rtimerAdaptation%delapsedCPU, 2))//' %')
    call output_line('  Time for error estimation     : '//&
                     trim(sys_sdL(rappDescriptor%rtimerErrorEstimation%delapsedCPU, 2))//'  '//&
                     trim(sys_sdL(dfraction*rappDescriptor%rtimerErrorEstimation%delapsedCPU, 2))//' %')
    call output_line('  Time for triangulation        : '//&
                     trim(sys_sdL(rappDescriptor%rtimerTriangulation%delapsedCPU, 2))//'  '//&
                     trim(sys_sdL(dfraction*rappDescriptor%rtimerTriangulation%delapsedCPU, 2))//' %')
    call output_line('  Time for coefficient assembly : '//&
                     trim(sys_sdL(rappDescriptor%rtimerAssemblyCoeff%delapsedCPU, 2))//'  '//&
                     trim(sys_sdL(dfraction*rappDescriptor%rtimerAssemblyCoeff%delapsedCPU, 2))//' %')
    call output_line('  Time for matrix assembly      : '//&
                     trim(sys_sdL(rappDescriptor%rtimerAssemblyMatrix%delapsedCPU, 2))//'  '//&
                     trim(sys_sdL(dfraction*rappDescriptor%rtimerAssemblyMatrix%delapsedCPU, 2))//' %')
    call output_line('  Time for vector assembly:       '//&
                     trim(sys_sdL(rappDescriptor%rtimerAssemblyVector%delapsedCPU, 2))//'  '//&
                     trim(sys_sdL(dfraction*rappDescriptor%rtimerAssemblyVector%delapsedCPU, 2))//' %')
    call output_line('  Time for pre-/post-processing : '//&
                     trim(sys_sdL(rappDescriptor%rtimerPrePostprocess%delapsedCPU, 2))//'  '//&
                     trim(sys_sdL(dfraction*rappDescriptor%rtimerPrePostprocess%delapsedCPU, 2))//' %')
    call output_separator (OU_SEP_MINUS)
    call output_line('  Time for total simulation     : '//&
                     trim(sys_sdL(dtotalTime, 2)))
    call output_lbrk()
  end subroutine codire_outputStatistics

  !*****************************************************************************

!<subroutine>

  subroutine codire_estimateTargetFuncError(rproblemLevel, rtimestep, rsolver,&
                                            rsolutionPrimal, rsolutionDual,&
                                            rcollection, rerror, derror, rrhs)

!<description>
    ! This subroutine estimates the error in the quantity of interest
!</description>

!<input>
    ! time-stepping algorithm
    type(t_timestep), intent(IN) :: rtimestep

    ! primal solution vector
    type(t_vectorBlock), intent(IN) :: rsolutionPrimal

    ! dual solution vector
    type(t_vectorBlock), intent(IN) :: rsolutionDual

    ! OPTIONAL: right-hand side vector
    type(t_vectorBlock), intent(IN), optional :: rrhs
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(INOUT) :: rsolver

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>

!<output>
    ! element-wise error distribution
    type(t_vectorScalar), intent(OUT) :: rerror

    ! global error
    real(DP), intent(OUT) :: derror
!</subroutine>


    ! local variables
    type(t_matrixScalar) :: rmatrix
    type(t_vectorBlock) :: rvector
    real(DP), dimension(:), pointer :: p_DsolutionDual, p_DnodalError, p_DlumpedMassMatrix
    integer :: i,convectionAFC,diffusionAFC
    integer :: lumpedMassMatrix,templateMatrix

    
    ! Create vector for nodal errors
    call lsysbl_createVectorBlock(rsolutionPrimal, rvector)

    ! Ok, this is a little bit tricky. We compute the standard residual vector
    ! for the zeroth iteration and switch off all types of stabilization.
    convectionAFC = collct_getvalue_int(rcollection, 'convectionAFC')
    diffusionAFC  = collct_getvalue_int(rcollection, 'diffusionAFC')

    call collct_setvalue_int(rcollection, 'convectionAFC', -abs(convectionAFC), .false.)
    call collct_setvalue_int(rcollection, 'diffusionAFC',  -abs(diffusionAFC), .false.)

    ! Calculate the standard Galerkin residual
    call codire_calcResidual(rproblemLevel, rtimestep, rsolver,&
                             rsolutionPrimal, rsolutionPrimal,&
                             rvector, rvector, 0, rcollection)
    
    ! Add contribution from the constant right-hand side
    if (present(rrhs)) call lsysbl_vectorLinearComb(rrhs, rvector, 1.0_DP, 1.0_DP)
    
    ! Ok, now we have to switch on all types of stabilization
    call collct_setvalue_int(rcollection, 'convectionAFC', convectionAFC, .false.)
    call collct_setvalue_int(rcollection, 'diffusionAFC',  diffusionAFC, .false.)
    
    ! Set pointers
    call lsysbl_getbase_double(rsolutionDual, p_DsolutionDual)
    call lsysbl_getbase_double(rvector, p_DnodalError)

    ! Compute/set pointer to the lumped mass matrix
    lumpedMassMatrix = collct_getvalue_int(rcollection, 'lumpedmassmatrix')
    if (lumpedMassMatrix > 0) then
      call lsyssc_getbase_double(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                 p_DlumpedMassMatrix)
    else
      templateMatrix = collct_getvalue_int(rcollection, 'templatematrix')
      call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(templateMatrix),&
                                  rmatrix, LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      call stdop_assembleSimpleMatrix(rmatrix, DER_FUNC, DER_FUNC) 
      call lsyssc_lumpMatrixScalar(rmatrix, LSYSSC_LUMP_DIAG)
      call lsyssc_getbase_double(rmatrix, p_DlumpedMassMatrix)
    end if
    
    ! Now we compute the nodal error contributions
    do i = 1, size(p_DsolutionDual,1)
      p_DnodalError(i) = abs(p_DnodalError(i)*p_DsolutionDual(i))/p_DlumpedMassMatrix(i)
    end do

    ! Compute the element-wise and the global error
    call lsyssc_createVector(rerror, rproblemLevel%rtriangulation%NEL, .false.)
    call pperr_scalar(rvector%RvectorBlock(1), PPERR_L1ERROR,&
                      derror, relementError=rerror)
 
    ! Release the nodal error
    call lsysbl_releaseVector(rvector)

    ! Clean temporal matrix if required
    if (lumpedMassMatrix .le. 0) call lsyssc_releaseMatrix(rmatrix)

  end subroutine codire_estimateTargetFuncError

  !*****************************************************************************

!<subroutine>

  subroutine codire_adaptTriangulation(rhadapt, rtriangulationSrc, rindicator,&
                                       rcollection, rtriangulationDest)

!<description>
    ! This subroutine performs h-adaptation for the given triangulation
!</description>

!<input>
    ! element-wise indicator
    type(t_vectorScalar), intent(IN) :: rindicator
!</input>

!<inputoutput>
    ! source triangulation structure
    type(t_triangulation), intent(INOUT), target :: rtriangulationSrc

    ! adaptation structure
    type(t_hadapt), intent(INOUT) :: rhadapt

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>

!<output>
    ! OPTIONAL: destination triangulation structure
    ! If it is not given, the source triangulation is updated
    type(t_triangulation), intent(OUT), optional, target :: rtriangulationDest
!</output>
!</subroutine>
    

    ! local variables
    type(t_triangulation), pointer :: p_rtriangulation
    integer, dimension(1) :: Ivalue = 0


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
      call codire_hadaptCallback1D(rcollection, HADAPT_OPR_INITCALLBACK, Ivalue, Ivalue)
      call hadapt_performAdaptation(rhadapt, rindicator, rcollection, codire_hadaptCallback1D)
      call codire_hadaptCallback1D(rcollection, -999, Ivalue, Ivalue)
      
    case (NDIM2D)
      call codire_hadaptCallback2D(rcollection, HADAPT_OPR_INITCALLBACK, Ivalue, Ivalue)
      call hadapt_performAdaptation(rhadapt, rindicator, rcollection, codire_hadaptCallback2D)
      call codire_hadaptCallback2D(rcollection, -999, Ivalue, Ivalue)

    case (NDIM3D)
      call codire_hadaptCallback3D(rcollection, HADAPT_OPR_INITCALLBACK, Ivalue, Ivalue)
      call hadapt_performAdaptation(rhadapt, rindicator, rcollection, codire_hadaptCallback3D)
      call codire_hadaptCallback3D(rcollection, -999, Ivalue, Ivalue)

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

  end subroutine codire_adaptTriangulation

  !*****************************************************************************

!<subroutine>

    subroutine codire_solveTransientPrimal(rappDescriptor, rparlist, ssectionName,&
                                           rbdrCond, rproblem, rtimestep, rsolver,&
                                           rsolution, rcollection)

!<description>
    ! This subroutine solves the transient primal flow problem
    !
    !  $$\frac{\partial u}{\partial t}+\nabla\cdot{\bf f}(u)=s(u)$$
    !
    ! for a scalar quantity $u$ in the domain $\Omega$.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(IN) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(IN) :: ssectionName

    ! boundary condition structure
    type(t_boundaryCondition), intent(IN) :: rbdrCond
!</input>

!<inputoutput>
    ! application descriptor
    type(t_codire), intent(INOUT) :: rappDescriptor

    ! problem structure
    type(t_problem), intent(INOUT) :: rproblem

    ! time-stepping structure
    type(t_timestep), intent(INOUT) :: rtimestep

    ! solver struchture
    type(t_solver), intent(INOUT), target :: rsolver

    ! primal solution vector
    type(t_vectorBlock), intent(INOUT) :: rsolution

    ! collection structure
    type(t_collection), intent(INOUT) :: rcollection    
!</inputoutput>
!</subroutine>

    ! Pointer to the current multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! Vector for the linear target functional and the right-hand side
    type(t_vectorBlock) :: rvectorBlock

    ! Vector for the element-wise error distribution
    type(t_vectorScalar) :: rvectorScalar
    
    ! Structure for h-adaptation
    type(t_hadapt) :: rhadapt

    ! Structure for the sparsity pattern
    type(t_graph) :: rgraph

    ! section names
    character(LEN=SYS_STRLEN) :: sadaptivityName
    character(LEN=SYS_STRLEN) :: soutputName

    ! local variables
    real(dp) :: derror,dstepUCD,dtimeUCD,dstepAdapt,dtimeAdapt
    integer :: templateMatrix
    integer :: nlmin


    ! Set pointer to maximum problem level
    p_rproblemLevel => rproblem%p_rproblemLevelMax

    ! Initialize timer for intermediate UCD exporter
    dtimeUCD = rtimestep%dinitialTime
    call parlst_getvalue_string(rparlist, ssectionName, 'output', soutputName)
    call parlst_getvalue_double(rparlist, trim(soutputName), 'dstepUCD', dstepUCD, 0.0_DP)


    !---------------------------------------------------------------------------
    ! Initialize the h-adaptation structure
    !---------------------------------------------------------------------------

    call parlst_getvalue_string(rparlist, ssectionName, 'adaptivity', sadaptivityName, '')
    if (trim(adjustl(sadaptivityName)) .ne. '') then

      call parlst_getvalue_double(rparlist, trim(sadaptivityName), 'dstepAdapt', dstepAdapt, 0.0_DP)
      call parlst_getvalue_double(rparlist, trim(sadaptivityName), 'dtimeAdapt', dtimeAdapt, 0.0_DP)

      if (dstepAdapt > 0.0_DP) then

        ! Initialize adaptation structure from parameter list
        call hadapt_initFromParameterlist(rhadapt, rparlist, sadaptivityName)

        ! Generate a dynamic graph for the sparsity pattern and attach
        ! it to the collection structure which is used to pass this type
        ! to the callback function for h-adaptation
        templateMatrix = collct_getvalue_int(rcollection, 'templateMatrix')
        call grph_createGraphFromMatrix(p_rproblemLevel%Rmatrix(templateMatrix), rgraph)
        call collct_setvalue_graph(rcollection, 'sparsitypattern', rgraph, .true.)
        
        ! Attach the primal solution vector to the collection structure
        call collct_setvalue_vec(rcollection, 'solutionvector', rsolution, .true.)

      end if   ! dstepAdapt > 0

    else

      dstepAdapt = 0.0_DP

    end if

    !---------------------------------------------------------------------------
    ! Initialize the solver for the primal problem
    !---------------------------------------------------------------------------

    ! Calculate the velocity field
    nlmin = solver_getMinimumMultigridlevel(rsolver)
    call codire_calcVelocityField(rappDescriptor, p_rproblemLevel,&
                                  rtimestep%dtime, rcollection, nlmin)
    
    ! Attach the boundary condition to the solver structure
    call solver_setBoundaryCondition(rsolver, rbdrCond, .true.)

    ! Set collection to primal problem mode
    call collct_setvalue_int(rcollection, 'primaldual', 1, .true.)

    !---------------------------------------------------------------------------
    ! Infinite time stepping loop
    !---------------------------------------------------------------------------

    timeloop: do
      
      ! Check for user interaction
      call codire_UserInterface
      
      !-------------------------------------------------------------------------
      ! Advance solution in time
      !-------------------------------------------------------------------------
      
      ! Start time measurement for solution procedure
      call stat_startTimer(rappDescriptor%rtimerSolution, STAT_TIMERSHORT)
      
      ! Check if right-hand side vector exists
      if (rappDescriptor%irhstype > 0) then

        ! Initialize right-hand side vector
        call codire_initRHS(rappDescriptor, p_rproblemLevel,&
                            0.0_DP, rvectorBlock)

        ! What time-stepping scheme should be used?
        select case(rtimestep%ctimestepType)
          
        case (SV_RK_SCHEME)
          
          ! Adopt explicit Runge-Kutta scheme
          call tstep_performRKStep(p_rproblemLevel, rtimestep,&
                                   rsolver, rsolution, codire_calcRHS,&
                                   codire_setBoundary, rcollection, rvectorBlock)
          
        case (SV_THETA_SCHEME)
          
          ! Adopt two-level theta-scheme
          call tstep_performThetaStep(p_rproblemLevel, rtimestep,&
                                      rsolver, rsolution, codire_calcResidual,&
                                      codire_calcJacobian, codire_applyJacobian,&
                                      codire_setBoundary, rcollection, rvectorBlock)
          
        case DEFAULT
          call output_line('Unsupported time-stepping algorithm!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'codire_solveTransient')
          call sys_halt()
        end select

        ! Release right-hand side vector
        call lsysbl_releaseVector(rvectorBlock)
        
      else

        ! What time-stepping scheme should be used?
        select case(rtimestep%ctimestepType)
          
        case (SV_RK_SCHEME)
          
          ! Adopt explicit Runge-Kutta scheme
          call tstep_performRKStep(p_rproblemLevel, rtimestep,&
                                   rsolver, rsolution, codire_calcRHS,&
                                   codire_setBoundary, rcollection)
          
        case (SV_THETA_SCHEME)
          
          ! Adopt two-level theta-scheme
          call tstep_performThetaStep(p_rproblemLevel, rtimestep,&
                                      rsolver, rsolution, codire_calcResidual,&
                                      codire_calcJacobian, codire_applyJacobian,&
                                      codire_setBoundary, rcollection)
          
        case DEFAULT
          call output_line('Unsupported time-stepping algorithm!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'codire_solveTransient')
          call sys_halt()
        end select

      end if

      ! Stop time measurement for solution procedure
      call stat_stopTimer(rappDescriptor%rtimerSolution)
      
      ! Reached final time, then exit the infinite time loop?
      if (rtimestep%dTime .ge. rtimestep%dfinalTime) exit timeloop


      !-------------------------------------------------------------------------
      ! Post-process intermediate solution
      !-------------------------------------------------------------------------
      
      if (dstepUCD .gt. 0._DP .and. rtimestep%dTime .ge. dtimeUCD) then

        ! Set time for next intermediate solution export
        dtimeUCD = dtimeUCD + dstepUCD

        ! Start time measurement for post-processing
        call stat_startTimer(rappDescriptor%rtimerPrepostProcess, STAT_TIMERSHORT)

        ! Export the intermediate solution
        call codire_outputSolution(rparlist, ssectionName, p_rproblemLevel,&
                                   rsolution, dtime=rtimestep%dTime)

        ! Stop time measurement for post-processing
        call stat_stopTimer(rappDescriptor%rtimerPrepostProcess)

      end if


      if (dstepAdapt .gt. 0.0_DP .and. rtimestep%dTime .ge. dtimeAdapt) then
      
        ! Set time for next adaptation step
        dtimeAdapt = dtimeAdapt + dstepAdapt

        !-----------------------------------------------------------------------
        ! Perform recovery-based error estimation
        !-----------------------------------------------------------------------
        
        ! Start time measurement for error estimation
        call stat_startTimer(rappDescriptor%rtimerErrorEstimation, STAT_TIMERSHORT)
        
        ! HERE WE NEED TO IMPLEMENT THE RECOVERY-BASED INDICATOR
        
        ! Stop time measurement for error estimation
        call stat_stopTimer(rappDescriptor%rtimerErrorEstimation)
        
        
        !-------------------------------------------------------------------------
        ! Perform h-adaptation
        !-------------------------------------------------------------------------
        
        ! Start time measurement for mesh adaptation
        call stat_startTimer(rappDescriptor%rtimerAdaptation, STAT_TIMERSHORT)
        
        ! Set the names of the template matrix and the solution vector
        rcollection%SquickAccess(1) = 'sparsitypattern'
        rcollection%SquickAccess(2) = 'solutionvector'
        
        ! Attach the primal solution vector to the collection structure
        call collct_setvalue_vec(rcollection, 'solutionvector', rsolution, .true.)
        
        ! Perform h-adaptation and update the triangulation structure
        call codire_adaptTriangulation(rhadapt, p_rproblemLevel%rtriangulation,&
                                       rvectorScalar, rcollection)
        
        ! Update the template matrix according to the sparsity pattern
        call grph_generateMatrix(rgraph, p_rproblemLevel%Rmatrix(templateMatrix))
        
        ! Resize the solution vector accordingly
        call lsysbl_resizeVectorBlock(rsolution,&
            p_rproblemLevel%Rmatrix(templateMatrix)%NEQ, .false.)
        
        ! Release element-wise error distribution
        call lsyssc_releaseVector(rvectorScalar)
        
        ! Stop time measurement for mesh adaptation
        call stat_stopTimer(rappDescriptor%rtimerAdaptation)

      
        !-------------------------------------------------------------------------
        ! Re-generate the discretization and coefficient matrices
        !-------------------------------------------------------------------------
        
        ! Start time measurement for generation of the triangulation
        call stat_startTimer(rappDescriptor%rtimerTriangulation, STAT_TIMERSHORT)
        
        ! Generate standard mesh from raw mesh
        call tria_initStandardMeshFromRaw(p_rproblemLevel%rtriangulation, rproblem%rboundary)
        
        ! Stop time measurement for generation of the triangulation
        call stat_stopTimer(rappDescriptor%rtimerTriangulation)
        
        
        ! Start time measurement for generation of constant coefficient matrices
        call stat_startTimer(rappDescriptor%rtimerAssemblyCoeff, STAT_TIMERSHORT)
        
        ! Re-initialize all constant coefficient matrices
        call codire_initProblemLevel(rappDescriptor, p_rproblemLevel, rcollection)
        
        ! Prepare internal data arrays of the solver structure
        call flagship_updateSolverMatrix(p_rproblemLevel, rsolver,&
                                         1, SYSTEM_INTERLEAVEFORMAT, UPDMAT_ALL)
        call solver_updateStructure(rsolver)
        
        ! Stop time measurement for generation of constant coefficient matrices
        call stat_stopTimer(rappDescriptor%rtimerAssemblyCoeff)

      end if
      
    end do timeloop

    
    ! Release adaptation structure
    if (dstepAdapt > 0.0_DP) then
      call hadapt_releaseAdaptation(rhadapt)
      call grph_releaseGraph(rgraph)
    end if

  end subroutine codire_solveTransientPrimal

  !*****************************************************************************

!<subroutine>

  subroutine codire_solvePseudoTransientPrimal(rappDescriptor, rparlist, ssectionName,&
                                               rbdrCond, rproblem, rtimestep, rsolver,&
                                               rsolution, rcollection)
!<description>
    ! This subroutine solves the pseudo-transient primal flow problem
    !
    ! $$\frac{\partial u}{\partial \tau}+\nabla\cdot{\bf f}(u)=s(u)$$
    !
    ! for a scalar quantity $u$ in the domain $\Omega$.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(IN) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(IN) :: ssectionName

    ! boundary condition structure
    type(t_boundaryCondition), intent(IN) :: rbdrCond
!</input>

!<inputoutput>
    ! application descriptor
    type(t_codire), intent(INOUT) :: rappDescriptor

    ! problem structure
    type(t_problem), intent(INOUT) :: rproblem

    ! time-stepping structure
    type(t_timestep), intent(INOUT) :: rtimestep

    ! solver struchture
    type(t_solver), intent(INOUT), target :: rsolver

    ! primal solution vector
    type(t_vectorBlock), intent(INOUT) :: rsolution

    ! collection structure
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! Pointer to the current multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! Vector for the right-hand side
    type(t_vectorBlock) :: rvectorBlock

    ! Vector for the element-wise error distribution
    type(t_vectorScalar) :: rvectorScalar

    ! Structure for h-adaptation
    type(t_hadapt) :: rhadapt

    ! Structure for the sparsity pattern
    type(t_graph) :: rgraph

    ! section names
    character(LEN=SYS_STRLEN) :: sadaptivityName

    ! local variables
    real(dp) :: derror
    integer :: templateMatrix
    integer :: nlmin, iadapt, nadapt


    ! Set pointer to maximum problem level
    p_rproblemLevel => rproblem%p_rproblemLevelMax

    !---------------------------------------------------------------------------
    ! Initialize the h-adaptation structure
    !---------------------------------------------------------------------------

    call parlst_getvalue_string(rparlist, ssectionName, 'adaptivity', sadaptivityName, '')
    if (trim(adjustl(sadaptivityName)) .ne. '') then

      call parlst_getvalue_int(rparlist, trim(sadaptivityName), 'nadapt', nadapt, 0)

      if (nadapt > 0) then

        ! Initialize adaptation structure from parameter list
        call hadapt_initFromParameterlist(rhadapt, rparlist, sadaptivityName)

        ! Generate a dynamic graph for the sparsity pattern and attach
        ! it to the collection structure which is used to pass this type
        ! to the callback function for h-adaptation
        templateMatrix = collct_getvalue_int(rcollection, 'templateMatrix')
        call grph_createGraphFromMatrix(p_rproblemLevel%Rmatrix(templateMatrix), rgraph)
        call collct_setvalue_graph(rcollection, 'sparsitypattern', rgraph, .true.)
        
        ! Attach the primal solution vector to the collection structure
        call collct_setvalue_vec(rcollection, 'solutionvector', rsolution, .true.)

      end if   ! nadapt > 0

    else   

      ! No h-adaptation
      nadapt = 0

    end if

    
    ! Adaptation loop
    adaptloop: do iadapt = 0, nadapt

      !-------------------------------------------------------------------------
      ! Compute the steady-state solution
      !-------------------------------------------------------------------------
      
      ! Start time measurement for solution procedure
      call stat_startTimer(rappDescriptor%rtimerSolution, STAT_TIMERSHORT)

      ! Calculate the velocity field
      nlmin = solver_getMinimumMultigridlevel(rsolver)
      call codire_calcVelocityField(rappDescriptor, p_rproblemLevel,&
                                    rtimestep%dtime, rcollection, nlmin)

      ! Attach the boundary condition to the solver structure
      call solver_setBoundaryCondition(rsolver, rbdrCond, .true.)

      ! Set collection to primal problem mode
      call collct_setvalue_int(rcollection, 'primaldual', 1, .true.)

      ! Reset the time-stepping algorithm
      call solver_resetTimestep(rtimestep, .false.)
      call solver_resetSolver(rsolver, .false.)

      ! Check if right-hand side vector exists
      if (rappDescriptor%irhstype > 0) then

        ! Initialize right-hand side vector
        call codire_initRHS(rappDescriptor, p_rproblemLevel,&
                            0.0_DP, rvectorBlock)

        ! Solve the primal problem with non-zero right-hand side
        call tstep_performPseudoStepping(p_rproblemLevel, rtimestep, rsolver,&
                                         rsolution, codire_calcRHS, codire_calcResidual,&
                                         codire_calcJacobian, codire_applyJacobian,&
                                         codire_setBoundary, rcollection, rvectorBlock)

        ! Release right-hand side vector
        call lsysbl_releaseVector(rvectorBlock)
        
      else
      
        ! Solve the primal problem with zero right-hand side
        call tstep_performPseudoStepping(p_rproblemLevel, rtimestep, rsolver,&
                                         rsolution, codire_calcRHS, codire_calcResidual,&
                                         codire_calcJacobian, codire_applyJacobian,&
                                         codire_setBoundary, rcollection)
  
      end if

      ! Stop time measurement for solution procedure
      call stat_stopTimer(rappDescriptor%rtimerSolution)
      
      
      !-------------------------------------------------------------------------
      ! Perform recovery-based error estimation
      !-------------------------------------------------------------------------

      ! Start time measurement for error estimation
      call stat_startTimer(rappDescriptor%rtimerErrorEstimation, STAT_TIMERSHORT)

      ! HERE WE NEED TO IMPLEMENT THE RECOVERY-BASED INDICATOR

      ! Stop time measurement for error estimation
      call stat_stopTimer(rappDescriptor%rtimerErrorEstimation)


      !-------------------------------------------------------------------------
      ! Perform h-adaptation
      !-------------------------------------------------------------------------

      ! Start time measurement for mesh adaptation
      call stat_startTimer(rappDescriptor%rtimerAdaptation, STAT_TIMERSHORT)
      
      ! Set the names of the template matrix and the solution vector
      rcollection%SquickAccess(1) = 'sparsitypattern'
      rcollection%SquickAccess(2) = 'solutionvector'
      
      ! Attach the primal solution vector to the collection structure
      call collct_setvalue_vec(rcollection, 'solutionvector', rsolution, .true.)
      
      ! Perform h-adaptation and update the triangulation structure
      call codire_adaptTriangulation(rhadapt, p_rproblemLevel%rtriangulation,&
                                     rvectorScalar, rcollection)

      ! Update the template matrix according to the sparsity pattern
      call grph_generateMatrix(rgraph, p_rproblemLevel%Rmatrix(templateMatrix))

      ! Resize the solution vector accordingly
      call lsysbl_resizeVectorBlock(rsolution,&
          p_rproblemLevel%Rmatrix(templateMatrix)%NEQ, .false.)
      
      ! Release element-wise error distribution
      call lsyssc_releaseVector(rvectorScalar)

      ! Stop time measurement for mesh adaptation
      call stat_stopTimer(rappDescriptor%rtimerAdaptation)

      
      !-------------------------------------------------------------------------
      ! Re-generate the discretization and coefficient matrices
      !-------------------------------------------------------------------------
      
      ! Start time measurement for generation of the triangulation
      call stat_startTimer(rappDescriptor%rtimerTriangulation, STAT_TIMERSHORT)

      ! Generate standard mesh from raw mesh
      call tria_initStandardMeshFromRaw(p_rproblemLevel%rtriangulation, rproblem%rboundary)

      ! Stop time measurement for generation of the triangulation
      call stat_stopTimer(rappDescriptor%rtimerTriangulation)

      
      ! Start time measurement for generation of constant coefficient matrices
      call stat_startTimer(rappDescriptor%rtimerAssemblyCoeff, STAT_TIMERSHORT)

      ! Re-initialize all constant coefficient matrices
      call codire_initProblemLevel(rappDescriptor, p_rproblemLevel, rcollection)

      ! Prepare internal data arrays of the solver structure
      call flagship_updateSolverMatrix(p_rproblemLevel, rsolver,&
                                       1, SYSTEM_INTERLEAVEFORMAT, UPDMAT_ALL)
      call solver_updateStructure(rsolver)
      
      ! Stop time measurement for generation of constant coefficient matrices
      call stat_stopTimer(rappDescriptor%rtimerAssemblyCoeff)
      
    end do adaptloop
    

    ! Release adaptation structure
    if (nadapt > 0) then
      call hadapt_releaseAdaptation(rhadapt)
      call grph_releaseGraph(rgraph)
    end if    
    
  end subroutine codire_solvePseudoTransientPrimal

  !*****************************************************************************

!<subroutine>

  subroutine codire_solveSteadyStatePrimal(rappDescriptor, rparlist, ssectionName,&
                                           rbdrCond, rproblem, rtimestep, rsolver,&
                                           rsolution, rcollection)

!<description>
    ! This subroutine solves the steady-state primal flow problem
    !
    ! $$\nabla\cdot{\bf f}(u)=s(u)$$
    !
    ! for a scalar quantity $u$ in the domain $\Omega$.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(IN) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(IN) :: ssectionName

    ! boundary condition structure
    type(t_boundaryCondition), intent(IN) :: rbdrCond
!</input>

!<inputoutput>
    ! application descriptor
    type(t_codire), intent(INOUT) :: rappDescriptor

    ! problem structure
    type(t_problem), intent(INOUT) :: rproblem

    ! time-stepping structure
    type(t_timestep), intent(INOUT) :: rtimestep

    ! solver struchture
    type(t_solver), intent(INOUT), target :: rsolver

    ! primal solution vector
    type(t_vectorBlock), intent(INOUT) :: rsolution

    ! collection structure
    type(t_collection), intent(INOUT) :: rcollection    
!</inputoutput>
!</subroutine>

    ! Pointer to the current multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! Vector for the linear target functional and the right-hand side
    type(t_vectorBlock) :: rvectorBlock

    ! Vector for the element-wise error distribution
    type(t_vectorScalar) :: rvectorScalar
    
    ! Structure for h-adaptation
    type(t_hadapt) :: rhadapt

    ! Structure for the sparsity pattern
    type(t_graph) :: rgraph

    ! section names
    character(LEN=SYS_STRLEN) :: sadaptivityName

    ! local variables
    real(dp) :: derror
    integer :: templateMatrix
    integer :: nlmin, iadapt, nadapt

    
    ! Adjust time stepping scheme
    rtimestep%ctimestepType = SV_THETA_SCHEME
    rtimestep%dinitialTime  = 0.0_DP
    rtimestep%dinitialStep  = 1.0_DP
    rtimestep%dfinalTime    = 1.0_DP
    rtimestep%theta         = 1.0_DP

    ! Set pointer to maximum problem level
    p_rproblemLevel => rproblem%p_rproblemLevelMax

    
    !---------------------------------------------------------------------------
    ! Initialize the h-adaptation structure
    !---------------------------------------------------------------------------

    call parlst_getvalue_string(rparlist, ssectionName, 'adaptivity', sadaptivityName, '')
    if (trim(adjustl(sadaptivityName)) .ne. '') then

      call parlst_getvalue_int(rparlist, trim(sadaptivityName), 'nadapt', nadapt, 0)

      if (nadapt > 0) then

        ! Initialize adaptation structure from parameter list
        call hadapt_initFromParameterlist(rhadapt, rparlist, sadaptivityName)
        
        ! Generate a dynamic graph for the sparsity pattern and attach
        ! it to the collection structure which is used to pass this type
        ! to the callback function for h-adaptation
        templateMatrix = collct_getvalue_int(rcollection, 'templateMatrix')
        call grph_createGraphFromMatrix(p_rproblemLevel%Rmatrix(templateMatrix), rgraph)
        call collct_setvalue_graph(rcollection, 'sparsitypattern', rgraph, .true.)
        
        ! Attach the primal solution vector to the collection structure
        call collct_setvalue_vec(rcollection, 'solutionvector', rsolution, .true.)

      end if   ! nadapt > 0

    else   

      ! No h-adaptation
      nadapt = 0

    end if


    adaptloop: do iadapt = 0, nadapt

      !-------------------------------------------------------------------------
      ! Compute steady-state solution for the primal problem
      !-------------------------------------------------------------------------
      
      ! Start time measurement for solution procedure
      call stat_startTimer(rappDescriptor%rtimerSolution, STAT_TIMERSHORT)
      
      ! Calculate the velocity field
      nlmin = solver_getMinimumMultigridlevel(rsolver)
      call codire_calcVelocityField(rappDescriptor, p_rproblemLevel,&
                                    rtimestep%dtime, rcollection, nlmin)
      
      ! Attach the boundary condition to the solver structure
      call solver_setBoundaryCondition(rsolver, rbdrCond, .true.)

      ! Set collection to primal problem mode
      call collct_setvalue_int(rcollection, 'primaldual', 1, .true.)
      
      ! Reset the time-stepping algorithm
      call solver_resetTimestep(rtimestep, .false.)
      call solver_resetSolver(rsolver, .false.)

      ! Check if right-hand side vector exists
      if (rappDescriptor%irhstype > 0) then

        ! Initialize right-hand side vector
        call codire_initRHS(rappDescriptor, p_rproblemLevel,&
                            0.0_DP, rvectorBlock)

        ! Solve the primal problem with non-zero right-hand side
        call tstep_performThetaStep(p_rproblemLevel, rtimestep, rsolver,&
                                    rsolution, codire_calcResidual,&
                                    codire_calcJacobian, codire_applyJacobian,&
                                    codire_setBoundary, rcollection, rvectorBlock)
        
        ! Release right-hand side vector
        call lsysbl_releaseVector(rvectorBlock)

      else

        ! Solve the primal problem with zero right-hand side
        call tstep_performThetaStep(p_rproblemLevel, rtimestep, rsolver,&
                                    rsolution, codire_calcResidual,&
                                    codire_calcJacobian, codire_applyJacobian,&
                                    codire_setBoundary, rcollection)

      end if

      ! Stop time measurement for solution procedure
      call stat_stopTimer(rappDescriptor%rtimerSolution)


      if (iadapt .eq. nadapt) exit adaptloop

      !-------------------------------------------------------------------------
      ! Perform recovery-based error estimation
      !-------------------------------------------------------------------------

      ! Start time measurement for error estimation
      call stat_startTimer(rappDescriptor%rtimerErrorEstimation, STAT_TIMERSHORT)

      ! HERE WE NEED TO IMPLEMENT THE RECOVERY-BASED INDICATOR

      ! Stop time measurement for error estimation
      call stat_stopTimer(rappDescriptor%rtimerErrorEstimation)


      !-------------------------------------------------------------------------
      ! Perform h-adaptation
      !-------------------------------------------------------------------------

      ! Start time measurement for mesh adaptation
      call stat_startTimer(rappDescriptor%rtimerAdaptation, STAT_TIMERSHORT)
      
      ! Set the names of the template matrix and the solution vector
      rcollection%SquickAccess(1) = 'sparsitypattern'
      rcollection%SquickAccess(2) = 'solutionvector'
      
      ! Attach the primal solution vector to the collection structure
      call collct_setvalue_vec(rcollection, 'solutionvector', rsolution, .true.)
      
      ! Perform h-adaptation and update the triangulation structure
      call codire_adaptTriangulation(rhadapt, p_rproblemLevel%rtriangulation,&
                                     rvectorScalar, rcollection)

      ! Update the template matrix according to the sparsity pattern
      call grph_generateMatrix(rgraph, p_rproblemLevel%Rmatrix(templateMatrix))

      ! Resize the solution vector accordingly
      call lsysbl_resizeVectorBlock(rsolution,&
          p_rproblemLevel%Rmatrix(templateMatrix)%NEQ, .false.)
      
      ! Release element-wise error distribution
      call lsyssc_releaseVector(rvectorScalar)

      ! Stop time measurement for mesh adaptation
      call stat_stopTimer(rappDescriptor%rtimerAdaptation)

      
      !-------------------------------------------------------------------------
      ! Re-generate the discretization and coefficient matrices
      !-------------------------------------------------------------------------
      
      ! Start time measurement for generation of the triangulation
      call stat_startTimer(rappDescriptor%rtimerTriangulation, STAT_TIMERSHORT)

      ! Generate standard mesh from raw mesh
      call tria_initStandardMeshFromRaw(p_rproblemLevel%rtriangulation, rproblem%rboundary)

      ! Stop time measurement for generation of the triangulation
      call stat_stopTimer(rappDescriptor%rtimerTriangulation)

      
      ! Start time measurement for generation of constant coefficient matrices
      call stat_startTimer(rappDescriptor%rtimerAssemblyCoeff, STAT_TIMERSHORT)

      ! Re-initialize all constant coefficient matrices
      call codire_initProblemLevel(rappDescriptor, p_rproblemLevel, rcollection)

      ! Prepare internal data arrays of the solver structure
      call flagship_updateSolverMatrix(p_rproblemLevel, rsolver,&
                                       1, SYSTEM_INTERLEAVEFORMAT, UPDMAT_ALL)
      call solver_updateStructure(rsolver)
      
      ! Stop time measurement for generation of constant coefficient matrices
      call stat_stopTimer(rappDescriptor%rtimerAssemblyCoeff)

    end do adaptloop


    ! Release adaptation structure
    if (nadapt > 0) then
      call hadapt_releaseAdaptation(rhadapt)
      call grph_releaseGraph(rgraph)
    end if

  end subroutine codire_solveSteadyStatePrimal

  !*****************************************************************************

!<subroutine>

  subroutine codire_solveSteadyStatePrimalDual(rappDescriptor, rparlist, ssectionName,&
                                               rbdrCondPrimal, rbdrCondDual, rproblem,&
                                               rtimestep, rsolver, rsolutionPrimal,&
                                               rsolutionDual, rcollection)

!<description>
    ! This subroutine solves the steady-state primal flow problem
    !
    ! $$\nabla\cdot{\bf f}(u)=s(u)$$
    !
    ! for a scalar quantity $u$ in the domain $\Omega$ both for the
    ! primal and the dual formulation and performs goal-oriented
    ! mesh adaptation if the optional parameter rhadapt is present.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(IN) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(IN) :: ssectionName

    ! boundary condition structure for the primal problem
    type(t_boundaryCondition), intent(IN) :: rbdrCondPrimal

    ! boundary condition structure for the dual problem
    type(t_boundaryCondition), intent(IN) :: rbdrCondDual
!</input>

!<inputoutput>
    ! application descriptor
    type(t_codire), intent(INOUT) :: rappDescriptor

    ! problem structure
    type(t_problem), intent(INOUT) :: rproblem

    ! time-stepping structure
    type(t_timestep), intent(INOUT) :: rtimestep

    ! solver struchture
    type(t_solver), intent(INOUT), target :: rsolver

    ! primal solution vector
    type(t_vectorBlock), intent(INOUT) :: rsolutionPrimal

    ! dual solution vector
    type(t_vectorBlock), intent(INOUT) :: rsolutionDual

    ! collection structure
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! Pointer to the current multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! Vector for the linear target functional and the right-hand side
    type(t_vectorBlock) :: rvectorBlock

    ! Vector for the element-wise error distribution
    type(t_vectorScalar) :: rvectorScalar
    
    ! Structure for h-adaptation
    type(t_hadapt) :: rhadapt

    ! Structure for the sparsity pattern
    type(t_graph) :: rgraph

    ! section names
    character(LEN=SYS_STRLEN) :: sadaptivityName

    ! local variables
    real(dp) :: derror
    integer :: templateMatrix
    integer :: nlmin, iadapt, nadapt


    ! Adjust time stepping scheme
    rtimestep%ctimestepType = SV_THETA_SCHEME
    rtimestep%dinitialTime  = 0.0_DP
    rtimestep%dinitialStep  = 1.0_DP
    rtimestep%dfinalTime    = 1.0_DP
    rtimestep%theta         = 1.0_DP

    ! Set pointer to maximum problem level
    p_rproblemLevel => rproblem%p_rproblemLevelMax


    !---------------------------------------------------------------------------
    ! Initialize the h-adaptation structure
    !---------------------------------------------------------------------------

    call parlst_getvalue_string(rparlist, ssectionName, 'adaptivity', sadaptivityName, '')
    if (trim(adjustl(sadaptivityName)) .ne. '') then

      call parlst_getvalue_int(rparlist, trim(sadaptivityName), 'nadapt', nadapt, 0)

      if (nadapt > 0) then
        
        ! Initialize adaptation structure from parameter list
        call hadapt_initFromParameterlist(rhadapt, rparlist, sadaptivityName)

        ! Generate a dynamic graph for the sparsity pattern and attach
        ! it to the collection structure which is used to pass this type
        ! to the callback function for h-adaptation
        templateMatrix = collct_getvalue_int(rcollection, 'templateMatrix')
        call grph_createGraphFromMatrix(p_rproblemLevel%Rmatrix(templateMatrix), rgraph)
        call collct_setvalue_graph(rcollection, 'sparsitypattern', rgraph, .true.)
        
        ! Attach the primal solution vector to the collection structure
        call collct_setvalue_vec(rcollection, 'solutionvector', rsolutionPrimal, .true.)

      end if   ! nadapt > 0

    else   

      ! No h-adaptation
      nadapt = 0

    end if

    
    adaptloop: do iadapt = 0, nadapt

      !-------------------------------------------------------------------------
      ! Compute steady-state solution for the primal problem
      !-------------------------------------------------------------------------
      
      ! Start time measurement for solution procedure
      call stat_startTimer(rappDescriptor%rtimerSolution, STAT_TIMERSHORT)
      
      ! Calculate the velocity field
      nlmin = solver_getMinimumMultigridlevel(rsolver)
      call codire_calcVelocityField(rappDescriptor, p_rproblemLevel,&
                                    rtimestep%dtime, rcollection, nlmin)
      
      ! Attach the boundary condition to the solver structure
      call solver_setBoundaryCondition(rsolver, rbdrCondPrimal, .true.)

      ! Set collection to primal problem mode
      call collct_setvalue_int(rcollection, 'primaldual', 1, .true.)
      
      ! Reset the time-stepping algorithm
      call solver_resetTimestep(rtimestep, .false.)
      call solver_resetSolver(rsolver, .false.)

      ! Check if right-hand side vector exists
      if (rappDescriptor%irhstype > 0) then

        ! Initialize right-hand side vector
        call codire_initRHS(rappDescriptor, p_rproblemLevel,&
                            0.0_DP, rvectorBlock)

        ! Solve the primal problem with non-zero right-hand side
        call tstep_performThetaStep(p_rproblemLevel, rtimestep, rsolver,&
                                    rsolutionPrimal, codire_calcResidual,&
                                    codire_calcJacobian, codire_applyJacobian,&
                                    codire_setBoundary, rcollection, rvectorBlock)
        
        ! Release right-hand side vector
        call lsysbl_releaseVector(rvectorBlock)

      else

        ! Solve the primal problem with zero right-hand side
        call tstep_performThetaStep(p_rproblemLevel, rtimestep, rsolver,&
                                    rsolutionPrimal, codire_calcResidual,&
                                    codire_calcJacobian, codire_applyJacobian,&
                                    codire_setBoundary, rcollection)

      end if

      ! Stop time measurement for solution procedure
      call stat_stopTimer(rappDescriptor%rtimerSolution)
      
      
      !-------------------------------------------------------------------------
      ! Compute the right-hand side for the dual problem
      !-------------------------------------------------------------------------

      ! Start time measurement for error estimation
      call stat_startTimer(rappDescriptor%rtimerErrorEstimation, STAT_TIMERSHORT)

      ! Initialize target functional
      call codire_initTargetFunc(rappDescriptor, p_rproblemLevel,&
                                 1.0_DP, rvectorBlock)
      
      ! Stop time measurement for error estimation
      call stat_stopTimer(rappDescriptor%rtimerErrorEstimation)


      !-------------------------------------------------------------------------
      ! Compute steady-state solution for the dual problem
      !-------------------------------------------------------------------------

      ! Start time measurement for solution procedure
      call stat_startTimer(rappDescriptor%rtimerSolution, STAT_TIMERSHORT)
      
      ! Calculate the velocity field
      nlmin = solver_getMinimumMultigridlevel(rsolver)
      call codire_calcVelocityField(rappDescriptor, p_rproblemLevel,&
                                    rtimestep%dtime, rcollection, nlmin)
      
      ! Attach the boundary condition to the solver structure
      call solver_setBoundaryCondition(rsolver, rbdrCondDual, .true.)
      
      ! Set collection to primal problem mode
      call collct_setvalue_int(rcollection, 'primaldual', 2, .true.)
      
      ! Create dual solution vector initialized by zeros
      call lsysbl_releaseVector(rsolutionDual)
      call lsysbl_createVectorBlock(rsolutionPrimal, rsolutionDual, .true.)
      
      ! Reset the time-stepping and solver algorithms
      call solver_resetTimestep(rtimestep, .false.)
      call solver_resetSolver(rsolver, .false.)
      
      ! Solve the dual problem
      call tstep_performThetaStep(p_rproblemLevel, rtimestep, rsolver, rsolutionDual,&
                                  codire_calcResidual, codire_calcJacobian,&
                                  codire_applyJacobian, codire_setBoundary,&
                                  rcollection, rvectorBlock)
    
      ! Release discretized target functional
      call lsysbl_releaseVector(rvectorBlock)

      ! Stop time measurement for solution procedure
      call stat_stopTimer(rappDescriptor%rtimerSolution)

      
      if (iadapt .eq. nadapt) exit adaptloop

      !-------------------------------------------------------------------------
      ! Perform goal-oriented error estimation
      !-------------------------------------------------------------------------

      ! Start time measurement for error estimation
      call stat_startTimer(rappDescriptor%rtimerErrorEstimation, STAT_TIMERSHORT)
      
      ! Calculate the velocity field
      nlmin = solver_getMinimumMultigridlevel(rsolver)
      call codire_calcVelocityField(rappDescriptor, p_rproblemLevel,&
                                    rtimestep%dtime, rcollection, nlmin)
      
      ! Attach the boundary condition to the solver structure
      call solver_setBoundaryCondition(rsolver, rbdrCondPrimal, .true.)
      
      ! Set collection to primal problem mode
      call collct_setvalue_int(rcollection, 'primaldual', 1, .true.)

      ! Check if right-hand side vector exists
      if (rappDescriptor%irhstype > 0) then

        ! Initialize right-hand side vector
        call codire_initRHS(rappDescriptor, p_rproblemLevel,&
                            0.0_DP, rvectorBlock)

        ! Compute the error in the quantity of interest
        call codire_estimateTargetFuncError(p_rproblemLevel, rtimestep, rsolver,&
                                            rsolutionPrimal, rsolutionDual,&
                                            rcollection, rvectorScalar, derror,&
                                            rvectorBlock)

        ! Release right-hand side vector
        call lsysbl_releaseVector(rvectorBlock)
        
      else

        ! Compute the error in the quantity of interest
        call codire_estimateTargetFuncError(p_rproblemLevel, rtimestep, rsolver,&
                                            rsolutionPrimal, rsolutionDual,&
                                            rcollection, rvectorScalar, derror)

      end if

      ! Stop time measurement for error estimation
      call stat_stopTimer(rappDescriptor%rtimerErrorEstimation)


      !-------------------------------------------------------------------------
      ! Perform h-adaptation
      !-------------------------------------------------------------------------

      ! Start time measurement for mesh adaptation
      call stat_startTimer(rappDescriptor%rtimerAdaptation, STAT_TIMERSHORT)
      
      ! Set the names of the template matrix and the solution vector
      rcollection%SquickAccess(1) = 'sparsitypattern'
      rcollection%SquickAccess(2) = 'solutionvector'
      
      ! Attach the primal solution vector to the collection structure
      call collct_setvalue_vec(rcollection, 'solutionvector', rsolutionPrimal, .true.)
      
      ! Perform h-adaptation and update the triangulation structure
      call codire_adaptTriangulation(rhadapt, p_rproblemLevel%rtriangulation,&
                                     rvectorScalar, rcollection)

      ! Update the template matrix according to the sparsity pattern
      call grph_generateMatrix(rgraph, p_rproblemLevel%Rmatrix(templateMatrix))

      ! Resize the solution vector accordingly
      call lsysbl_resizeVectorBlock(rsolutionPrimal,&
          p_rproblemLevel%Rmatrix(templateMatrix)%NEQ, .false.)
      
      ! Release element-wise error distribution
      call lsyssc_releaseVector(rvectorScalar)

      ! Stop time measurement for mesh adaptation
      call stat_stopTimer(rappDescriptor%rtimerAdaptation)

      
      !-------------------------------------------------------------------------
      ! Re-generate the discretization and coefficient matrices
      !-------------------------------------------------------------------------
      
      ! Start time measurement for generation of the triangulation
      call stat_startTimer(rappDescriptor%rtimerTriangulation, STAT_TIMERSHORT)

      ! Generate standard mesh from raw mesh
      call tria_initStandardMeshFromRaw(p_rproblemLevel%rtriangulation, rproblem%rboundary)

      ! Stop time measurement for generation of the triangulation
      call stat_stopTimer(rappDescriptor%rtimerTriangulation)

      
      ! Start time measurement for generation of constant coefficient matrices
      call stat_startTimer(rappDescriptor%rtimerAssemblyCoeff, STAT_TIMERSHORT)

      ! Re-initialize all constant coefficient matrices
      call codire_initProblemLevel(rappDescriptor, p_rproblemLevel, rcollection)

      ! Prepare internal data arrays of the solver structure
      call flagship_updateSolverMatrix(p_rproblemLevel, rsolver,&
                                       1, SYSTEM_INTERLEAVEFORMAT, UPDMAT_ALL)
      call solver_updateStructure(rsolver)
      
      ! Stop time measurement for generation of constant coefficient matrices
      call stat_stopTimer(rappDescriptor%rtimerAssemblyCoeff)

    end do adaptloop


    ! Release adaptation structure
    if (nadapt > 0) then
      call hadapt_releaseAdaptation(rhadapt)
      call grph_releaseGraph(rgraph)
    end if

  end subroutine codire_solveSteadyStatePrimalDual

  !*****************************************************************************
  ! AUXILIARY ROUTINES
  !*****************************************************************************

!<subroutine>

  subroutine codire_parseCmdlArguments(rparlist)

!<description>
    ! This subroutine parses the commandline arguments and modifies the
    ! parameter values in the global parameter list.
!</description>

!<inputoutput>
    ! parameter list
    type(t_parlist), intent(INOUT) :: rparlist
!</inputoutput>
!</subroutine>

    ! local variables
    character(LEN=SYS_STRLEN) :: cbuffer
    integer :: iarg, narg
    
    iarg = 1; narg = command_argument_count()
    
    cmdarg: do
      ! Retrieve next command line argument
      call get_command_argument(iarg,cbuffer)
      select case(trim(adjustl(cbuffer)))
        
      case ('-A','--adaptivity')
        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        call parlst_setvalue(rparlist, '', 'adaptivity', trim(adjustl(cbuffer)))

      case ('-B','--benchmark')
        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        call parlst_setvalue(rparlist, '', 'benchmark', trim(adjustl(cbuffer)))
       
      case ('-DC','--dualconv')
        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        call parlst_setvalue(rparlist, '', 'dualconv', trim(adjustl(cbuffer)))
        
      case ('-DD','--dualdiff')
        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        call parlst_setvalue(rparlist, '', 'dualdiff', trim(adjustl(cbuffer)))

      case ('-E','--errorest')
        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        call parlst_setvalue(rparlist, '', 'errorest', trim(adjustl(cbuffer)))
        
      case ('-I','--io')
        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        call parlst_setvalue(rparlist, '', 'inputoutput', trim(adjustl(cbuffer)))
        
      case ('-PC','--primalconv')
        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        call parlst_setvalue(rparlist, '', 'primalconv', trim(adjustl(cbuffer)))
        
      case ('-PD','--primaldiff')
        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        call parlst_setvalue(rparlist, '', 'primaldiff', trim(adjustl(cbuffer)))

      case ('-S','--solver')
        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        call parlst_setvalue(rparlist, '', 'solver', trim(adjustl(cbuffer)))
        
      case ('-T','--timestep')
        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        call parlst_setvalue(rparlist, '', 'timestep', trim(adjustl(cbuffer)))
        
      case DEFAULT
        iarg = iarg+1
        if (iarg .ge. narg) exit cmdarg
      end select
    end do cmdarg

  end subroutine codire_parseCmdlArguments

  !*****************************************************************************

!<subroutine>

  subroutine codire_UserInterface()

!<description>
    ! This subroutine enables the user to interact with the simulation.
!</description>

!</subroutine>
    
!!$    ! local variables
!!$    integer, external :: signal_SIGINT
!!$    
!!$    ! Perform intermediate output
!!$    if (signal_SIGINT(-1) > 0 ) call codire_postprocess(rsolution)
    
  end subroutine codire_UserInterface

end module codire_application
