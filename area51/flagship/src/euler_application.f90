!##############################################################################
!# ****************************************************************************
!# <name> euler_application </name>
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
!# and the tuple/triple of fluxes ${\bf F}=(F^1,F^2,F^3)$
!# for each coordinate direction
!#
!#   $${\bf F}=[\rho{\bf v},
!#              \rho{\bf v}\otimes{\bf v}+o{\mathcal I},
!#              \rho H{\bf v}]^T$$
!#
!# in the two- or three-dimensional domain $\Omega$.
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
!# 1.) euler
!#     -> The main routine of the application called from the main problem
!#
!#
!# The following auxiliary routines are available:
!#
!# 1.) euler_parseCmdlArguments
!#     -> Parses the list of commandline arguments and overwrites
!#        parameter values from the parameter files
!#
!# 2.) euler_initApplication
!#     -> Initializes the application descriptor by the parameter
!#        settings given by the parameter list
!#
!# 3.) euler_initCollection
!#     -> Initializes the collection structure based on the parameter
!#        settings given by the application descriptor
!#
!# 4.) euler_initSolvers
!#     -> Initializes the solve structures from the parameter list
!#
!# 5.) euler_initProblem
!#     -> Initializes the global problem structure based on the
!#        parameter settings given by the parameter list
!#
!# 6.) euler_initProblemLevel
!#     -> Initializes the individual problem levels based in the
!#        parameter settings given by the application descriptor
!#
!# 7.) euler_initAllProblemLevels
!#     -> Initializes all problem levels attached to the global
!#        problem structure based on the parameter settings
!#        given by the parameter list
!#
!# 8.) euler_initSolution
!#     -> Initializes the solution vector based on the parameter
!#        settings given by the parameter list
!#
!# 9.) euler_outputSolution
!#     -> Outputs the solution vector to file in UCD format
!#
!# 10.) euler_outputStatistics
!#      -> Outputs the application statitics
!#
!# 11.) euler_estimateRecoveryError
!#      -> Estimates the solution error using recovery techniques
!#
!# 12.) euler_adaptTriangulation
!#      -> Performs h-adaptation for the given triangulation
!#
!# 13.) euler_solveTransientPrimal
!#      -> Solves the primal formulation of the time-dependent 
!#         compressible Euler equations.
!# 
!# </purpose>
!##############################################################################

module euler_application

  use afcstabilisation
  use bilinearformevaluation
  use boundaryfilter
  use collection
  use euler_basic
  use euler_callback
  use euler_callback1d
  use euler_callback2d
  use euler_callback3d
  use flagship_basic
  use fparser
  use fsystem
  use genoutput
  use graph
  use groupfemsystem
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
  use solveraux
  use spatialdiscretisation
  use statistics
  use stdoperators
  use storage
  use thermodynamics
  use timestep
  use timestepaux
  use ucd

  implicit none

  private
  public :: euler
  public :: euler_initApplication
  public :: euler_initCollection
  public :: euler_initSolvers
  public :: euler_initProblem
  public :: euler_initProblemLevel
  public :: euler_initAllProblemLevels
  public :: euler_initSolution
  public :: euler_outputSolution
  public :: euler_outputStatistics
  public :: euler_estimateRecoveryError
  public :: euler_adaptTriangulation
  public :: euler_solveTransientPrimal
  

contains

  !*****************************************************************************

!<subroutine>

  subroutine euler(rparlist)

!<description>
    ! This is the main application for the compressible Euler
    ! equations.  It is a so-called driver routine which can be used
    ! to start a standalone Euler simulation.
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
    ! for the compressible Euler flow application
    type(t_euler) :: rappDescriptor

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

    ! Parameter file and section names
    character(LEN=SYS_STRLEN) :: sindatfileName
    character(LEN=SYS_STRLEN) :: sbdrcondName
    character(LEN=SYS_STRLEN) :: algorithm

    ! local variables
    integer :: systemMatrix
    

    ! Start total time measurement
    call stat_clearTimer(rtimerTotal)
    call stat_startTimer(rtimerTotal)
    
    !---------------------------------------------------------------------------
    ! Pre-processing
    !---------------------------------------------------------------------------

    ! Overwrite global configuration from command line arguments. After
    ! this subroutine has been called, the parameter list remains unchanged.
    call euler_parseCmdlArguments(rparlist)

    ! Initialize global collection structure
    call collct_init(rcollection)

    ! Initialize the application descriptor
    call euler_initApplication(rparlist, 'euler', rappDescriptor)

    ! Start time measurement for pre-processing
    call stat_startTimer(rappDescriptor%rtimerPrepostProcess, STAT_TIMERSHORT)

    ! Initialize the global collection
    call euler_initCollection(rappDescriptor, rparlist, 'euler', rcollection)

    ! Initialize the solver structures
    call euler_initSolvers(rparlist, 'euler', rtimestep, rsolver)

    ! Initialize the abstract problem structure
    call euler_initProblem(rparlist, 'euler',&
                           solver_getMinimumMultigridlevel(rsolver),&
                           solver_getMaximumMultigridlevel(rsolver),&
                           rproblem, rcollection)

    ! Initialize the individual problem levels
    call euler_initAllProblemLevels(rappDescriptor, rproblem, rcollection)

    ! Prepare internal data arrays of the solver structure
    systemMatrix = collct_getvalue_int(rcollection, 'systemMatrix') 
    call flagship_updateSolverMatrix(rproblem%p_rproblemLevelMax, rsolver,&
                                     systemMatrix, rappDescriptor%isystemFormat, UPDMAT_ALL)
    call solver_updateStructure(rsolver)

    ! Stop time measurement for pre-processing
    call stat_stopTimer(rappDescriptor%rtimerPrePostprocess)

    
    !---------------------------------------------------------------------------
    ! Solution algorithm
    !---------------------------------------------------------------------------

    if (rtimestep%dfinalTime > 0) then
      
      ! Get global configuration from parameter list
      call parlst_getvalue_string(rparlist, 'euler', 'algorithm', algorithm)
      call parlst_getvalue_string(rparlist, 'euler', 'indatfile', sindatfileName)
      
      ! The boundary condition for the primal problem is required for all 
      ! solution strategies so initialize it from the parameter file
      call parlst_getvalue_string(rparlist, 'euler', 'sprimalbdrcondname', sbdrcondName)
      call bdrf_readBoundaryCondition(rbdrCondPrimal, sindatfileName,&
                                      '['//trim(sbdrcondName)//']', rappDescriptor%ndimension)
      
      ! What solution algorithm should be applied?
      select case(trim(algorithm))

      case ('transient_primal')
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Solve the primal formulation for the time-dependent problem
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call euler_solveTransientPrimal(rappDescriptor, rparlist, 'euler',&
                                        rbdrCondPrimal, rproblem, rtimestep,&
                                        rsolver, rsolutionPrimal, rcollection)
        call euler_outputSolution(rparlist, 'euler', rproblem%p_rproblemLevelMax,&
                                  rsolutionPrimal, dtime=rtimestep%dTime)

        
      case DEFAULT
        call output_line(trim(algorithm)//' is not a valid solution algorithm!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'euler')
        call sys_halt()
      end select
    end if

    !---------------------------------------------------------------------------
    ! Post-processing
    !---------------------------------------------------------------------------

    ! Start time measurement for pre-processing
    call stat_startTimer(rappDescriptor%rtimerPrepostProcess, STAT_TIMERSHORT)
    
    ! Release time-stepping
    call tstep_releaseTimestep(rtimestep)
    
    ! Release solver
    call solver_releaseSolver(rsolver)
    
    ! Release problem structure
    call problem_releaseProblem(rproblem)

    ! Release boundary conditions
    call bdrf_release(rbdrCondPrimal)
    call bdrf_release(rbdrCondDual)

    ! Release vectors
    call lsysbl_releaseVector(rsolutionPrimal)
    call lsysbl_releaseVector(rsolutionDual)

    ! Release collection
    call collct_done(rcollection)

    ! Stop time measurement for pre-processing
    call stat_stopTimer(rappDescriptor%rtimerPrePostprocess)

    ! Stop time measurement for total time measurement
    call stat_stopTimer(rtimerTotal)

    ! Output statistics
    call euler_outputStatistics(rappDescriptor, rtimerTotal)

  end subroutine euler

  !*****************************************************************************

!<subroutine>

  subroutine euler_initApplication(rparlist, ssectionName, rappDescriptor)

!<description>
    ! This subroutine initializes the application descriptor by the
    ! parameter settings given by the parameter list
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(IN) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(IN) :: ssectionName
!</input>

!<output>
    ! application descriptor
    type(t_euler), intent(OUT) :: rappDescriptor
!</output>
!</subroutine>

    ! section names
    character(LEN=SYS_STRLEN) :: sindatfileName

    
    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist, ssectionName, 'indatfile', sindatfileName)

    ! Read in all constants and predefined expressions from the parameter file
    call fparser_parseFileForKeyword(sindatfileName, 'defconst', FPAR_CONSTANT)
    call fparser_parseFileForKeyword(sindatfileName, 'defexpr',  FPAR_EXPRESSION)

    ! Get application specifig parameters from the parameterlist
    call parlst_getvalue_int(rparlist, ssectionName,&
                             'ndimension', rappDescriptor%ndimension)
    call parlst_getvalue_int(rparlist, ssectionName,& 
                             'imasstype', rappDescriptor%imasstype)
    call parlst_getvalue_int(rparlist, ssectionName,& 
                             'imassantidiffusiontype', rappDescriptor%imassantidiffusiontype)
    call parlst_getvalue_int(rparlist, ssectionName,&
                             'idissipationtype', rappDescriptor%idissipationtype)
    call parlst_getvalue_int(rparlist, ssectionName,&
                             'isystemCoupling', rappDescriptor%isystemCoupling)
    call parlst_getvalue_int(rparlist, ssectionName,&
                             'isystemPrecond', rappDescriptor%isystemPrecond)
    call parlst_getvalue_int(rparlist, ssectionName,&
                             'celement', rappDescriptor%celement)
    call parlst_getvalue_int(rparlist, ssectionName,&
                             'imatrixformat', rappDescriptor%imatrixFormat)
    call parlst_getvalue_int(rparlist, ssectionName,&
                             'isystemformat', rappDescriptor%isystemFormat)

  end subroutine euler_initApplication

  !*****************************************************************************

!<subroutine>

  subroutine euler_initCollection(rappDescriptor, rparlist, ssectionName, rcollection)

!<description>
    ! This subroutine initializes the collection based on
    ! the parameter settings of the application descriptor
!</description>

!<input>
    ! application descriptor
    type(t_euler), intent(IN) :: rappDescriptor
    
    ! parameter list
    type(t_parlist), intent(IN) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(IN) :: ssectionName
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
    call collct_setvalue_int(rcollection, 'imassantidiffusiontype',&
                             rappDescriptor%imassantidiffusiontype, .true.)
    call collct_setvalue_int(rcollection, 'idissipationtype',&
                             rappDescriptor%idissipationtype, .true.)
    call collct_setvalue_int(rcollection, 'isystemCoupling',&
                             rappDescriptor%isystemCoupling, .true.)
    call collct_setvalue_int(rcollection, 'isystemPrecond',&
                             rappDescriptor%isystemPrecond, .true.)
    call collct_setvalue_int(rcollection, 'celement',&
                             rappDescriptor%celement, .true.)
    call collct_setvalue_int(rcollection, 'imatrixformat',&
                             rappDescriptor%imatrixFormat, .true.)
    call collct_setvalue_int(rcollection, 'isystemformat',&
                             rappDescriptor%isystemFormat, .true.)

    ! Add timer structures
    call collct_setvalue_timer(rcollection, 'timerSolution',&
                               rappDescriptor%rtimerSolution, .true.)
    call collct_setvalue_timer(rcollection, 'timerLinearSolver',&
                               rappDescriptor%rtimerLinearSolver, .true.)
    call collct_setvalue_timer(rcollection, 'timerNonlinearSolver',&
                               rappDescriptor%rtimerNonlinearSolver, .true.)
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
    call parlst_getvalue_int(rparlist, ssectionName, 'Discretisation', i)
    call collct_setvalue_int(rcollection, 'Discretisation', i, .true.)

    call parlst_getvalue_int(rparlist, ssectionName, 'InviscidAFC', i)
    call collct_setvalue_int(rcollection, 'InviscidAFC', i, .true.)

    call parlst_getvalue_int(rparlist, ssectionName, 'SystemMatrix', i)
    call collct_setvalue_int(rcollection, 'SystemMatrix', i, .true.)

    call parlst_getvalue_int(rparlist, ssectionName, 'JacobianMatrix', i)
    call collct_setvalue_int(rcollection, 'JacobianMatrix', i, .true.)

    call parlst_getvalue_int(rparlist, ssectionName, 'TemplateMatrix', i)
    call collct_setvalue_int(rcollection, 'TemplateMatrix', i, .true.)

    if ((rappDescriptor%imasstype .ne. MASS_ZERO) .or.&
        (rappDescriptor%imassantidiffusiontype .ne. MASS_ZERO)) then

      call parlst_getvalue_int(rparlist, ssectionName, 'ConsistentMassMatrix', i)
      call collct_setvalue_int(rcollection, 'ConsistentMassMatrix', i, .true.)

      call parlst_getvalue_int(rparlist, ssectionName, 'LumpedMassMatrix', i)
      call collct_setvalue_int(rcollection, 'LumpedMassMatrix', i, .true.)

    else

      call collct_setvalue_int(rcollection, 'ConsistentMassMatrix', 0, .true.)
      call collct_setvalue_int(rcollection, 'LumpedMassMatrix',     0, .true.)

    end if

    select case(rappDescriptor%ndimension)
    case (NDIM1D)
      call parlst_getvalue_int(rparlist, ssectionName, 'CoeffMatrix_CX', i)
      call collct_setvalue_int(rcollection, 'CoeffMatrix_CX', i, .true.)
      call collct_setvalue_int(rcollection, 'CoeffMatrix_CY', 0, .true.)
      call collct_setvalue_int(rcollection, 'CoeffMatrix_CZ', 0, .true.)
    case (NDIM2D)
      call parlst_getvalue_int(rparlist, ssectionName, 'CoeffMatrix_CX', i)
      call collct_setvalue_int(rcollection, 'CoeffMatrix_CX', i, .true.)
      call parlst_getvalue_int(rparlist, ssectionName, 'CoeffMatrix_CY', i)
      call collct_setvalue_int(rcollection, 'CoeffMatrix_CY', i, .true.)
      call collct_setvalue_int(rcollection, 'CoeffMatrix_CZ', 0, .true.)
    case (NDIM3D)
      call parlst_getvalue_int(rparlist, ssectionName, 'CoeffMatrix_CX', i)
      call collct_setvalue_int(rcollection, 'CoeffMatrix_CX', i, .true.)
      call parlst_getvalue_int(rparlist, ssectionName, 'CoeffMatrix_CY', i)
      call collct_setvalue_int(rcollection, 'CoeffMatrix_CY', i, .true.)
      call parlst_getvalue_int(rparlist, ssectionName, 'CoeffMatrix_CZ', i)
      call collct_setvalue_int(rcollection, 'CoeffMatrix_CZ', i, .true.)
    case DEFAULT
      call collct_setvalue_int(rcollection, 'CoeffMatrix_CX',  0, .true.)
      call collct_setvalue_int(rcollection, 'CoeffMatrix_CY',  0, .true.)
      call collct_setvalue_int(rcollection, 'CoeffMatrix_CZ',  0, .true.)
    end select
  
  end subroutine euler_initCollection

  !*****************************************************************************

!<subroutine>

  subroutine euler_initSolvers(rparlist, ssectionName, rtimestep, rsolver)

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
    
  end subroutine euler_initSolvers

  !*****************************************************************************

!<subroutine>

  subroutine euler_initProblem(rparlist, ssectionName, nlmin, nlmax,&
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
    character(LEN=SYS_STRLEN) :: sinviscidName

    ! abstract problem descriptor
    type(t_problemDescriptor) :: rproblemDescriptor

    ! pointer to the problem level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! local variables
    integer :: inviscidAFC
    integer :: iconvToTria
    

    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist, ssectionName, 'inviscid', sinviscidName)
    call parlst_getvalue_string(rparlist, ssectionName, 'trifile', rproblemDescriptor%trifile)
    call parlst_getvalue_string(rparlist, ssectionName, 'prmfile', rproblemDescriptor%prmfile, '')
    call parlst_getvalue_int(rparlist, ssectionName, 'ndimension', rproblemDescriptor%ndimension)
    
    ! Set additional problem descriptor
    rproblemDescriptor%ndiscretisation = 1   ! one discretisation
    rproblemDescriptor%nafcstab        = 1   ! for inviscid fluxes
    rproblemDescriptor%nlmin           = nlmin
    rproblemDescriptor%nlmax           = nlmax
    rproblemDescriptor%nmatrixScalar   = rproblemDescriptor%ndimension + 5
    rproblemDescriptor%nmatrixBlock    = 2   ! for system matrix and Jacobian
    rproblemDescriptor%nvectorScalar   = 0
    rproblemDescriptor%nvectorBlock    = 0

    ! Check if quadrilaterals should be converted to triangles
    call parlst_getvalue_int(rparlist, ssectionName, 'iconvtotria', iconvToTria, 0)
    if (iconvToTria .ne. 0) then
      rproblemDescriptor%iproblemSpec = rproblemDescriptor%iproblemSpec &
                                      + PROBDESC_MSPEC_CONVTRIANGLES
    end if

    ! Initialize problem structure
    call problem_initProblem(rproblemDescriptor, rproblem)

    
    ! Initialize the stabilisation structure
    inviscidAFC = collct_getvalue_int(rcollection, 'inviscidAFC')
    
    ! Loop over all problem levels
    p_rproblemLevel => rproblem%p_rproblemLevelMax
    do while(associated(p_rproblemLevel))

      if (inviscidAFC > 0) then
        call afcstab_initFromParameterlist(rparlist, sinviscidName,&
                                           p_rproblemLevel%Rafcstab(inviscidAFC))
      end if
      
      ! Switch to next coarser level
      p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
    end do

  end subroutine euler_initProblem
  
  !*****************************************************************************

!<subroutine>

  subroutine euler_initProblemLevel(rappDescriptor, rproblemLevel, rcollection)

!<description>
    ! This subroutine initielizes the individual problem level. It
    ! generates the discretization, the template matrix and the
    ! coefficient matrices as duplicates of the template matrix.
!</description>

!<input>
    ! application descriptor
    type(t_euler), intent(IN) :: rappDescriptor
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT), target :: rproblemLevel

    ! collection structure
    type(t_collection), intent(INOUT) :: rcollection
!</output>
!</subroutine>

    ! local variables
    integer :: templateMatrix
    integer :: systemMatrix
    integer :: jacobianMatrix
    integer :: consistentMassMatrix
    integer :: lumpedMassMatrix
    integer :: coeffMatrix_CX
    integer :: coeffMatrix_CY
    integer :: coeffMatrix_CZ
    integer :: inviscidAFC
    integer :: discretisation
    integer :: ivar,jvar
    
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    type(t_triangulation) , pointer :: p_rtriangulation
    type(t_boundary) , pointer :: p_rboundary


    ! Retrieve application specific parameters from the collection
    templateMatrix       = collct_getvalue_int(rcollection, 'templatematrix')
    systemMatrix         = collct_getvalue_int(rcollection, 'systemmatrix')
    jacobianMatrix       = collct_getvalue_int(rcollection, 'jacobianmatrix')
    consistentMassMatrix = collct_getvalue_int(rcollection, 'consistentmassmatrix')
    lumpedMassMatrix     = collct_getvalue_int(rcollection, 'lumpedmassmatrix')
    coeffMatrix_CX       = collct_getvalue_int(rcollection, 'coeffmatrix_cx')
    coeffMatrix_CY       = collct_getvalue_int(rcollection, 'coeffmatrix_cy')
    coeffMatrix_CZ       = collct_getvalue_int(rcollection, 'coeffmatrix_cz')
    inviscidAFC          = collct_getvalue_int(rcollection, 'inviscidAFC')
    discretisation       = collct_getvalue_int(rcollection, 'discretisation')

    ! Set pointers
    p_rtriangulation  => rproblemLevel%rtriangulation
    p_rboundary       => rproblemLevel%p_rproblem%rboundary


    ! Create discretisation structure
    if (discretisation > 0) then
    
      ! Initialize the discretization structure
      p_rdiscretisation => rproblemLevel%Rdiscretisation(discretisation)
      if (p_rdiscretisation%ndimension .eq. 0) then
        select case(rappDescriptor%isystemformat)
        case (SYSTEM_INTERLEAVEFORMAT)
          call spdiscr_initBlockDiscr(p_rdiscretisation, 1, rproblemLevel%rtriangulation)
        case (SYSTEM_BLOCKFORMAT)
          call spdiscr_initBlockDiscr(p_rdiscretisation,&
                                      euler_getNVAR(rappDescriptor), p_rtriangulation)
        case DEFAULT
          call output_line('Unsupported system format!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'euler_initProblemLevel')
          call sys_halt()
        end select
      end if

      ! Get spatial dimension
      select case(p_rtriangulation%ndim)
      case (NDIM1D)
        select case(rappDescriptor%celement)
        case (-1,1,11)
          ! P1=Q1 finite elements
          call spdiscr_initDiscr_simple(p_rdiscretisation%RspatialDiscr(1), &
                                        EL_E001_1D, SPDISC_CUB_AUTOMATIC,&
                                        p_rtriangulation, p_rboundary)
          ! P2=Q2 finite elements
          call spdiscr_initDiscr_simple(p_rdiscretisation%RspatialDiscr(1), &
                                        EL_E002_1D, SPDISC_CUB_AUTOMATIC,&
                                        p_rtriangulation, p_rboundary)
        case DEFAULT
          call output_line('Unsupproted element type!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'euler_initProblemLevel')
          call sys_halt()
        end select
        
      case (NDIM2D)
        select case(rappDescriptor%celement)
        case (1)
          ! P1 finite elements
          call spdiscr_initDiscr_simple(p_rdiscretisation%RspatialDiscr(1), &
                                        EL_E001, SPDISC_CUB_AUTOMATIC,&
                                        p_rtriangulation, p_rboundary)
        case (2)
          ! P2 finite elements
          call spdiscr_initDiscr_simple(p_rdiscretisation%RspatialDiscr(1), &
                                        EL_E002, SPDISC_CUB_AUTOMATIC,&
                                        p_rtriangulation, p_rboundary)
        case (11)
          ! Q1 finite elements
          call spdiscr_initDiscr_simple(p_rdiscretisation%RspatialDiscr(1), &
                                        EL_E011, SPDISC_CUB_AUTOMATIC,&
                                        p_rtriangulation, p_rboundary)
        case (12)
          ! Q2 finite elements
          call spdiscr_initDiscr_simple(p_rdiscretisation%RspatialDiscr(1), &
                                        EL_E013, SPDISC_CUB_AUTOMATIC,&
                                        p_rtriangulation, p_rboundary)
        case (-1)
          ! mixed P1/Q1 finite elements
          call spdiscr_initDiscr_triquad(p_rdiscretisation%RspatialDiscr(1), &
                                         EL_E001, EL_E011, SPDISC_CUB_AUTOMATIC,&
                                         SPDISC_CUB_AUTOMATIC,&
                                         p_rtriangulation, p_rboundary)
        case (-2)
          ! mixed P2/Q2 finite elements
          call spdiscr_initDiscr_triquad(p_rdiscretisation%RspatialDiscr(1), &
                                         EL_E002, EL_E013, SPDISC_CUB_AUTOMATIC,&
                                         SPDISC_CUB_AUTOMATIC,&
                                         p_rtriangulation, p_rboundary)
        case DEFAULT
          call output_line('Unsupproted element type!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'euler_initProblemLevel')
          call sys_halt()
        end select
        
      case (NDIM3D)
        select case(rappDescriptor%celement)
        case (1)
          ! P1 finite elements
          call spdiscr_initDiscr_simple(p_rdiscretisation%RspatialDiscr(1), &
                                        EL_E001_3D, SPDISC_CUB_AUTOMATIC,&
                                        p_rtriangulation, p_rboundary)
        case (11)
          ! Q1 finite elements
          call spdiscr_initDiscr_simple(p_rdiscretisation%RspatialDiscr(1), &
                                        EL_E010_3D, SPDISC_CUB_AUTOMATIC,&
                                        p_rtriangulation, p_rboundary)
        case DEFAULT
          call output_line('Unsupproted element type!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'euler_initProblemLevel')
          call sys_halt()
        end select
        
      case DEFAULT
        call output_line('Invalid number of spatial dimensions',&
                         OU_CLASS_ERROR,OU_MODE_STD,'euler_initProblemLevel')
        call sys_halt()
      end select
    
      ! Duplicate scalar discretisation structure for block matrix format
      if (rappDescriptor%isystemFormat .eq. SYSTEM_BLOCKFORMAT) then
        do ivar = 2, euler_getNVAR(rappDescriptor)
          call spdiscr_duplicateDiscrSc(&
              p_rdiscretisation%RspatialDiscr(1),&
              p_rdiscretisation%RspatialDiscr(ivar), .true.)
        end do
      end if

    end if
    

    ! If the template matrix has no structure data then generate the
    ! finite element matrix sparsity structure based on the spatial
    ! descretization and store it as the template matrix. Otherwise we
    ! assume that the template matrix has been generated externally.
    if (.not.lsyssc_hasMatrixStructure(rproblemLevel%Rmatrix(templateMatrix))) then
      call bilf_createMatrixStructure(p_rdiscretisation%RspatialDiscr(1),&
                                      rappDescriptor%imatrixFormat,&
                                      rproblemLevel%Rmatrix(templateMatrix))
    end if


    ! Create system matrix
    if (systemMatrix > 0) then
      select case(rappDescriptor%isystemFormat)

      case (SYSTEM_INTERLEAVEFORMAT)
        ! The global operator is stored as an interleave matrix with
        ! NVAR components. However, the row and column structure of
        ! the template matrix can be adopted without modification
        if (lsyssc_hasMatrixStructure(rproblemLevel%Rmatrix(systemMatrix))) then

          ! Release pseudo block matrix
          call lsysbl_releaseMatrix(rproblemLevel%RmatrixBlock(systemMatrix))

          ! Resize scalar matrix
          call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(systemMatrix),&
                                   rproblemLevel%Rmatrix(templateMatrix)%NEQ,&
                                   rproblemLevel%Rmatrix(templateMatrix)%NCOLS,&
                                   rproblemLevel%Rmatrix(templateMatrix)%NA,&
                                   .false., .false., bforce=.true.)
          
        else   ! System matrix has no structure

          call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(templateMatrix),&
                                      rproblemLevel%Rmatrix(systemMatrix),&
                                      LSYSSC_DUP_SHARE, LSYSSC_DUP_REMOVE)

          ! Set number of variables per node
          rproblemLevel%Rmatrix(systemMatrix)%NVAR = euler_getNVAR(rappDescriptor)
          
          ! What matrix format should be used?
          select case(rappDescriptor%imatrixFormat)
          case (LSYSSC_MATRIX7)
            rproblemLevel%Rmatrix(systemMatrix)%cmatrixFormat = LSYSSC_MATRIX7INTL
            
          case (LSYSSC_MATRIX9)
            rproblemLevel%Rmatrix(systemMatrix)%cmatrixFormat = LSYSSC_MATRIX9INTL
            
          case DEFAULT
            call output_line('Unsupported matrix format!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'euler_initProblemLevel')
            call sys_halt()
          end select
          
          ! What kind of global operator should be adopted?
          select case(rappDescriptor%isystemCoupling)
          case (SYSTEM_SEGREGATED)
            rproblemLevel%Rmatrix(systemMatrix)%cinterleavematrixFormat = LSYSSC_MATRIXD
            
          case (SYSTEM_ALLCOUPLED)
            rproblemLevel%Rmatrix(systemMatrix)%cinterleavematrixFormat = LSYSSC_MATRIX1
            
          case DEFAULT
            call output_line('Unsupported interleave matrix format!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'euler_initProblemLevel')
            call sys_halt()
          end select
        
          ! Create global operator physically
          call lsyssc_allocEmptyMatrix(rproblemLevel%Rmatrix(systemMatrix),&
                                       LSYSSC_SETM_UNDEFINED)
        end if

        ! Create pseudo block matrix from global operator
        call lsysbl_createMatFromScalar(rproblemLevel%Rmatrix(systemMatrix),&
                                        rproblemLevel%RmatrixBlock(systemMatrix),&
                                        p_rdiscretisation)


      case (SYSTEM_BLOCKFORMAT)
        ! The global operator is stored as a block matrix with
        ! NVARxNVAR blocks made up from scalar matrices

        if ((rproblemLevel%RmatrixBlock(systemMatrix)%nblocksPerRow .ne. 0) .and.&
            (rproblemLevel%RmatrixBlock(systemMatrix)%nblocksPerCol .ne. 0)) then

          ! What kind of global operator should be adopted?
          select case(rappDescriptor%isystemCoupling)
            
          case (SYSTEM_SEGREGATED)
            ! Create only NVAR diagonal blocks
            do ivar = 1, euler_getNVAR(rappDescriptor)
              call lsyssc_resizeMatrix(&
                  rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                  rproblemLevel%Rmatrix(templateMatrix), .false., .false., .true.)
            end do

          case (SYSTEM_ALLCOUPLED)
            ! Create all NVAR x NVAR blocks
            do ivar = 1, euler_getNVAR(rappDescriptor)
              do jvar = 1, euler_getNVAR(rappDescriptor)
                call lsyssc_resizeMatrix(&
                    rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,jvar),&
                    rproblemLevel%Rmatrix(templateMatrix), .false., .false., .true.)
              end do
            end do

          case DEFAULT
            call output_line('Unsupported block matrix format!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'euler_initProblemLevel')
            call sys_halt()
          end select

        else   ! System matrix has no structure

          ! Create empty NVARxNVAR block matrix directly
          call lsysbl_createEmptyMatrix(rproblemLevel%RmatrixBlock(systemMatrix),&
                                        euler_getNVAR(rappDescriptor),&
                                        euler_getNVAR(rappDescriptor))

          ! Specify matrix as 'group matrix'
          rproblemLevel%RmatrixBlock(systemMatrix)%imatrixSpec = LSYSBS_MSPEC_GROUPMATRIX
          
          ! What kind of global operator should be adopted?
          select case(rappDescriptor%isystemCoupling)
            
          case (SYSTEM_SEGREGATED)
            ! Create only NVAR diagonal blocks
            do ivar = 1, euler_getNVAR(rappDescriptor)
              call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(templateMatrix),&
                                          rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                                          LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
            end do
          
          case (SYSTEM_ALLCOUPLED)
            ! Create all NVAR x NVAR blocks
            do ivar = 1, euler_getNVAR(rappDescriptor)
              do jvar = 1, euler_getNVAR(rappDescriptor)
                call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(templateMatrix),&
                                            rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,jvar),&
                                            LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
              end do
            end do
            
          case DEFAULT
            call output_line('Unsupported block matrix format!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'euler_initProblemLevel')
            call sys_halt()
          end select

        end if
        
        ! Update internal structure of block matrix
        call lsysbl_updateMatStrucInfo(rproblemLevel%RmatrixBlock(systemMatrix))
        
      case DEFAULT
        call output_line('Unsupported system format!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'euler_initProblemLevel')
        call sys_halt()
      end select
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
    if (inviscidAFC > 0) then
      if (rproblemLevel%Rafcstab(inviscidAFC)%iSpec .eq. AFCSTAB_UNDEFINED) then
        call gfsys_initStabilisation(rproblemLevel%RmatrixBlock(systemMatrix),&
                                     rproblemLevel%Rafcstab(inviscidAFC))
      else
        call afcstab_resizeStabilisation(rproblemLevel%Rafcstab(inviscidAFC),&
                                         rproblemLevel%Rmatrix(templateMatrix))
        rproblemLevel%Rafcstab(inviscidAFC)%iSpec =&
            iand(rproblemLevel%Rafcstab(inviscidAFC)%iSpec,&
            not(AFCSTAB_SUBDIAGONALEDGES))
      end if
    end if
    
  end subroutine euler_initProblemLevel

  !*****************************************************************************

!<subroutine>

  subroutine euler_initAllProblemLevels(rappDescriptor, rproblem, rcollection)

!<description>
    ! This subroutine initializes the all problem levels attached to
    ! the global problem structure. It generates the discretization,
    ! the template matrix and the coefficient matrices as duplicates
    ! of the template matrix.
!</description>

!<input>
    ! application descriptor
    type(t_euler), intent(IN) :: rappDescriptor
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
      call euler_initProblemLevel(rappDescriptor, p_rproblemLevel, rcollection)
      
      ! Switch to next coarser level
      p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
    end do

  end subroutine euler_initAllProblemLevels

  !*****************************************************************************

!<subroutine>

  subroutine euler_initSolution(rparlist, ssectionName, rproblemLevel, dtime, rvector)

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

!<inputoutput>
    ! solution vector
    type(t_vectorBlock), intent(INOUT) :: rvector
!</inputoutput>
!</subroutine>

    ! section names
    character(LEN=SYS_STRLEN) :: sindatfileName
    character(LEN=SYS_STRLEN) :: ssolutionName

    ! local variables
    type(t_fparser) :: rfparser
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Ddata
    real(DP), dimension(NDIM3D+1) :: Dvalue
    integer :: isolutiontype, icomp, iblock, ivar, nvar, ieq, neq, ndim

    ! symbolic variable names
    character(LEN=*), dimension(4), parameter ::&
                      cvariables = (/ (/'x'/), (/'y'/), (/'z'/), (/'t'/) /)

    
    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist, ssectionName, 'indatfile', sindatfileName)
    call parlst_getvalue_string(rparlist, ssectionName, 'ssolutionname', ssolutionName)
    call parlst_getvalue_int(rparlist, ssectionName, 'isolutiontype', isolutiontype)

    
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
      
      ! Initialization
      Dvalue = 0.0_DP
      icomp  = 0
      ndim   = size(p_DvertexCoords, 1)
      
      ! Loop over all blocks of the global solution vector
      do iblock = 1, rvector%nblocks

        ! Set pointer to data array
        call lsyssc_getbase_double(rvector%RvectorBlock(iblock), p_Ddata)

        ! Initialization for scalar subvectir
        neq  = rvector%RvectorBlock(iblock)%NEQ
        nvar = rvector%RvectorBlock(iblock)%NVAR

        ! Loop over all equations of the scalar subvector
        do ieq = 1, neq

          ! Set coordinates and evalution time
          Dvalue(1:ndim)   = p_DvertexCoords(:,ieq)
          Dvalue(NDIM3D+1) = dtime

          ! Loop over all variables of the solution vector
          do ivar = 1, nvar
            
            call fparser_evalFunction(rfparser, icomp+ivar, Dvalue, p_Ddata((ieq-1)*nvar+ivar))
            
          end do   ! ivar

        end do   ! ieq

        ! Increase component number for next subvector
        icomp = icomp+nvar

      end do   ! iblock
            
      call fparser_release(rfparser)


    case DEFAULT
      call output_line('Invalid type of solution profile!',&
                       OU_CLASS_ERROR, OU_MODE_STD, 'euler_initSolution')
      call sys_halt()
    end select

  end subroutine euler_initSolution

  !*****************************************************************************

!<subroutine>

  subroutine euler_outputSolution(rparlist, ssectionName, rproblemLevel,&
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
    type(t_ucdExport) :: rexport
    type(t_vectorScalar) :: rvector1, rvector2, rvector3
    real(DP), dimension(:), pointer :: p_Dsolution, p_Ddata1, p_Ddata2, p_Ddata3
    integer :: iformatUCD, isystemFormat, isize, ndim


    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist, ssectionName, 'output', soutputName)
    call parlst_getvalue_string(rparlist, trim(soutputName), 'ucdsolution', ucdsolution)
    call parlst_getvalue_int(rparlist, trim(soutputName), 'iformatucd', iformatUCD)
    call parlst_getvalue_int(rparlist, ssectionName, 'isystemformat', isystemformat)

    ! Initialize the UCD exporter
    call flagship_initUCDexport(rproblemLevel, ucdsolution,&
                                iformatUCD, rexport, ifilenumber)

    ! Increase filenumber by one
    ifilenumber = ifilenumber+1

    ! Set simulation time
    if (present(dtime)) call ucd_setSimulationTime(rexport, dtime)

    ! Set pointers
    call lsysbl_getbase_double(rsolutionPrimal, p_Dsolution)
    isize = size(p_Dsolution)/euler_getNVAR(rproblemLevel)
    ndim  = rproblemLevel%rtriangulation%ndim

    ! Create auxiliary vectors
    select case(ndim)
    case (NDIM1D)
      call lsyssc_createVector(rvector1, isize, .false.)
      call lsyssc_getbase_double(rvector1, p_Ddata1)
      
    case (NDIM2D)
      call lsyssc_createVector(rvector1, isize, .false.)
      call lsyssc_createVector(rvector2, isize, .false.)
      call lsyssc_getbase_double(rvector1, p_Ddata1)
      call lsyssc_getbase_double(rvector2, p_Ddata2)

    case (NDIM3D)
      call lsyssc_createVector(rvector1, isize, .false.)
      call lsyssc_createVector(rvector2, isize, .false.)
      call lsyssc_createVector(rvector3, isize, .false.)
      call lsyssc_getbase_double(rvector1, p_Ddata1)
      call lsyssc_getbase_double(rvector2, p_Ddata2)
      call lsyssc_getbase_double(rvector3, p_Ddata3)

    case DEFAULT
      call output_line('Invalid number of spatial dimensions',&
                       OU_CLASS_ERROR,OU_MODE_STD,'euler_outputSolution')
      call sys_halt()
    end select


    select case(isystemFormat)
    case(SYSTEM_INTERLEAVEFORMAT)
      
      select case(ndim)
      case (NDIM1D)
        call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR1D, 'velocity_x', p_Dsolution, p_Ddata1)
        call ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1)

        call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR1D, 'density', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'density', UCD_VAR_STANDARD, p_Ddata1)

        call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR1D, 'energy', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'energy', UCD_VAR_STANDARD, p_Ddata1)
        
        call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR1D, 'pressure', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'pressure', UCD_VAR_STANDARD, p_Ddata1)

        call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR1D, 'machnumber', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'machnumber', UCD_VAR_STANDARD, p_Ddata1)
        
        
      case (NDIM2D)
        call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR2D, 'velocity_x', p_Dsolution, p_Ddata1)
        call euler_getVarInterleaveFormat(rvector2%NEQ, NVAR2D, 'velocity_y', p_Dsolution, p_Ddata2)
        call ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1, p_Ddata2)

        call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR2D, 'density', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'density', UCD_VAR_STANDARD, p_Ddata1)

        call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR2D, 'energy', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'energy', UCD_VAR_STANDARD, p_Ddata1)
        
        call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR2D, 'pressure', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'pressure', UCD_VAR_STANDARD, p_Ddata1)

        call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR2D, 'machnumber', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'machnumber', UCD_VAR_STANDARD, p_Ddata1)

      case (NDIM3D)
        call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR3D, 'velocity_x', p_Dsolution, p_Ddata1)
        call euler_getVarInterleaveFormat(rvector2%NEQ, NVAR3D, 'velocity_y', p_Dsolution, p_Ddata2)
        call euler_getVarInterleaveFormat(rvector3%NEQ, NVAR3D, 'velocity_z', p_Dsolution, p_Ddata3)
        call ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1, p_Ddata2, p_Ddata3)
        
        call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR3D, 'density', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'density', UCD_VAR_STANDARD, p_Ddata1)

        call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR3D, 'energy', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'energy', UCD_VAR_STANDARD, p_Ddata1)
        
        call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR3D, 'pressure', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'pressure', UCD_VAR_STANDARD, p_Ddata1)

        call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR3D, 'machnumber', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'machnumber', UCD_VAR_STANDARD, p_Ddata1)

      end select
      
      
    case (SYSTEM_BLOCKFORMAT)

      select case(ndim)
      case (NDIM1D)
        call euler_getVarBlockFormat(rvector1%NEQ, NVAR1D, 'velocity_x', p_Dsolution, p_Ddata1)
        call ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1)

        call euler_getVarBlockFormat(rvector1%NEQ, NVAR1D, 'density', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'density', UCD_VAR_STANDARD, p_Ddata1)

        call euler_getVarBlockFormat(rvector1%NEQ, NVAR1D, 'energy', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'energy', UCD_VAR_STANDARD, p_Ddata1)
        
        call euler_getVarBlockFormat(rvector1%NEQ, NVAR1D, 'pressure', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'pressure', UCD_VAR_STANDARD, p_Ddata1)

        call euler_getVarBlockFormat(rvector1%NEQ, NVAR1D, 'machnumber', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'machnumber', UCD_VAR_STANDARD, p_Ddata1)
        
        
      case (NDIM2D)
        call euler_getVarBlockFormat(rvector1%NEQ, NVAR2D, 'velocity_x', p_Dsolution, p_Ddata1)
        call euler_getVarBlockFormat(rvector2%NEQ, NVAR2D, 'velocity_y', p_Dsolution, p_Ddata2)
        call ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1, p_Ddata2)

        call euler_getVarBlockFormat(rvector1%NEQ, NVAR2D, 'density', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'density', UCD_VAR_STANDARD, p_Ddata1)

        call euler_getVarBlockFormat(rvector1%NEQ, NVAR2D, 'energy', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'energy', UCD_VAR_STANDARD, p_Ddata1)
        
        call euler_getVarBlockFormat(rvector1%NEQ, NVAR2D, 'pressure', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'pressure', UCD_VAR_STANDARD, p_Ddata1)

        call euler_getVarBlockFormat(rvector1%NEQ, NVAR2D, 'machnumber', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'machnumber', UCD_VAR_STANDARD, p_Ddata1)

      case (NDIM3D)
        call euler_getVarBlockFormat(rvector1%NEQ, NVAR3D, 'velocity_x', p_Dsolution, p_Ddata1)
        call euler_getVarBlockFormat(rvector2%NEQ, NVAR3D, 'velocity_y', p_Dsolution, p_Ddata2)
        call euler_getVarBlockFormat(rvector3%NEQ, NVAR3D, 'velocity_z', p_Dsolution, p_Ddata3)
        call ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1, p_Ddata2, p_Ddata3)
        
        call euler_getVarBlockFormat(rvector1%NEQ, NVAR3D, 'density', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'density', UCD_VAR_STANDARD, p_Ddata1)

        call euler_getVarBlockFormat(rvector1%NEQ, NVAR3D, 'energy', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'energy', UCD_VAR_STANDARD, p_Ddata1)
        
        call euler_getVarBlockFormat(rvector1%NEQ, NVAR3D, 'pressure', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'pressure', UCD_VAR_STANDARD, p_Ddata1)

        call euler_getVarBlockFormat(rvector1%NEQ, NVAR3D, 'machnumber', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'machnumber', UCD_VAR_STANDARD, p_Ddata1)

      end select
      
    end select

    ! Release temporal memory
    call lsyssc_releaseVector(rvector1)
    call lsyssc_releaseVector(rvector2)
    call lsyssc_releaseVector(rvector3)

    ! Write UCD file
    call ucd_write  (rexport)
    call ucd_release(rexport)

  end subroutine euler_outputSolution

  !*****************************************************************************

!<subroutine>

  subroutine euler_outputStatistics(rappDescriptor, rtimerTotal)

!<description>
    ! This subroutine output application statistics
!</description>

!<input>
    ! application descriptor
    type(t_euler), intent(INOUT) :: rappDescriptor

    ! timer for total time measurement
    type(t_timer), intent(IN) :: rtimerTotal
!</input>
!</subroutine>

    ! local variable
    type(t_timer) :: rtimerSolution
    real(DP) :: dtotalTime, dfraction

    call output_lbrk()
    call output_line('Time measurement:')
    call output_line('-----------------')
    
    rtimerSolution = rappDescriptor%rtimerSolution
    call stat_subTimers(rappDescriptor%rtimerAssemblyMatrix, rtimerSolution)
    call stat_subTimers(rappDescriptor%rtimerAssemblyVector, rtimerSolution)

    dtotalTime = max(rtimerTotal%delapsedCPU, rtimerTotal%delapsedReal)
    dfraction  = 100.0_DP/dtotalTime

    call output_line('Time for computing solution   : '//&
                     trim(adjustl(sys_sdE(rtimerSolution%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*rtimerSolution%delapsedCPU, 5)))//' %')
    call output_line('Time for mesh adaptivity      : '//&
                     trim(adjustl(sys_sdE(rappDescriptor%rtimerAdaptation%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*rappDescriptor%rtimerAdaptation%delapsedCPU, 5)))//' %')
    call output_line('Time for error estimation     : '//&
                     trim(adjustl(sys_sdE(rappDescriptor%rtimerErrorEstimation%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*rappDescriptor%rtimerErrorEstimation%delapsedCPU, 5)))//' %')
    call output_line('Time for triangulation        : '//&
                     trim(adjustl(sys_sdE(rappDescriptor%rtimerTriangulation%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*rappDescriptor%rtimerTriangulation%delapsedCPU, 5)))//' %')
    call output_line('Time for coefficient assembly : '//&
                     trim(adjustl(sys_sdE(rappDescriptor%rtimerAssemblyCoeff%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*rappDescriptor%rtimerAssemblyCoeff%delapsedCPU, 5)))//' %')
    call output_line('Time for matrix assembly      : '//&
                     trim(adjustl(sys_sdE(rappDescriptor%rtimerAssemblyMatrix%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*rappDescriptor%rtimerAssemblyMatrix%delapsedCPU, 5)))//' %')
    call output_line('Time for vector assembly:       '//&
                     trim(adjustl(sys_sdE(rappDescriptor%rtimerAssemblyVector%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*rappDescriptor%rtimerAssemblyVector%delapsedCPU, 5)))//' %')
    call output_line('Time for pre-/post-processing : '//&
                     trim(adjustl(sys_sdE(rappDescriptor%rtimerPrePostprocess%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*rappDescriptor%rtimerPrePostprocess%delapsedCPU, 5)))//' %')
    call output_lbrk()
    call output_line('Time for total simulation     : '//&
                     trim(adjustl(sys_sdE(dtotalTime, 5))))
    call output_lbrk()

  end subroutine euler_outputStatistics

  !*****************************************************************************

!<subroutine>

  subroutine euler_estimateRecoveryError(rparlist, ssectionName, rproblemLevel,&
                                          rsolution, dtime, rerror, derror)

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
    type(t_parlist), intent(IN) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(IN) :: ssectionName

    ! solution vector
    type(t_vectorBlock), intent(IN) :: rsolution

    ! problem level structure
    type(t_problemLevel), intent(IN) :: rproblemLevel

    ! simulation time
    real(DP), intent(IN) :: dtime
!</input>

!<output>
    ! element-wise error distribution
    type(t_vectorScalar), intent(OUT) :: rerror

    ! global error
    real(DP), intent(OUT) :: derror
!</output>
!</subroutine>

    ! section names
    character(LEN=SYS_STRLEN) :: sindatfileName
    character(LEN=SYS_STRLEN) :: serrorestimatorName
    character(LEN=SYS_STRLEN) :: sexactsolutionName

    ! local variables
    type(t_vectorScalar) :: rvectorScalar, rvectorTmp
    real(DP), dimension(:), pointer :: p_Ddata
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:,:), pointer :: p_IneighboursAtElement
    logical, dimension(:), pointer :: p_BisactiveElement
    character(LEN=SYS_STRLEN) :: serrorvariable
    real(DP) :: dnoiseFilter, dabsFilter, dvalue,&
                dprotectLayerTolerance, derrorTmp
    integer :: ierrorEstimator, igridindicator, iexactsolutiontype
    integer :: iprotectLayer, nprotectLayers, ierrorVariable, nerrorVariables
    integer :: h_BisactiveElement
    
    ! symbolic variable names
    character(LEN=*), dimension(4), parameter ::&
                      cvariables = (/ (/'x'/), (/'y'/), (/'z'/), (/'t'/) /)


    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist, ssectionName, 'indatfile', sindatfileName)
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

    nerrorVariables = parlst_querysubstrings(rparlist, trim(serrorestimatorName), 'serrorvariable')

    ! Loop over all error variables
    do ierrorVariable = 1, nerrorVariables
      
      ! Get name of error variable
      call parlst_getvalue_string(rparlist, trim(serrorestimatorName), 'serrorvariable',&
                                  serrorVariable, isubString=ierrorVariable)

      ! Extract scalar variable from vector of conservative variables
      call euler_getVariable(rsolution, serrorVariable, rvectorScalar)

      ! What type of error estimator are we?
      select case(ierrorEstimator)
        
      case (ERREST_L2PROJECTION)
        call lsyssc_createVector(rvectorTmp, rproblemLevel%rtriangulation%NEL, .false.)
        call ppgrd_calcGradientError(rvectorScalar, derrorTmp,&
                                     PPGRD_INTERPOL, 0, rvectorTmp)

      case (ERREST_SPR_VERTEX)
        call lsyssc_createVector(rvectorTmp, rproblemLevel%rtriangulation%NEL, .false.)
        call ppgrd_calcGradientError(rvectorScalar, derrorTmp,&
                                     PPGRD_ZZTECHNIQUE, PPGRD_NODEPATCH, rvectorTmp)

      case (ERREST_SPR_ELEMENT)
        call lsyssc_createVector(rvectorTmp, rproblemLevel%rtriangulation%NEL, .false.)
        call ppgrd_calcGradientError(rvectorScalar, derrorTmp,&
                                     PPGRD_ZZTECHNIQUE, PPGRD_ELEMPATCH, rvectorTmp)

      case (ERREST_SPR_FACE)
        call lsyssc_createVector(rvectorTmp, rproblemLevel%rtriangulation%NEL, .false.)
        call ppgrd_calcGradientError(rvectorScalar, derrorTmp,&
                                     PPGRD_ZZTECHNIQUE, PPGRD_FACEPATCH, rvectorTmp)
      
      case (ERREST_LIMAVR)
        call lsyssc_createVector(rvectorTmp, rproblemLevel%rtriangulation%NEL, .false.)
        call ppgrd_calcGradientError(rvectorScalar, derrorTmp,&
                                     PPGRD_LATECHNIQUE, 0, rvectorTmp)
      
      case (ERREST_SECONDDIFF)
        call parlst_getvalue_double(rparlist, trim(serrorestimatorName), 'dnoiseFilter', dnoiseFilter)
        call parlst_getvalue_double(rparlist, trim(serrorestimatorName), 'dabsFilter', dabsFilter)
        call ppind_secondDifference(rvectorScalar, dnoiseFilter, dabsFilter, rvectorTmp)
        
        ! This is no error estimator
        derrorTmp = 1.0
        
      case DEFAULT
        call output_line('Invalid type of error estimator!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'euler_estimateRecoveryError')
        call sys_halt()
      end select

      
      ! Compute the global and element-wise error
      if (ierrorVariable .eq. 1) then

        ! Initialize the global error
        derror = derrorTmp
        
        ! Initialize the element-wise error
        call lsyssc_copyVector(rvectorTmp, rerror)
      else
        
        ! Update the global error
        derror = derror + derrorTmp

        ! Update the element-wise error
        call lsyssc_vectorLinearComb(rvectorTmp, rerror, 1.0_DP, 1.0_DP)
      end if

      ! Release scalar variable and temporal error
      call lsyssc_releaseVector(rvectorScalar)
      call lsyssc_releaseVector(rvectorTmp)
    end do
    
    ! Scale the global and element-wise error by the number of error variables
    dvalue = 1.0_DP/real(nerrorVariables, DP)
    derror = derror*dvalue
    call lsyssc_scaleVector(rerror, dvalue)


    !---------------------------------------------------------------------------
    ! Compute the effectivity index
    !---------------------------------------------------------------------------
    
    call output_lbrk()
    call output_line('Error Analysis')
    call output_line('--------------')
    call output_line('estimated H1-error: '//trim(sys_sdEP(derror,15,6)))
    

    !---------------------------------------------------------------------------
    ! Apply the adaptation strategy
    !---------------------------------------------------------------------------

    ! Well at the moment we use the estimator 'as is'

    
    !---------------------------------------------------------------------------
    ! Calculate protection layers
    !---------------------------------------------------------------------------
    
    if (nprotectLayers > 0) then

      ! Create auxiliary memory
      h_BisactiveElement = ST_NOHANDLE
      call storage_new('codire_estimateRecoveryError',' BisactiveElement',&
                       rproblemLevel%rtriangulation%NEL, ST_LOGICAL,&
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
        call doProtectionLayerUniform(p_IverticesAtElement, p_IneighboursAtElement,&
                                      rproblemLevel%rtriangulation%NEL,&
                                      dprotectLayerTolerance, p_Ddata, p_BisActiveElement)
      end do

      ! Release memory
      call storage_free(h_BisactiveElement)
      
    end if

  contains
    
    ! Here, the real working routines follow.
    
    !**************************************************************
    ! Compute one uniformly distributed protection layer

    subroutine doProtectionLayerUniform(IverticesAtElement, IneighboursAtElement, NEL,&
                                        dthreshold, Ddata, BisactiveElement)

      integer, dimension(:,:), intent(IN) :: IverticesAtElement
      integer, dimension(:,:), intent(IN) :: IneighboursAtElement     
      real(DP), intent(IN) :: dthreshold
      integer, intent(IN) :: NEL
      
      real(DP), dimension(:), intent(INOUT) :: Ddata
      logical, dimension(:), intent(INOUT) :: BisactiveElement
      
      
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

  end subroutine euler_estimateRecoveryError

  !*****************************************************************************

!<subroutine>

  subroutine euler_adaptTriangulation(rhadapt, rtriangulationSrc, rindicator,&
                                      rcollection, rtriangulationDest)

!<description>
    ! This subroutine performs h-adaptation for the given triangulation
!</description>

!<inputoutput>
    ! adaptation structure
    type(t_hadapt), intent(INOUT) :: rhadapt

    ! source triangulation structure
    type(t_triangulation), intent(INOUT), target :: rtriangulationSrc
    
    ! element-wise indicator
    type(t_vectorScalar), intent(INOUT) :: rindicator

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
    integer :: isystemFormat


    ! Check if adaptation structure has been prepared
    if (rhadapt%iSpec .eq. HADAPT_HAS_PARAMETERS) then
      
      ! Initialize adaptation structure from triangulation
      call hadapt_initFromTriangulation(rhadapt, rtriangulationSrc)
     
    else
     
      ! Refresh adaptation structure
      call hadapt_refreshAdaptation(rhadapt, rtriangulationSrc)
      
    end if

    isystemFormat = collct_getvalue_int(rcollection, 'isystemformat')

    select case (isystemFormat)

    case (SYSTEM_INTERLEAVEFORMAT)

      ! How many spatial dimensions do we have?
      select case(rtriangulationSrc%ndim)
      case (NDIM1D)
        call euler_hadaptCallbackScalar1D(rcollection, HADAPT_OPR_INITCALLBACK, Ivalue, Ivalue)
        call hadapt_performAdaptation(rhadapt, rindicator, rcollection, euler_hadaptCallbackScalar1D)
        
      case (NDIM2D)
        call euler_hadaptCallbackScalar2D(rcollection, HADAPT_OPR_INITCALLBACK, Ivalue, Ivalue)
        call hadapt_performAdaptation(rhadapt, rindicator, rcollection, euler_hadaptCallbackScalar2D)
        
      case (NDIM3D)
        call euler_hadaptCallbackScalar3D(rcollection, HADAPT_OPR_INITCALLBACK, Ivalue, Ivalue)
        call hadapt_performAdaptation(rhadapt, rindicator, rcollection, euler_hadaptCallbackScalar3D)
      end select

    case (SYSTEM_BLOCKFORMAT)

      ! How many spatial dimensions do we have?
      select case(rtriangulationSrc%ndim)
      case (NDIM1D)
        call euler_hadaptCallbackBlock1D(rcollection, HADAPT_OPR_INITCALLBACK, Ivalue, Ivalue)
        call hadapt_performAdaptation(rhadapt, rindicator, rcollection, euler_hadaptCallbackBlock1D)
        
      case (NDIM2D)
        call euler_hadaptCallbackBlock2D(rcollection, HADAPT_OPR_INITCALLBACK, Ivalue, Ivalue)
        call hadapt_performAdaptation(rhadapt, rindicator, rcollection, euler_hadaptCallbackBlock2D)
        
      case (NDIM3D)
        call euler_hadaptCallbackBlock3D(rcollection, HADAPT_OPR_INITCALLBACK, Ivalue, Ivalue)
        call hadapt_performAdaptation(rhadapt, rindicator, rcollection, euler_hadaptCallbackBlock3D)
      end select

    case DEFAULT

      ! How many spatial dimensions do we have?
      select case(rtriangulationSrc%ndim)
      case (NDIM1D)
        call hadapt_performAdaptation(rhadapt, rindicator)
        
      case (NDIM2D)
        call hadapt_performAdaptation(rhadapt, rindicator)
        
      case (NDIM3D)
        call hadapt_performAdaptation(rhadapt, rindicator)
      end select

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

  end subroutine euler_adaptTriangulation

  !*****************************************************************************

!<subroutine>

    subroutine euler_solveTransientPrimal(rappDescriptor, rparlist, ssectionName,&
                                          rbdrCond, rproblem, rtimestep, rsolver,&
                                          rsolution, rcollection)

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
    ! parameter list
    type(t_parlist), intent(IN) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(IN) :: ssectionName

    ! boundary condition structure
    type(t_boundaryCondition), intent(IN) :: rbdrCond
!</input>

!<inputoutput>
    ! application descriptor
    type(t_euler), intent(INOUT) :: rappDescriptor

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

    ! Pointer to the multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! Pointer to the discretisation structure
    type(t_blockDiscretisation), pointer :: p_rdiscretisation

    ! Vector for the element-wise error distribution
    type(t_vectorScalar) :: relementError
    
    ! Structure for h-adaptation
    type(t_hadapt) :: rhadapt

    ! Structure for the sparsity pattern
    type(t_graph) :: rgraph

    ! section names
    character(LEN=SYS_STRLEN) :: sadaptivityName
    character(LEN=SYS_STRLEN) :: soutputName

    ! local variables
    real(dp) :: derror, dstepUCD, dtimeUCD, dstepAdapt, dtimeAdapt
    integer :: templateMatrix, systemMatrix, discretisation
    integer :: isize, ipreadapt, npreadapt
    integer, external :: signal_SIGINT


    ! Start time measurement for pre-processing
    call stat_startTimer(rappDescriptor%rtimerPrePostprocess, STAT_TIMERSHORT)

    ! Get position of discretisation structure
    call parlst_getvalue_int(rparlist, ssectionName, 'discretisation', discretisation)

    ! Set pointer to maximum problem level
    p_rproblemLevel => rproblem%p_rproblemLevelMax
    p_rdiscretisation => p_rproblemLevel%Rdiscretisation(discretisation)

    ! Create the solution vector
    call lsysbl_createVectorBlock(p_rdiscretisation, rsolution, .false., ST_DOUBLE)
    if (p_rdiscretisation%ncomponents .ne. euler_getNVAR(p_rproblemLevel)) then
      rsolution%RvectorBlock(1)%NVAR = euler_getNVAR(p_rproblemLevel)
      isize = rsolution%NEQ*euler_getNVAR(p_rproblemLevel)
      call lsysbl_resizeVectorBlock(rsolution, isize, .false., .false.)
    end if
    
    ! Initialize the solution vector and impose boundary conditions explicitly
    call euler_initSolution(rparlist, ssectionName, p_rproblemLevel, 0.0_DP, rsolution)
    select case(rappDescriptor%ndimension)
    case (NDIM1D)
      call bdrf_filterVectorExplicit(rbdrCond, p_rproblemLevel%rtriangulation,&
                                     rsolution, 0.0_DP, rproblem%rboundary,&
                                     euler_calcBoundaryvalues1d)
    case (NDIM2D)
      call bdrf_filterVectorExplicit(rbdrCond, p_rproblemLevel%rtriangulation,&
                                     rsolution, 0.0_DP, rproblem%rboundary,&
                                     euler_calcBoundaryvalues2d)
    case (NDIM3D)
      call bdrf_filterVectorExplicit(rbdrCond, p_rproblemLevel%rtriangulation,&
                                     rsolution, 0.0_DP, rproblem%rboundary,&
                                     euler_calcBoundaryvalues3d)
    end select

    ! Initialize timer for intermediate UCD exporter
    dtimeUCD = rtimestep%dinitialTime
    call parlst_getvalue_string(rparlist, ssectionName, 'output', soutputName)
    call parlst_getvalue_double(rparlist, trim(soutputName), 'dstepUCD', dstepUCD, 0.0_DP)


    !---------------------------------------------------------------------------
    ! Initialize the h-adaptation structure and perform pre-adaptation
    !---------------------------------------------------------------------------
    
    call parlst_getvalue_string(rparlist, ssectionName, 'adaptivity', sadaptivityName, '')
    if (trim(adjustl(sadaptivityName)) .ne. '') then

      call parlst_getvalue_double(rparlist, trim(sadaptivityName), 'dstepAdapt', dstepAdapt)
      call parlst_getvalue_double(rparlist, trim(sadaptivityName), 'dtimeAdapt', dtimeAdapt)
      call parlst_getvalue_int(rparlist, trim(sadaptivityName), 'npreadapt', npreadapt)


      if ((dstepAdapt > 0.0_DP) .or. (npreadapt > 0)) then

        ! Initialize adaptation structure from parameter list
        call hadapt_initFromParameterlist(rhadapt, rparlist, sadaptivityName)

        ! Generate a dynamic graph for the sparsity pattern and attach
        ! it to the collection structure which is used to pass this
        ! type to the callback function for h-adaptation
        templateMatrix = collct_getvalue_int(rcollection, 'templateMatrix')
        call grph_createGraphFromMatrix(p_rproblemLevel%Rmatrix(templateMatrix), rgraph)
        call collct_setvalue_graph(rcollection, 'sparsitypattern', rgraph, .true.)
        
        ! Attach the primal solution vector to the collection structure
        call collct_setvalue_vec(rcollection, 'solutionvector', rsolution, .true.)


        ! Perform pre-adaptation?
        if (npreadapt > 0) then

          ! Set the names of the template matrix and the solution vector
          rcollection%SquickAccess(1) = 'sparsitypattern'
          rcollection%SquickAccess(2) = 'solutionvector'
          
          ! Perform number of pre-adaptation steps
          do ipreadapt = 1, npreadapt

            ! Compute the error estimator using recovery techniques
            call euler_estimateRecoveryError(rparlist, ssectionname, p_rproblemLevel,&
                                             rsolution, 0.0_DP, relementError, derror)

            ! Perform h-adaptation and update the triangulation structure
            call euler_adaptTriangulation(rhadapt, p_rproblemLevel%rtriangulation,&
                                          relementError, rcollection)
            
            ! Release element-wise error distribution
            call lsyssc_releaseVector(relementError)

            ! Generate standard mesh from raw mesh
            call tria_initStandardMeshFromRaw(p_rproblemLevel%rtriangulation, rproblem%rboundary)
            
            ! Update the template matrix according to the sparsity pattern
            templateMatrix = collct_getvalue_int(rcollection, 'templateMatrix')
            call grph_generateMatrix(rgraph, p_rproblemLevel%Rmatrix(templateMatrix))
            
            ! Re-initialize all constant coefficient matrices
            call euler_initProblemLevel(rappDescriptor, p_rproblemLevel, rcollection)

            ! Resize the solution vector accordingly
            systemMatrix = collct_getvalue_int(rcollection, 'systemMatrix')
            call lsysbl_resizeVecBlockIndMat(p_rproblemLevel%RmatrixBlock(systemMatrix),&
                                             rsolution, .false., .true.)

            ! Re-generate the initial solution vector and impose boundary conditions explicitly
            call euler_initSolution(rparlist, ssectionname, p_rproblemLevel, 0.0_DP, rsolution)
            select case(rappDescriptor%ndimension)
            case (NDIM1D)
              call bdrf_filterVectorExplicit(rbdrCond, p_rproblemLevel%rtriangulation,&
                                             rsolution, 0.0_DP, rproblem%rboundary,&
                                             euler_calcBoundaryvalues1d)
            case (NDIM2D)
              call bdrf_filterVectorExplicit(rbdrCond, p_rproblemLevel%rtriangulation,&
                                             rsolution, 0.0_DP, rproblem%rboundary,&
                                             euler_calcBoundaryvalues2d)
            case (NDIM3D)
              call bdrf_filterVectorExplicit(rbdrCond, p_rproblemLevel%rtriangulation,&
                                             rsolution, 0.0_DP, rproblem%rboundary,&
                                             euler_calcBoundaryvalues3d)
            end select
          end do

          ! Prepare internal data arrays of the solver structure
          systemMatrix = collct_getvalue_int(rcollection, 'systemMatrix')
          call flagship_updateSolverMatrix(p_rproblemLevel, rsolver, systemMatrix,&
                                           rappDescriptor%isystemFormat, UPDMAT_ALL)
          call solver_updateStructure(rsolver)

        end if   ! npreadapt > 0
        
      end if   ! dstepAdapt > 0
      
    else
      
      dstepAdapt = 0.0_DP
      
    end if
    
    ! Attach the boundary condition to the solver structure
    call solver_setBoundaryCondition(rsolver, rbdrCond, .true.)

    ! Set collection to primal problem mode
    call collct_setvalue_int(rcollection, 'primaldual', 1, .true.)

    ! Stop time measurement for pre-processing
    call stat_stopTimer(rappDescriptor%rtimerPrePostprocess)

    
    !---------------------------------------------------------------------------
    ! Infinite time stepping loop
    !---------------------------------------------------------------------------

    timeloop: do
      
      ! Check for user interaction
      if (signal_SIGINT(-1) > 0 )&
      call euler_outputSolution(rparlist, ssectionName, p_rproblemLevel,&
                                rsolution, dtime=rtimestep%dTime)
      
      !-------------------------------------------------------------------------
      ! Advance solution in time
      !-------------------------------------------------------------------------
      
      ! Start time measurement for solution procedure
      call stat_startTimer(rappDescriptor%rtimerSolution, STAT_TIMERSHORT)
      
      ! What time-stepping scheme should be used?
      select case(rtimestep%ctimestepType)
        
      case (TSTEP_RK_SCHEME)
        
        ! Adopt explicit Runge-Kutta scheme
        call tstep_performRKStep(p_rproblemLevel, rtimestep, rsolver,&
                                 rsolution, euler_nlsolverCallback, rcollection)
        
      case (TSTEP_THETA_SCHEME)
        
        ! Adopt two-level theta-scheme
        call tstep_performThetaStep(p_rproblemLevel, rtimestep, rsolver,&
                                    rsolution, euler_nlsolverCallback, rcollection)
        
      case DEFAULT
        call output_line('Unsupported time-stepping algorithm!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'euler_solveTransient')
        call sys_halt()
      end select
      
      ! Stop time measurement for solution procedure
      call stat_stopTimer(rappDescriptor%rtimerSolution)
      
      ! Reached final time, then exit the infinite time loop?
      if (rtimestep%dTime .ge. rtimestep%dfinalTime) exit timeloop


      !-------------------------------------------------------------------------
      ! Post-process intermediate solution
      !-------------------------------------------------------------------------
      
      if ((dstepUCD .gt. 0.0_DP) .and. (rtimestep%dTime .ge. dtimeUCD)) then

        ! Set time for next intermediate solution export
        dtimeUCD = dtimeUCD + dstepUCD

        ! Start time measurement for post-processing
        call stat_startTimer(rappDescriptor%rtimerPrepostProcess, STAT_TIMERSHORT)

        ! Export the intermediate solution
        call euler_outputSolution(rparlist, ssectionName, p_rproblemLevel,&
                                  rsolution, dtime=rtimestep%dTime)

        ! Stop time measurement for post-processing
        call stat_stopTimer(rappDescriptor%rtimerPrepostProcess)
        
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
        call stat_startTimer(rappDescriptor%rtimerErrorEstimation, STAT_TIMERSHORT)
        
        ! Compute the error estimator using recovery techniques
        call euler_estimateRecoveryError(rparlist, ssectionname, p_rproblemLevel,&
                                         rsolution, rtimestep%dTime, relementError, derror)
        
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
        call euler_adaptTriangulation(rhadapt, p_rproblemLevel%rtriangulation,&
                                       relementError, rcollection)
        
        ! Release element-wise error distribution
        call lsyssc_releaseVector(relementError)

        ! Update the template matrix according to the sparsity pattern
        call grph_generateMatrix(rgraph, p_rproblemLevel%Rmatrix(templateMatrix))
        
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
        call euler_initProblemLevel(rappDescriptor, p_rproblemLevel, rcollection)
        
        ! Resize the solution vector accordingly
        systemMatrix = collct_getvalue_int(rcollection, 'systemMatrix')
        call lsysbl_resizeVecBlockIndMat(p_rproblemLevel%RmatrixBlock(systemMatrix),&
                                         rsolution, .false., .true.)

        ! Prepare internal data arrays of the solver structure
        call flagship_updateSolverMatrix(p_rproblemLevel, rsolver, systemMatrix,&
                                         rappDescriptor%isystemFormat, UPDMAT_ALL)
        call solver_updateStructure(rsolver)
        
        ! Stop time measurement for generation of constant coefficient matrices
        call stat_stopTimer(rappDescriptor%rtimerAssemblyCoeff)

      end if

    end do timeloop

    
    ! Release adaptation structure
    if ((dstepAdapt > 0.0_DP) .or. (npreadapt > 0)) then
      call hadapt_releaseAdaptation(rhadapt)
      call grph_releaseGraph(rgraph)
    end if
    
  end subroutine euler_solveTransientPrimal


  !*****************************************************************************
  ! AUXILIARY ROUTINES
  !*****************************************************************************

!<subroutine>

  subroutine euler_parseCmdlArguments(rparlist)

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

  end subroutine euler_parseCmdlArguments

end module euler_application
