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
!#     -> The application's main routine called from the main problem
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
!# 6.) euler_initProblemLevels
!#     -> Initializes the individual problem levels based in the
!#        parameter settings given by the application descriptor
!#
!# 7.) euler_initSolution
!#     -> Initializes the solution vector based on the parameter
!#        settings given by the parameter list
!#
!# 8.) euler_outputSolution
!#     -> Outputs the solution vector to file in UCD format
!#
!# 11.) euler_outputStatistics
!#      -> Outputs the application statitics

!# 
!# </purpose>
!##############################################################################

module euler_application

  use afcstabilisation
  use boundaryfilter
  use collection
  use euler_basic
  use euler_callback
  use euler_callback1d
  use euler_callback2d
  use euler_callback3d
  use flagship_basic
  use fparser
  use linearformevaluation
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use problem
  use solver
  use statistics
  use stdoperators
  use storage
  use thermodynamics
  use ucd

  implicit none

  private
  public :: euler

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
    character(LEN=SYS_STRLEN) :: sinputoutputName
    character(LEN=SYS_STRLEN) :: sbenchmarkName
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
    call euler_parseCmdlArguments(rparlist)

    ! Initialize global collection structure
    call collct_init(rcollection)

    ! Initialize the application descriptor
    call euler_initApplication(rparlist, 'euler', rappDescriptor)

    ! Start time measurement for pre-processing
    call stat_startTimer(rappDescriptor%rtimerPrepostProcess, STAT_TIMERSHORT)

    ! Initialize the global collection
    call euler_initCollection(rappDescriptor, rcollection)

    ! Initialize the solver structures
    call euler_initSolvers(rparlist, 'euler', rtimestep, rsolver)

    ! Initialize the abstract problem structure
    call euler_initProblem(rparlist, 'euler',&
                           solver_getMinimumMultigridlevel(rsolver),&
                           solver_getMaximumMultigridlevel(rsolver),&
                           rproblem, rcollection)

    ! Initialize the individual problem levels
    call euler_initProblemLevels(rappDescriptor, rproblem, rcollection)

    ! Initialize the primal solution vector
    call euler_initSolution(rparlist, 'euler', rproblem%p_rproblemLevelMax,&
                            0.0_DP, rsolutionPrimal)

    ! Prepare internal data arrays of the solver structure
    call flagship_updateSolverMatrix(rproblem%p_rproblemLevelMax, rsolver,&
                                     1, rappDescriptor%isystemFormat, UPDMAT_ALL)
    call solver_updateStructure(rsolver)

    ! Stop time measurement for pre-processing
    call stat_stopTimer(rappDescriptor%rtimerPrePostprocess)

    
    !---------------------------------------------------------------------------
    ! Solution algorithm
    !---------------------------------------------------------------------------

    if (rtimestep%dfinalTime > 0) then
      
      ! Get global configuration from parameter list
      call parlst_getvalue_string(rparlist, 'euler', "benchmark", sbenchmarkName)
      call parlst_getvalue_string(rparlist, trim(sbenchmarkName), "algorithm", algorithm)
      call parlst_getvalue_string(rparlist, 'euler', "inputoutput", sinputoutputName)
      call parlst_getvalue_string(rparlist, trim(sinputoutputName), "indatfile", sindatfileName)
      
      ! The boundary condition for the primal problem is required for all 
      ! solution strategies so initialize it from the parameter file
      call parlst_getvalue_string(rparlist, trim(sinputoutputName),&
                                  "sprimalbdrcondname", sbdrcondName)
      call bdrf_readBoundaryCondition(rbdrCondPrimal, sindatfileName,&
                                      '['//trim(sbdrcondName)//']', rappDescriptor%ndimension)

      ! Impose primal boundary conditions explicitely
      select case(rappDescriptor%ndimension)
      case (NDIM1D)
        call bdrf_filterVectorExplicit(rbdrCondPrimal,&
                                       rproblem%p_rproblemLevelMax%rtriangulation,&
                                       rsolutionPrimal, 0.0_DP, rproblem%rboundary,&
                                       euler_calcBoundaryvalues1d)
      case (NDIM2D)
        call bdrf_filterVectorExplicit(rbdrCondPrimal,&
                                       rproblem%p_rproblemLevelMax%rtriangulation,&
                                       rsolutionPrimal, 0.0_DP, rproblem%rboundary,&
                                       euler_calcBoundaryvalues2d)
      case (NDIM3D)
        call bdrf_filterVectorExplicit(rbdrCondPrimal,&
                                       rproblem%p_rproblemLevelMax%rtriangulation,&
                                       rsolutionPrimal, 0.0_DP, rproblem%rboundary,&
                                       euler_calcBoundaryvalues3d)
      end select

      
      ! What solution algorithm should be applied?
      select case(trim(algorithm))

      
      case DEFAULT
        call output_line(trim(algorithm)//' is not a valid solution algorithm!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'euler_start')
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
    
    ! Release boundary conditions
    call bdrf_release(rbdrCondPrimal)
    call bdrf_release(rbdrCondDual)

    ! Release vectors
    call lsysbl_releaseVector(rsolutionPrimal)
    call lsysbl_releaseVector(rsolutionDual)

    ! Release function parser
    call fparser_done()

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
    type(t_euler), intent(OUT) :: rappDescriptor
!</output>
!</subroutine>

    ! section names
    character(LEN=SYS_STRLEN) :: sbenchmarkName
    character(LEN=SYS_STRLEN) :: sinputoutputName
    character(LEN=SYS_STRLEN) :: sindatfileName

    
    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist, ssectionName, "benchmark",   sbenchmarkName)
    call parlst_getvalue_string(rparlist, ssectionName, "inputoutput", sinputoutputName)
    call parlst_getvalue_string(rparlist, trim(sinputoutputName), "indatfile", sindatfileName)

    ! Initialize function parser
    call fparser_init()
    call fparser_parseFileForKeyword(sindatfileName, 'defconst', FPAR_CONSTANT)
    call fparser_parseFileForKeyword(sindatfileName, 'defexpr',  FPAR_EXPRESSION)

    ! Get application specifig parameters from the parameterlist
    call parlst_getvalue_int(rparlist, trim(sbenchmarkName),&
                             "ndimension", rappDescriptor%ndimension)
    call parlst_getvalue_int(rparlist, trim(sbenchmarkName),&
                             "imasstype", rappDescriptor%imasstype)
    call parlst_getvalue_int(rparlist, trim(sbenchmarkName),&
                             "imassantidiffusiontype", rappDescriptor%imassantidiffusiontype)
    call parlst_getvalue_int(rparlist, trim(sbenchmarkName),&
                             "idissipationtype", rappDescriptor%idissipationtype)
    call parlst_getvalue_int(rparlist, trim(sbenchmarkName),&
                             "icoupled", rappDescriptor%icoupled)
    call parlst_getvalue_int(rparlist, trim(sbenchmarkName),&
                             "iprecond", rappDescriptor%iprecond)
    call parlst_getvalue_int(rparlist, trim(sbenchmarkName),&
                             "ieltype", rappDescriptor%ieltype)
    call parlst_getvalue_int(rparlist, trim(sbenchmarkName),&
                             "imatrixformat", rappDescriptor%imatrixFormat)
    call parlst_getvalue_int(rparlist, trim(sbenchmarkName),&
                             "isystemformat", rappDescriptor%isystemFormat)

  end subroutine euler_initApplication

  !*****************************************************************************

!<subroutine>

  subroutine euler_initCollection(rappDescriptor, rcollection)

!<description>
    ! This subroutine initializes the collection based on
    ! the parameter settings of the application descriptor
!</description>

!<input>
    ! application descriptor
    type(t_euler), intent(IN) :: rappDescriptor
!</input>

!<inputoutput>
    ! collection structure
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! Add parameter settings from application descriptor
    call collct_setvalue_int(rcollection, 'ndimension',&
                             rappDescriptor%ndimension, .true.)
    call collct_setvalue_int(rcollection, 'imasstype',&
                             rappDescriptor%imasstype, .true.)
    call collct_setvalue_int(rcollection, 'imassantidiffusiontype',&
                             rappDescriptor%imassantidiffusiontype, .true.)
    call collct_setvalue_int(rcollection, 'idissipationtype',&
                             rappDescriptor%idissipationtype, .true.)
    call collct_setvalue_int(rcollection, 'icoupled',&
                             rappDescriptor%icoupled, .true.)
    call collct_setvalue_int(rcollection, 'iprecond',&
                             rappDescriptor%iprecond, .true.)
    call collct_setvalue_int(rcollection, 'ieltype',&
                             rappDescriptor%ieltype, .true.)
    call collct_setvalue_int(rcollection, 'imatrixformat',&
                             rappDescriptor%imatrixFormat, .true.)
    call collct_setvalue_int(rcollection, 'isystemformat',&
                             rappDescriptor%isystemFormat, .true.)

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
    call collct_setvalue_int(rcollection, 'inviscidAFC', 1, .true.)

    call collct_setvalue_int(rcollection, 'systemmatrix',    1, .true.)
    call collct_setvalue_int(rcollection, 'jacobianmatrix',  2, .true.)
    call collct_setvalue_int(rcollection, 'templatematrix',  3, .true.)

    if ((rappDescriptor%imasstype .ne. MASS_ZERO) .or.&
        (rappDescriptor%imassantidiffusiontype .ne. MASS_ZERO)) then
      call collct_setvalue_int(rcollection, 'consistentmassmatrix', 4, .true.)
      call collct_setvalue_int(rcollection, 'lumpedmassmatrix',     5, .true.)
    else
      call collct_setvalue_int(rcollection, 'consistentmassmatrix', 0, .true.)
      call collct_setvalue_int(rcollection, 'lumpedmassmatrix',     0, .true.)
    end if

    select case(rappDescriptor%ndimension)
    case (NDIM1D)
      call collct_setvalue_int(rcollection, 'coeffmatrix_cx', 6, .true.)
      call collct_setvalue_int(rcollection, 'coeffmatrix_cy', 0, .true.)
      call collct_setvalue_int(rcollection, 'coeffmatrix_cz', 0, .true.)
    case (NDIM2D)
      call collct_setvalue_int(rcollection, 'coeffmatrix_cx', 6, .true.)
      call collct_setvalue_int(rcollection, 'coeffmatrix_cy', 7, .true.)
      call collct_setvalue_int(rcollection, 'coeffmatrix_cz', 0, .true.)
    case (NDIM3D)
      call collct_setvalue_int(rcollection, 'coeffmatrix_cx', 6, .true.)
      call collct_setvalue_int(rcollection, 'coeffmatrix_cy', 7, .true.)
      call collct_setvalue_int(rcollection, 'coeffmatrix_cz', 8, .true.)
    case DEFAULT
      call collct_setvalue_int(rcollection, 'coeffmatrix_cx', 0, .true.)
      call collct_setvalue_int(rcollection, 'coeffmatrix_cy', 0, .true.)
      call collct_setvalue_int(rcollection, 'coeffmatrix_cz', 0, .true.)
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
    call parlst_getvalue_string(rparlist, ssectionName, "timestep", stimestepName)
    call parlst_getvalue_string(rparlist, ssectionName, "solver",   ssolverName)

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
    character(LEN=SYS_STRLEN) :: sbenchmarkName
    character(LEN=SYS_STRLEN) :: sinputoutputName
    character(LEN=SYS_STRLEN) :: sinviscidName

    ! abstract problem descriptor
    type(t_problemDescriptor) :: rproblemDescriptor

    ! pointer to the problem level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! local variables
    integer :: inviscidAFC
    integer :: iconvToTria
    

    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist, ssectionName, "benchmark", sbenchmarkName)
    call parlst_getvalue_string(rparlist, ssectionName, "inputoutput", sinputoutputName)
    call parlst_getvalue_string(rparlist, ssectionName, "primalinviscid",  sinviscidName)
    
    call parlst_getvalue_int(rparlist, trim(adjustl(sbenchmarkName)),&
                             "ndimension", rproblemDescriptor%ndimension)
    call parlst_getvalue_string(rparlist, trim(adjustl(sinputoutputName)),&
                                "trifile", rproblemDescriptor%trifile, '')
    call parlst_getvalue_string(rparlist, trim(adjustl(sinputoutputName)),&
                                "prmfile", rproblemDescriptor%prmfile, '')
    
    ! Set additional problem descriptor
    rproblemDescriptor%nafcstab      = 1   ! for inviscid fluxes
    rproblemDescriptor%nlmin         = nlmin
    rproblemDescriptor%nlmax         = nlmax
    rproblemDescriptor%nmatrixScalar = rproblemDescriptor%ndimension + 5
    rproblemDescriptor%nmatrixBlock  = 2
    rproblemDescriptor%nvectorScalar = 0
    rproblemDescriptor%nvectorBlock  = 0

    ! Check if quadrilaterals should be converted to triangles
    call parlst_getvalue_int(rparlist, trim(adjustl(sbenchmarkName)),&
                             "iconvtotria", iconvToTria, 0)
    if (iconvToTria .eq. 1)&
        rproblemDescriptor%iproblemSpec = rproblemDescriptor%iproblemSpec &
                                        + PROBDESC_MSPEC_CONVTRIANGLES   

    ! Initialize problem structure
    call problem_initProblem(rproblemDescriptor, rproblem)

    
    ! Initialize the stabilisation structure
    inviscidAFC = collct_getvalue_int(rcollection, 'inviscidAFC')
    
    ! loop over all problem levels
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

  subroutine euler_initProblemLevels(rappDescriptor, rproblem, rcollection)

!<description>
    ! This subroutine initielizes the individual problem levels, that is,
    ! it generates the discretization, the template matrix and the required
    ! coefficient matrices as duplicates of the template matrix.
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
!</output>
!</subroutine>

    ! pointer to the problem level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! local variables
    integer :: templateMatrix
    integer :: systemMatrix
    integer :: jacobianMatrix
    integer :: consistentMassMatrix
    integer :: lumpedMassMatrix
    integer :: coeffMatrix_CX
    integer :: coeffMatrix_CY
    integer :: coeffMatrix_CZ
    integer :: ivar,jvar
    

    ! retrieve application specific parameters from the collection
    templateMatrix       = collct_getvalue_int(rcollection, 'templatematrix')
    systemMatrix         = collct_getvalue_int(rcollection, 'systemmatrix')
    jacobianMatrix       = collct_getvalue_int(rcollection, 'jacobianmatrix')
    consistentMassMatrix = collct_getvalue_int(rcollection, 'consistentmassmatrix')
    lumpedMassMatrix     = collct_getvalue_int(rcollection, 'lumpedmassmatrix')
    coeffMatrix_CX       = collct_getvalue_int(rcollection, 'coeffmatrix_cx')
    coeffMatrix_CY       = collct_getvalue_int(rcollection, 'coeffmatrix_cy')
    coeffMatrix_CZ       = collct_getvalue_int(rcollection, 'coeffmatrix_cz')

    ! loop over all problem levels
    p_rproblemLevel => rproblem%p_rproblemLevelMax
    do while(associated(p_rproblemLevel))

      ! Initialize the discretization structure
      select case(rappDescriptor%isystemformat)
      case (SYSTEM_INTERLEAVEFORMAT)
        call spdiscr_initBlockDiscr(p_rproblemLevel%rdiscretisation, 1,&
                                    p_rproblemLevel%rtriangulation)
      case (SYSTEM_BLOCKFORMAT)
        call spdiscr_initBlockDiscr(p_rproblemLevel%rdiscretisation,&
                                    euler_getNVAR(rappDescriptor),&
                                    p_rproblemLevel%rtriangulation)
      case DEFAULT
        call output_line('Unsupported system format!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'euler_initProblemLevels')
        call sys_halt()
      end select

      ! Get spatial dimension
      select case(p_rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        select case(rappDescriptor%ieltype)
        case (-1,1,11)
          ! P1=Q1 finite elements
          call spdiscr_initDiscr_simple(p_rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                        EL_E001_1D, SPDISC_CUB_AUTOMATIC,&
                                        p_rproblemLevel%rtriangulation,&
                                        p_rproblemLevel%p_rproblem%rboundary)
          ! P2=Q2 finite elements
          call spdiscr_initDiscr_simple(p_rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                        EL_E002_1D, SPDISC_CUB_AUTOMATIC,&
                                        p_rproblemLevel%rtriangulation,&
                                        p_rproblemLevel%p_rproblem%rboundary)
        case DEFAULT
          call output_line('Unsupproted element type!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'euler_initProblemLevels')
          call sys_halt()
        end select

      case (NDIM2D)
        select case(rappDescriptor%ieltype)
        case (1)
          ! P1 finite elements
          call spdiscr_initDiscr_simple(p_rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                        EL_E001, SPDISC_CUB_AUTOMATIC,&
                                        p_rproblemLevel%rtriangulation,&
                                        p_rproblemLevel%p_rproblem%rboundary)
        case (2)
          ! P2 finite elements
          call spdiscr_initDiscr_simple(p_rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                        EL_E002, SPDISC_CUB_AUTOMATIC,&
                                        p_rproblemLevel%rtriangulation,&
                                        p_rproblemLevel%p_rproblem%rboundary)
        case (11)
          ! Q1 finite elements
          call spdiscr_initDiscr_simple(p_rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                        EL_E011, SPDISC_CUB_AUTOMATIC,&
                                        p_rproblemLevel%rtriangulation,&
                                        p_rproblemLevel%p_rproblem%rboundary)
        case (12)
          ! Q2 finite elements
          call spdiscr_initDiscr_simple(p_rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                        EL_E013, SPDISC_CUB_AUTOMATIC,&
                                        p_rproblemLevel%rtriangulation,&
                                        p_rproblemLevel%p_rproblem%rboundary)
        case (-1)
          ! mixed P1/Q1 finite elements
          call spdiscr_initDiscr_triquad(p_rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                         EL_E001, EL_E011, SPDISC_CUB_AUTOMATIC,&
                                         SPDISC_CUB_AUTOMATIC,&
                                         p_rproblemLevel%rtriangulation,&
                                         p_rproblemLevel%p_rproblem%rboundary)
        case (-2)
          ! mixed P2/Q2 finite elements
          call spdiscr_initDiscr_triquad(p_rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                         EL_E002, EL_E013, SPDISC_CUB_AUTOMATIC,&
                                         SPDISC_CUB_AUTOMATIC,&
                                         p_rproblemLevel%rtriangulation,&
                                         p_rproblemLevel%p_rproblem%rboundary)
        case DEFAULT
          call output_line('Unsupproted element type!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'euler_initProblemLevels')
          call sys_halt()
        end select

      case (NDIM3D)
        select case(rappDescriptor%ieltype)
        case (1)
          ! P1 finite elements
          call spdiscr_initDiscr_simple(p_rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                        EL_E001_3D, SPDISC_CUB_AUTOMATIC,&
                                        p_rproblemLevel%rtriangulation,&
                                        p_rproblemLevel%p_rproblem%rboundary)
        case (11)
          ! Q1 finite elements
          call spdiscr_initDiscr_simple(p_rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                        EL_E010_3D, SPDISC_CUB_AUTOMATIC,&
                                        p_rproblemLevel%rtriangulation,&
                                        p_rproblemLevel%p_rproblem%rboundary)
        case DEFAULT
          call output_line('Unsupproted element type!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'euler_initProblemLevels')
          call sys_halt()
        end select
        
      case DEFAULT
        call output_line('Invalid number of spatial dimensions',&
                         OU_CLASS_ERROR,OU_MODE_STD,'euler_initProblemLevels')
        call sys_halt()
      end select

      ! Duplicate scalar discretisation structure for block matrix format
      if (rappDescriptor%isystemFormat .eq. SYSTEM_BLOCKFORMAT) then
        do ivar = 2, euler_getNVAR(rappDescriptor)
          call spdiscr_duplicateDiscrSc(&
              p_rproblemLevel%rdiscretisation%RspatialDiscr(1),&
              p_rproblemLevel%rdiscretisation%RspatialDiscr(ivar), .true.)
        end do
      end if
      
      ! Generate the finite element matrix sparsity structure based on the
      ! spatial descretization and store it as the template matrix
      call bilf_createMatrixStructure(p_rproblemLevel%rdiscretisation%RspatialDiscr(1),&
                                      rappDescriptor%imatrixFormat,&
                                      p_rproblemLevel%Rmatrix(templateMatrix))

      ! Create consistent (and lumped) mass matrix as duplicate of the template matrix
      if (consistentMassMatrix > 0) then
        call lsyssc_duplicateMatrix(p_rproblemLevel%Rmatrix(templateMatrix),&
                                    p_rproblemLevel%Rmatrix(consistentMassMatrix),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
        call stdop_assembleSimpleMatrix(p_rproblemLevel%Rmatrix(consistentMassMatrix),&
                                        DER_FUNC, DER_FUNC) 
        if (lumpedMassMatrix > 0) then
          call lsyssc_duplicateMatrix(p_rproblemLevel%Rmatrix(consistentMassMatrix),&
                                      p_rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                      LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)
          call lsyssc_lumpMatrixScalar(p_rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                       LSYSSC_LUMP_DIAG)
        end if
      elseif (lumpedMassMatrix > 0) then
        call lsyssc_duplicateMatrix(p_rproblemLevel%Rmatrix(templateMatrix),&
                                    p_rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
        call stdop_assembleSimpleMatrix(p_rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                        DER_FUNC, DER_FUNC) 
        call lsyssc_lumpMatrixScalar(p_rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                     LSYSSC_LUMP_DIAG)
      end if

      ! Create coefficient matrix (phi, dphi/dx) duplicate of the template matrix
      if (coeffMatrix_CX > 0) then
        call lsyssc_duplicateMatrix(p_rproblemLevel%Rmatrix(templateMatrix),&
                                    p_rproblemLevel%Rmatrix(coeffMatrix_CX),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
        call stdop_assembleSimpleMatrix(p_rproblemLevel%Rmatrix(coeffMatrix_CX),&
                                        DER_DERIV3D_X, DER_FUNC)
      end if

      ! Create coefficient matrix (phi, dphi/dy) duplicate of the template matrix
      if (coeffMatrix_CY > 0) then
        call lsyssc_duplicateMatrix(p_rproblemLevel%Rmatrix(templateMatrix),&
                                    p_rproblemLevel%Rmatrix(coeffMatrix_CY),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
        call stdop_assembleSimpleMatrix(p_rproblemLevel%Rmatrix(coeffMatrix_CY),&
                                        DER_DERIV3D_Y, DER_FUNC)
      end if

      ! Create coefficient matrix (phi, dphi/dz) duplicate of the template matrix
      if (coeffMatrix_CZ > 0) then
        call lsyssc_duplicateMatrix(p_rproblemLevel%Rmatrix(templateMatrix),&
                                    p_rproblemLevel%Rmatrix(coeffMatrix_CZ),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
        call stdop_assembleSimpleMatrix(p_rproblemLevel%Rmatrix(coeffMatrix_CZ),&
                                        DER_DERIV3D_Z, DER_FUNC)
      end if      

      ! Create system matrix
      if (systemMatrix > 0) then
        select case(rappDescriptor%isystemFormat)
        case (SYSTEM_INTERLEAVEFORMAT)
          ! The global operator is stored as an interleave matrix with
          ! NVAR components. However, the row and column structure of
          ! the template matrix can be adopted
          call lsyssc_duplicateMatrix(p_rproblemLevel%Rmatrix(templateMatrix),&
                                      p_rproblemLevel%Rmatrix(systemMatrix),&
                                      LSYSSC_DUP_SHARE, LSYSSC_DUP_REMOVE)
          p_rproblemLevel%Rmatrix(systemMatrix)%NVAR = euler_getNVAR(rappDescriptor)

          ! What matrix format should be used?
          select case(rappDescriptor%imatrixFormat)
          case (LSYSSC_MATRIX7)
            p_rproblemLevel%Rmatrix(systemMatrix)%cmatrixFormat = LSYSSC_MATRIX7INTL
            
          case (LSYSSC_MATRIX9)
            p_rproblemLevel%Rmatrix(systemMatrix)%cmatrixFormat = LSYSSC_MATRIX9INTL
            
          case DEFAULT
            call output_line('Unsupported matrix format!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'euler_initProblemLevels')
            call sys_halt()
          end select

          ! What kind of global operator should be adopted?
          select case(rappDescriptor%icoupled)
          case (FLOW_SEGREGATED)
            p_rproblemLevel%Rmatrix(systemMatrix)%cinterleavematrixFormat = LSYSSC_MATRIXD
            
          case (FLOW_ALLCOUPLED)
            p_rproblemLevel%Rmatrix(systemMatrix)%cinterleavematrixFormat = LSYSSC_MATRIX1
            
          case DEFAULT
            call output_line('Unsupported interleave matrix format!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'euler_initProblemLevels')
            call sys_halt()
          end select
          
          ! Create global operator physically
          call lsyssc_allocEmptyMatrix(p_rproblemLevel%Rmatrix(systemMatrix),&
                                       LSYSSC_SETM_UNDEFINED)
          
          ! Create pseudo block matrix from global operator
          call lsysbl_createMatFromScalar(p_rproblemLevel%Rmatrix(systemMatrix),&
                                          p_rproblemLevel%RmatrixBlock(systemMatrix),&
                                          p_rproblemLevel%rdiscretisation)

        case (SYSTEM_BLOCKFORMAT)
          ! The global operator is stored as a block matrix with NVARxNVAR
          ! blocks. Thus, create empty NVARxNVAR block matrix directly
          call lsysbl_createEmptyMatrix(p_rproblemLevel%RmatrixBlock(systemMatrix),&
                                        euler_getNVAR(rappDescriptor))
          p_rproblemLevel%RmatrixBlock(systemMatrix)%imatrixSpec = LSYSBS_MSPEC_GROUPMATRIX

          ! What kind of global operator should be adopted?
          select case(rappDescriptor%icoupled)
          case (FLOW_SEGREGATED)
            ! Create only NVAR diagonal blocks
            do ivar = 1, euler_getNVAR(rappDescriptor)
              call lsyssc_duplicateMatrix(p_rproblemLevel%Rmatrix(templateMatrix),&
                                          p_rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                                          LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
            end do

          case (FLOW_ALLCOUPLED)
            ! Create all NVAR x NVAR blocks
            do ivar = 1, euler_getNVAR(rappDescriptor)
              do jvar = 1, euler_getNVAR(rappDescriptor)
                call lsyssc_duplicateMatrix(p_rproblemLevel%Rmatrix(templateMatrix),&
                                            p_rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,jvar),&
                                            LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
              end do
            end do
          
          case DEFAULT
            call output_line('Unsupported block matrix format!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'euler_initProblemLevels')
            call sys_halt()
          end select

          ! Update internal structure of block matrix
          call lsysbl_updateMatStrucInfo(p_rproblemLevel%RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Unsupported system format!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'euler_initProblemLevels')
          call sys_halt()
        end select
      end if     

      ! Switch to next coarser level
      p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
    end do

  end subroutine euler_initProblemLevels

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
    character(LEN=SYS_STRLEN) :: sbenchmarkName
    character(LEN=SYS_STRLEN) :: sinputoutputName
    character(LEN=SYS_STRLEN) :: sindatfileName
    character(LEN=SYS_STRLEN) :: ssolutionName

    ! local variables
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    type(t_fparser) :: rfparser
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Ddata
    real(DP), dimension(NDIM3D+1) :: Dvalue
    integer :: isolutiontype,isystemformat
    integer :: isize,iblock,ivar,nvar,ieq,neq,ndim

    ! symbolic variable names
    character(LEN=*), dimension(4), parameter ::&
                      cvariables = (/ (/'x'/), (/'y'/), (/'z'/), (/'t'/) /)

    
    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist, ssectionName, "benchmark",   sbenchmarkName)
    call parlst_getvalue_string(rparlist, ssectionName, "inputoutput", sinputoutputName)
    call parlst_getvalue_string(rparlist, trim(sinputoutputName), "indatfile", sindatfileName)
    call parlst_getvalue_string(rparlist, trim(sinputoutputName), "ssolutionname", ssolutionName)
    call parlst_getvalue_int(rparlist, trim(sbenchmarkName), "isolutiontype", isolutiontype)
    call parlst_getvalue_int(rparlist, trim(sbenchmarkName), "isystemformat", isystemformat)
    
    ! Create new solution vector based on the spatial discretisation
    p_rdiscretisation => rproblemLevel%rdiscretisation
    
    select case(isystemformat)
    case (SYSTEM_INTERLEAVEFORMAT)
      call lsysbl_createVectorBlock(p_rdiscretisation, rvector, .false., ST_DOUBLE)
      
      rvector%RvectorBlock(1)%NVAR = euler_getNVAR(rproblemLevel)
      isize = rvector%NEQ*euler_getNVAR(rproblemLevel)
      call lsysbl_resizeVectorBlock(rvector, isize, .false., .false.)
      
    case (SYSTEM_BLOCKFORMAT)
      call lsysbl_createVectorBlock(p_rdiscretisation, rvector, .false., ST_DOUBLE)

    case DEFAULT
      call output_line('Unsupported system format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'euler_initSolution')
      call sys_halt()
    end select
      
    
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
      
      select case(isystemformat)
      case (SYSTEM_INTERLEAVEFORMAT)
        ! Get number of equations of first scalar subvector
        neq  = rvector%RvectorBlock(1)%NEQ
        nvar = rvector%RvectorBlock(1)%NVAR
        call lsyssc_getbase_double(rvector%RvectorBlock(1), p_Ddata)

        ! Loop over all equations of scalar subvector
        do ieq = 1, neq
          Dvalue(1:ndim)   = p_DvertexCoords(:,ieq)
          Dvalue(NDIM3D+1) = dtime
          
          ! Loop over all variables of the solution vector
          do ivar = 1, nvar
            call fparser_evalFunction(rfparser, ivar, Dvalue, p_Ddata((ieq-1)*nvar+ivar))
          end do
        end do
        
      case (SYSTEM_BLOCKFORMAT)
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
        
      end select
      
      call fparser_release(rfparser)


    case DEFAULT
      call output_line('Invalid type of solution profile!',&
                       OU_CLASS_ERROR, OU_MODE_STD, 'euler_initSolution')
      call sys_halt()
    end select

  end subroutine euler_initSolution

  !*****************************************************************************

!<subroutine>

  subroutine euler_outputSolution(rparlist, ssectionName, rproblem,&
                                  rsolutionPrimal, rsolutionDual, dtime)

!<description>
    ! This subroutine exports the solution vector to file in UCD format
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(IN) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(IN) :: ssectionName

    ! problem  structure
    type(t_problem), intent(IN) :: rproblem

    ! solution vector for primal problem
    type(t_vectorBlock), intent(IN) :: rsolutionPrimal

    ! OPTIONAL: solution vector for dual problem
    type(t_vectorBlock), intent(IN), optional :: rsolutionDual

    ! OPTIONAL: simulation time
    real(DP), intent(IN), optional :: dtime
!</input>
!</subroutine>

    ! section names
    character(LEN=SYS_STRLEN) :: sbenchmarkName
    character(LEN=SYS_STRLEN) :: sinputoutputName
    character(LEN=SYS_STRLEN) :: ucdsolution

    ! persistent variable
    integer, save :: ifilenumber = 1
    
    ! local variables
    type(t_ucdExport) :: rexport
    type(t_vectorScalar) :: rvector1,rvector2,rvector3
    real(DP), dimension(:), pointer :: p_Dsolution,p_Ddata1,p_Ddata2,p_Ddata3
    integer :: iformatUCD,isystemFormat,isize,ndim


    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist, ssectionName, "inputoutput", sinputoutputName)
    call parlst_getvalue_string(rparlist, ssectionName, "benchmark", sbenchmarkName)
    call parlst_getvalue_string(rparlist, trim(sinputoutputName), "ucdsolution", ucdsolution)
    call parlst_getvalue_int(rparlist, trim(sinputoutputName), "iformatucd", iformatUCD)
    call parlst_getvalue_int(rparlist, trim(sbenchmarkName), "isystemformat", isystemformat)

    ! Initialize the UCD exporter
    call flagship_initUCDexport(rproblem%p_rproblemLevelMax, ucdsolution,&
                                iformatUCD, rexport, ifilenumber)

    ! Increase filenumber by one
    ifilenumber = ifilenumber+1

    ! Set simulation time
    if (present(dtime)) call ucd_setSimulationTime(rexport, dtime)

    ! Set pointers
    call lsysbl_getbase_double(rsolutionPrimal, p_Dsolution)
    isize = size(p_Dsolution)/euler_getNVAR(rproblem%p_rproblemLevelMax)
    ndim  = rproblem%p_rproblemLevelMax%rdiscretisation%ndimension

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


    select case(isystemformat)
    case(SYSTEM_INTERLEAVEFORMAT)
      
      select case(ndim)
      case (NDIM1D)
        call getVarInterleaveformat(rvector1%NEQ, NVAR1D, 'velocity_x', p_Dsolution, p_Ddata1)
        call ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1)

        call getVarInterleaveformat(rvector1%NEQ, NVAR1D, 'density', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'density', UCD_VAR_STANDARD, p_Ddata1)

        call getVarInterleaveformat(rvector1%NEQ, NVAR1D, 'energy', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'energy', UCD_VAR_STANDARD, p_Ddata1)
        
        call getVarInterleaveformat(rvector1%NEQ, NVAR1D, 'pressure', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'pressure', UCD_VAR_STANDARD, p_Ddata1)

        call getVarInterleaveformat(rvector1%NEQ, NVAR1D, 'machnumber', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'machnumber', UCD_VAR_STANDARD, p_Ddata1)
        
        
      case (NDIM2D)
        call getVarInterleaveformat(rvector1%NEQ, NVAR2D, 'velocity_x', p_Dsolution, p_Ddata1)
        call getVarInterleaveformat(rvector2%NEQ, NVAR2D, 'velocity_y', p_Dsolution, p_Ddata2)
        call ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1, p_Ddata2)

        call getVarInterleaveformat(rvector1%NEQ, NVAR2D, 'density', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'density', UCD_VAR_STANDARD, p_Ddata1)

        call getVarInterleaveformat(rvector1%NEQ, NVAR2D, 'energy', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'energy', UCD_VAR_STANDARD, p_Ddata1)
        
        call getVarInterleaveformat(rvector1%NEQ, NVAR2D, 'pressure', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'pressure', UCD_VAR_STANDARD, p_Ddata1)

        call getVarInterleaveformat(rvector1%NEQ, NVAR2D, 'machnumber', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'machnumber', UCD_VAR_STANDARD, p_Ddata1)

      case (NDIM3D)
        call getVarInterleaveformat(rvector1%NEQ, NVAR3D, 'velocity_x', p_Dsolution, p_Ddata1)
        call getVarInterleaveformat(rvector2%NEQ, NVAR3D, 'velocity_y', p_Dsolution, p_Ddata2)
        call getVarInterleaveformat(rvector3%NEQ, NVAR3D, 'velocity_z', p_Dsolution, p_Ddata3)
        call ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1, p_Ddata2, p_Ddata3)
        
        call getVarInterleaveformat(rvector1%NEQ, NVAR3D, 'density', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'density', UCD_VAR_STANDARD, p_Ddata1)

        call getVarInterleaveformat(rvector1%NEQ, NVAR3D, 'energy', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'energy', UCD_VAR_STANDARD, p_Ddata1)
        
        call getVarInterleaveformat(rvector1%NEQ, NVAR3D, 'pressure', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'pressure', UCD_VAR_STANDARD, p_Ddata1)

        call getVarInterleaveformat(rvector1%NEQ, NVAR3D, 'machnumber', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'machnumber', UCD_VAR_STANDARD, p_Ddata1)

      end select
      
      
    case (SYSTEM_BLOCKFORMAT)

      select case(ndim)
      case (NDIM1D)
        call getVarBlockformat(rvector1%NEQ, NVAR1D, 'velocity_x', p_Dsolution, p_Ddata1)
        call ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1)

        call getVarBlockformat(rvector1%NEQ, NVAR1D, 'density', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'density', UCD_VAR_STANDARD, p_Ddata1)

        call getVarBlockformat(rvector1%NEQ, NVAR1D, 'energy', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'energy', UCD_VAR_STANDARD, p_Ddata1)
        
        call getVarBlockformat(rvector1%NEQ, NVAR1D, 'pressure', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'pressure', UCD_VAR_STANDARD, p_Ddata1)

        call getVarBlockformat(rvector1%NEQ, NVAR1D, 'machnumber', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'machnumber', UCD_VAR_STANDARD, p_Ddata1)
        
        
      case (NDIM2D)
        call getVarBlockformat(rvector1%NEQ, NVAR2D, 'velocity_x', p_Dsolution, p_Ddata1)
        call getVarBlockformat(rvector2%NEQ, NVAR2D, 'velocity_y', p_Dsolution, p_Ddata2)
        call ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1, p_Ddata2)

        call getVarBlockformat(rvector1%NEQ, NVAR2D, 'density', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'density', UCD_VAR_STANDARD, p_Ddata1)

        call getVarBlockformat(rvector1%NEQ, NVAR2D, 'energy', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'energy', UCD_VAR_STANDARD, p_Ddata1)
        
        call getVarBlockformat(rvector1%NEQ, NVAR2D, 'pressure', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'pressure', UCD_VAR_STANDARD, p_Ddata1)

        call getVarBlockformat(rvector1%NEQ, NVAR2D, 'machnumber', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'machnumber', UCD_VAR_STANDARD, p_Ddata1)

      case (NDIM3D)
        call getVarBlockformat(rvector1%NEQ, NVAR3D, 'velocity_x', p_Dsolution, p_Ddata1)
        call getVarBlockformat(rvector2%NEQ, NVAR3D, 'velocity_y', p_Dsolution, p_Ddata2)
        call getVarBlockformat(rvector3%NEQ, NVAR3D, 'velocity_z', p_Dsolution, p_Ddata3)
        call ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1, p_Ddata2, p_Ddata3)
        
        call getVarBlockformat(rvector1%NEQ, NVAR3D, 'density', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'density', UCD_VAR_STANDARD, p_Ddata1)

        call getVarBlockformat(rvector1%NEQ, NVAR3D, 'energy', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'energy', UCD_VAR_STANDARD, p_Ddata1)
        
        call getVarBlockformat(rvector1%NEQ, NVAR3D, 'pressure', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'pressure', UCD_VAR_STANDARD, p_Ddata1)

        call getVarBlockformat(rvector1%NEQ, NVAR3D, 'machnumber', p_Dsolution, p_Ddata1)
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

  contains

    ! Here, the working routines follow
    
    !*************************************************************
    ! Extract the variables from the global solution 
    ! vector which is stored in interleave format

    pure subroutine getVarInterleaveformat(neq, nvar, cvariable, Ddata, Dvalue)
      
      integer, intent(IN) :: neq,nvar
      character(LEN=*), intent(IN) :: cvariable
      real(DP), dimension(nvar,neq), intent(IN) :: Ddata
      real(DP), dimension(:), intent(OUT) :: Dvalue

      ! local variables
      real(DP) :: p
      integer :: ieq,ivar

      select case (cvariable)
      case ('density')
        do ieq = 1, neq
          Dvalue(ieq) = Ddata(1, ieq)
        end do

      case ('velocity_x')
        do ieq = 1, neq
          Dvalue(ieq) = Ddata(2, ieq)/Ddata(1, ieq)
        end do

      case ('velocity_y')
        do ieq = 1, neq
          Dvalue(ieq) = Ddata(3, ieq)/Ddata(1, ieq)
        end do

      case ('velocity_z')
        do ieq = 1, neq
          Dvalue(ieq) = Ddata(4, ieq)/Ddata(1, ieq)
        end do

      case ('energy')
        do ieq = 1, neq
          Dvalue(ieq) = Ddata(nvar, ieq)/Ddata(1, ieq)
        end do
        
      case ('pressure')
        select case (nvar)
        case (NVAR1D)
          do ieq = 1, neq
            Dvalue(ieq) = thdyn_pressure(GAMMA_AIR,&
                Ddata(nvar, ieq)/Ddata(1, ieq), Ddata(1, ieq),&
                Ddata(2, ieq)/Ddata(1, ieq))
          end do

        case (NVAR2D)
          do ieq = 1, neq
            Dvalue(ieq) = thdyn_pressure(GAMMA_AIR,&
                Ddata(nvar, ieq)/Ddata(1, ieq), Ddata(1, ieq),&
                Ddata(2, ieq)/Ddata(1, ieq),&
                Ddata(3, ieq)/Ddata(1, ieq))
          end do

        case (NVAR3D)
          do ieq = 1, neq
            Dvalue(ieq) = thdyn_pressure(GAMMA_AIR,&
                Ddata(nvar, ieq)/Ddata(1, ieq), Ddata(1, ieq),&
                Ddata(2, ieq)/Ddata(1, ieq),&
                Ddata(3, ieq)/Ddata(1, ieq),&
                Ddata(4, ieq)/Ddata(1, ieq))
          end do
        end select

      case ('machnumber')
        select case (nvar)
        case (NVAR1D)
          do ieq = 1, neq
            p = thdyn_pressure(GAMMA_AIR,&
                Ddata(nvar, ieq)/Ddata(1, ieq), Ddata(1, ieq),&
                Ddata(2, ieq)/Ddata(1, ieq))

            Dvalue(ieq) = thdyn_Machnumber(GAMMA_AIR,&
                p, Ddata(1, ieq),&
                Ddata(2, ieq)/Ddata(1, ieq))
          end do

        case (NVAR2D)
          do ieq = 1, neq
            p = thdyn_pressure(GAMMA_AIR,&
                Ddata(nvar, ieq)/Ddata(1, ieq), Ddata(1, ieq),&
                Ddata(2, ieq)/Ddata(1, ieq),&
                Ddata(3, ieq)/Ddata(1, ieq))

            Dvalue(ieq) = thdyn_Machnumber(GAMMA_AIR,&
                p, Ddata(1, ieq),&
                Ddata(2, ieq)/Ddata(1, ieq),&
                Ddata(3, ieq)/Ddata(1, ieq))
          end do
          
        case (NVAR3D)
          do ieq = 1, neq
            p = thdyn_pressure(GAMMA_AIR,&
                Ddata(nvar, ieq)/Ddata(1, ieq), Ddata(1, ieq),&
                Ddata(2, ieq)/Ddata(1, ieq),&
                Ddata(3, ieq)/Ddata(1, ieq),&
                Ddata(4, ieq)/Ddata(1, ieq))

            Dvalue(ieq) = thdyn_Machnumber(GAMMA_AIR,&
                p, Ddata(1, ieq),&
                Ddata(2, ieq)/Ddata(1, ieq),&
                Ddata(3, ieq)/Ddata(1, ieq),&
                Ddata(4, ieq)/Ddata(1, ieq))
          end do
        end select
      end select
      
    end subroutine getVarInterleaveformat

    !*************************************************************
    ! Extract the variables from the global solution 
    ! vector which is stored in block format

    pure subroutine getVarBlockformat(neq, nvar, cvariable, Ddata, Dvalue)
      
      integer, intent(IN) :: neq,nvar
      character(LEN=*), intent(IN) :: cvariable
      real(DP), dimension(neq,nvar), intent(IN) :: Ddata
      real(DP), dimension(:), intent(OUT) :: Dvalue

      ! local variables
      real(DP) :: p
      integer :: ieq,ivar

      select case (cvariable)
      case ('density')
        do ieq = 1, neq
          Dvalue(ieq) = Ddata(ieq, 1)
        end do

      case ('velocity_x')
        do ieq = 1, neq
          Dvalue(ieq) = Ddata(ieq, 2)/Ddata(ieq, 1)
        end do

      case ('velocity_y')
        do ieq = 1, neq
          Dvalue(ieq) = Ddata(ieq, 3)/Ddata(ieq, 1)
        end do

      case ('velocity_z')
        do ieq = 1, neq
          Dvalue(ieq) = Ddata(ieq, 4)/Ddata(ieq, 1)
        end do

      case ('energy')
        do ieq = 1, neq
          Dvalue(ieq) = Ddata(ieq, nvar)/Ddata(ieq, 1)
        end do
        
      case ('pressure')
        select case (nvar)
        case (NVAR1D)
          do ieq = 1, neq
            Dvalue(ieq) = thdyn_pressure(GAMMA_AIR,&
                Ddata(ieq, nvar)/Ddata(ieq, 1), Ddata(ieq, 1),&
                Ddata(ieq, 2)/Ddata(ieq, 1))
          end do

        case (NVAR2D)
          do ieq = 1, neq
            Dvalue(ieq) = thdyn_pressure(GAMMA_AIR,&
                Ddata(ieq, nvar)/Ddata(ieq, 1), Ddata(ieq, 1),&
                Ddata(ieq, 2)/Ddata(ieq, 1),&
                Ddata(ieq, 3)/Ddata(ieq, 1))
          end do

        case (NVAR3D)
          do ieq = 1, neq
            Dvalue(ieq) = thdyn_pressure(GAMMA_AIR,&
                Ddata(ieq, nvar)/Ddata(ieq, 1), Ddata(ieq, 1),&
                Ddata(ieq, 2)/Ddata(ieq, 1),&
                Ddata(ieq, 3)/Ddata(ieq, 1),&
                Ddata(ieq, 4)/Ddata(ieq, 1))
          end do
        end select

      case ('machnumber')
        select case (nvar)
        case (NVAR1D)
          do ieq = 1, neq
            p = thdyn_pressure(GAMMA_AIR,&
                Ddata(ieq, nvar)/Ddata(ieq, 1), Ddata(ieq, 1),&
                Ddata(ieq, 2)/Ddata(ieq, 1))

            Dvalue(ieq) = thdyn_Machnumber(GAMMA_AIR,&
                p, Ddata(ieq, 1),&
                Ddata(ieq, 2)/Ddata(ieq, 1))
          end do

        case (NVAR2D)
          do ieq = 1, neq
            p = thdyn_pressure(GAMMA_AIR,&
                Ddata(ieq, nvar)/Ddata(ieq, 1), Ddata(ieq, 1),&
                Ddata(ieq, 2)/Ddata(ieq, 1),&
                Ddata(ieq, 3)/Ddata(ieq, 1))

            Dvalue(ieq) = thdyn_Machnumber(GAMMA_AIR,&
                p, Ddata(ieq, 1),&
                Ddata(ieq, 2)/Ddata(ieq, 1),&
                Ddata(ieq, 3)/Ddata(ieq, 1))
          end do
          
        case (NVAR3D)
          do ieq = 1, neq
            p = thdyn_pressure(GAMMA_AIR,&
                Ddata(ieq, nvar)/Ddata(ieq, 1), Ddata(ieq, 1),&
                Ddata(ieq, 2)/Ddata(ieq, 1),&
                Ddata(ieq, 3)/Ddata(ieq, 1),&
                Ddata(ieq, 4)/Ddata(ieq, 1))

            Dvalue(ieq) = thdyn_Machnumber(GAMMA_AIR,&
                p, Ddata(ieq, 1),&
                Ddata(ieq, 2)/Ddata(ieq, 1),&
                Ddata(ieq, 3)/Ddata(ieq, 1),&
                Ddata(ieq, 4)/Ddata(ieq, 1))
          end do
        end select
      end select
      
    end subroutine getVarBlockformat

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
  end subroutine euler_outputStatistics

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
        call parlst_setvalue(rparlist, '', "adaptivity", trim(adjustl(cbuffer)))

      case ('-B','--benchmark')
        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        call parlst_setvalue(rparlist, '', "benchmark", trim(adjustl(cbuffer)))
       
      case ('-DC','--dualconv')
        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        call parlst_setvalue(rparlist, '', "dualconv", trim(adjustl(cbuffer)))
        
      case ('-DD','--dualdiff')
        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        call parlst_setvalue(rparlist, '', "dualdiff", trim(adjustl(cbuffer)))

      case ('-E','--errorest')
        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        call parlst_setvalue(rparlist, '', "errorest", trim(adjustl(cbuffer)))
        
      case ('-I','--io')
        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        call parlst_setvalue(rparlist, '', "inputoutput", trim(adjustl(cbuffer)))
        
      case ('-PC','--primalconv')
        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        call parlst_setvalue(rparlist, '', "primalconv", trim(adjustl(cbuffer)))
        
      case ('-PD','--primaldiff')
        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        call parlst_setvalue(rparlist, '', "primaldiff", trim(adjustl(cbuffer)))

      case ('-S','--solver')
        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        call parlst_setvalue(rparlist, '', "solver", trim(adjustl(cbuffer)))
        
      case ('-T','--timestep')
        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        call parlst_setvalue(rparlist, '', "timestep", trim(adjustl(cbuffer)))
        
      case DEFAULT
        iarg = iarg+1
        if (iarg .ge. narg) exit cmdarg
      end select
    end do cmdarg

  end subroutine euler_parseCmdlArguments

end module euler_application
