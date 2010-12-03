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
!# 2.) hydro_initSolvers
!#     -> Initializes the solve structures from the parameter list.
!#
!# 3.) hydro_initProblem
!#     -> Initializes the global problem structure based on the
!#        parameter settings given by the parameter list. This routine
!#        is quite universal, that is, it prepares the internal
!#        structure of the global problem and generates a linked
!#        list of problem levels used in the multigrid hierarchy.
!#
!# 4.) hydro_initProblemLevel
!#     -> Initializes the individual problem level based on the
!#        parameter settings given by the parameter list.
!#        This routine is called repeatedly by the global
!#        initialisation routine hydro_initAllProblemLevels.
!#
!# 5.) hydro_initAllProblemLevels
!#     -> Initializes ALL problem levels attached to the global
!#        problem structure based on the parameter settings
!#        given by the parameter list.
!#
!# 6.) hydro_initSolution
!#     -> Initializes the solution vector based on the parameter
!#        settings given by the parameter list
!#
!# 7.) hydro_outputSolution
!#     -> Outputs the solution vector to file in UCD format
!#
!# 8.) hydro_outputStatistics
!#     -> Outputs the application statitics
!#
!# 9.) hydro_estimateRecoveryError
!#      -> Estimates the solution error using recovery techniques
!#
!# 10.) hydro_adaptTriangulation
!#      -> Performs h-adaptation for the given triangulation
!#
!# 11.) hydro_solveTransientPrimal
!#      -> Solves the primal formulation of the time-dependent
!#         compressible Euler equations.
!#
!# 12.) hydro_projectSolution
!#      -> Performs conservative projection of the solution from
!#         a given FE-space to another FE-space
!#
!# The following auxiliary routines are available:
!#
!# 1.) hydro_parseCmdlArguments
!#     -> Parses the list of commandline arguments and overwrites
!#        parameter values from the parameter files
!#
!# </purpose>
!##############################################################################

module hydro_application

  use afcstabilisation
  use bilinearformevaluation
  use boundary
  use boundarycondaux
  use boundaryfilter
  use collection
  use cubature
  use derivatives
  use element
  use hydro_basic
  use hydro_callback
  use hydro_callback1d
  use hydro_callback2d
  use hydro_callback3d
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
  use scalarpde
  use solveraux
  use spatialdiscretisation
  use statistics
  use stdoperators
  use storage
  use timestep
  use timestepaux
  use ucd

  implicit none

  private
  public :: hydro_app
  public :: hydro_adaptTriangulation
  public :: hydro_adjustParameterlist
  public :: hydro_estimateRecoveryError
  public :: hydro_initAllProblemLevels
  public :: hydro_initProblem
  public :: hydro_initProblemLevel
  public :: hydro_initSolution
  public :: hydro_initSolvers
  public :: hydro_outputSolution
  public :: hydro_outputStatistics
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

    ! Parameter file and section names
    character(LEN=SYS_STRLEN) :: sindatfileName
    character(LEN=SYS_STRLEN) :: sbdrcondName
    character(LEN=SYS_STRLEN) :: algorithm

    ! local variables
    integer :: isystemFormat, systemMatrix, ndimension


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
    call hydro_parseCmdlArguments(rparlist)
    call hydro_adjustParameterlist(rparlist, ssectionName)

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
    call hydro_initSolvers(rparlist, ssectionName, rtimestep, rsolver)

    ! Initialize the abstract problem structure
    call hydro_initProblem(rparlist, ssectionName,&
        solver_getMinimumMultigridlevel(rsolver),&
        solver_getMaximumMultigridlevel(rsolver), rproblem)

    ! Initialize the individual problem levels
    call hydro_initAllProblemLevels(rparlist, ssectionName,&
        rproblem, rcollection)

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

    if (rtimestep%dfinalTime .gt. rtimestep%dinitialTime) then

      ! Get global configuration from parameter list
      call parlst_getvalue_string(rparlist,&
          ssectionName, 'algorithm', algorithm)
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'ndimension', ndimension)
      call parlst_getvalue_string(rparlist,&
          ssectionName, 'sprimalbdrcondname', sbdrcondName)
      
      ! The boundary condition for the primal problem is required for
      ! all solution strategies so initialize it from the parameter file
      call bdrc_readBoundaryCondition(rbdrCondPrimal,&
          sindatfileName, '['//trim(sbdrcondName)//']',&
          ndimension, hydro_parseBoundaryCondition)

      ! What solution algorithm should be applied?
      if (trim(algorithm) .eq. 'transient_primal') then
        
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Solve the primal formulation for the time-dependent problem
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call hydro_solveTransientPrimal(rparlist, ssectionName,&
            rbdrCondPrimal, rproblem, rtimestep, rsolver,&
            rsolutionPrimal, rcollection)

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
    call hydro_outputStatistics(rtimerTotal, rcollection)

    ! Release collection
    call collct_done(rcollection)

  end subroutine hydro_app

  !*****************************************************************************

!<subroutine>

  subroutine hydro_initSolvers(rparlist, ssectionName,&
      rtimestep, rsolver)

!<description>
    ! This subroutine initializes the time-stepping structure and
    ! the top-level solver structure from the parameter list
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
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'timestep', stimestepName)
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'solver',   ssolverName)

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

  end subroutine hydro_initSolvers

  !*****************************************************************************

!<subroutine>

  subroutine hydro_initProblem(rparlist, ssectionName,&
      nlmin, nlmax, rproblem)

!<description>
    ! This subroutine initializes the abstract problem structure
    ! based on the parameters settings given by the parameter list
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
    character(LEN=SYS_STRLEN) :: sinviscid

    ! abstract problem descriptor
    type(t_problemDescriptor) :: rproblemDescriptor

    ! pointer to the problem level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! local variables
    integer :: inviscidAFC, iconvToTria


    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'trifile', rproblemDescriptor%trifile)
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'prmfile', rproblemDescriptor%prmfile, '')
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'ndimension', rproblemDescriptor%ndimension)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'iconvtotria', iconvToTria, 0)
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'inviscid', sinviscid)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'inviscidAFC', inviscidAFC)

    ! Set additional problem descriptor
    rproblemDescriptor%ndiscretisation = 1   ! one discretisation
    rproblemDescriptor%nafcstab        = 1   ! stabilisation for inviscid fluxes
    rproblemDescriptor%nlmin           = nlmin
    rproblemDescriptor%nlmax           = nlmax
    rproblemDescriptor%nmatrixScalar   = rproblemDescriptor%ndimension + 5
    rproblemDescriptor%nmatrixBlock    = 2   ! system matrix and Jacobian
    rproblemDescriptor%nvectorScalar   = 0
    rproblemDescriptor%nvectorBlock    = 0

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

      if (inviscidAFC > 0) then
        call afcstab_initFromParameterlist(rparlist, sinviscid,&
            p_rproblemLevel%Rafcstab(inviscidAFC))
      end if

      ! Switch to next coarser level
      p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
    end do

  end subroutine hydro_initProblem

  !*****************************************************************************

!<subroutine>

  subroutine hydro_initProblemLevel(rparlist, ssectionName,&
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
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    type(t_triangulation) , pointer :: p_rtriangulation
    type(t_boundary) , pointer :: p_rboundary
    character(len=SYS_STRLEN) :: slimitingvariable
    integer :: templateMatrix
    integer :: systemMatrix
    integer :: jacobianMatrix
    integer :: consistentMassMatrix
    integer :: lumpedMassMatrix
    integer :: coeffMatrix_CX
    integer :: coeffMatrix_CY
    integer :: coeffMatrix_CZ
    integer :: coeffMatrix_CXX
    integer :: coeffMatrix_CYY
    integer :: coeffMatrix_CZZ
    integer :: coeffMatrix_CXY
    integer :: coeffMatrix_CXZ
    integer :: coeffMatrix_CYZ
    integer :: inviscidAFC
    integer :: discretisation
    integer :: celement
    integer :: isystemFormat
    integer :: isystemCoupling
    integer :: imatrixFormat

    integer :: i,j,ivar,jvar,ivariable,nvariable,nvartransformed,nmatrices
    integer :: nsumcubRefBilForm,nsumcubRefLinForm,nsumcubRefEval

    ! Retrieve application specific parameters from the collection
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'templatematrix', templateMatrix)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'systemmatrix', systemMatrix)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'jacobianmatrix', jacobianMatrix)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'consistentmassmatrix', consistentMassMatrix)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'lumpedmassmatrix', lumpedMassMatrix)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_cx', coeffMatrix_CX)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_cy', coeffMatrix_CY)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_cz', coeffMatrix_CZ)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_cxx', coeffMatrix_CXX)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_cyy', coeffMatrix_CYY)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_czz', coeffMatrix_CZZ)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_cxy', coeffMatrix_CXY)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_cxz', coeffMatrix_CXZ)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_cyz', coeffMatrix_CYZ)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'inviscidAFC', inviscidAFC)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'discretisation', discretisation)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'imatrixFormat', imatrixFormat)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'isystemFormat', isystemFormat)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'isystemCoupling', isystemCoupling)
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
        select case(isystemFormat)
        case (SYSTEM_INTERLEAVEFORMAT)
          call spdiscr_initBlockDiscr(p_rdiscretisation, 1,&
              rproblemLevel%rtriangulation)

        case (SYSTEM_BLOCKFORMAT)
          call spdiscr_initBlockDiscr(p_rdiscretisation,&
              hydro_getNVAR(rproblemLevel), p_rtriangulation)

        case DEFAULT
          call output_line('Unsupported system format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_initProblemLevel')
          call sys_halt()
        end select
      end if

      ! Get spatial dimension
      select case(p_rtriangulation%ndim)
      case (NDIM1D)
        select case(celement)
        case (-1,1,11)
          ! P1=Q1 finite elements
          call spdiscr_initDiscr_simple(&
              p_rdiscretisation%RspatialDiscr(1),&
              EL_E001_1D, SPDISC_CUB_AUTOMATIC,&
              p_rtriangulation, p_rboundary)

        case (-2,2,12)
          ! P2=Q2 finite elements
          call spdiscr_initDiscr_simple(&
              p_rdiscretisation%RspatialDiscr(1),&
              EL_E002_1D, SPDISC_CUB_AUTOMATIC,&
              p_rtriangulation, p_rboundary)

        case DEFAULT
          call output_line('Unsupproted element type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_initProblemLevel')
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
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_initProblemLevel')
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
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_initProblemLevel')
          call sys_halt()
        end select

      case DEFAULT
        call output_line('Invalid number of spatial dimensions',&
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_initProblemLevel')
        call sys_halt()
      end select

      ! Duplicate scalar discretisation structure for block matrix format
      if (isystemFormat .eq. SYSTEM_BLOCKFORMAT) then
        do ivar = 2, hydro_getNVAR(rproblemLevel)
          call spdiscr_duplicateDiscrSc(&
              p_rdiscretisation%RspatialDiscr(1),&
              p_rdiscretisation%RspatialDiscr(ivar), .true.)
        end do
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

    end if


    ! If the template matrix has no structure data then generate the
    ! finite element matrix sparsity structure based on the spatial
    ! descretisation and store it as the template matrix. Otherwise we
    ! assume that the template matrix has been generated externally.
    if (.not.lsyssc_hasMatrixStructure(rproblemLevel%Rmatrix(templateMatrix))) then
      call bilf_createMatrixStructure(&
          p_rdiscretisation%RspatialDiscr(1), imatrixFormat,&
          rproblemLevel%Rmatrix(templateMatrix))

    end if


    ! Create system matrix
    if (systemMatrix > 0) then
      select case(isystemFormat)

      case (SYSTEM_INTERLEAVEFORMAT)
        ! The global operator is stored as an interleave matrix with
        ! NVAR components. However, the row and column structure of
        ! the template matrix can be adopted without modification
        if (lsyssc_hasMatrixStructure(rproblemLevel%Rmatrix(systemMatrix))) then

          ! Release pseudo block matrix
          call lsysbl_releaseMatrix(rproblemLevel%RmatrixBlock(systemMatrix))

          ! Resize scalar matrix
          call lsyssc_resizeMatrix(&
              rproblemLevel%Rmatrix(systemMatrix),&
              rproblemLevel%Rmatrix(templateMatrix)%NEQ,&
              rproblemLevel%Rmatrix(templateMatrix)%NCOLS,&
              rproblemLevel%rmatrix(templateMatrix)%NA,&
              .false., .false., bforce=.true.)

        else   ! System matrix has no structure

          call lsyssc_duplicateMatrix(&
              rproblemLevel%Rmatrix(templateMatrix),&
              rproblemLevel%Rmatrix(systemMatrix),&
              LSYSSC_DUP_SHARE, LSYSSC_DUP_REMOVE)

          ! Set number of variables per node
          rproblemLevel%Rmatrix(systemMatrix)%NVAR = hydro_getNVAR(rproblemLevel)

          ! What matrix format should be used?
          select case(imatrixFormat)
          case (LSYSSC_MATRIX7)
            rproblemLevel%Rmatrix(systemMatrix)%cmatrixFormat = LSYSSC_MATRIX7INTL

          case (LSYSSC_MATRIX9)
            rproblemLevel%Rmatrix(systemMatrix)%cmatrixFormat = LSYSSC_MATRIX9INTL

          case DEFAULT
            call output_line('Unsupported matrix format!',&
                OU_CLASS_ERROR,OU_MODE_STD,'hydro_initProblemLevel')
            call sys_halt()
          end select

          ! What kind of global operator should be adopted?
          select case(isystemCoupling)
          case (SYSTEM_SEGREGATED)
            rproblemLevel%Rmatrix(systemMatrix)%cinterleavematrixFormat = LSYSSC_MATRIXD

          case (SYSTEM_ALLCOUPLED)
            rproblemLevel%Rmatrix(systemMatrix)%cinterleavematrixFormat = LSYSSC_MATRIX1

          case DEFAULT
            call output_line('Unsupported interleave matrix format!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'hydro_initProblemLevel')
            call sys_halt()
          end select

          ! Create global operator physically
          call lsyssc_allocEmptyMatrix(&
              rproblemLevel%Rmatrix(systemMatrix), LSYSSC_SETM_UNDEFINED)

        end if

        ! Create pseudo block matrix from global operator
        call lsysbl_createMatFromScalar(&
            rproblemLevel%Rmatrix(systemMatrix),&
            rproblemLevel%RmatrixBlock(systemMatrix), p_rdiscretisation)



      case (SYSTEM_BLOCKFORMAT)
        ! The global operator is stored as a block matrix with
        ! NVARxNVAR blocks made up from scalar matrices

        if ((rproblemLevel%RmatrixBlock(systemMatrix)%nblocksPerRow .ne. 0) .and.&
            (rproblemLevel%RmatrixBlock(systemMatrix)%nblocksPerCol .ne. 0)) then

          ! What kind of global operator should be adopted?
          select case(isystemCoupling)

          case (SYSTEM_SEGREGATED)
            ! Create only NVAR diagonal blocks
            do ivar = 1, hydro_getNVAR(rproblemLevel)
              call lsyssc_resizeMatrix(&
                  rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                  rproblemLevel%Rmatrix(templateMatrix), .false., .false., .true.)
            end do

          case (SYSTEM_ALLCOUPLED)
            ! Create all NVAR x NVAR blocks
            do ivar = 1, hydro_getNVAR(rproblemLevel)
              do jvar = 1, hydro_getNVAR(rproblemLevel)
                call lsyssc_resizeMatrix(&
                    rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar ,jvar),&
                    rproblemLevel%Rmatrix(templateMatrix), .false., .false., .true.)
              end do
            end do

          case DEFAULT
            call output_line('Unsupported block matrix format!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'hydro_initProblemLevel')
            call sys_halt()
          end select

        else   ! System matrix has no structure

          ! Create empty NVARxNVAR block matrix directly
          call lsysbl_createEmptyMatrix(&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              hydro_getNVAR(rproblemLevel),&
              hydro_getNVAR(rproblemLevel))

          ! Specify matrix as 'group matrix'
          rproblemLevel%RmatrixBlock(systemMatrix)%imatrixSpec = LSYSBS_MSPEC_GROUPMATRIX

          ! What kind of global operator should be adopted?
          select case(isystemCoupling)

          case (SYSTEM_SEGREGATED)
            ! Create only NVAR diagonal blocks
            do ivar = 1, hydro_getNVAR(rproblemLevel)
              call lsyssc_duplicateMatrix(&
                  rproblemLevel%Rmatrix(templateMatrix),&
                  rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                  LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
            end do

          case (SYSTEM_ALLCOUPLED)
            ! Create all NVAR x NVAR blocks
            do ivar = 1, hydro_getNVAR(rproblemLevel)
              do jvar = 1, hydro_getNVAR(rproblemLevel)
                call lsyssc_duplicateMatrix(&
                    rproblemLevel%Rmatrix(templateMatrix),&
                    rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,jvar),&
                    LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
              end do
            end do

          case DEFAULT
            call output_line('Unsupported block matrix format!',&
                OU_CLASS_ERROR,OU_MODE_STD,'hydro_initProblemLevel')
            call sys_halt()
          end select

        end if

        ! Update internal structure of block matrix
        call lsysbl_updateMatStrucInfo(rproblemLevel%RmatrixBlock(systemMatrix))

      case DEFAULT
        call output_line('Unsupported system format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_initProblemLevel')
        call sys_halt()
      end select
    end if


    ! Create consistent (and lumped) mass matrix as duplicate of the template matrix
    if (consistentMassMatrix > 0) then
      if (lsyssc_isMatrixStructureShared(&
          rproblemLevel%Rmatrix(consistentMassMatrix),&
          rproblemLevel%Rmatrix(templateMatrix))) then

        call lsyssc_resizeMatrix(&
            rproblemLevel%Rmatrix(consistentMassMatrix),&
            rproblemLevel%Rmatrix(templateMatrix), .false., .false., .true.)

      else
        call lsyssc_duplicateMatrix(&
            rproblemLevel%Rmatrix(templateMatrix),&
            rproblemLevel%Rmatrix(consistentMassMatrix),&
            LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

      end if
      call stdop_assembleSimpleMatrix(&
          rproblemLevel%Rmatrix(consistentMassMatrix),&
          DER_FUNC, DER_FUNC)
      if (lumpedMassMatrix > 0) then
        call lsyssc_duplicateMatrix(&
            rproblemLevel%Rmatrix(consistentMassMatrix),&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)

        call lsyssc_lumpMatrixScalar(&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            LSYSSC_LUMP_DIAG)

      end if
    elseif (lumpedMassMatrix > 0) then
      call lsyssc_duplicateMatrix(&
          rproblemLevel%Rmatrix(templateMatrix),&
          rproblemLevel%Rmatrix(lumpedMassMatrix),&
          LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

      call stdop_assembleSimpleMatrix(&
          rproblemLevel%Rmatrix(lumpedMassMatrix),&
          DER_FUNC, DER_FUNC)

      call lsyssc_lumpMatrixScalar(&
          rproblemLevel%Rmatrix(lumpedMassMatrix),&
          LSYSSC_LUMP_DIAG)

    end if


    ! Create coefficient matrix (phi, dphi/dx) as duplicate of the
    ! template matrix
    if (coeffMatrix_CX > 0) then
      if (lsyssc_isMatrixStructureShared(&
          rproblemLevel%Rmatrix(coeffMatrix_CX),&
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
          rproblemLevel%Rmatrix(coeffMatrix_CX),&
          DER_DERIV3D_X, DER_FUNC)
    end if

    
    ! Create coefficient matrix (phi, dphi/dy) as duplicate of the
    ! template matrix
    if (coeffMatrix_CY > 0) then
      if (lsyssc_isMatrixStructureShared(&
          rproblemLevel%Rmatrix(coeffMatrix_CY),&
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
          rproblemLevel%Rmatrix(coeffMatrix_CY),&
          DER_DERIV3D_Y, DER_FUNC)
    end if


    ! Create coefficient matrix (phi, dphi/dz) as duplicate of the
    ! template matrix
    if (coeffMatrix_CZ > 0) then
      if (lsyssc_isMatrixStructureShared(&
          rproblemLevel%Rmatrix(coeffMatrix_CZ),&
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
          rproblemLevel%Rmatrix(coeffMatrix_CZ),&
          DER_DERIV3D_Z, DER_FUNC)
    end if

    
    ! Create coefficient matrix (dphi/dx, dphi/dx) as duplicate of the
    ! template matrix
    if (coeffMatrix_CXX > 0) then
      if (lsyssc_isMatrixStructureShared(&
          rproblemLevel%Rmatrix(coeffMatrix_CXX),&
          rproblemLevel%Rmatrix(templateMatrix))) then

        call lsyssc_resizeMatrix(&
            rproblemLevel%Rmatrix(coeffMatrix_CXX),&
            rproblemLevel%Rmatrix(templateMatrix),&
            .false., .false., .true.)

      else
        call lsyssc_duplicateMatrix(&
            rproblemLevel%Rmatrix(templateMatrix),&
            rproblemLevel%Rmatrix(coeffMatrix_CXX),&
            LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

      end if
      call stdop_assembleSimpleMatrix(&
          rproblemLevel%Rmatrix(coeffMatrix_CXX),&
          DER_DERIV3D_X, DER_DERIV3D_X)
    end if

    
    ! Create coefficient matrix (dphi/dy, dphi/dy) as duplicate of the
    ! template matrix
    if (coeffMatrix_CYY > 0) then
      if (lsyssc_isMatrixStructureShared(&
          rproblemLevel%Rmatrix(coeffMatrix_CYY),&
          rproblemLevel%Rmatrix(templateMatrix))) then

        call lsyssc_resizeMatrix(&
            rproblemLevel%Rmatrix(coeffMatrix_CYY),&
            rproblemLevel%Rmatrix(templateMatrix),&
            .false., .false., .true.)

      else
        call lsyssc_duplicateMatrix(&
            rproblemLevel%Rmatrix(templateMatrix),&
            rproblemLevel%Rmatrix(coeffMatrix_CYY),&
            LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

      end if
      call stdop_assembleSimpleMatrix(&
          rproblemLevel%Rmatrix(coeffMatrix_CYY),&
          DER_DERIV3D_Y, DER_DERIV3D_Y)
    end if


    ! Create coefficient matrix (dphi/dz, dphi/dz) as duplicate of the
    ! template matrix
    if (coeffMatrix_CZZ > 0) then
      if (lsyssc_isMatrixStructureShared(&
          rproblemLevel%Rmatrix(coeffMatrix_CZZ),&
          rproblemLevel%Rmatrix(templateMatrix))) then

        call lsyssc_resizeMatrix(&
            rproblemLevel%Rmatrix(coeffMatrix_CZZ),&
            rproblemLevel%Rmatrix(templateMatrix),&
            .false., .false., .true.)

      else
        call lsyssc_duplicateMatrix(&
            rproblemLevel%Rmatrix(templateMatrix),&
            rproblemLevel%Rmatrix(coeffMatrix_CZZ),&
            LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

      end if
      call stdop_assembleSimpleMatrix(&
          rproblemLevel%Rmatrix(coeffMatrix_CZZ),&
          DER_DERIV3D_Z, DER_DERIV3D_Z)
    end if


    ! Create coefficient matrix (dphi/dx, dphi/dy) as duplicate of the
    ! template matrix
    if (coeffMatrix_CXY > 0) then
      if (lsyssc_isMatrixStructureShared(&
          rproblemLevel%Rmatrix(coeffMatrix_CXY),&
          rproblemLevel%Rmatrix(templateMatrix))) then

        call lsyssc_resizeMatrix(&
            rproblemLevel%Rmatrix(coeffMatrix_CXY),&
            rproblemLevel%Rmatrix(templateMatrix),&
            .false., .false., .true.)

      else
        call lsyssc_duplicateMatrix(&
            rproblemLevel%Rmatrix(templateMatrix),&
            rproblemLevel%Rmatrix(coeffMatrix_CXY),&
            LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

      end if
      call stdop_assembleSimpleMatrix(&
          rproblemLevel%Rmatrix(coeffMatrix_CXY),&
          DER_DERIV3D_X, DER_DERIV3D_Y)
    end if

    
    ! Create coefficient matrix (dphi/dx, dphi/dz) as duplicate of the
    ! template matrix
    if (coeffMatrix_CXZ > 0) then
      if (lsyssc_isMatrixStructureShared(&
          rproblemLevel%Rmatrix(coeffMatrix_CXZ),&
          rproblemLevel%Rmatrix(templateMatrix))) then

        call lsyssc_resizeMatrix(&
            rproblemLevel%Rmatrix(coeffMatrix_CXZ),&
            rproblemLevel%Rmatrix(templateMatrix),&
            .false., .false., .true.)

      else
        call lsyssc_duplicateMatrix(&
            rproblemLevel%Rmatrix(templateMatrix),&
            rproblemLevel%Rmatrix(coeffMatrix_CXZ),&
            LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

      end if
      call stdop_assembleSimpleMatrix(&
          rproblemLevel%Rmatrix(coeffMatrix_CXZ),&
          DER_DERIV3D_X, DER_DERIV3D_Z)
    end if
    

    ! Create coefficient matrix (dphi/dy, dphi/dz) as duplicate of the
    ! template matrix
    if (coeffMatrix_CYZ > 0) then
      if (lsyssc_isMatrixStructureShared(&
          rproblemLevel%Rmatrix(coeffMatrix_CYZ),&
          rproblemLevel%Rmatrix(templateMatrix))) then

        call lsyssc_resizeMatrix(&
            rproblemLevel%Rmatrix(coeffMatrix_CYZ),&
            rproblemLevel%Rmatrix(templateMatrix),&
            .false., .false., .true.)

      else
        call lsyssc_duplicateMatrix(&
            rproblemLevel%Rmatrix(templateMatrix),&
            rproblemLevel%Rmatrix(coeffMatrix_CYZ),&
            LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

      end if
      call stdop_assembleSimpleMatrix(&
          rproblemLevel%Rmatrix(coeffMatrix_CYZ),&
          DER_DERIV3D_Y, DER_DERIV3D_Z)
    end if


    ! Resize stabilisation structures if necessary and remove the
    ! indicator for the subdiagonal edge structure. If they are
    ! needed, then they are re-generated on-the-fly.
    if (inviscidAFC > 0) then
      if (rproblemLevel%Rafcstab(inviscidAFC)%istabilisationSpec&
          .eq. AFCSTAB_UNDEFINED) then

        ! Get number of expressions for limiting variables
        nvariable = max(1,&
            parlst_querysubstrings(rparlist,&
            ssectionName, 'slimitingvariable'))
        
        ! Initialise number of limiting variables
        nvartransformed = 1

        ! Determine maximum number of limiting variables in a single set
        do ivariable = 1, nvariable
          call parlst_getvalue_string(rparlist,&
              ssectionName, 'slimitingvariable',&
              slimitingvariable, isubstring=ivariable)
          nvartransformed = max(nvartransformed,&
              hydro_getNVARtransformed(rproblemLevel, slimitingvariable))
        end do

        ! Initialise stabilisation structure
        call gfsys_initStabilisation(&
            rproblemLevel%RmatrixBlock(systemMatrix),&
            rproblemLevel%Rafcstab(inviscidAFC),&
            nvartransformed, p_rdiscretisation)

        ! Compute number of matrices to by copied
        nmatrices = 0
        if (coeffMatrix_CX   > 0) nmatrices = nmatrices+1
        if (coeffMatrix_CY   > 0) nmatrices = nmatrices+1
        if (coeffMatrix_CZ   > 0) nmatrices = nmatrices+1
        if (coeffMatrix_CXX  > 0) nmatrices = nmatrices+1
        if (coeffMatrix_CYY  > 0) nmatrices = nmatrices+1
        if (coeffMatrix_CZZ  > 0) nmatrices = nmatrices+1
        if (coeffMatrix_CXY  > 0) nmatrices = nmatrices+1
        if (coeffMatrix_CXZ  > 0) nmatrices = nmatrices+1
        if (coeffMatrix_CYZ  > 0) nmatrices = nmatrices+1

        ! Initialise memory for constant coefficient matrices
        call afcstab_initMatrixCoeffs(rproblemLevel%Rafcstab(inviscidAFC), nmatrices)

      else
        ! Resize stabilisation structure
        call afcstab_resizeStabilisation(&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rproblemLevel%Rmatrix(templateMatrix))

        rproblemLevel%Rafcstab(inviscidAFC)%istabilisationSpec =&
            iand(rproblemLevel%Rafcstab(inviscidAFC)%istabilisationSpec,&
            not(AFCSTAB_HAS_OFFDIAGONALEDGES))
      end if
      
      ! Copy constant coefficient matrices to stabilisation structure
      nmatrices = 0
      if (coeffMatrix_CX > 0) then
        nmatrices = nmatrices+1
        call afcstab_CopyMatrixCoeffs(&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX), (/nmatrices/))
      end if
      if (coeffMatrix_CY > 0) then
        nmatrices = nmatrices+1
        call afcstab_CopyMatrixCoeffs(&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rproblemLevel%Rmatrix(coeffMatrix_CY:coeffMatrix_CY), (/nmatrices/))
      end if
      if (coeffMatrix_CZ > 0) then
        nmatrices = nmatrices+1
        call afcstab_CopyMatrixCoeffs(&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rproblemLevel%Rmatrix(coeffMatrix_CZ:coeffMatrix_CZ), (/nmatrices/))
      end if
      if (coeffMatrix_CXX > 0) then
        nmatrices = nmatrices+1
        call afcstab_CopyMatrixCoeffs(&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rproblemLevel%Rmatrix(coeffMatrix_CXX:coeffMatrix_CXX), (/nmatrices/))
      end if
      if (coeffMatrix_CYY > 0) then
        nmatrices = nmatrices+1
        call afcstab_CopyMatrixCoeffs(&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rproblemLevel%Rmatrix(coeffMatrix_CYY:coeffMatrix_CYY), (/nmatrices/))
      end if
      if (coeffMatrix_CZZ > 0) then
        nmatrices = nmatrices+1
        call afcstab_CopyMatrixCoeffs(&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rproblemLevel%Rmatrix(coeffMatrix_CZZ:coeffMatrix_CZZ), (/nmatrices/))
      end if
      if (coeffMatrix_CXY > 0) then
        nmatrices = nmatrices+1
        call afcstab_CopyMatrixCoeffs(&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rproblemLevel%Rmatrix(coeffMatrix_CXY:coeffMatrix_CXY), (/nmatrices/))
      end if
      if (coeffMatrix_CXZ > 0) then
        nmatrices = nmatrices+1
        call afcstab_CopyMatrixCoeffs(&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rproblemLevel%Rmatrix(coeffMatrix_CXZ:coeffMatrix_CXZ), (/nmatrices/))
      end if
      if (coeffMatrix_CYZ > 0) then
        nmatrices = nmatrices+1
        call afcstab_CopyMatrixCoeffs(&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rproblemLevel%Rmatrix(coeffMatrix_CYZ:coeffMatrix_CYZ), (/nmatrices/))
      end if
      
    end if

  end subroutine hydro_initProblemLevel

  !*****************************************************************************

!<subroutine>

  subroutine hydro_initAllProblemLevels(rparlist, ssectionName,&
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
      call hydro_initProblemLevel(rparlist, ssectionName,&
          p_rproblemLevel, rcollection)

      ! Switch to next coarser level
      p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
    end do

  end subroutine hydro_initAllProblemLevels

  !*****************************************************************************

!<subroutine>

  subroutine hydro_initSolution(rparlist, ssectionName,&
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
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_afcstab) :: rafcstab
    type(t_linearForm) :: rform
    type(t_vectorBlock) :: rvectorBlock, rvectorHigh, rvectorAux
    type(t_matrixScalar), target :: rlumpedMassMatrix, rconsistentMassMatrix
    type(t_matrixScalar), pointer :: p_rlumpedMassMatrix, p_rConsistentMassMatrix
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
    type(t_fparser), pointer :: p_rfparser
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Ddata
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP), dimension(:), pointer :: Dnorm0
    real(DP) :: depsAbsSolution, depsRelSolution, dnorm
    character(len=SYS_STRLEN), dimension(:), pointer :: SsolutionFailsafeVariables
    character(LEN=SYS_STRLEN) :: ssolutionName
    integer :: isolutiontype, nexpression, nsolutionfailsafe
    integer :: icomp, iblock, ivar, nvar, ieq, neq, ndim, iter
    integer :: lumpedMassMatrix, consistentMassMatrix, systemMatrix
    integer :: nmaxIterationsSolution, ivariable, nvariable

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

      call output_line('Initialisation if solution by graymap image is not yet supported!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'hydro_initSolution')
      call sys_halt()


    case (SOLUTION_ANALYTIC_POINTVALUE)

      !-------------------------------------------------------------------------
      ! Initialize the nodal values by the data of an analytical expression
      !-------------------------------------------------------------------------

      ! Initialize total number of expressions
      nexpression = 0

      ! Compute total number of expressions
      do iblock = 1, rvector%nblocks
        nexpression = nexpression + rvector%RvectorBlock(iblock)%NVAR
      end do

      ! Check if array of solution names is available
      if (parlst_querysubstrings(rparlist, ssectionName,&
          'ssolutionname') .lt. nexpression) then
        call output_line('Invalid number of expressions!',&
            OU_CLASS_ERROR, OU_MODE_STD, 'hydro_initSolution')
        call sys_halt()
      end if

      ! Get function parser from collection structure
      p_rfparser => collct_getvalue_pars(rcollection, 'rfparser')

      ! Set pointers
      call storage_getbase_double2D(&
          rproblemLevel%rtriangulation%h_DvertexCoords, p_DvertexCoords)

      ! Get number of spatial dimensions
      ndim = rproblemLevel%rtriangulation%ndim

      ! Initialize variable values
      Dvalue = 0.0_DP
      nexpression = 0

      ! Loop over all blocks of the global solution vector
      do iblock = 1, rvector%nblocks

        ! Set pointer to data array
        call lsyssc_getbase_double(rvector%RvectorBlock(iblock), p_Ddata)

        ! Initialisation for scalar subvector
        neq  = rvector%RvectorBlock(iblock)%NEQ
        nvar = rvector%RvectorBlock(iblock)%NVAR

        ! Loop over all equations of the scalar subvector
        do ieq = 1, neq

          ! Set coordinates and evalution time
          Dvalue(1:ndim)   = p_DvertexCoords(:,ieq)
          Dvalue(NDIM3D+1) = dtime

          ! Loop over all variables of the solution vector
          do ivar = 1, nvar

            ! Get the function name of the component used for evaluating the initial solution.
            call parlst_getvalue_string(rparlist, ssectionName,&
                'ssolutionName', ssolutionName, isubstring=nexpression+ivar)

            ! Get the number of the component used for evaluating the initial solution
            icomp = fparser_getFunctionNumber(p_rfparser, ssolutionname)

            ! Evaluate the function parser
            call fparser_evalFunction(p_rfparser, icomp,&
                Dvalue, p_Ddata((ieq-1)*nvar+ivar))

          end do   ! ivar
        end do   ! ieq

        ! Increase number of processed expressions
        nexpression = nexpression + nvar

      end do   ! iblock
   

    case (SOLUTION_ANALYTIC_L2_CONSISTENT,&
          SOLUTION_ANALYTIC_L2_LUMPED)

      !-------------------------------------------------------------------------
      ! Initialize the FE-function by the L2-projection of the analytical data
      !-------------------------------------------------------------------------

      ! Initialize total number of expressions
      nexpression = 0

      ! Compute total number of expressions
      do iblock = 1, rvector%nblocks
        nexpression = nexpression + rvector%RvectorBlock(iblock)%NVAR
      end do

      ! Check if array of solution names is available
      if (parlst_querysubstrings(rparlist, ssectionName,&
          'ssolutionname') .lt. nexpression) then
        call output_line('Invalid number of expressions!',&
            OU_CLASS_ERROR, OU_MODE_STD, 'hydro_initSolution')
        call sys_halt()
      end if

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
      
      ! Initialize number of expressions
      nexpression = 0

      ! Loop over all blocks of the global solution vector
      do iblock = 1, rvector%nblocks
        
        ! Set pointer to spatial discretisation
        p_rspatialDiscr => rvector%RvectorBlock(iblock)%p_rspatialDiscr
        
        ! Scalar vectors in interleaved format have to be treated differently
        if (rvector%RvectorBlock(iblock)%NVAR .eq. 1) then

          ! Get the function name of the component used for evaluating the initial solution.
          call parlst_getvalue_string(rparlist, ssectionName,&
              'ssolutionName', ssolutionName, isubstring=nexpression+1)

          ! Set the number of the component used for evaluating the initial solution
          rcollection%IquickAccess(1) = fparser_getFunctionNumber(p_rfparser, ssolutionname)
          
          ! Assemble the linear form for the scalar subvector
          call linf_buildVectorScalar2(rform, .true.,&
              rvector%RvectorBlock(iblock), hydro_coeffVectorAnalytic, rcollection)

          ! Increase number of processed expressions
          nexpression = nexpression + 1
        else

          ! Convert scalar vector in interleaved format to true block vector
          call lsysbl_convertScalarBlockVector(&
              rvector%RvectorBlock(iblock), rvectorBlock)

          ! Loop over all blocks
          do ivar = 1, rvectorBlock%nblocks
            
            ! Get the function name of the component used for evaluating the initial solution.
            call parlst_getvalue_string(rparlist, ssectionName,&
                'ssolutionName', ssolutionName, isubstring=nexpression+ivar)
            
            ! Set the number of the component used for evaluating the initial solution
            rcollection%IquickAccess(1) = fparser_getFunctionNumber(p_rfparser, ssolutionname)

            ! Assemble the linear form for the scalar subvector
            call linf_buildVectorScalar2(rform, .true.,&
                rvectorBlock%RvectorBlock(ivar), hydro_coeffVectorAnalytic,&
                rcollection)
          end do

          ! Convert block vector back to scalar vector in interleaved format
          call lsysbl_convertBlockScalarVector(&
              rvectorBlock,rvector%RvectorBlock(iblock))
          
          ! Increase number of processed expressions
          nexpression = nexpression + rvectorBlock%nblocks
          
          ! Release temporal block vector
          call lsysbl_releaseVector(rvectorBlock)
        end if

      end do
      
      ! Store norm of load vector (if required)
      if (isolutionType .eq. SOLUTION_ANALYTIC_L2_CONSISTENT) then
        allocate(Dnorm0(rvector%nblocks))
        do iblock = 1, rvector%nblocks
          Dnorm0(iblock) = lsyssc_vectorNorm(&
              rvector%RvectorBlock(iblock), LINALG_NORML2)
        end do
      end if

      ! Compute the lumped L2-projection
      do iblock = 1, rvector%nblocks
        call lsyssc_invertedDiagMatVec(p_rlumpedMassMatrix,&
            rvector%RvectorBlock(iblock), 1.0_DP, rvector%RvectorBlock(iblock))
      end do

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
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'nsolutionfailsafe', nsolutionfailsafe, 0)

        ! Compute auxiliary vectors for high-order solution and increment
        call lsysbl_duplicateVector(rvector, rvectorHigh,&
            LSYSSC_DUP_TEMPLATE, LSYSSC_DUP_COPY)
        call lsysbl_duplicateVector(rvector, rvectorAux,&
            LSYSSC_DUP_TEMPLATE, LSYSSC_DUP_EMPTY)
        
        ! Compute the consistent L2-projection by Richardson iteration
        do iblock = 1, rvector%nblocks
          richardson: do iter = 1, nmaxIterationsSolution
            ! Compute the increment for each scalar subvector
            call lsyssc_scalarMatVec(p_rconsistentMassMatrix,&
                rvectorHigh%RvectorBlock(iblock),&
                rvectorAux%RvectorBlock(iblock), 1.0_DP, 0.0_DP)
            call lsyssc_invertedDiagMatVec(p_rlumpedMassMatrix,&
                rvectorAux%RvectorBlock(iblock), 1.0_DP,&
                rvectorAux%RvectorBlock(iblock))
            call lsyssc_vectorLinearComb(rvector%RvectorBlock(iblock),&
                rvectorAux%RvectorBlock(iblock), 1.0_DP, -1.0_DP)

            ! Update the scalar subvector of thesolution
            call lsyssc_vectorLinearComb(rvectorAux%RvectorBlock(iblock),&
                rvectorHigh%RvectorBlock(iblock), 1.0_DP, 1.0_DP)
            
            ! Check for convergence
            dnorm = lsyssc_vectorNorm(&
                rvectorAux%RvectorBlock(iblock), LINALG_NORML2)
            if ((dnorm .le. depsAbsSolution) .or.&
                (dnorm .le. depsRelSolution*Dnorm0(iblock))) exit richardson
          end do richardson
        end do
        
        ! Initialise stabilisation structure by hand
        rafcstab%istabilisationSpec= AFCSTAB_UNDEFINED
        rafcstab%bprelimiting = .false.
        rafcstab%ctypeAFCstabilisation = AFCSTAB_FEMFCT_MASS
        call gfsys_initStabilisation(&
            rproblemLevel%RmatrixBlock(systemMatrix), rafcstab)
        call afcstab_generateVerticesAtEdge(&
            p_rconsistentMassMatrix, rafcstab)
        
        ! Compute the raw antidiffusive mass fluxes. Note that we may supply any
        ! callback function for assembling the antidiffusive fluxes since it 
        ! will not be used for assembling antidiffusive mass fluxes !!!
        call gfsys_buildFluxFCT(rafcstab, rvectorHigh,&
            hydro_calcFluxFCTScDiss1d_sim, 0.0_DP, 0.0_DP, 1.0_DP, .true.,&
            p_rconsistentMassMatrix, rcollection=rcollection)
        
        ! Attach section name to collection structure
        rcollection%SquickAccess(1) = ssectionName

        if (nsolutionfailsafe .gt. 0) then

          ! Get number of failsafe variables
          nvariable = max(1,&
              parlst_querysubstrings(rparlist,&
              ssectionName, 'ssolutionfailsafevariable'))
          
          ! Allocate character array that stores all failsafe variable names
          allocate(SsolutionFailsafeVariables(nvariable))
          
          ! Initialize character array with failsafe variable names
          do ivariable = 1, nvariable
            call parlst_getvalue_string(rparlist,&
                ssectionName, 'ssolutionfailsafevariable',&
                SsolutionFailsafevariables(ivariable), isubstring=ivariable)
          end do

          ! Compute and apply FEM-FCT correction
          call hydro_calcCorrectionFCT(rproblemLevel, rvector, 1.0_DP,&
              .false., AFCSTAB_FCTALGO_STANDARD-AFCSTAB_FCTALGO_CORRECT,&
              rvector, rcollection, rafcstab, 'ssolutionconstrainvariable')
          
          ! Apply failsafe flux correction
          call afcstab_failsafeLimiting(rafcstab, p_rlumpedMassMatrix,&
              SsolutionFailsafeVariables, 1.0_DP,&
              nsolutionfailsafe, hydro_getVariable, rvector)
          
          ! Deallocate temporal memory
          deallocate(SsolutionFailsafeVariables)

        else
          
          ! Compute and apply FEM-FCT correction
          call hydro_calcCorrectionFCT(rproblemLevel, rvector, 1.0_DP,&
              .false., AFCSTAB_FCTALGO_STANDARD+AFCSTAB_FCTALGO_SCALEBYMASS,&
              rvector, rcollection, rafcstab, 'ssolutionconstrainvariable')
        end if
        
        ! Release stabilisation structure
        call afcstab_releaseStabilisation(rafcstab)

        ! Release auxiliary vectors
        call lsysbl_releaseVector(rvectorHigh)
        call lsysbl_releaseVector(rvectorAux)

        ! Release temporal memory
        deallocate(Dnorm0)
      end if

      ! Release temporal matrices (if any)
      call lsyssc_releaseMatrix(rconsistentMassMatrix)
      call lsyssc_releaseMatrix(rlumpedMassMatrix)
        

    case DEFAULT
      call output_line('Invalid type of solution profile!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'hydro_initSolution')
      call sys_halt()
    end select
    
  end subroutine hydro_initSolution

  !*****************************************************************************

!<subroutine>

  subroutine hydro_outputSolution(rparlist, ssectionName,&
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
    type(t_vectorScalar) :: rvector1, rvector2, rvector3
    real(DP), dimension(:), pointer :: p_Dsolution, p_Ddata1, p_Ddata2, p_Ddata3
    character(len=SYS_NAMELEN) :: cvariable
    integer :: iformatUCD, isystemFormat, isize, ndim, nvar, ivariable, nvariable


    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'output', soutputName)
    call parlst_getvalue_string(rparlist,&
        trim(soutputName), 'sucdsolution', sucdsolution)
    call parlst_getvalue_int(rparlist,&
        trim(soutputName), 'iformatucd', iformatUCD)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'isystemformat', isystemformat)

    ! Initialize the UCD exporter
    call flagship_initUCDexport(rproblemLevel, sucdsolution,&
        iformatUCD, rexport, ifilenumber)

    ! Increase filenumber by one
    ifilenumber = ifilenumber+1

    ! Set simulation time
    if (present(dtime)) call ucd_setSimulationTime(rexport, dtime)

    ! Add primal solution vector
    if (present(rsolutionPrimal)) then
      
      ! Set pointers
      call lsysbl_getbase_double(rsolutionPrimal, p_Dsolution)
      nvar  = hydro_getNVAR(rproblemLevel)
      isize = size(p_Dsolution)/nvar
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
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_outputSolution')
        call sys_halt()
      end select

      ! Get number of variables to be written
      nvariable = max(1,&
          parlst_querysubstrings(rparlist,&
          trim(soutputName), 'sucdvariable'))
      
      select case(isystemFormat)
      case(SYSTEM_INTERLEAVEFORMAT)
        
        ! Loop over all variables
        do ivariable = 1, nvariable
          
          ! Get variable name
          call parlst_getvalue_string(rparlist, trim(soutputName),&
              'sucdvariable', cvariable, isubstring=ivariable)
          
          if (trim(cvariable) .eq. 'velocity') then
            
            ! Special treatment of velocity vector
            select case(ndim)
            case (NDIM1D)
              call hydro_getVarInterleaveFormat(rvector1%NEQ, NVAR1D,&
                  'velocity_x', p_Dsolution, p_Ddata1)
              call ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1)

            case (NDIM2D)
              call hydro_getVarInterleaveFormat(rvector1%NEQ, NVAR2D,&
                  'velocity_x', p_Dsolution, p_Ddata1)
              call hydro_getVarInterleaveFormat(rvector2%NEQ, NVAR2D,&
                  'velocity_y', p_Dsolution, p_Ddata2)
              call ucd_addVarVertBasedVec(rexport, 'velocity',&
                  p_Ddata1, p_Ddata2)

            case (NDIM3D)
              call hydro_getVarInterleaveFormat(rvector1%NEQ, NVAR3D,&
                  'velocity_x', p_Dsolution, p_Ddata1)
              call hydro_getVarInterleaveFormat(rvector2%NEQ, NVAR3D,&
                  'velocity_y', p_Dsolution, p_Ddata2)
              call hydro_getVarInterleaveFormat(rvector3%NEQ, NVAR3D,&
                  'velocity_z', p_Dsolution, p_Ddata3)
              call ucd_addVarVertBasedVec(rexport, 'velocity',&
                  p_Ddata1, p_Ddata2, p_Ddata3)
            end select

          elseif (trim(cvariable) .eq. 'momentum') then
            
            ! Special treatment of momentum vector
            select case(ndim)
            case (NDIM1D)
              call hydro_getVarInterleaveFormat(rvector1%NEQ, NVAR1D,&
                  'momentum_x', p_Dsolution, p_Ddata1)
              call ucd_addVarVertBasedVec(rexport, 'momentum', p_Ddata1)

            case (NDIM2D)
              call hydro_getVarInterleaveFormat(rvector1%NEQ, NVAR2D,&
                  'momentum_x', p_Dsolution, p_Ddata1)
              call hydro_getVarInterleaveFormat(rvector2%NEQ, NVAR2D,&
                  'momentum_y', p_Dsolution, p_Ddata2)
              call ucd_addVarVertBasedVec(rexport, 'momentum',&
                  p_Ddata1, p_Ddata2)

            case (NDIM3D)
              call hydro_getVarInterleaveFormat(rvector1%NEQ, NVAR3D,&
                  'momentum_x', p_Dsolution, p_Ddata1)
              call hydro_getVarInterleaveFormat(rvector2%NEQ, NVAR3D,&
                  'momentum_y', p_Dsolution, p_Ddata2)
              call hydro_getVarInterleaveFormat(rvector3%NEQ, NVAR3D,&
                  'momentum_z', p_Dsolution, p_Ddata3)
              call ucd_addVarVertBasedVec(rexport, 'momentum',&
                  p_Ddata1, p_Ddata2, p_Ddata3)
            end select

          else

            ! Standard treatment for scalar quantity
            call hydro_getVarInterleaveFormat(rvector1%NEQ,  nvar,&
                cvariable, p_Dsolution, p_Ddata1)
            call ucd_addVariableVertexBased (rexport, cvariable,&
                UCD_VAR_STANDARD, p_Ddata1)
            
          end if
        end do
        
      case (SYSTEM_BLOCKFORMAT)

        ! Loop over all variables
        do ivariable = 1, nvariable
          
          ! Get variable name
          call parlst_getvalue_string(rparlist, trim(soutputName),&
              'sucdvariable', cvariable, isubstring=ivariable)
          
          if (trim(cvariable) .eq. 'velocity') then

            ! Special treatment of velocity vector
            select case(ndim)
            case (NDIM1D)
              call hydro_getVarBlockFormat(rvector1%NEQ, NVAR1D,&
                  'velocity_x', p_Dsolution, p_Ddata1)
              call ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1)

            case (NDIM2D)
              call hydro_getVarBlockFormat(rvector1%NEQ, NVAR2D,&
                  'velocity_x', p_Dsolution, p_Ddata1)
              call hydro_getVarBlockFormat(rvector2%NEQ, NVAR2D,&
                  'velocity_y', p_Dsolution, p_Ddata2)
              call ucd_addVarVertBasedVec(rexport, 'velocity',&
                  p_Ddata1, p_Ddata2)

            case (NDIM3D)
              call hydro_getVarBlockFormat(rvector1%NEQ, NVAR3D,&
                  'velocity_x', p_Dsolution, p_Ddata1)
              call hydro_getVarBlockFormat(rvector2%NEQ, NVAR3D,&
                  'velocity_y', p_Dsolution, p_Ddata2)
              call hydro_getVarBlockFormat(rvector3%NEQ, NVAR3D,&
                  'velocity_z', p_Dsolution, p_Ddata3)
              call ucd_addVarVertBasedVec(rexport, 'velocity',&
                  p_Ddata1, p_Ddata2, p_Ddata3)
            end select
            
          elseif (trim(cvariable) .eq. 'momentum') then

            ! Special treatment of momentum vector
            select case(ndim)
            case (NDIM1D)
              call hydro_getVarBlockFormat(rvector1%NEQ, NVAR1D,&
                  'momentum_x', p_Dsolution, p_Ddata1)
              call ucd_addVarVertBasedVec(rexport, 'momentum', p_Ddata1)

            case (NDIM2D)
              call hydro_getVarBlockFormat(rvector1%NEQ, NVAR2D,&
                  'momentum_x', p_Dsolution, p_Ddata1)
              call hydro_getVarBlockFormat(rvector2%NEQ, NVAR2D,&
                  'momentum_y', p_Dsolution, p_Ddata2)
              call ucd_addVarVertBasedVec(rexport, 'momentum',&
                  p_Ddata1, p_Ddata2)

            case (NDIM3D)
              call hydro_getVarBlockFormat(rvector1%NEQ, NVAR3D,&
                  'momentum_x', p_Dsolution, p_Ddata1)
              call hydro_getVarBlockFormat(rvector2%NEQ, NVAR3D,&
                  'momentum_y', p_Dsolution, p_Ddata2)
              call hydro_getVarBlockFormat(rvector3%NEQ, NVAR3D,&
                  'momentum_z', p_Dsolution, p_Ddata3)
              call ucd_addVarVertBasedVec(rexport, 'momentum',&
                  p_Ddata1, p_Ddata2, p_Ddata3)
            end select

          else
            
            ! Standard treatment for scalar quantity
            call hydro_getVarBlockFormat(rvector1%NEQ, nvar,&
                cvariable, p_Dsolution, p_Ddata1)
            call ucd_addVariableVertexBased (rexport, cvariable,&
                UCD_VAR_STANDARD, p_Ddata1)
            
          end if
        end do

      case default
        call output_line('Invalid system format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_outputSolution')
        call sys_halt()
      end select
      
      ! Release temporal memory
      call lsyssc_releaseVector(rvector1)
      call lsyssc_releaseVector(rvector2)
      call lsyssc_releaseVector(rvector3)

    end if

    ! Write UCD file
    call ucd_write(rexport)
    call ucd_release(rexport)

  end subroutine hydro_outputSolution

  !*****************************************************************************

!<subroutine>

  subroutine hydro_outputStatistics(rtimerTotal, rcollection)

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
    call output_line('Time for vector assembly:       '//&
                     trim(adjustl(sys_sdE(p_rtimerAssemblyVector%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerAssemblyVector%delapsedCPU, 5)))//' %')
    call output_line('Time for pre-/post-processing : '//&
                     trim(adjustl(sys_sdE(p_rtimerPrePostprocess%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerPrePostprocess%delapsedCPU, 5)))//' %')
    call output_lbrk()
    call output_line('Time for total simulation     : '//&
                     trim(adjustl(sys_sdE(dtotalTime, 5))))
    call output_lbrk()

  end subroutine hydro_outputStatistics

  !*****************************************************************************

!<subroutine>

  subroutine hydro_estimateRecoveryError(rparlist, ssectionName,&
      rproblemLevel, rsolution, dtime, rerror, derror)

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

!<output>
    ! element-wise error distribution
    type(t_vectorScalar), intent(out) :: rerror

    ! global error
    real(DP), intent(out) :: derror
!</output>
!</subroutine>

    ! section names
    character(LEN=SYS_STRLEN) :: sindatfileName
    character(LEN=SYS_STRLEN) :: serrorestimatorName
    character(LEN=SYS_STRLEN) :: sexactsolutionName

    ! local variables
    type(t_vectorScalar) :: rvectorScalar, rvectorTmp
    real(DP), dimension(:), pointer :: p_Ddata, p_DdataTmp
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:,:), pointer :: p_IneighboursAtElement
    logical, dimension(:), pointer :: p_BisactiveElement
    character(LEN=SYS_STRLEN) :: serrorvariable
    real(DP) :: dnoiseFilter, dabsFilter, dvalue,&
                dprotectLayerTolerance, derrorTmp
    integer :: ierrorEstimator, igridindicator, iexactsolutiontype
    integer :: iprotectLayer, nprotectLayers, ierrorVariable, nerrorVariables
    integer :: h_BisactiveElement


    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'indatfile', sindatfileName)
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'errorestimator', serrorestimatorName)
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'sexactsolutionname', sexactsolutionName, '')
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'iexactsolutiontype', iexactsolutiontype, 0)
    call parlst_getvalue_int(rparlist,&
        trim(serrorestimatorName), 'ierrorestimator', ierrorestimator)
    call parlst_getvalue_int(rparlist,&
        trim(serrorestimatorName), 'igridindicator', igridindicator)
    call parlst_getvalue_int(rparlist,&
        trim(serrorestimatorName), 'nprotectLayers', nprotectLayers, 0)
    call parlst_getvalue_double(rparlist,&
        trim(serrorestimatorName), 'dprotectLayerTolerance',&
        dprotectLayerTolerance, 0.0_DP)


    !---------------------------------------------------------------------------
    ! Perform recovery-based error estimation
    !---------------------------------------------------------------------------

    nerrorVariables = parlst_querysubstrings(rparlist,&
        trim(serrorestimatorName), 'serrorvariable')

    ! Loop over all error variables
    do ierrorVariable = 1, nerrorVariables

      ! Get name of error variable
      call parlst_getvalue_string(rparlist, trim(serrorestimatorName),&
          'serrorvariable', serrorVariable, isubString=ierrorVariable)

      ! Extract scalar variable from vector of conservative variables
      call hydro_getVariable(rsolution, serrorVariable, rvectorScalar)

      ! What type of error estimator are we?
      select case(ierrorEstimator)

      case (ERREST_L2PROJECTION)
        call lsyssc_createVector(rvectorTmp,&
            rproblemLevel%rtriangulation%NEL, .false.)
        call ppgrd_calcGradientError(rvectorScalar, derrorTmp,&
            PPGRD_INTERPOL, 0, rvectorTmp)

      case (ERREST_SPR_VERTEX)
        call lsyssc_createVector(rvectorTmp,&
            rproblemLevel%rtriangulation%NEL, .false.)
        call ppgrd_calcGradientError(rvectorScalar, derrorTmp,&
            PPGRD_ZZTECHNIQUE, PPGRD_NODEPATCH, rvectorTmp)

      case (ERREST_SPR_ELEMENT)
        call lsyssc_createVector(rvectorTmp,&
            rproblemLevel%rtriangulation%NEL, .false.)
        call ppgrd_calcGradientError(rvectorScalar, derrorTmp,&
            PPGRD_ZZTECHNIQUE, PPGRD_ELEMPATCH, rvectorTmp)

      case (ERREST_SPR_FACE)
        call lsyssc_createVector(rvectorTmp,&
            rproblemLevel%rtriangulation%NEL, .false.)
        call ppgrd_calcGradientError(rvectorScalar, derrorTmp,&
            PPGRD_ZZTECHNIQUE, PPGRD_FACEPATCH, rvectorTmp)

      case (ERREST_LIMAVR)
        call lsyssc_createVector(rvectorTmp,&
            rproblemLevel%rtriangulation%NEL, .false.)
        call ppgrd_calcGradientError(rvectorScalar, derrorTmp,&
            PPGRD_LATECHNIQUE, 0, rvectorTmp)

      case (ERREST_SECONDDIFF)
        call parlst_getvalue_double(rparlist,&
            trim(serrorestimatorName), 'dnoiseFilter', dnoiseFilter)
        call parlst_getvalue_double(rparlist,&
            trim(serrorestimatorName), 'dabsFilter', dabsFilter)
        call ppind_secondDifference(rvectorScalar, dnoiseFilter,&
            dabsFilter, rvectorTmp)

        ! This is no error estimator
        derrorTmp = 1.0

      case DEFAULT
        call output_line('Invalid type of error estimator!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_estimateRecoveryError')
        call sys_halt()
      end select


      ! Compute the root-mean square value
      call lsyssc_getbase_double(rvectorTmp, p_Ddata)
      dvalue = sqrt(sum(p_Ddata**2)/real(rvectorTmp%NEQ, DP))
      if (abs(dvalue) .gt. SYS_EPSREAL) then
        dvalue = 1.0_DP/dvalue
        call lsyssc_scaleVector(rvectorTmp, dvalue)
      end if

      ! Compute the global and element-wise error
      if (ierrorVariable .eq. 1) then

        ! Initialize the global error
        derror = derrorTmp

        ! Initialize the element-wise error
        call lsyssc_copyVector(rvectorTmp, rerror)
      else

        ! Update the global error
        derror = max(derror, derrorTmp)

        ! Update the element-wise error
        call lsyssc_getbase_double(rvectorTmp, p_DdataTmp)
        call lsyssc_getbase_double(rerror, p_Ddata)
        p_Ddata = max(p_Ddata, p_DdataTmp)
!!$        call lsyssc_vectorLinearComb(rvectorTmp, rerror, 1.0_DP, 1.0_DP)
      end if

      ! Release scalar variable and temporal error
      call lsyssc_releaseVector(rvectorScalar)
      call lsyssc_releaseVector(rvectorTmp)
    end do

!!$    ! Scale the global and element-wise error by the number of error variables
!!$    if (nerrorVariables .gt. 1) then
!!$      dvalue = 1.0_DP/real(nerrorVariables, DP)
!!$      derror = derror*dvalue
!!$      call lsyssc_scaleVector(rerror, dvalue)
!!$    end if


    !---------------------------------------------------------------------------
    ! Calculate protection layers
    !---------------------------------------------------------------------------

    if (nprotectLayers > 0) then

      ! Create auxiliary memory
      h_BisactiveElement = ST_NOHANDLE
      call storage_new('transp_estimateRecoveryError',' BisactiveElement',&
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

  end subroutine hydro_estimateRecoveryError

  !*****************************************************************************

!<subroutine>

  subroutine hydro_adaptTriangulation(rparlist, ssectionName, rhadapt ,&
      rtriangulationSrc, rindicator, rcollection, rtriangulationDest)

!<description>
    ! This subroutine performs h-adaptation for the given triangulation
!</description>

!<inputoutput>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName

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
    integer :: isystemFormat


    ! Check if adaptation structure has been prepared
    if (rhadapt%iSpec .eq. HADAPT_HAS_PARAMETERS) then

      ! Initialize adaptation structure from triangulation
      call hadapt_initFromTriangulation(rhadapt, rtriangulationSrc)

    else

      ! Refresh adaptation structure
      call hadapt_refreshAdaptation(rhadapt, rtriangulationSrc)

    end if

    ! Get parameters from parameter list
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'isystemformat', isystemFormat)

    ! What type of system format are we?
    select case (isystemFormat)

    case (SYSTEM_INTERLEAVEFORMAT)

      ! How many spatial dimensions do we have?
      select case(rtriangulationSrc%ndim)
      case (NDIM1D)
        call hydro_hadaptCallbackScalar1D(&
            HADAPT_OPR_INITCALLBACK, rcollection)
        call hadapt_performAdaptation(rhadapt, rindicator,&
            rcollection, hydro_hadaptCallbackScalar1D)

      case (NDIM2D)
        call hydro_hadaptCallbackScalar2D(&
            HADAPT_OPR_INITCALLBACK, rcollection)
        call hadapt_performAdaptation(rhadapt, rindicator,&
            rcollection, hydro_hadaptCallbackScalar2D)

      case (NDIM3D)
        call hydro_hadaptCallbackScalar3D(&
            HADAPT_OPR_INITCALLBACK, rcollection)
        call hadapt_performAdaptation(rhadapt, rindicator,&
            rcollection, hydro_hadaptCallbackScalar3D)
      end select

    case (SYSTEM_BLOCKFORMAT)

      ! How many spatial dimensions do we have?
      select case(rtriangulationSrc%ndim)
      case (NDIM1D)
        call hydro_hadaptCallbackBlock1D(&
            HADAPT_OPR_INITCALLBACK, rcollection)
        call hadapt_performAdaptation(rhadapt, rindicator,&
            rcollection, hydro_hadaptCallbackBlock1D)

      case (NDIM2D)
        call hydro_hadaptCallbackBlock2D(&
            HADAPT_OPR_INITCALLBACK, rcollection)
        call hadapt_performAdaptation(rhadapt, rindicator,&
            rcollection, hydro_hadaptCallbackBlock2D)

      case (NDIM3D)
        call hydro_hadaptCallbackBlock3D(&
            HADAPT_OPR_INITCALLBACK, rcollection)
        call hadapt_performAdaptation(rhadapt, rindicator,&
            rcollection, hydro_hadaptCallbackBlock3D)
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

  end subroutine hydro_adaptTriangulation

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

    ! primal solution vector
    type(t_vectorBlock), intent(inout), target :: rsolution

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
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

    ! Timer structures
    type(t_timer), pointer :: p_rtimerPrePostprocess
    type(t_timer), pointer :: p_rtimerSolution
    type(t_timer), pointer :: p_rtimerErrorEstimation
    type(t_timer), pointer :: p_rtimerAdaptation
    type(t_timer), pointer :: p_rtimerTriangulation
    type(t_timer), pointer :: p_rtimerAssemblyCoeff

    ! local variables
    type(t_ucdExport) :: rimport
    character(LEN=SYS_STRLEN) :: sadaptivityName
    character(LEN=SYS_STRLEN) :: soutputName
    character(LEN=SYS_STRLEN) :: sucdimport
    real(dp) :: derror, dstepUCD, dtimeUCD, dstepAdapt, dtimeAdapt
    integer :: templateMatrix, systemMatrix, isystemFormat, discretisation
    integer :: isize, ipreadapt, npreadapt, ndimension
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
        ssectionName, 'ndimension', ndimension)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'discretisation', discretisation)

    ! Set pointer to maximum problem level
    p_rproblemLevel   => rproblem%p_rproblemLevelMax
    p_rdiscretisation => p_rproblemLevel%Rdiscretisation(discretisation)

    ! Create the solution vector
    call lsysbl_createVectorBlock(p_rdiscretisation,&
        rsolution, .false., ST_DOUBLE)
    if (p_rdiscretisation%ncomponents .ne.&
        hydro_getNVAR(p_rproblemLevel)) then
      rsolution%RvectorBlock(1)%NVAR = hydro_getNVAR(p_rproblemLevel)
      isize = rsolution%NEQ*hydro_getNVAR(p_rproblemLevel)
      call lsysbl_resizeVectorBlock(rsolution, isize, .false., .false.)
    end if

    ! Initialize the solution vector and impose boundary conditions
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
            call hydro_estimateRecoveryError(rparlist, ssectionname,&
                p_rproblemLevel, rsolution, rtimestep%dinitialTime,&
                relementError, derror)

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

            ! Generate standard mesh from raw mesh
            call tria_initStandardMeshFromRaw(&
                p_rproblemLevel%rtriangulation, rproblem%rboundary)

            ! Update the template matrix according to the sparsity pattern
            call grph_generateMatrix(rgraph,&
                p_rproblemLevel%Rmatrix(templateMatrix))

            ! Re-initialize all constant coefficient matrices
            call hydro_initProblemLevel(rparlist,&
                ssectionName, p_rproblemLevel, rcollection)

            ! Resize the solution vector accordingly
            call parlst_getvalue_int(rparlist,&
                ssectionName, 'systemMatrix', systemMatrix)
            call lsysbl_resizeVecBlockIndMat(&
                p_rproblemLevel%RmatrixBlock(systemMatrix),&
                rsolution, .false., .true.)

            ! Re-generate the initial solution vector and impose
            ! boundary conditions explicitly
            call hydro_initSolution(rparlist, ssectionname,&
                p_rproblemLevel, rtimestep%dinitialTime, rsolution,&
                rcollection)

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
          call parlst_getvalue_int(rparlist,&
              ssectionName, 'systemMatrix', systemMatrix)
          call parlst_getvalue_int(rparlist,&
              ssectionName, 'isystemFormat', isystemFormat)
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

      ! Prepare quick access arrays of the collection
      rcollection%SquickAccess(1) = ssectionName

      ! What time-stepping scheme should be used?
      select case(rtimestep%ctimestepType)

      case (TSTEP_RK_SCHEME)

        ! Adopt explicit Runge-Kutta scheme
        call tstep_performRKStep(p_rproblemLevel, rtimestep, rsolver,&
            rsolution, hydro_nlsolverCallback, rcollection)

      case (TSTEP_THETA_SCHEME)

        ! Adopt two-level theta-scheme
        call tstep_performThetaStep(p_rproblemLevel, rtimestep,&
            rsolver, rsolution, hydro_nlsolverCallback, rcollection)

      case DEFAULT
        call output_line('Unsupported time-stepping algorithm!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_solveTransientPrimal')
        call sys_halt()
      end select

      ! Perform linearised FEM-FCT post-processing
      call hydro_calcLinearisedFCT(rbdrCond, p_rproblemLevel,&
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
        call hydro_estimateRecoveryError(rparlist, ssectionname,&
            p_rproblemLevel, rsolution, rtimestep%dTime,&
            relementError, derror)

        ! Stop time measurement for error estimation
        call stat_stopTimer(p_rtimerErrorEstimation)


        !-------------------------------------------------------------------------
        ! Perform h-adaptation
        !-------------------------------------------------------------------------

        ! Start time measurement for mesh adaptation
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
        call grph_generateMatrix(rgraph,&
            p_rproblemLevel%Rmatrix(templateMatrix))

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


        ! Start time measurement for generation of constant
        ! coefficient matrices
        call stat_startTimer(p_rtimerAssemblyCoeff, STAT_TIMERSHORT)

        ! Re-initialize all constant coefficient matrices
        call hydro_initProblemLevel(rparlist, ssectionName,&
            p_rproblemLevel, rcollection)

        ! Resize the solution vector accordingly
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'systemmatrix', systemMatrix)
        call lsysbl_resizeVecBlockIndMat(&
            p_rproblemLevel%RmatrixBlock(systemMatrix),&
            rsolution, .false., .true.)

        ! Prepare internal data arrays of the solver structure
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'isystemformat', isystemFormat)
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

  end subroutine hydro_solveTransientPrimal

  !*****************************************************************************

!<subroutine>

  subroutine hydro_projectSolution(rsourceVector, rdestVector)
    
!<description>
    ! This subroutine performs conservative projection of the given solution
    ! stored in rsourceVector to another FE-space and stores the result in
    ! rdestVector. An FCT algorithm is used ensure monotonicity preservation.
!</description>

!<input>
    ! Source vector
    type(t_vectorBlock), intent(inout), target :: rsourceVector
!</input>

!<inputoutput>
    ! Destination vector
    type(t_vectorBlock), intent(inout) :: rdestVector
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_linearForm) :: rform
    type(t_collection) :: rcollection
    type(t_matrixScalar) :: rmatrix1,rmatrix2
    type(t_vectorScalar) :: rvector
    real(DP), dimension(:), pointer :: p_Ddata
    real(DP) :: dmass1, dmass2
    integer :: iblock
    
    ! Set up the linear form
    rform%itermCount = 1
    rform%Idescriptors(1) = DER_FUNC
    
    ! Set up the collection structure
    call collct_init(rcollection)
    rcollection%IquickAccess(1) = SYSTEM_BLOCKFORMAT
    rcollection%p_rvectorQuickAccess1 => rsourceVector
    
    ! Assemble the linear form for destination vector
    do iblock = 1, rdestVector%nblocks
      
      ! Create sparsity pattern of mass matrices
      call bilf_createMatrixStructure(&
          rsourceVector%p_rblockDiscr%RspatialDiscr(1),&
          LSYSSC_MATRIX9, rmatrix1)
      call bilf_createMatrixStructure(&
          rdestVector%p_rblockDiscr%RspatialDiscr(1),&
          LSYSSC_MATRIX9, rmatrix2)

      ! Create mass matrices
      call stdop_assembleSimpleMatrix(rmatrix1, DER_FUNC, DER_FUNC, 1.0_DP, .true.)
      call stdop_assembleSimpleMatrix(rmatrix2, DER_FUNC, DER_FUNC, 1.0_DP, .true.)
      
      ! Compute the lumped mass matrices
      call lsyssc_lumpMatrixScalar(rmatrix1, LSYSSC_LUMP_DIAG)
      call lsyssc_lumpMatrixScalar(rmatrix2, LSYSSC_LUMP_DIAG)

      ! Set the number of the scalar subvector to the collection structure
      rcollection%IquickAccess(2) = iblock
      
      ! Assemble the linear form for the scalar subvector
      call linf_buildVectorScalar2(rform, .true.,&
          rdestVector%RvectorBlock(iblock), hydro_coeffVectorFE, rcollection)

      ! Compute the lumped L2-projection
      call lsyssc_invertedDiagMatVec(rmatrix2, rdestVector%RvectorBlock(iblock),&
          1.0_DP, rdestVector%RvectorBlock(iblock))
      
      ! Compute density-mass
      call lsyssc_duplicateVector(rsourceVector%RvectorBlock(iblock), rvector,&
          LSYSSC_DUP_TEMPLATE, LSYSSC_DUP_EMPTY)
      call lsyssc_scalarMatVec(rmatrix1, rsourceVector%RvectorBlock(iblock),&
          rvector, 1.0_DP, 0.0_DP)
      call lsyssc_getbase_double(rvector, p_Ddata)
      dmass1 = sum(p_Ddata)
      call lsyssc_releaseVector(rvector)
      
      ! Compute density-mass
      call lsyssc_duplicateVector(rdestVector%RvectorBlock(iblock), rvector,&
          LSYSSC_DUP_TEMPLATE, LSYSSC_DUP_EMPTY)
      call lsyssc_scalarMatVec(rmatrix2, rdestVector%RvectorBlock(iblock),&
          rvector, 1.0_DP, 0.0_DP)
      call lsyssc_getbase_double(rvector, p_Ddata)
      dmass2 = sum(p_Ddata)
      call lsyssc_releaseVector(rvector)
      
      ! Release matrices
      call lsyssc_releaseMatrix(rmatrix1)
      call lsyssc_releaseMatrix(rmatrix2)

    end do
    
    ! Release the collection structure
    call collct_done(rcollection)

  end subroutine hydro_projectSolution

  !*****************************************************************************
  ! AUXILIARY ROUTINES
  !*****************************************************************************

!<subroutine>

  subroutine hydro_parseCmdlArguments(rparlist)

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

  end subroutine hydro_parseCmdlArguments

  !*****************************************************************************

!<subroutine>

  subroutine hydro_adjustParameterlist(rparlist, ssectionName)

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

    ! Check which coefficient matrices for inviscid part need to be build
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

  end subroutine hydro_adjustParameterlist

end module hydro_application
