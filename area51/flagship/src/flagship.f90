!##############################################################################
!# ****************************************************************************
!# <name> afc </name>
!# ****************************************************************************
!#
!# <purpose>
!# This application solves the generic time-dependent conservation law
!#
!#   $$\frac{\partial u}{\partial t}+\nabla\cdot{\bf f}(u)=0$$
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
!# </purpose>
!##############################################################################

program afc

  use collection
  use externalstorage
  use fparser
  use fsystem
  use genoutput
  use hadaptaux1d
  use hadaptaux2d
  use hadaptaux3d
  use hadaptivity
  use linearalgebra
  use linearformevaluation
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use signal
  use statistics
  use storage
  use ucd

  use afc_adaptation
  use afc_basic
  use afc_callback
  use afc_estimation
  use afc_init
  use afc_postprocessing
  use afc_problem
  use boundaryfilter
  use errorestimation
  use problem
  use solver
  use timestep

  implicit none

!<globals>

  !****************************************************************************
  ! Global structures required for this application

  ! Collection structure
  type(t_collection) :: rcollection

  ! Structure for primal problem
  type(t_problem) :: rproblem
  
  ! Time-stepping structure
  type(t_timestep) :: rtimestep

  ! Global solver structure
  type(t_solver), target :: rsolver

  ! Global (primal + dual) solution vectors
  type(t_vectorBlock) :: rprimalSolution
  type(t_vectorBlock) :: rdualSolution

  ! Global (primal + dual) boundary conditions
  type(t_boundaryCondition), dimension(2) :: RboundaryCondition

  ! Global right-hand side vector
  type(t_vectorBlock) :: rf

  ! Structure for grid adaptivity
  type(t_hadapt) :: rhadapt

  ! Indicator for grid adaptation
  type(t_vectorScalar) :: rgridIndicator

  ! Structure for error estimation
  type(t_errorEstimator) :: rerrorEstimator


  !****************************************************************************
  ! Parameter settings required for this application
  
  ! name of the parameter file containing the problem configuration
  character(LEN=SYS_STRLEN) :: parmfile = ''

  ! section name for benchmark configuration
  character(LEN=SYS_STRLEN) :: sbenchmarkName

  ! section name for file input/output
  character(LEN=SYS_STRLEN) :: sinputoutputName

  ! section name for the top-level solver
  character(LEN=SYS_STRLEN) :: ssolverName

  ! section name for time-stepping scheme
  character(LEN=SYS_STRLEN) :: stimestepName
  
  ! section name for adaptivity
  character(LEN=SYS_STRLEN) :: sadaptivityName
  
  ! section name for error estimator
  character(LEN=SYS_STRLEN) :: serrorestName
  
  
  ! Interruption by user
  integer, external :: signal_SIGINT
  integer, external :: signal_SIGQUIT
 
!</globals>


  ! Initialization
  call afc_initialize()

  ! Solution procedure
  if (rtimestep%dfinalTime > 0) then
    
    select case(iflowtype)
      
    case (FLOW_TRANSIENT)
      call afc_solveTransient()

    case (FLOW_PSEUDOTRANSIENT)
      call afc_solvePseudoTransient()
      
    case (FLOW_STEADYSTATE)
!!$ call afc_solveSteadystate()
      call afc_solvePrimalDual()
      
    case DEFAULT
      call output_line('Unsupported flow type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'afc')
      call sys_halt()
    end select
  end if


  ! Postprocessing
  call afc_postprocess
  
  ! Stop timer for total time
  call stat_stopTimer(rtimer_total)

  ! Statistics
  call afc_statistics()

  ! Finalization
  call afc_finalize()
  
  ! That's it.


contains
 
  !*****************************************************************************

!<subroutine>

  subroutine afc_initialize()

!<description>
    ! This subroutine performs all initialization tasks for this application.
    ! In particular, the storage subsystem is initialized and the command line
    ! parameters and those from the parameter file are parsed. The global data
    ! structures are generated and filled with data.
!</description>
!</subroutine>

    ! Entries in the diffusion operator
    real(DP), dimension(:,:), allocatable :: Ddiffusion
    
    ! Rotation for anisotropic diffusion
    real(DP) :: ddiffusionRotation

    ! Type of finite elements
    integer :: ieltype

    ! Type of matrix format
    integer :: imatrixFormat

    ! Convert grid to triangular grid
    integer :: iconvToTria
    
    ! Minimum/maximum multigrid levels
    integer :: nlmin, nlmax

    ! Local parameters
    character(LEN=SYS_STRLEN) :: trifile,prmfile,indatfile,pgmfile,tmpdir

    ! Spatial dimension
    integer :: ndimension

    ! Error variable
    logical :: berror


    ! Initialize Feat subsystem
    call system_init()
    sys_haltmode = SYS_HALT_THROWFPE
    
    ! Initialize storage subsystem
    call storage_init(300, 100)  
    
    ! Clear all timers
    call stat_clearTimer(rtimer_total)
    call stat_clearTimer(rtimer_solution)
    call stat_clearTimer(rtimer_adaptivity)
    call stat_clearTimer(rtimer_errorestimation)
    call stat_clearTimer(rtimer_triangulation)
    call stat_clearTimer(rtimer_assembly_coeff)
    call stat_clearTimer(rtimer_assembly_matrix)
    call stat_clearTimer(rtimer_assembly_resrhs)
    call stat_clearTimer(rtimer_prepostprocess)
    
    ! Start total time measurement
    call stat_startTimer(rtimer_total)

    ! Start time measurement for pre-processing
    call stat_startTimer(rtimer_prepostprocess, STAT_TIMERSHORT)
    
    ! Initialize collection
    call collct_init(rcollection)
    
    ! Initialize signal handler for SIGINT and SIGQUIT
    call fsignal(SIGINT, signal_SIGINT)
    call fsignal(SIGQUIT, signal_SIGQUIT)
    
    ! Initialize data input files
    call afc_UserInterface(0)
    
    ! Get parameters from file
    call parlst_init(rparlist)
    call parlst_readfromfile(rparlist, trim(adjustl(parmfile)))
    call afc_getParameters(rparlist)
    
    ! Initialize external storage
    call exstor_init (300, 100)
    call parlst_getvalue_string(rparlist, trim(adjustl(sinputoutputName)),&
                                "tmpdir", tmpdir, ' ')
    if (trim(adjustl(tmpdir)) .ne. '')&
        call exstor_attachDirectory(trim(adjustl(tmpdir)))

    ! Initialize h-adaptivity
    if (trim(adjustl(sadaptivityName)) .ne. '')&
        call hadapt_initFromParameterlist(rhadapt, rparlist, sadaptivityName)

    ! Initialize error estimator
    if (trim(adjustl(serrorestName)) .ne. '')&
        call errest_initErrorEstimator(rerrorEstimator, rparlist, serrorestName)
    
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

    ! Get global minimum/maximum multigrid level
    nlmin = solver_getNLMIN(rsolver)
    nlmax = solver_getNLMAX(rsolver)

    ! Initialize problem structure
    call parlst_getvalue_string(rparlist, trim(adjustl(sinputoutputName)),&
                                "indatfile", indatfile)
    call parlst_getvalue_string(rparlist, trim(adjustl(sinputoutputName)),&
                                  "trifile", trifile, '')
    call parlst_getvalue_string(rparlist, trim(adjustl(sinputoutputName)),&
                                "prmfile", prmfile, '')
    call parlst_getvalue_int(rparlist, trim(adjustl(sbenchmarkName)),&
                             "iconvToTria", iconvToTria, 0)
    call parlst_getvalue_int(rparlist, '', "ndimension", ndimension)
    
    call afc_initProblem(rproblem, trifile, prmfile, indatfile,&
                         nlmin, nlmax, iconvToTria, ndimension)
    
    call bdrf_readBoundaryCondition(RboundaryCondition(CDEQ_BOUNDARY_PRIMAL),&
                                    indatfile, '[boundary_conditions_primal]',&
                                    ndimension)
    call bdrf_readBoundaryCondition(RboundaryCondition(CDEQ_BOUNDARY_DUAL),&
                                    indatfile, '[boundary_conditions_dual]',&
                                    ndimension, berror)
    
    ! Initialize global solver structure. Note that this function call is required to
    ! prepare all subnodes of the entire solver structure prior to attaching matrices.
    call solver_updateStructure(rsolver)
    call solver_setBoundaryCondition(rsolver,&
                                     RboundaryCondition(CDEQ_BOUNDARY_PRIMAL), .true.)
    
    ! Initialize function parser
    call fparser_init()
    call afc_initParser(indatfile)
    
    ! Get matrix format and element type
    call parlst_getvalue_int(rparlist, trim(adjustl(sbenchmarkName)),&
                             "imatrixformat", imatrixFormat, 9)
    call parlst_getvalue_int(rparlist, trim(adjustl(sbenchmarkName)),&
                             "ieltype", ieltype)

    ! What dimension are we?
    select case(ndimension)

    case(NDIM1D)
      ! Allocate 1x1 diffusion matrix
      allocate(Ddiffusion(1,1))

      ! Initialize the diffusion matrix
      call parlst_getvalue_double(rparlist, trim(adjustl(sbenchmarkName)),&
                                  "ddiffusionA11", Ddiffusion(1,1), 0._DP)

      call afc_initDiffusion1D(Ddiffusion)

      ! Deallocate temporal array
      deallocate(Ddiffusion)

      ! Initialize constant coefficient operators in 1D
      call afc_initConstOperators1D(RboundaryCondition(CDEQ_BOUNDARY_PRIMAL),&
                                    rproblem, ieltype, imatrixFormat, nlmin, nlmax)

      ! Initialize global solution from file
      call afc_initSolutionFromParser(RboundaryCondition(CDEQ_BOUNDARY_PRIMAL),&
                                      rproblem%p_rproblemLevelMax, rprimalSolution,&
                                      indatfile, rtimestep%dinitialTime)
      
    case (NDIM2D)
      ! Allocate 2x2 diffusion matrix
      allocate(Ddiffusion(2,2))

      ! Initialize the diffusion matrix
      call parlst_getvalue_double(rparlist, trim(adjustl(sbenchmarkName)),&
                                  "ddiffusionA11", Ddiffusion(1,1), 0._DP)
      call parlst_getvalue_double(rparlist, trim(adjustl(sbenchmarkName)),&
                                  "ddiffusionA12", Ddiffusion(1,2), 0._DP)
      call parlst_getvalue_double(rparlist, trim(adjustl(sbenchmarkName)),&
                                  "ddiffusionA21", Ddiffusion(2,1), 0._DP)
      call parlst_getvalue_double(rparlist, trim(adjustl(sbenchmarkName)),&
                                  "ddiffusionA22", Ddiffusion(2,2), 0._DP)
      call parlst_getvalue_double(rparlist, trim(adjustl(sbenchmarkName)),&
                                  "ddiffusionRotation", ddiffusionRotation, 0._DP)

      call afc_initDiffusion2D(Ddiffusion, ddiffusionRotation)

      ! Deallocate temporal array
      deallocate(Ddiffusion)
      
      ! Initialize constant coefficient operators in 2D
      call afc_initConstOperators2D(RboundaryCondition(CDEQ_BOUNDARY_PRIMAL),&
                                    rproblem, ieltype, imatrixFormat, nlmin, nlmax)

      ! Initialize global solution from file or graymap image
      call parlst_getvalue_string(rparlist, trim(adjustl(sinputoutputName)),&
                                  "pgmfile", pgmfile, ' ')
      if (trim(adjustl(pgmfile)) .eq. '') then
        call afc_initSolutionFromParser(RboundaryCondition(CDEQ_BOUNDARY_PRIMAL),&
                                        rproblem%p_rproblemLevelMax, rprimalSolution,&
                                        indatfile, rtimestep%dinitialTime)
      else
        call afc_initSolutionFromImage(RboundaryCondition(CDEQ_BOUNDARY_PRIMAL),&
                                       rproblem%p_rproblemLevelMax, rprimalSolution, pgmfile)
      end if

    case (NDIM3D)
      ! Allocate 3x3 diffusion matrix
      allocate(Ddiffusion(3,3))

      ! Initialize the diffusion matrix
      call parlst_getvalue_double(rparlist, trim(adjustl(sbenchmarkName)),&
                                  "ddiffusionA11", Ddiffusion(1,1), 0._DP)
      call parlst_getvalue_double(rparlist, trim(adjustl(sbenchmarkName)),&
                                  "ddiffusionA12", Ddiffusion(1,2), 0._DP)
      call parlst_getvalue_double(rparlist, trim(adjustl(sbenchmarkName)),&
                                  "ddiffusionA13", Ddiffusion(1,3), 0._DP)
      call parlst_getvalue_double(rparlist, trim(adjustl(sbenchmarkName)),&
                                  "ddiffusionA21", Ddiffusion(2,1), 0._DP)
      call parlst_getvalue_double(rparlist, trim(adjustl(sbenchmarkName)),&
                                  "ddiffusionA22", Ddiffusion(2,2), 0._DP)
      call parlst_getvalue_double(rparlist, trim(adjustl(sbenchmarkName)),&
                                  "ddiffusionA23", Ddiffusion(2,3), 0._DP)
      call parlst_getvalue_double(rparlist, trim(adjustl(sbenchmarkName)),&
                                  "ddiffusionA31", Ddiffusion(3,1), 0._DP)
      call parlst_getvalue_double(rparlist, trim(adjustl(sbenchmarkName)),&
                                  "ddiffusionA32", Ddiffusion(3,2), 0._DP)
      call parlst_getvalue_double(rparlist, trim(adjustl(sbenchmarkName)),&
                                  "ddiffusionA33", Ddiffusion(3,3), 0._DP)
      call parlst_getvalue_double(rparlist, trim(adjustl(sbenchmarkName)),&
                                  "ddiffusionRotation", ddiffusionRotation, 0._DP)

      call afc_initDiffusion3D(Ddiffusion, ddiffusionRotation)

      ! Deallocate temporal array
      deallocate(Ddiffusion)

      ! Initialize constant coefficient operators in 3D
      call afc_initConstOperators3D(RboundaryCondition(CDEQ_BOUNDARY_PRIMAL),&
                                    rproblem, ieltype, imatrixFormat, nlmin, nlmax)

      ! Initialize global solution from file
      call afc_initSolutionFromParser(RboundaryCondition(CDEQ_BOUNDARY_PRIMAL),&
                                      rproblem%p_rproblemLevelMax, rprimalSolution,&
                                      indatfile, rtimestep%dinitialTime)
      
    case default
      call output_line('Unsupported spatial dimension!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'afc_initialize')
      call sys_halt()
    end select
    
    ! Initialize global right-hand side
    call afc_initRHS(RboundaryCondition(CDEQ_BOUNDARY_PRIMAL),&
                     rproblem%p_rproblemLevelMax, rf, indatfile)
    
    ! Initialize velocity vector
    call afc_initVelocity(rproblem, indatfile, nlmin, nlmax)
    
    ! Initialize global system operator. Actually, this is a dummy call to assemble
    ! all matrix structures for all possible multigrid levels throughout the complete
    ! solver structure. In case of constant velocity, no more matrix assemblies and/or
    ! updates of the solver structure need to be performed.
    call afc_calcPrecond(rproblem%p_rproblemLevelMax, rtimestep, rsolver,&
                         rprimalSolution, 1, nlmin, nlmax)
    call solver_updateStructure(rsolver)
    
    ! Stop time measurement for pre-processing
    call stat_stopTimer(rtimer_prepostprocess)

  end subroutine afc_initialize
  
  !*****************************************************************************

!<subroutine>

  subroutine afc_solveTransient()

!<description>
    ! This subroutine solves the transient flow problem
    !
    !  $$\frac{\partial u}{\partial t}+\nabla\cdot{\bf f}(u)=0$$
    !
    ! for a scalar quantity $u$ in the domain $\Omega$.
!</description>
!</subroutine>

    ! Name of greymap file
    character(LEN=SYS_STRLEN) :: pgmfile

    ! Name of bechmark file
    character(LEN=SYS_STRLEN) :: indatfile

    ! Name of UCD output file for error estimator
    character(LEN=SYS_STRLEN) :: errorfile

    ! File names for solution output
    character(LEN=SYS_STRLEN) :: gmvfile, avsfile, vtkfile

    ! Time for solution output and grid adaptation
    real(DP) :: ttimegmv, ttimeavs, ttimevtk, ttimeadapt

    ! Time step for solution output and grid adaptation
    real(DP) :: tstepgmv, tstepavs, tstepvtk, tstepadapt

    ! Number of pre-adaptation steps
    integer :: npreadapt

    ! Type of UCD output format
    integer :: ierrorUCD

    ! Type of finite element
    integer :: ieltype

    ! local variables
    integer :: ilev,ipreadapt


    ! Retrieve parameters
    call parlst_getvalue_int(rparlist, trim(adjustl(sbenchmarkName)),&
                             "ieltype", ieltype)
    call parlst_getvalue_int(rparlist, trim(adjustl(sbenchmarkName)),&
                             "npreadapt", npreadapt, 0)
    call parlst_getvalue_int(rparlist, trim(adjustl(sinputoutputName)),&
                             "ierrorUCD", ierrorUCD, 1)
    call parlst_getvalue_string(rparlist, trim(adjustl(sinputoutputName)),&
                                "indatfile", indatfile)
    call parlst_getvalue_string(rparlist, trim(adjustl(sinputoutputName)),&
                                "pgmfile", pgmfile, ' ')
    call parlst_getvalue_string(rparlist, trim(adjustl(sinputoutputName)),&
                                "errorfile", errorfile, ' ')
    call parlst_getvalue_string(rparlist, trim(adjustl(sinputoutputName)),&
                                "gmvfile", gmvfile, ' ')
    call parlst_getvalue_string(rparlist, trim(adjustl(sinputoutputName)),&
                                "avsfile", avsfile, ' ')
    call parlst_getvalue_string(rparlist, trim(adjustl(sinputoutputName)),&
                                "vtkfile", vtkfile, ' ')
    call parlst_getvalue_double(rparlist, trim(adjustl(sinputoutputName)),&
                                "tstepgmv", tstepgmv, 0.0_DP)
    call parlst_getvalue_double(rparlist, trim(adjustl(sinputoutputName)),&
                                "ttimegmv", ttimegmv, 0.0_DP)
    call parlst_getvalue_double(rparlist, trim(adjustl(sinputoutputName)),&
                                "tstepavs", tstepavs, 0.0_DP)
    call parlst_getvalue_double(rparlist, trim(adjustl(sinputoutputName)),&
                                "ttimeavs", ttimeavs, 0.0_DP)
    call parlst_getvalue_double(rparlist, trim(adjustl(sinputoutputName)),&
                                "tstepvtk", tstepvtk, 0.0_DP)
    call parlst_getvalue_double(rparlist, trim(adjustl(sinputoutputName)),&
                                "ttimevtk", ttimevtk, 0.0_DP)


    ! Check if adaptivity and error estimator are available
    if (trim(adjustl(sadaptivityName)) .ne. '' .or.&
        trim(adjustl(serrorestName))   .ne. '') then

      ! Retrieve parameters for grid adaptation
      call parlst_getvalue_double(rparlist, trim(adjustl(sadaptivityName)),&
                                  "ttimeadapt", ttimeadapt, 0.0_DP)
      call parlst_getvalue_double(rparlist, trim(adjustl(sadaptivityName)),&
                                  "tstepadapt", tstepadapt, 0.0_DP)
      
      ! Perform required number of pre-adaptation steps
      do ipreadapt = 1, npreadapt
      
        ! Perform error estimation
        call afc_prepareErrorEstimator(rproblem%p_rproblemLevelMax, rprimalSolution, rerrorEstimator)
        call afc_performErrorEstimation(rerrorEstimator, rgridIndicator,&
                                        rhadapt%drefinementTolerance)
        call errest_clearErrorEstimator(rerrorEstimator)
      
        ! Output grid indicator (if required)
        if (trim(errorfile) .ne. '') then
          call afc_outputVectorScalar(rproblem%p_rproblemLevelMax, rgridIndicator,&
                                      rtimestep%dTime, errorfile, ierrorUCD, .false.)
        end if
      
        ! What dimension are we?
        select case(rproblem%p_rproblemLevelMax%rtriangulation%ndim)
        case (NDIM1D)
          ! Perform grid adaptation in 1D
          call afc_performGridAdaptation(rcollection, rhadapt, rproblem%p_rproblemLevelMax,&
                                         rprimalSolution, rgridIndicator, fcb_hadaptCallback1D)
          
          call lsyssc_releaseVector(rgridIndicator)
          call afc_reinitConstOperators1D(rcollection, rhadapt, rproblem%p_rproblemLevelMax,&
                                          rprimalSolution, ieltype, fcb_hadaptCallback1D)

          ! Re-initialize the solution
          call lsysbl_releaseVector(rprimalSolution)
          call afc_initSolutionFromParser(RboundaryCondition(CDEQ_BOUNDARY_PRIMAL),&
                                          rproblem%p_rproblemLevelMax, rprimalSolution,&
                                          indatfile, rtimestep%dinitialTime)

        case (NDIM2D)
          ! Perform grid adaptation in 2D
          call afc_performGridAdaptation(rcollection, rhadapt, rproblem%p_rproblemLevelMax,&
                                         rprimalSolution, rgridIndicator, fcb_hadaptCallback2D)
          call lsyssc_releaseVector(rgridIndicator)
          call afc_reinitConstOperators2D(rcollection, rhadapt, rproblem%p_rproblemLevelMax,&
                                          rprimalSolution, ieltype, fcb_hadaptCallback2D)

          ! Re-initialize the solution
          call lsysbl_releaseVector(rprimalSolution)
          if (trim(adjustl(pgmfile)) .eq. '') then
            call afc_initSolutionFromParser(RboundaryCondition(CDEQ_BOUNDARY_PRIMAL),&
                                            rproblem%p_rproblemLevelMax, rprimalSolution,&
                                            indatfile, rtimestep%dinitialTime)
          else
            call afc_initSolutionFromImage(RboundaryCondition(CDEQ_BOUNDARY_PRIMAL),&
                                           rproblem%p_rproblemLevelMax, rprimalSolution, pgmfile)
          end if

        case (NDIM3D)
          ! Perform grid adaptation in 3D
          call afc_performGridAdaptation(rcollection, rhadapt, rproblem%p_rproblemLevelMax,&
                                         rprimalSolution, rgridIndicator, fcb_hadaptCallback3D)
          call lsyssc_releaseVector(rgridIndicator)
          call afc_reinitConstOperators3D(rcollection, rhadapt, rproblem%p_rproblemLevelMax,&
                                          rprimalSolution, ieltype, fcb_hadaptCallback3D)
          
          ! Re-initialize the solution
          call lsysbl_releaseVector(rprimalSolution)
          call afc_initSolutionFromParser(RboundaryCondition(CDEQ_BOUNDARY_PRIMAL),&
                                          rproblem%p_rproblemLevelMax, rprimalSolution,&
                                          indatfile, rtimestep%dinitialTime)
        end select
            
        ! Re-initialize global right-hand side
        call lsysbl_releaseVector(rf)
        call afc_initRHS(RboundaryCondition(CDEQ_BOUNDARY_PRIMAL),&
                         rproblem%p_rproblemLevelMax, rf, indatfile)

        ! We must enforce velocity update so that all internal
        ! structure are regenerated for the current level
        bvelocityUpdate = .true.

        ! Update global system operator and solver structure
        ilev = rproblem%p_rproblemLevelMax%ilev
        call afc_calcPrecond(rproblem%p_rproblemLevelMax, rtimestep, rsolver,&
                             rprimalSolution, 1, ilev, ilev)
        call solver_updateStructure(rsolver)
      end do
    else
      tstepadapt = 0._DP
      ttimeadapt = 0._DP
    end if

    ! Initialize time
    if (tstepadapt > 0) ttimeadapt = max(ttimeadapt, tstepadapt, rtimestep%dTime)
    if (tstepgmv   > 0) ttimegmv   = max(ttimegmv,   tstepgmv,   rtimestep%dTime)
    if (tstepavs   > 0) ttimeavs   = max(ttimeavs,   tstepavs,   rtimestep%dTime)
    if (tstepvtk   > 0) ttimevtk   = max(ttimevtk,   tstepvtk,   rtimestep%dTime)

    ! Infinite time stepping loop
    timeloop: do

      ! Check for user interaction
      call afc_UserInterface(1)

      !-------------------------------------------------------------------------
      ! Advance solution in time
      !-------------------------------------------------------------------------
      
      ! Start time measurement for solution procedure
      call stat_startTimer(rtimer_solution, STAT_TIMERSHORT)
      
      ! What time-stepping scheme should be used?
      select case(rtimestep%ctimestepType)

      case (SV_RK_SCHEME)
        
        ! Adopt explicit Runge-Kutta scheme
        if (irhstype .eq. 0) then
          call timestep_performRKStep(rproblem%p_rproblemLevelMax, rtimestep, rsolver,&
                                      rprimalSolution, 1,&
                                      fcb_calcRHS, fcb_setBoundary)
        else
          call timestep_performRKStep(rproblem%p_rproblemLevelMax, rtimestep, rsolver,&
                                      rprimalSolution, 1,&
                                      fcb_calcRHS, fcb_setBoundary, rf)
        end if
        
        
      case (SV_THETA_SCHEME)
        
        ! Adopt two-level theta-scheme
        if (irhstype .eq. 0) then
          call timestep_performThetaStep(rproblem%p_rproblemLevelMax, rtimestep, rsolver,&
                                         rprimalSolution, 1,&
                                         fcb_calcResidual, fcb_calcJacobian,&
                                         fcb_applyJacobian, fcb_setBoundary)
        else
          call timestep_performThetaStep(rproblem%p_rproblemLevelMax, rtimestep, rsolver,&
                                         rprimalSolution, 1,&
                                         fcb_calcResidual, fcb_calcJacobian,&
                                         fcb_applyJacobian, fcb_setBoundary, rf)
        end if
        
        
      case DEFAULT
        call output_line('Unsupported time-stepping algorithm!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'afc_solveTransient')
        call sys_halt()
      end select

      ! Stop time measurement for solution procedure
      call stat_stopTimer(rtimer_solution)

      
      ! Reached final time?
      if (rtimestep%dTime .ge. rtimestep%dfinalTime) exit timeloop


      !-------------------------------------------------------------------------
      ! Postprocess intermediate solution
      !-------------------------------------------------------------------------
      
      ! GMV output
      if (tstepgmv > 0._DP .and. rtimestep%dTime .ge. ttimegmv) then
        call afc_outputSolution(rproblem%p_rproblemLevelMax, rprimalSolution,&
                                rtimestep%dTime, gmvfile, 1)
        ttimegmv = ttimegmv+tstepgmv
      end if
      
      ! AVS output
      if (tstepavs > 0._DP .and. rtimestep%dTime .ge. ttimeavs) then
        call afc_outputSolution(rproblem%p_rproblemLevelMax, rprimalSolution,&
                                rtimestep%dTime, avsfile, 2)
        ttimeavs = ttimeavs+tstepavs
      end if
      
      ! VTK output
      if (tstepvtk > 0._DP .and. rtimestep%dTime .ge. ttimevtk) then
        call afc_outputSolution(rproblem%p_rproblemLevelMax, rprimalSolution,&
                                rtimestep%dTime, vtkfile, 3)
        ttimevtk = ttimevtk+tstepvtk
      end if
      
      
      !-------------------------------------------------------------------------
      ! Adaptive grid refinement
      !-------------------------------------------------------------------------

      if (tstepadapt > 0 .and. rtimestep%dTime .ge. ttimeadapt) then
        
        ! Perform error estimation
        call afc_prepareErrorEstimator(rproblem%p_rproblemLevelMax, rprimalSolution, rerrorEstimator)
        call afc_performErrorEstimation(rerrorEstimator, rgridIndicator,&
                                        rhadapt%drefinementTolerance)
        call errest_clearErrorEstimator(rerrorEstimator)          
        
        ! Output grid indicator (if required)
        if (trim(errorfile) .ne. '') then
          call afc_outputVectorScalar(rproblem%p_rproblemLevelMax, rgridIndicator,&
                                      rtimestep%dTime, errorfile, ierrorUCD, .false.)
        end if
        
        ! What dimension are we?
        select case(rproblem%p_rproblemLevelMax%rtriangulation%ndim)
        case (NDIM1D)
          ! Perform grid adaptation in 1D
          call afc_performGridAdaptation(rcollection, rhadapt, rproblem%p_rproblemLevelMax,&
                                         rprimalSolution, rgridIndicator, fcb_hadaptCallback1D)
          call lsyssc_releaseVector(rgridIndicator)
          call afc_reinitConstOperators1D(rcollection, rhadapt, rproblem%p_rproblemLevelMax,&
                                          rprimalSolution, ieltype, fcb_hadaptCallback1D)

        case (NDIM2D)
          ! Perform grid adaptation in 2D
          call afc_performGridAdaptation(rcollection, rhadapt, rproblem%p_rproblemLevelMax,&
                                         rprimalSolution, rgridIndicator, fcb_hadaptCallback2D)
          call lsyssc_releaseVector(rgridIndicator)
          call afc_reinitConstOperators2D(rcollection, rhadapt, rproblem%p_rproblemLevelMax,&
                                          rprimalSolution, ieltype, fcb_hadaptCallback2D)

        case (NDIM3D)
          ! Perform grid adaptation in 3D
          call afc_performGridAdaptation(rcollection, rhadapt, rproblem%p_rproblemLevelMax,&
                                         rprimalSolution, rgridIndicator, fcb_hadaptCallback3D)
          call lsyssc_releaseVector(rgridIndicator)
          call afc_reinitConstOperators3D(rcollection, rhadapt, rproblem%p_rproblemLevelMax,&
                                          rprimalSolution, ieltype, fcb_hadaptCallback3D)
        end select
        
        ! We must enforce velocity update so that all internal
        ! structure are regenerated for the current level
        bvelocityUpdate = .true.

        ! Update global system operator and solver structure
        ilev = rproblem%p_rproblemLevelMax%ilev
        call afc_calcPrecond(rproblem%p_rproblemLevelMax, rtimestep, rsolver,&
                             rprimalSolution, 1, ilev, ilev)
        call solver_updateStructure(rsolver)
          
        ttimeadapt = ttimeadapt+tstepadapt
      end if

    end do timeloop

  end subroutine afc_solveTransient

  !*****************************************************************************

!<subroutine>

  subroutine afc_solvePseudoTransient()

!<description>
    ! This subroutine solves the pseudo-transient flow problem
    !
    ! $$\frac{\partial u}{\partial \tau}+\nabla\cdot{\bf f}(u)=0$$
    !
    ! for a scalar quantity $u$ in the domain $\Omega$.
!</description>
!</subroutine>

    ! Type of finite element
    integer :: ieltype
    
    ! File names for solution error output
    character(LEN=SYS_STRLEN) :: errorfile

    ! Type of UCD output format
    integer :: ierrorUCD

    ! Number of adaptation steps
    integer :: nadapt

    ! local variables
    integer :: ilev,iadapt


    ! Retrieve parameters
    call parlst_getvalue_int(rparlist, trim(adjustl(sbenchmarkName)),&
                             "ieltype", ieltype)
    call parlst_getvalue_int(rparlist, trim(adjustl(sbenchmarkName)),&
                             "nadapt", nadapt, 0)    
    call parlst_getvalue_string(rparlist, trim(adjustl(sinputoutputName)),&
                                "errorfile", errorfile, ' ')
    call parlst_getvalue_int(rparlist, trim(adjustl(sinputoutputName)),&
                             "ierrorUCD", ierrorUCD, 1)
    
    
    ! Adaptation loop
    adaptloop: do iadapt = 0, nadapt
      
      ! Reset time stepping algorithm
      call solver_resetTimestep(rtimestep, .true.)
      
      ! Infinite time stepping loop
      timeloop: do
        
        ! Check for user interaction
        call afc_UserInterface(1)
        
        !-----------------------------------------------------------------------
        ! Advance solution in time
        !-----------------------------------------------------------------------
        
        ! Start time measurement for solution procedure
        call stat_startTimer(rtimer_solution, STAT_TIMERSHORT)
        
        ! What time-stepping scheme should be used?
        select case(rtimestep%ctimestepType)
          
        case (SV_RK_SCHEME)
          
          ! Adopt explicit Runge-Kutta scheme
          if (irhstype .eq. 0) then
            call timestep_performRKStep(rproblem%p_rproblemLevelMax, rtimestep, rsolver,&
                                        rprimalSolution, 1,&
                                        fcb_calcRHS, fcb_setBoundary)
          else
            call timestep_performRKStep(rproblem%p_rproblemLevelMax, rtimestep, rsolver,&
                                        rprimalSolution, 1,&
                                        fcb_calcRHS, fcb_setBoundary, rf)
          end if
          
          
        case (SV_THETA_SCHEME)
          
          ! Adopt two-level theta-scheme
          if (irhstype .eq. 0) then
            call timestep_performThetaStep(rproblem%p_rproblemLevelMax, rtimestep, rsolver,&
                                           rprimalSolution, 1,&
                                           fcb_calcResidual, fcb_calcJacobian,&
                                           fcb_applyJacobian, fcb_setBoundary)
          else
            call timestep_performThetaStep(rproblem%p_rproblemLevelMax, rtimestep, rsolver,&
                                           rprimalSolution, 1,&
                                           fcb_calcResidual, fcb_calcJacobian,&
                                           fcb_applyJacobian, fcb_setBoundary, rf)
          end if
          
          
        case DEFAULT
          call output_line('Unsupported time-stepping algorithm!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'afc_solvePseudoTransient')
          call sys_halt()
        end select
        
        ! Stop time measurement for solution procedure
        call stat_stopTimer(rtimer_solution)
        
        
        ! Reached final time?
        if (rtimestep%dTime .ge. rtimestep%dfinalTime) exit timeloop
        
        ! Reached steady state limit?
        if (rtimestep%depsSteady > 0.0_DP) then

          ! Check if steady-state residual exceeds tolerance
          if ((rsolver%dfinalDefect < rsolver%dinitialDefect) .and.&
              ((iadapt .eq. nadapt) .and.&
               (rsolver%dinitialDefect .le.&
                rtimestep%dStep*rtimestep%depsSteady) .or.&
               (rsolver%dinitialDefect .le.&
                rtimestep%dStep*sqrt(rtimestep%depsSteady)))) exit timeloop
        end if
        
      end do timeloop
      
      !-------------------------------------------------------------------------
      ! Postprocess intermediate solution
      !-------------------------------------------------------------------------

      call afc_postprocess()


      !-------------------------------------------------------------------------
      ! Adaptive grid refinement
      !-------------------------------------------------------------------------
      
      if (trim(adjustl(sadaptivityName)) .eq. '' .or.&
                                  nadapt .eq. 0 ) exit adaptloop

      ! Perform error estimation
      call afc_prepareErrorEstimator(rproblem%p_rproblemLevelMax, rprimalSolution, rerrorEstimator)
      call afc_performErrorEstimation(rerrorEstimator, rgridIndicator,&
                                      rhadapt%drefinementTolerance)
      call errest_clearErrorEstimator(rerrorEstimator)          
      
      ! Output grid indicator (if required)
      if (trim(errorfile) .ne. '') then
        call afc_outputVectorScalar(rproblem%p_rproblemLevelMax, rgridIndicator,&
                                    rtimestep%dTime, errorfile, ierrorUCD, .false.)
      end if
        
      ! What dimension are we?
      select case(rproblem%p_rproblemLevelMax%rtriangulation%ndim)
      case (NDIM1D)
        ! Perform grid adaptation in 1D
        call afc_performGridAdaptation(rcollection, rhadapt, rproblem%p_rproblemLevelMax,&
                                       rprimalSolution, rgridIndicator, fcb_hadaptCallback1D)
        call lsyssc_releaseVector(rgridIndicator)
        call afc_reinitConstOperators1D(rcollection, rhadapt, rproblem%p_rproblemLevelMax,&
                                        rprimalSolution, ieltype, fcb_hadaptCallback1D)

      case (NDIM2D)
        ! Perform grid adaptation in 2D
        call afc_performGridAdaptation(rcollection, rhadapt, rproblem%p_rproblemLevelMax,&
                                       rprimalSolution, rgridIndicator, fcb_hadaptCallback2D)
        call lsyssc_releaseVector(rgridIndicator)
        call afc_reinitConstOperators2D(rcollection, rhadapt, rproblem%p_rproblemLevelMax,&
                                        rprimalSolution, ieltype, fcb_hadaptCallback2D)

      case (NDIM3D)
        ! Perform grid adaptation in 3D
        call afc_performGridAdaptation(rcollection, rhadapt, rproblem%p_rproblemLevelMax,&
                                       rprimalSolution, rgridIndicator, fcb_hadaptCallback3D)
        call lsyssc_releaseVector(rgridIndicator)
        call afc_reinitConstOperators3D(rcollection, rhadapt, rproblem%p_rproblemLevelMax,&
                                        rprimalSolution, ieltype, fcb_hadaptCallback3D)
      end select

      ! We must enforce velocity update so that all internal
      ! structure are regenerated for the current level
      bvelocityUpdate = .true.
      
      ! Update global system operator and solver structure
      ilev = rproblem%p_rproblemLevelMax%ilev
      call afc_calcPrecond(rproblem%p_rproblemLevelMax, rtimestep, rsolver,&
                           rprimalSolution, 1, ilev, ilev)
      call solver_updateStructure(rsolver)

    end do adaptloop

  end subroutine afc_solvePseudoTransient

  !*****************************************************************************

!<subroutine>

  subroutine afc_solveSteadystate()

!<description>
    ! This subroutine solves the steady-state flow problem
    !
    ! $$\nabla\cdot{\bf f}(u)=0$$
    !
    ! for a scalar quantity $u$ in the domain $\Omega$.
!</description>
!</subroutine>

    ! Pointer to the full multigrid solver
    type(t_solverMultigrid), pointer :: p_rsolverFMG

    ! Pointer to solver
    type(t_solver), pointer :: p_rsolver

    ! Pointer to the current multigrid level
    type(t_problemLevel), pointer :: rproblemLevel

    ! Vector for fine grid solution
    type(t_vectorBlock) :: rsolutionFine
    
    ! Name of graymap file
    character(LEN=SYS_STRLEN) :: pgmfile

    ! Name of bechmark file
    character(LEN=SYS_STRLEN) :: indatfile

    ! Type of finite element
    integer :: ieltype
    
    ! File names for solution error output
    character(LEN=SYS_STRLEN) :: errorfile

    ! Type of UCD output format
    integer :: ierrorUCD

    ! Number of adaptation steps
    integer :: nadapt

    ! local variables
    integer :: ilev,iadapt

    
    ! Retrieve parameters
    call parlst_getvalue_int(rparlist, trim(adjustl(sbenchmarkName)),&
                             "ieltype", ieltype)
    call parlst_getvalue_int(rparlist, trim(adjustl(sbenchmarkName)),&
                             "nadapt", nadapt, 0)    
    call parlst_getvalue_string(rparlist, trim(adjustl(sinputoutputName)),&
                                "errorfile", errorfile, ' ')
    call parlst_getvalue_int(rparlist, trim(adjustl(sinputoutputName)),&
                             "ierrorUCD", ierrorUCD, 1)
    
    ! Adjust time stepping scheme
    rtimestep%ctimestepType = SV_THETA_SCHEME
    rtimestep%dinitialTime  = 0.0_DP
    rtimestep%dfinalTime    = 1.0_DP
    rtimestep%dinitialStep  = 1.0_DP
    rtimestep%dStep         = 1.0_DP
    rtimestep%dStep1        = 1.0_DP
    rtimestep%theta         = 1.0_DP

    ! Set pointer
    p_rsolver => rsolver
    
    ! Walk down the solver structure until applicable solver is reached
    do while (associated(p_rsolver))
      if ((p_rsolver%csolverType .eq. SV_FMG)) exit
      p_rsolver => solver_getNextSolver(p_rsolver)
    end do
    
    if (.not. associated(p_rsolver)) then
      call output_line('Unsupported/invalid solver type!',&
                       OU_CLASS_ERROR, OU_MODE_STD, 'afc_solveSteadystate')
      call sys_halt()
    end if
    
    ! Set pointers
    p_rsolverFMG => p_rsolver%p_solverMultigrid
    p_rsolver    => p_rsolverFMG%p_solverCoarsegrid

    ! Check if there is only one multigrid level
    if (p_rsolverFMG%nlmin .eq. p_rsolverFMG%nlmax) then

      !-------------------------------------------------------------------------
      ! Compute steady-state solution on single grid
      !-------------------------------------------------------------------------

      ! Check if solution matches the required multigrid level
      if (rproblem%p_rproblemLevelMax%ilev .ne. p_rsolverFMG%nlmax) then
        call output_line('Solution does not match the required level!',&
                         OU_CLASS_ERROR, OU_MODE_STD, 'afc_solveSteadystate')
        call sys_halt()
      end if

      ! Start time measurement for solution procedure
      call stat_startTimer(rtimer_solution, STAT_TIMERSHORT)
      
      if (irhstype .eq. 0) then
        call timestep_performThetaStep(rproblem%p_rproblemLevelMax, rtimestep, p_rsolver,&
                                       rprimalSolution, 1,&
                                       fcb_calcResidual, fcb_calcJacobian,&
                                       fcb_applyJacobian, fcb_setBoundary)
      else
        call timestep_performThetaStep(rproblem%p_rproblemLevelMax, rtimestep, p_rsolver,&
                                       rprimalSolution, 1,&
                                       fcb_calcResidual, fcb_calcJacobian,&
                                       fcb_applyJacobian, fcb_setBoundary, rf)
      end if
      
      ! Stop time measurement for solution procedure
      call stat_stopTimer(rtimer_solution)

      ! Set pointer
      rproblemLevel => rproblem%p_rproblemLevelMax

    else

      !-------------------------------------------------------------------------
      ! Compute steady-state solution on sequence of nested grids
      !-------------------------------------------------------------------------

      ! Retrieve local variables
      call parlst_getvalue_string(rparlist, trim(adjustl(sinputoutputName)),&
                                  "indatfile", indatfile)
      call parlst_getvalue_string(rparlist, trim(adjustl(sinputoutputName)),&
                                  "pgmfile", pgmfile, ' ')

      ! Set pointer
      rproblemLevel => rproblem%p_rproblemLevelMin
      
      !-------------------------------------------------------------------------
      ! Phase 1: Choose initial solution on coarsest level.
      !-------------------------------------------------------------------------
      phase1: do
        if (associated(rproblemLevel)) then
        
          ! Are we already on coarse grid level?
          if (rproblemLevel%ilev .eq. p_rsolverFMG%nlmin) then
            exit phase1
          else           
            ! Proceed with next finer level
            rproblemLevel => rproblemLevel%p_rproblemLevelFine
          end if
        else
          call output_line('Could not reach coarse grid level!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'afc_solveSteadystate')
          call sys_halt()
        end if
      end do phase1

      ! Re-initialize the solutionon coarsest grid
      call lsysbl_releaseVector(rprimalSolution)
      if (trim(adjustl(pgmfile)) .eq. '') then
        call afc_initSolutionFromParser(RboundaryCondition(CDEQ_BOUNDARY_PRIMAL),&
                                        rproblemLevel, rprimalSolution, indatfile,&
                                        rtimestep%dinitialTime)
      else
        call afc_initSolutionFromImage(RboundaryCondition(CDEQ_BOUNDARY_PRIMAL),&
                                       rproblemLevel, rprimalSolution, pgmfile)
      end if
      
      ! Re-initialize global right-hand side on coarsest grid
      call lsysbl_releaseVector(rf)
      call afc_initRHS(RboundaryCondition(CDEQ_BOUNDARY_PRIMAL), rproblemLevel, rf, indatfile)


      !-------------------------------------------------------------------------
      ! Phase 2: Solve solution on current grid and perform prolongation
      !          to next finer grid once the solution converged
      !-------------------------------------------------------------------------
      phase2: do
        if (associated(rproblemLevel)) then

          ! Reset solver and timestep to original values
          call solver_resetSolver(p_rsolver, .true.)
          call solver_resetTimestep(rtimestep, .true.)

          ! Adjust the solver hierarchy and update the structure
          call solver_adjustHierarchy(p_rsolver, nlmaxOpt=rproblemLevel%ilev)
          call solver_updateStructure(p_rsolver)
          
          ! We must enforce velocity update so that all internal
          ! structure are regenerated for the current level
          bvelocityUpdate = .true.

          ! Calculate the primal preconditioner for the current level
          call afc_calcPrecond(rproblemLevel, rtimestep, p_rsolver,&
                               rprimalSolution, 1, rproblemLevel%ilev, rproblemLevel%ilev)
          call solver_updateStructure(p_rsolver)

          ! Perform one step of the fully implicit theta-scheme
          if (irhstype .eq. 0) then
            call timestep_performThetaStep(rproblemLevel, rtimestep, p_rsolver,&
                                           rprimalSolution, 1,&
                                           fcb_calcResidual, fcb_calcJacobian,&
                                           fcb_applyJacobian, fcb_setBoundary)

          else
            call timestep_performThetaStep(rproblemLevel, rtimestep, p_rsolver,&
                                           rprimalSolution, 1,&
                                           fcb_calcResidual, fcb_calcJacobian,&
                                           fcb_applyJacobian, fcb_setBoundary, rf)
          end if

          ! Are we already on finest grid level?
          if (rproblemLevel%ilev .eq. p_rsolverFMG%nlmax) then
            exit phase2
          else           
            ! Proceed with next finer level
            rproblemLevel => rproblemLevel%p_rproblemLevelFine
          end if
          
          ! Create new solution vector on the fine grid
          call lsysbl_createVecBlockByDiscr(rproblemLevel%rdiscretisation,&
                                            rsolutionFine, .true., ST_DOUBLE)

          ! Prolongate solution from coarse grid to fine grid
          call solver_prolongationBlock(rproblemLevel%p_rproblemLevelCoarse%rtriangulation,&
                                        rproblemLevel%rtriangulation,&
                                        rprimalSolution, rsolutionFine)

          ! Swap solution vectors: rsolution <-> rsolutionFine
          call lsysbl_swapVectors(rprimalSolution, rsolutionFine)

          ! Release axuiliary solution vector
          call lsysbl_releaseVector(rsolutionFine)

        else
          call output_line('Could not reach fine grid level!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'afc_solveSteadystate')
          call sys_halt()
        end if
      end do phase2
    end if

    
    ! Adaptation loop
    adaptloop: do iadapt = 1, nadapt

      
      !-------------------------------------------------------------------------
      ! Postprocess intermediate solution
      !-------------------------------------------------------------------------
      
      call afc_postprocess()
    
      
      !-------------------------------------------------------------------------
      ! Adaptive grid refinement
      !-------------------------------------------------------------------------
      
      if (trim(adjustl(sadaptivityName)) .eq. '' .or.&
                                  nadapt .eq. 0 ) exit adaptloop

      ! Perform error estimation
      call afc_prepareErrorEstimator(rproblemLevel, rprimalSolution, rerrorEstimator)
      call afc_performErrorEstimation(rerrorEstimator, rgridIndicator,&
                                      rhadapt%drefinementTolerance)
      call errest_clearErrorEstimator(rerrorEstimator)          
      
      ! Output grid indicator (if required)
      if (trim(errorfile) .ne. '') then
        call afc_outputVectorScalar(rproblemLevel, rgridIndicator,&
                                    rtimestep%dTime, errorfile, ierrorUCD, .false.)
      end if
        
      ! What dimension are we?
      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        ! Perform grid adaptation in 1D
        call afc_performGridAdaptation(rcollection, rhadapt, rproblemLevel,&
                                       rprimalSolution, rgridIndicator, fcb_hadaptCallback1D)
        call lsyssc_releaseVector(rgridIndicator)
        call afc_reinitConstOperators1D(rcollection, rhadapt, rproblemLevel,&
                                        rprimalSolution, ieltype, fcb_hadaptCallback1D)
        
      case (NDIM2D)
        ! Perform grid adaptation in 2D
        call afc_performGridAdaptation(rcollection, rhadapt, rproblemLevel,&
                                       rprimalSolution, rgridIndicator, fcb_hadaptCallback2D)
        call lsyssc_releaseVector(rgridIndicator)
        call afc_reinitConstOperators2D(rcollection, rhadapt, rproblemLevel,&
                                        rprimalSolution, ieltype, fcb_hadaptCallback2D)
        
      case (NDIM3D)
        ! Perform grid adaptation in 3D
        call afc_performGridAdaptation(rcollection, rhadapt, rproblemLevel,&
                                       rprimalSolution, rgridIndicator, fcb_hadaptCallback3D)
        call lsyssc_releaseVector(rgridIndicator)
        call afc_reinitConstOperators3D(rcollection, rhadapt, rproblemLevel,&
                                        rprimalSolution, ieltype, fcb_hadaptCallback3D)
      end select

      ! We must enforce velocity update so that all internal
      ! structure are regenerated for the current level
      bvelocityUpdate = .true.

      ! Update global system operator and solver structure
      call afc_calcPrecond(rproblemLevel, rtimestep, p_rsolver,&
                           rprimalSolution, 1, rproblemLevel%ilev, rproblemLevel%ilev)
      call solver_updateStructure(p_rsolver)           
      
      ! Reset solver and timestep to original values
      call solver_resetSolver(p_rsolver, .false.)
      call solver_resetTimestep(rtimestep, .false.)

      ! Perform one step of the fully implicit theta-scheme
      if (irhstype .eq. 0) then
        call timestep_performThetaStep(rproblemLevel, rtimestep, p_rsolver,&
                                       rprimalSolution, 1,&
                                       fcb_calcResidual, fcb_calcJacobian,&
                                       fcb_applyJacobian, fcb_setBoundary)

      else
        call timestep_performThetaStep(rproblemLevel, rtimestep, p_rsolver,&
                                       rprimalSolution, 1,&
                                       fcb_calcResidual, fcb_calcJacobian,&
                                       fcb_applyJacobian, fcb_setBoundary, rf)
      end if

    end do adaptloop

  end subroutine afc_solveSteadystate
  
  !*****************************************************************************

!<subroutine>

  subroutine afc_solvePrimalDual()

!<description>
    ! This subroutine solves the primal and dual problem
    ! for the steady-state flow problem
    !
    ! $$\nabla\cdot{\bf f}(u)=0$$
    !
    ! for a scalar quantity $u$ in the domain $\Omega$.
!</description>
!</subroutine>

    ! Pointer to the full multigrid solver
    type(t_solverMultigrid), pointer :: p_rsolverFMG

    ! Pointer to solver
    type(t_solver), pointer :: p_rsolver
    
    ! Vector for discretized target functional
    type(t_vectorBlock) :: rtargetFunc

    ! Scalar vector for discretized target functional
    type(t_vectorScalar) :: rvector
    type(t_vectorBlock) :: rvectorBlock
    
    ! Matrix for transposed coefficient matrices
    type(t_matrixScalar) :: rmatrix

    ! Linear form to evaluate the target functional
    type(t_linearForm) :: rform

    ! Name of bechmark file
    character(LEN=SYS_STRLEN) :: indatfile

    ! Type of finite element
    integer :: ieltype
    
    ! File names for solution error output
    character(LEN=SYS_STRLEN) :: errorfile

    ! Type of UCD output format
    integer :: ierrorUCD

    ! Number of adaptation steps
    integer :: nadapt

    ! local variables
    real(DP), dimension(:,:), pointer :: p_Dcoords
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer :: ilev,iadapt,ive,nve,iel,ivt
    real(DP), dimension(2,3) :: Dcoords3
    real(DP), dimension(2,4) :: Dcoords4

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata1, p_Ddata2


    ! Retrieve parameters
    call parlst_getvalue_int(rparlist, trim(adjustl(sbenchmarkName)),&
                             "ieltype", ieltype)
    call parlst_getvalue_int(rparlist, trim(adjustl(sbenchmarkName)),&
                             "nadapt", nadapt, 0)    
    call parlst_getvalue_string(rparlist, trim(adjustl(sinputoutputName)),&
                                "errorfile", errorfile, ' ')
    call parlst_getvalue_int(rparlist, trim(adjustl(sinputoutputName)),&
                             "ierrorUCD", ierrorUCD, 1)

    ! Adjust time stepping scheme
    rtimestep%ctimestepType = SV_THETA_SCHEME
    rtimestep%dinitialTime  = 0.0_DP
    rtimestep%dfinalTime    = 1.0_DP
    rtimestep%dinitialStep  = 1.0_DP
    rtimestep%dStep         = 1.0_DP
    rtimestep%dStep1        = 1.0_DP
    rtimestep%theta         = 1.0_DP

    ! Set pointer
    p_rsolver => rsolver
    
    ! Walk down the solver structure until applicable solver is reached
    do while (associated(p_rsolver))
      if ((p_rsolver%csolverType .eq. SV_FMG)) exit
      p_rsolver => solver_getNextSolver(p_rsolver)
    end do
    
    if (.not. associated(p_rsolver)) then
      call output_line('Unsupported/invalid solver type!',&
                       OU_CLASS_ERROR, OU_MODE_STD, 'afc_solvePrimalDual')
      call sys_halt()
    end if
    
    ! Set pointers
    p_rsolverFMG => p_rsolver%p_solverMultigrid
    p_rsolver    => p_rsolverFMG%p_solverCoarsegrid

    ! Check if there is only one multigrid level
    if (p_rsolverFMG%nlmin .eq. p_rsolverFMG%nlmax) then

      ! Adaptation loop
      adaptloop: do iadapt = 1, nadapt
        
        !-------------------------------------------------------------------------
        ! Compute primal steady-state solution on single grid
        !-------------------------------------------------------------------------
        
        ! Check if solution matches the required multigrid level
        if (rproblem%p_rproblemLevelMax%ilev .ne. p_rsolverFMG%nlmax) then
          call output_line('Solution does not match the required level!',&
                           OU_CLASS_ERROR, OU_MODE_STD, 'afc_solvePrimalDual')
          call sys_halt()
        end if
        
        ! Start time measurement for solution procedure
        call stat_startTimer(rtimer_solution, STAT_TIMERSHORT)
        
        call timestep_performThetaStep(rproblem%p_rproblemLevelMax, rtimestep, p_rsolver,&
                                       rprimalSolution, 1,&
                                       fcb_calcResidual, fcb_calcJacobian,&
                                       fcb_applyJacobian, fcb_setBoundary)
        
        ! Stop time measurement for solution procedure
        call stat_stopTimer(rtimer_solution)
        
        !-------------------------------------------------------------------------
        ! Compute dual steady-state solution on single grid
        !-------------------------------------------------------------------------
        
        ! Create dual solution vector
        call lsysbl_releaseVector(rdualSolution)
        call lsysbl_duplicateVector(rprimalSolution, rdualSolution,&
                                    LSYSSC_DUP_TEMPLATE, LSYSSC_DUP_EMPTY)
        call lsysbl_clearVector(rdualSolution)

        ! Impose boundary conditions for the dual problem
        call solver_setBoundaryCondition(rsolver,&
                                         RboundaryCondition(CDEQ_BOUNDARY_DUAL), .true.)
        call solver_updateStructure(rsolver)

        ! We must enforce velocity update so that all internal
        ! structure are regenerated for the current level
        bvelocityUpdate = .true.
        
        ! Update global system operator and solver structure
        call afc_calcPrecond(rproblem%p_rproblemLevelMax, rtimestep, rsolver,&
                             rdualSolution, 2, rproblem%p_rproblemLevelMax%ilev,&
                             rproblem%p_rproblemLevelMax%ilev)
        call solver_updateStructure(rsolver)

        ! Reset solver and timestep to original values
        call solver_resetSolver(p_rsolver, .false.)
        call solver_resetTimestep(rtimestep, .false.)

        ! Set up the corresponding linear form (1,Phi_j):
        rform%itermCount      = 1
        rform%Idescriptors(1) = DER_FUNC
        
        ! Build the discretized target functional
        call linf_buildVectorScalar(rproblem%p_rproblemLevelMax%rdiscretisation%RspatialDiscr(1),&
                                    rform, .true., rvector, fcb_coeffTargetFunc)

        ! Convert the temporal scalar vector to a 1-block vector
        call lsysbl_convertVecFromScalar(rvector, rtargetFunc,&
                                         rproblem%p_rproblemLevelMax%rdiscretisation)
        
        ! Release temporal scalar vector
        call lsyssc_releaseVector(rvector)
        
        
        ! Start time measurement for solution procedure
        call stat_startTimer(rtimer_solution, STAT_TIMERSHORT)
        
        call timestep_performThetaStep(rproblem%p_rproblemLevelMax, rtimestep, p_rsolver,&
                                       rdualSolution, 2, &
                                       fcb_calcResidual, fcb_calcJacobian,&
                                       fcb_applyJacobian, fcb_setBoundary, rtargetFunc)
        
        ! Stop time measurement for solution procedure
        call stat_stopTimer(rtimer_solution)
        
        ! Release the discretized target functional
        call lsysbl_releaseVector(rtargetFunc)
        
        !-------------------------------------------------------------------------
        ! Compute Galerkin orthogonality error
        !-------------------------------------------------------------------------
        
        ! Compute Galerkin residual
        call afc_calcGalerkinResidual(rproblem%p_rproblemLevelMax, rprimalSolution, rvectorBlock, 1)
        
        call lsysbl_getbase_double(rvectorBlock, p_Ddata1)
        call lsysbl_getbase_double(rdualSolution, p_Ddata2)
        
        p_Ddata2 = abs(p_Ddata1*p_Ddata2)
        
        ! That's it
        call lsysbl_releaseVector(rvectorBlock)

        !-------------------------------------------------------------------------
        ! Postprocess intermediate solution
        !-------------------------------------------------------------------------
        
        call afc_postprocess()
        
        
        !-------------------------------------------------------------------------
        ! Adaptive grid refinement
        !-------------------------------------------------------------------------
        
        if (trim(adjustl(sadaptivityName)) .eq. '' .or.&
                                    nadapt .eq. 0 ) exit adaptloop

        call lsyssc_createVector(rgridIndicator, rproblem%p_rproblemLevelMax%rtriangulation%NEL, .true.)
        call lsyssc_getbase_double(rproblem%p_rproblemLevelMax%Rmatrix(CDEQ_MATRIX_ML), p_Ddata1)
        p_Ddata2 = p_Ddata2/p_Ddata1

        call lsyssc_getbase_double(rgridIndicator, p_Ddata1)

        call storage_getbase_int2D(rproblem%p_rproblemLevelMax%rtriangulation%h_IverticesAtElement,&
                                   p_IverticesAtElement)
        call storage_getbase_double2D(rproblem%p_rproblemLevelMax%rtriangulation%h_DvertexCoords,&
                                      p_Dcoords)

        ! Compute element contributions
        do iel = 1, rproblem%p_rproblemLevelMax%rtriangulation%NEL
          
          nve = tria_getNVE(p_IverticesAtElement, iel)
          do ive = 1, nve
            ivt = p_IverticesAtElement(ive, iel)
            p_Ddata1(iel) = p_Ddata1(iel) + p_Ddata2(ivt)
          end do

          p_Ddata1(iel) = p_Ddata1(iel)/nve
          
          if (nve .eq. 3) then
            Dcoords3 = p_Dcoords(:,p_IverticesAtElement(1:3,iel))
            p_Ddata1(iel) = p_Ddata1(iel) * gaux_getArea_tria2D(Dcoords3)
                                            
          else
            Dcoords4 = p_Dcoords(:,p_IverticesAtElement(:,iel))
            p_Ddata1(iel) = p_Ddata1(iel) * gaux_getArea_quad2D(Dcoords4)
                end if
        end do

        ! Output grid indicator (if required)
        if (trim(errorfile) .ne. '') then
          call afc_outputVectorScalar(rproblem%p_rproblemLevelMax, rgridIndicator,&
                                      rtimestep%dTime, errorfile, ierrorUCD, .false.)
        end if
        
        ! What dimension are we?
        select case(rproblem%p_rproblemLevelMax%rtriangulation%ndim)
        case (NDIM1D)
          ! Perform grid adaptation in 1D
          call afc_performGridAdaptation(rcollection, rhadapt, rproblem%p_rproblemLevelMax,&
                                         rprimalSolution, rgridIndicator, fcb_hadaptCallback1D)
          call lsyssc_releaseVector(rgridIndicator)
          call afc_reinitConstOperators1D(rcollection, rhadapt, rproblem%p_rproblemLevelMax,&
                                          rprimalSolution, ieltype, fcb_hadaptCallback1D)
          
        case (NDIM2D)
          ! Perform grid adaptation in 2D
          call afc_performGridAdaptation(rcollection, rhadapt, rproblem%p_rproblemLevelMax,&
                                         rprimalSolution, rgridIndicator, fcb_hadaptCallback2D)
          call lsyssc_releaseVector(rgridIndicator)
          call afc_reinitConstOperators2D(rcollection, rhadapt, rproblem%p_rproblemLevelMax,&
                                          rprimalSolution, ieltype, fcb_hadaptCallback2D)
          
        case (NDIM3D)
          ! Perform grid adaptation in 3D
          call afc_performGridAdaptation(rcollection, rhadapt, rproblem%p_rproblemLevelMax,&
                                         rprimalSolution, rgridIndicator, fcb_hadaptCallback3D)
          call lsyssc_releaseVector(rgridIndicator)
          call afc_reinitConstOperators3D(rcollection, rhadapt, rproblem%p_rproblemLevelMax,&
                                          rprimalSolution, ieltype, fcb_hadaptCallback3D)
        end select
        
        ! Impose boundary conditions for the primal problem
        call solver_setBoundaryCondition(rsolver,&
                                         RboundaryCondition(CDEQ_BOUNDARY_PRIMAL), .true.)
        call solver_updateStructure(rsolver)

        ! We must enforce velocity update so that all internal
        ! structure are regenerated for the current level
        bvelocityUpdate = .true.

        ! Update global system operator and solver structure
        call afc_calcPrecond(rproblem%p_rproblemLevelMax, rtimestep, p_rsolver,&
                             rprimalSolution, 1, rproblem%p_rproblemLevelMax%ilev,&
                             rproblem%p_rproblemLevelMax%ilev)
        call solver_updateStructure(p_rsolver)

        ! Reset solver and timestep to original values
        call solver_resetSolver(p_rsolver, .false.)
        call solver_resetTimestep(rtimestep, .false.)
        
      end do adaptloop

    else
      
      call output_line('Cannot compute steady-state solution on multigrids!',&
                       OU_CLASS_ERROR, OU_MODE_STD, 'afc_solvePrimalDual')
      call sys_halt()
     
    end if

  end subroutine afc_solvePrimalDual

  !*****************************************************************************

!<subroutine>

  subroutine afc_postprocess()

!<description>
    ! This subroutine performs all post-processing tasks
!</description>
!</subroutine>

    ! Time step for solution output
    real(DP) :: tstepgmv, tstepavs, tstepvtk

    ! File names for solution output
    character(LEN=SYS_STRLEN) :: gmvfile, avsfile, vtkfile

    ! local variables
    logical :: bcompatible
    
    
    ! Retrieve parameters
    call parlst_getvalue_string(rparlist, trim(adjustl(sinputoutputName)),&
                                "gmvfile", gmvfile, ' ')
    call parlst_getvalue_double(rparlist, trim(adjustl(sinputoutputName)),&
                                "tstepgmv", tstepgmv, 0.0_DP)
    call parlst_getvalue_string(rparlist, trim(adjustl(sinputoutputName)),&
                                "avsfile", avsfile, ' ')
    call parlst_getvalue_double(rparlist, trim(adjustl(sinputoutputName)),&
                                "tstepavs", tstepavs, 0.0_DP)
    call parlst_getvalue_string(rparlist, trim(adjustl(sinputoutputName)),&
                                "vtkfile", vtkfile, ' ')
    call parlst_getvalue_double(rparlist, trim(adjustl(sinputoutputName)),&
                                "tstepvtk", tstepvtk, 0.0_DP)
    
    ! Perform solution output in GMV format?
    if (tstepgmv .ge. 0._DP) then
      call lsysbl_isVectorCompatible(rprimalSolution,&
                                     rdualSolution, bcompatible)
      if (bcompatible) then
        call afc_outputSolution(rproblem%p_rproblemLevelMax, rprimalSolution,&
                                rtimestep%dTime, gmvfile, 1, rdualSolution)
      else
        call afc_outputSolution(rproblem%p_rproblemLevelMax, rprimalSolution,&
                                rtimestep%dTime, gmvfile, 1)
      end if
    end if
    
    ! Perform solution output in AVS format?
    if (tstepavs .ge. 0._DP) then
      call lsysbl_isVectorCompatible(rprimalSolution,&
                                     rdualSolution, bcompatible)
      if (bcompatible) then
        call afc_outputSolution(rproblem%p_rproblemLevelMax, rprimalSolution,&
                                rtimestep%dTime, avsfile, 2, rdualSolution)
      else
        call afc_outputSolution(rproblem%p_rproblemLevelMax, rprimalSolution,&
                                rtimestep%dTime, avsfile, 2)
      end if
    end if
    
    ! Perform solution output in VTK format?
    if (tstepvtk .ge. 0._DP) then
      call lsysbl_isVectorCompatible(rprimalSolution,&
                                     rdualSolution, bcompatible)
      if (bcompatible) then
        call afc_outputSolution(rproblem%p_rproblemLevelMax, rprimalSolution,&
                                rtimestep%dTime, vtkfile, 3, rdualSolution)
      else
        call afc_outputSolution(rproblem%p_rproblemLevelMax, rprimalSolution,&
                                rtimestep%dTime, vtkfile, 3)
      end if
    end if

  end subroutine afc_postprocess

  !*****************************************************************************

!<subroutine>

  subroutine afc_statistics()

!<description>
    ! This subroutine provides all statistics of the simulation
!</description>
!</subroutine>
    
    ! local variables
    real(DP) :: dtotalTime,dfraction
    character(LEN=SYS_STRLEN) :: indatfile

    call output_lbrk(nlbrk=5)
    call output_separator (OU_SEP_STAR)
    write(*,FMT='(24X,A)') '*** S T A T I S T I C S ***'
    call output_separator (OU_SEP_STAR)
    call output_lbrk(nlbrk=3)
    
    if (rtimestep%dTime > 0) then
      call solver_infoTimestep(rtimestep)
      call output_lbrk()

      call solver_infoSolver(rsolver)
      call output_lbrk()

      call hadapt_infoStatistics(rhadapt)
      call output_lbrk()
    end if

    ! Statistics about the solution error (if available)
    call parlst_getvalue_string(rparlist, trim(adjustl(sinputoutputName)),&
                                "indatfile", indatfile)
    call afc_calcSolutionError(rproblem%p_rproblemLevelMax, rprimalSolution, indatfile,&
                               rtimestep%dTime, rcollection)

    ! Statistics about the time measurement
    call output_lbrk(nlbrk=5)
    call output_separator (OU_SEP_STAR)
    write(*,FMT='(24X,A)') '*** TIME MEASUREMENT ***'
    call output_separator (OU_SEP_STAR)

    call stat_subTimers(rtimer_assembly_matrix, rtimer_solution)
    call stat_subTimers(rtimer_assembly_resrhs, rtimer_solution)

    dtotalTime = max(rtimer_total%delapsedCPU, rtimer_total%delapsedReal)
    dfraction  = 100.0_DP/dtotalTime

    call output_line('  Time for computing solution   : '//&
                     trim(sys_sdL(rtimer_solution%delapsedCPU, 2))//'  '//&
                     trim(sys_sdL(dfraction*rtimer_solution%delapsedCPU, 2))//' %')
    call output_line('  Time for mesh adaptivity      : '//&
                     trim(sys_sdL(rtimer_adaptivity%delapsedCPU, 2))//'  '//&
                     trim(sys_sdL(dfraction*rtimer_adaptivity%delapsedCPU, 2))//' %')
    call output_line('  Time for error estimation     : '//&
                     trim(sys_sdL(rtimer_errorestimation%delapsedCPU, 2))//'  '//&
                     trim(sys_sdL(dfraction*rtimer_errorestimation%delapsedCPU, 2))//' %')
    call output_line('  Time for triangulation        : '//&
                     trim(sys_sdL(rtimer_triangulation%delapsedCPU, 2))//'  '//&
                     trim(sys_sdL(dfraction*rtimer_triangulation%delapsedCPU, 2))//' %')
    call output_line('  Time for coefficient assembly : '//&
                     trim(sys_sdL(rtimer_assembly_coeff%delapsedCPU, 2))//'  '//&
                     trim(sys_sdL(dfraction*rtimer_assembly_coeff%delapsedCPU, 2))//' %')
    call output_line('  Time for matrix assembly      : '//&
                     trim(sys_sdL(rtimer_assembly_matrix%delapsedCPU, 2))//'  '//&
                     trim(sys_sdL(dfraction*rtimer_assembly_matrix%delapsedCPU, 2))//' %')
    call output_line('  Time for residual/rhs assembly: '//&
                     trim(sys_sdL(rtimer_assembly_resrhs%delapsedCPU, 2))//'  '//&
                     trim(sys_sdL(dfraction*rtimer_assembly_resrhs%delapsedCPU, 2))//' %')
    call output_line('  Time for pre-/post-processing : '//&
                     trim(sys_sdL(rtimer_prepostprocess%delapsedCPU, 2))//'  '//&
                     trim(sys_sdL(dfraction*rtimer_prepostprocess%delapsedCPU, 2))//' %')
    call output_separator (OU_SEP_MINUS)
    call output_line('  Time for total simulation     : '//&
                     trim(sys_sdL(dtotalTime, 2)))
    call output_lbrk()
  end subroutine afc_statistics

  !*****************************************************************************

!<subroutine>

  subroutine afc_finalize()

!<description>
    ! This subroutine performs all finalization tasks
!</description>
!</subroutine>

    ! Release collection
    call collct_done(rcollection)
    
    ! Release parameter list
    call parlst_done(rparlist)
    
    ! Release solvers
    call solver_releaseTimestep(rtimestep)
    call solver_releaseSolver(rsolver)
    
    ! Release problem structure
    call problem_releaseProblem(rproblem)
    
    ! Release boundary conditions
    call bdrf_release(RboundaryCondition(CDEQ_BOUNDARY_PRIMAL))
    call bdrf_release(RboundaryCondition(CDEQ_BOUNDARY_DUAL))
    
    ! Release vectors
    call lsysbl_releaseVector(rprimalSolution)
    call lsysbl_releaseVector(rdualSolution)
    call lsysbl_releaseVector(rf)
    
    ! Release grid adaptation structure
    call hadapt_releaseAdaptation(rhadapt)
    
    ! Release error estomator
    call errest_releaseErrorEstimator(rerrorEstimator)
    
    ! Release function parser for velocity
    call fparser_release(rvelocityParser)
    
    ! Release function parser for right-hand side
    call fparser_release(rrhsParser)
    
    ! Release function parser
    call fparser_done()
    
    ! Release data of callback routine
    call fcb_hadaptCallback1D(rcollection, HADAPT_OPR_DONECALLBACK,&
                              (/0/), (/0/))
    call fcb_hadaptCallback2D(rcollection, HADAPT_OPR_DONECALLBACK,&
                              (/0/), (/0/))
    call fcb_hadaptCallback3D(rcollection, HADAPT_OPR_DONECALLBACK,&
                              (/0/), (/0/))
    
    ! Release external storage
    call exstor_info(.true.)
    call exstor_done()
    call output_lbrk()

    ! Release storage
    call storage_info(.true.)
    call storage_done()
    call output_lbrk()

  end subroutine afc_finalize

  !*****************************************************************************

!<subroutine>

  subroutine afc_getParameters(rparlist)

!<description>
    ! This subroutine sets all global parameters from the parameter list
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(IN) :: rparlist
!</input>
!</subroutine>

    ! local variables
    integer :: iarg,narg
    character(LEN=SYS_STRLEN) :: cbuffer

    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist, '', "benchmark",   sbenchmarkName)
    call parlst_getvalue_string(rparlist, '', "inputoutput", sinputoutputName)
    call parlst_getvalue_string(rparlist, '', "timestep",    stimestepName)
    call parlst_getvalue_string(rparlist, '', "solver",      ssolverName)
    call parlst_getvalue_string(rparlist, '', "adaptivity",  sadaptivityName, ' ')
    call parlst_getvalue_string(rparlist, '', "errorest",    serrorestName ,' ')

    ! Overwrite global configuration from command line arguments
    iarg = 1; narg = command_argument_count()
    
    cmdarg: do
      ! Retrieve next command line argument
      call get_command_argument(iarg,cbuffer)
      select case(trim(adjustl(cbuffer)))
        
      case ('-A','--adaptivity')
        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        sadaptivityName = trim(adjustl(cbuffer))

      case ('-B','--benchmark')
        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        sbenchmarkName = trim(adjustl(cbuffer))
        
      case ('-E','--errorest')
        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        serrorestName = trim(adjustl(cbuffer))
        
      case ('-I','--io')
        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        sinputoutputName = trim(adjustl(cbuffer))
        
      case ('-S','--solver')
        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        ssolverName = trim(adjustl(cbuffer))
        
      case ('-T','--timestep')
        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        stimestepName = trim(adjustl(cbuffer))
        
      case DEFAULT
        iarg = iarg+1
        if (iarg .ge. narg) exit cmdarg
      end select
    end do cmdarg


    ! Get benchmark configuration
    call parlst_getvalue_int(rparlist, trim(adjustl(sbenchmarkName)),&
                             "idiffusiontype", idiffusiontype, 0)
    call parlst_getvalue_int(rparlist, trim(adjustl(sbenchmarkName)),&
                             "ivelocitytype", ivelocitytype)
    call parlst_getvalue_int(rparlist, trim(adjustl(sbenchmarkName)),&
                             "iflowtype", iflowtype)
    call parlst_getvalue_int(rparlist, trim(adjustl(sbenchmarkName)),&
                             "irhstype", irhstype, 0)

  end subroutine afc_getParameters

  !*****************************************************************************

!<subroutine>

  subroutine afc_UserInterface(imode)

!<description>
    ! This subroutine enables the user to interact with the
    ! simulation. Depending on the "mode" different tasks are
    ! performed.
!</description>

!<input>
    ! interaction mode
    integer, intent(IN) :: imode
!</input>
!</subroutine>

    include '.version'

    ! local variables
    integer, external :: signal_SIGINT
    integer :: narg
    character(LEN=SYS_STRLEN) :: cbuffer
    
    select case(imode)
    case (0)
      ! Welcome screen
      call output_lbrk()
      call output_separator(OU_SEP_STAR)
      write(*,FMT='(2X,A,T10,A,T30,A,T50,A,T72,A)') '***','AFC','Version '//VERSION,'Build '//BUILD,'***'
      call output_separator(OU_SEP_STAR)
      call output_line('  Algebraic Flux Correction for scalar conservation laws')
      call output_lbrk()
      call output_line('  Authors:  D. Kuzmin, M. Moeller (2004-2008)')
      call output_line('            Applied Mathematics, Dortmund University of Technology')
      call output_line('            Vogelpothsweg 87, 44227 Dortmund, Germany')
      call output_separator(OU_SEP_STAR)
      call getenv('HOST',cbuffer)
      call output_line('  Hostname:        '//trim(adjustl(cbuffer)))
      call getenv('USER',cbuffer)
      call output_line('  User:            '//trim(adjustl(cbuffer)))
      call getenv('HOSTTYPE',cbuffer)
      call output_line('  Hosttype:        '//trim(adjustl(cbuffer)))
      
      ! Get number of command line arguments
      narg = command_argument_count()
      if (narg < 1) then
        call output_separator(OU_SEP_STAR)
        call output_lbrk()
        call output_line('  Usage: afc [OPTIONS] FILE')
        call output_lbrk()
        call output_line('    FILE            name of parameter file')
        call output_lbrk()
        call output_line('    -A,--adaptivity section name for adaptivity')
        call output_line('    -B,--benchmark  section name for benchmark')
        call output_line('    -E,--errorest   section name for error estimator')
        call output_line('    -I,--io         section name for input/output')
        call output_line('    -S,--solver     section name for top-level solver')
        call output_line('    -T,--timestep   section name for time-stepping scheme')
        call output_lbrk()
        call sys_halt()
      end if
      
      ! Get parameter file from command line arguments (last argumant)
      call get_command_argument(narg,cbuffer)
      parmfile = trim(adjustl(cbuffer))
      call output_line('  Parameterfile:   '//trim(adjustl(cbuffer)))
      call output_separator(OU_SEP_STAR)
      call output_lbrk()

      
    case (1)
      ! Perform intermediate output
      if (signal_SIGINT(-1) > 0 ) call afc_postprocess()
      
    end select
  end subroutine afc_UserInterface
end program afc

!*****************************************************************************

!<function>

function signal_SIGINT(signum) result(sigcount)

  use fsystem
  use genoutput
  use signal

!<description>
  ! This subroutine performs signal handling for SIGINT. In essence,
  ! it counts the number if SIGINT's received and terminates if user
  ! sent SIGINT more than three times.
!</description>

!<input>
  integer, intent(IN) :: signum
!</input>

!<result>
  ! signal
  integer :: sigcount
!</result>
!</function>

  ! local variables
  integer, save :: icount = 0
  
  sigcount = icount
  
  if (signum .eq. -1) then
    
    ! Reset counter
    icount = 0

  elseif (signum .eq. SIGINT) then
    
    ! Increase counter
    icount = icount+1
    if (icount .ge. 3) then
      call output_line('Simulation terminated due to user interruption (SIGINT)')
      call sys_halt()
    end if

  end if
end function signal_SIGINT

!*****************************************************************************

!<function>

function signal_SIGQUIT(signum) result(sigcount)

  use fsystem
  use genoutput
  use signal

!<description>
  ! This subroutine performs signal handling for SIGQUIT.
!</description>

!<input>
  integer, intent(IN) :: signum
!</input>

!<result>
  ! signal
  integer :: sigcount
!</result>
!</function>

  call output_line('Simulation terminated due to user interruption (SIGQUIT)')
  stop
end function signal_SIGQUIT
