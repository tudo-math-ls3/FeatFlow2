!##############################################################################
!# ****************************************************************************
!# <name> zpinch_application </name>
!# ****************************************************************************
!#
!# <purpose>
!# This application solves the time-dependent magnetohydrodynamic equations
!# in the one-, two- or three-dimensional domain $\Omega$.
!#
!# The spatial discretization is perform by means of the algebraic
!# flux correction (AFC) paradigm by Kuzmin, Moeller and Turek. In
!# particular, high-resolution finite element schemes of TVD- and
!# FCT-type are available. For the temporal discretization, the 
!# two-level theta-scheme is employed, whereby $\theta\in(0,1]$.
!#
!# Dynamic mesh adaptation is based on the red-green strategy, whereby
!# mixed triangulations are supported. Error estimation is based on the
!# scalar tracer quantity required for evaluating the source term.
!#
!#
!# The following routines are available:
!#
!# 1.) zpinch_app
!#     -> The main routine of the application called from the main
!#        program. The routine gets all required information from the
!#        parameter list which needs to be initialized and filled in
!#        the main program. It then works black-box, that is, it
!#        determines the solution algorithm to be used and performs
!#        the simulation. The user should only have to modify this
!#        routine if another solution algorithm is implemented.
!#
!# 2.) zpinch_initSolvers
!#     -> Initializes the solve structures from the parameter list.
!#
!# 3.) zpinch_initProblem
!#     -> Initializes the global problem structure based on the
!#        parameter settings given by the parameter list. This routine
!#        is quite universal, that is, it prepares the internal
!#        structure of the global problem and generates a linked
!#        list of problem levels used in the multigrid hierarchy.
!#
!# 4.) zpinch_applySourceTerm
!#     -> Initializes the source term and applies it to the solution
!#
!# 5.) zpinch_calcAdaptationIndicator
!#     -> Calculates the element-wise indicator for refinement/-coarsening
!#        based on the detectation of shocks and contact discontinuities
!#
!# 6.) zpinch_outputSolution
!#     -> Outputs the solution vector to file in UCD format
!#
!# 7.) zpinch_adaptTriangulation
!#      -> Performs h-adaptation for the given triangulation
!#
!# 8.) zpinch_solveTransientPrimal
!#     -> Solves the primal formulation of the time-dependent 
!#        simplified MHD equations
!#
!#
!# The following auxiliary routines are available:
!#
!# 1.) zpinch_parseCmdlArguments
!#     -> Parses the list of commandline arguments and overwrites
!#        parameter values from the parameter files
!#
!# </purpose>
!##############################################################################

module zpinch_application

  use afcstabilisation
  use boundaryfilter
  use collection
  use euler_application
  use euler_basic
  use euler_callback
  use euler_callback1d
  use euler_callback2d
  use euler_callback3d
  use flagship_basic
  use fparser
  use genoutput
  use graph
  use hadaptaux
  use hadaptivity
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
  use transport_application
  use transport_basic
  use transport_callback
  use transport_callback1d
  use transport_callback2d
  use transport_callback3d
  use ucd
  use zpinch_callback
  use zpinch_callback2d

  implicit none

  private
  public :: zpinch_app

contains

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_app(rparlist)

!<description>
    ! This is the main application for the simplified MHD equations. It
    ! is a so-called driver routine which can be used to start a
    ! standalone MHD simulation.
!</description>

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
    type(t_boundaryCondition) :: rbdrCondEuler, rbdrCondTransport

    ! Problem structure which holds all internal data (vectors/matrices)
    type(t_problem) :: rproblem

    ! Time-stepping structures
    type(t_timestep) :: rtimestep
    
    ! Global solver structure and points
    type(t_solver) :: rsolver
    type(t_solver), pointer :: p_rsolver => null()

    ! Solution vectors and pointers
    type(t_vectorBlock), dimension(2), target :: rsolution
    type(t_vectorBlock), pointer :: p_rsolutionEuler, p_rsolutionTransport

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
    character(LEN=SYS_STRLEN) :: ssectionNameEuler
    character(LEN=SYS_STRLEN) :: ssectionNameTransport

    ! local variables
    integer :: isystemFormat, systemMatrix
    integer :: nlmin, nlmax, ndimension

    ! Start total time measurement
    call stat_startTimer(rtimerTotal)
    
    !---------------------------------------------------------------------------
    ! Pre-processing
    !---------------------------------------------------------------------------
    
    ! Start time measurement
    call stat_startTimer(rtimerPrePostprocess, STAT_TIMERSHORT)

    p_rsolutionEuler => rsolution(1)
    p_rsolutionTransport => rsolution(2)

    ! Retrieve section names of sub-applications
    call parlst_getvalue_string(rparlist,&
        'zpinch', 'application_euler', ssectionNameEuler)
    call parlst_getvalue_string(rparlist,&
        'zpinch', 'application_transport', ssectionNameTransport)

    ! Overwrite configuration from command line arguments. After this
    ! subroutine has been called, the parameter list remains unchanged
    ! unless the used updates some parameter values interactively.
    call zpinch_parseCmdlArguments(rparlist)
    call euler_adjustParameterlist(rparlist, ssectionNameEuler)
    call transp_adjustParameterlist(rparlist, ssectionNameTransport)

    
    ! Initialize global collection structure
    call collct_init(rcollection)

    ! Attach the parameter list and the timers 
    ! to the collection for the Euler model
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
    ! and functions from the parameter files
    call parlst_getvalue_string(rparlist,&
        'zpinch', 'indatfile', sindatfileName)
    call fparser_parseFileForKeyword(rfparser,&
        sindatfileName, 'defconst', FPAR_CONSTANT)
    call fparser_parseFileForKeyword(rfparser,&
        sindatfileName, 'defexpr', FPAR_EXPRESSION)
    call fparser_parseFileForKeyword(rfparser,&
        sindatfileName, 'deffunc', FPAR_FUNCTION)

    call parlst_getvalue_string(rparlist,&
        ssectionNameEuler, 'indatfile', sindatfileName)
    call fparser_parseFileForKeyword(rfparser,&
        sindatfileName, 'defconst', FPAR_CONSTANT)
    call fparser_parseFileForKeyword(rfparser,&
        sindatfileName, 'defexpr', FPAR_EXPRESSION)
    call fparser_parseFileForKeyword(rfparser,&
        sindatfileName, 'deffunc', FPAR_FUNCTION)

    call parlst_getvalue_string(rparlist,&
        ssectionNameTransport, 'indatfile', sindatfileName)
    call fparser_parseFileForKeyword(rfparser,&
        sindatfileName, 'defconst', FPAR_CONSTANT)
    call fparser_parseFileForKeyword(rfparser,&
        sindatfileName, 'defexpr', FPAR_EXPRESSION)
    call fparser_parseFileForKeyword(rfparser,&
        sindatfileName, 'deffunc', FPAR_FUNCTION)

    ! Attach the function parser to the collection
    call collct_setvalue_pars(rcollection,&
        'rfparser', rfparser, .true.)
    
    ! Initialize the solver structures
    call zpinch_initSolvers(rparlist, 'zpinch', rtimestep, rsolver)

    ! Initialize the abstract problem structure
    nlmin = solver_getMinimumMultigridlevel(rsolver)
    nlmax = solver_getMaximumMultigridlevel(rsolver)
    
    ! Initialize the abstract problem structure
    call zpinch_initProblem(rparlist, 'zpinch',&
        nlmin, nlmax, rproblem, rcollection)

    ! Initialize the individual problem levels
    call euler_initAllProblemLevels(rparlist,&
        ssectionNameEuler, rproblem, rcollection)
    call transp_initAllProblemLevels(rparlist,&
        ssectionNameTransport, rproblem, rcollection)

    ! Prepare internal data arrays of the solver structure for Euler model
    call parlst_getvalue_int(rparlist,&
        ssectionNameEuler, 'systemMatrix', systemMatrix)
    call parlst_getvalue_int(rparlist,&
        ssectionNameEuler, 'isystemFormat', isystemFormat)

    ! Get first solver subnode
    p_rsolver => solver_getNextSolver(rsolver, 1)
    call flagship_updateSolverMatrix(rproblem%p_rproblemLevelMax,&
        p_rsolver, systemMatrix, isystemFormat, UPDMAT_ALL)

    ! Prepare internal data arrays of the solver structure for transport model
    call parlst_getvalue_int(rparlist,&
        ssectionNameTransport, 'systemMatrix', systemMatrix)

    ! Get second solver subnode
    p_rsolver => solver_getNextSolver(rsolver, 2)
    call flagship_updateSolverMatrix(rproblem%p_rproblemLevelMax,&
        p_rsolver, systemMatrix, SYSTEM_INTERLEAVEFORMAT, UPDMAT_ALL)

    ! Update complete solver structure
    call solver_updateStructure(rsolver)
    
    ! Stop time measurement for pre-processing
    call stat_stopTimer(rtimerPrePostprocess)
    
    
    !---------------------------------------------------------------------------
    ! Solution algorithm
    !---------------------------------------------------------------------------

    if (rtimestep%dfinalTime .gt. rtimestep%dinitialTime) then

      ! Get global configuration from parameter list
      call parlst_getvalue_string(rparlist, 'zpinch',&
          'algorithm', algorithm)
      call parlst_getvalue_int(rparlist, 'zpinch',&
          'ndimension', ndimension)

      ! Initialize the boundary condition for the Euler model
      call parlst_getvalue_string(rparlist, ssectionNameEuler,&
          'sprimalbdrcondname', sbdrcondName)
      call parlst_getvalue_string(rparlist, ssectionNameEuler,&
          'indatfile', sindatfileName)

      ! The boundary condition for the primal problem is required for
      ! all solution strategies so initialize it from the parameter file
      call bdrf_readBoundaryCondition(rbdrCondEuler, sindatfileName,&
          '['//trim(sbdrcondName)//']', ndimension)
      
      ! Initialize the boundary condition for the transport model
      call parlst_getvalue_string(rparlist, ssectionNameTransport,&
          'sprimalbdrcondname', sbdrcondName)
      call parlst_getvalue_string(rparlist, ssectionNameTransport,&
          'indatfile', sindatfileName)

      ! The boundary condition for the primal problem is required for
      ! all solution strategies so initialize it from the parameter file
      call bdrf_readBoundaryCondition(rbdrCondTransport,&
          sindatfileName, '['//trim(sbdrcondName)//']', ndimension)

      ! What solution algorithm should be applied?
      select case(trim(algorithm))
        
      case ('transient_primal')
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Solve the primal formulation for the time-dependent problem
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call zpinch_solveTransientPrimal(rparlist, 'zpinch',&
            ssectionNameEuler, ssectionNameTransport, rbdrCondEuler,&
            rbdrCondTransport, rproblem, rtimestep, rsolver,&
            rsolution, rcollection)
        
        call zpinch_outputSolution(rparlist, 'zpinch',&
            rproblem%p_rproblemLevelMax, rsolution, rtimestep%dTime)
        
      case DEFAULT
        call output_line(trim(algorithm)//' is not a valid solution algorithm!',&
            OU_CLASS_ERROR,OU_MODE_STD,'zpinch_app')
        call sys_halt()
      end select

    else

      ! Just output the computational mesh and exit
      call zpinch_outputSolution(rparlist, 'zpinch',&
          rproblem%p_rproblemLevelMax)

    end if

    !---------------------------------------------------------------------------
    ! Post-processing
    !---------------------------------------------------------------------------
    
    ! Start time measurement for pre-processing
    call stat_startTimer(rtimerPrepostProcess, STAT_TIMERSHORT)
    
    ! Release time-stepping
    call tstep_releaseTimestep(rtimestep)
    
    ! Release solver
    call solver_releaseSolver(rsolver)

    ! Release problem structure
    call problem_releaseProblem(rproblem)

    ! Release boundary conditions
    call bdrf_release(rbdrCondEuler)
    call bdrf_release(rbdrCondTransport)
    
    ! Release vectors
    call lsysbl_releaseVector(p_rsolutionEuler)
    call lsysbl_releaseVector(p_rsolutionTransport)

    ! Release function parser
    call fparser_release(rfparser)

    ! Stop time measurement for pre-processing
    call stat_stopTimer(rtimerPrePostprocess)

    ! Stop time measurement for total time measurement
    call stat_stopTimer(rtimerTotal)

    ! Output statistics
    call euler_outputStatistics(rtimerTotal, rcollection)
    
    ! Release collection
    call collct_done(rcollection)
    
  end subroutine zpinch_app
  
  !*****************************************************************************

!<subroutine>

  subroutine zpinch_initSolvers(rparlist, ssectionName, rtimestep, rsolver)

!<description>
    ! This subroutine initializes the time-stepping structure and
    ! the top-level solver structure from the parameter list
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section names in parameter list
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
        
    
    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'timestep', stimestepName)
    
    ! Initialize time-stepping
    call tstep_createTimestep(rparlist, stimestepName, rtimestep, 2)
    

    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'solver', ssolverName)
    
    ! Initialize solver structures
    call solver_createSolver(rparlist, ssolverName, rsolver)
    call solver_adjustHierarchy(rsolver)
    call solver_updateStructure(rsolver)
    
  end subroutine zpinch_initSolvers

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_initProblem(rparlist, ssectionName,&
      nlmin, nlmax, rproblem, rcollection)

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

!<inputoutput>
    ! collection
    type(t_collection), intent(inout) :: rcollection
!</intputoutput>

!<output>
    ! problem structure
    type(t_problem), intent(out) :: rproblem
!</output>
!</subroutine>

    ! section names
    character(LEN=SYS_STRLEN) :: ssectionNameEuler
    character(LEN=SYS_STRLEN) :: ssectionNameTransport
    character(LEN=SYS_STRLEN) :: sconvection
    character(LEN=SYS_STRLEN) :: sinviscid

    ! abstract problem descriptor
    type(t_problemDescriptor) :: rproblemDescriptor

    ! pointer to the problem level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! local variables
    integer :: convectionAFC
    integer :: inviscidAFC
    integer :: iconvToTria
    

    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'application_euler', ssectionNameEuler)
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'application_transport', ssectionNameTransport)

    call parlst_getvalue_string(rparlist,&
        ssectionName, 'trifile', rproblemDescriptor%trifile)
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'prmfile', rproblemDescriptor%prmfile, '')
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'ndimension', rproblemDescriptor%ndimension)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'iconvtotria', iconvToTria, 0)

    call parlst_getvalue_string(rparlist,&
        ssectionNameTransport, 'convection', sconvection)
    call parlst_getvalue_int(rparlist,&
        ssectionNameTransport, 'convectionAFC', convectionAFC)
    call parlst_getvalue_string(rparlist,&
        ssectionNameEuler, 'inviscid', sinviscid)
    call parlst_getvalue_int(rparlist,&
        ssectionNameEuler, 'inviscidAFC', inviscidAFC)
    
    ! Set additional problem descriptor
    rproblemDescriptor%ndiscretisation = 2
    rproblemDescriptor%nafcstab        = 2   ! for inviscid and convective stabilization
    rproblemDescriptor%nlmin           = nlmin
    rproblemDescriptor%nlmax           = nlmax
    rproblemDescriptor%nmatrixScalar   = rproblemDescriptor%ndimension + 10
    rproblemDescriptor%nmatrixBlock    = 2
    rproblemDescriptor%nvectorScalar   = 0
    rproblemDescriptor%nvectorBlock    = 2   ! external velocity field

    ! Check if quadrilaterals should be converted to triangles
    if (iconvToTria .ne. 0)&
        rproblemDescriptor%iproblemSpec = rproblemDescriptor%iproblemSpec &
                                        + PROBDESC_MSPEC_CONVTRIANGLES   

    ! Initialize problem structure
    call problem_initProblem(rproblemDescriptor, rproblem)
    
    ! Loop over all problem levels
    p_rproblemLevel => rproblem%p_rproblemLevelMax
    do while(associated(p_rproblemLevel))

      if (convectionAFC > 0) then
        call afcstab_initFromParameterlist(rparlist, sconvection,&
            p_rproblemLevel%Rafcstab(convectionAFC))
      end if
      
      if (inviscidAFC > 0) then
        call afcstab_initFromParameterlist(rparlist, sinviscid,&
            p_rproblemLevel%Rafcstab(inviscidAFC))
      end if
      
      ! Switch to next coarser level
      p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
    end do

  end subroutine zpinch_initProblem

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_redistributeMassEuler(rparlist, ssectionName,&
      rproblemLevel, rsolution, rcollection)

!<description>
    ! This subroutine redistributes the mass so that the total mass
    ! is equal to the total mass of Bank and Shadid
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section names in parameter list
    character(LEN=*), intent(in) :: ssectionName

    ! problem level structure
    type(t_problemLevel), intent(in), target :: rproblemLevel
!</input>

!<inputoutput>
    ! solution vector
    type(t_vectorBlock), intent(inout) :: rsolution

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_triangulation), pointer :: p_rtriangulation
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_ML, p_u
    real(DP) :: ddensity1, ddensity2, ddensity3
    real(DP) :: dmass1, dmass2, dmass3, dradius
    integer :: lumpedMassMatrix, isystemFormat
    integer :: ivt

    ! Get global configuration from parameter list
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'lumpedmassmatrix', lumpedMassMatrix)
    call parlst_getvalue_int(rparlist,&
          ssectionName, 'isystemformat', isystemFormat)

    ! Compute total densities
    ddensity1 = SYS_PI
    ddensity2 = SYS_PI * (1.05_DP**2-1._DP) * 1e6_DP
    ddensity3 = SYS_PI * (1.5_DP**2-1.05_DP**2) * 5e-1_DP

    ! Set pointer to triangulation
    p_rtriangulation => rproblemLevel%rtriangulation

    ! Set pointer to coordinate vector
    call storage_getbase_double2D(&
        p_rtriangulation%h_DvertexCoords, p_DvertexCoords)

    ! Set pointer to lumped mass matrix
    call lsyssc_getbase_double(&
        rproblemLevel%Rmatrix(lumpedMassMatrix), p_ML)

    ! Set pointer to conservative variables
    call lsysbl_getbase_double(rsolution, p_u)

    ! Compute total masses
    dmass1 = 0.0_DP
    dmass2 = 0.0_DP
    dmass3 = 0.0_DP
    
    ! Loop over all vertices
    do ivt = 1, p_rtriangulation%NVT
      
      ! Compute distance to origin
      dradius = sqrt(p_DvertexCoords(1,ivt)**2 +&
                     p_DvertexCoords(2,ivt)**2)
      
      if (dradius .lt. 1.0_DP) then
        dmass1 = dmass1 + p_ML(ivt)
      elseif (dradius .gt. 1.05_DP) then
        dmass3 = dmass3 + p_ML(ivt)
      else
        dmass2 = dmass2 + p_ML(ivt)
      end if
    end do

    ddensity1 = ddensity1/dmass1
    ddensity2 = ddensity2/dmass2
    ddensity3 = ddensity3/dmass3

    select case(isystemFormat)
    case(SYSTEM_INTERLEAVEFORMAT)
      call doRedistributeInterleaveFormat(rproblemLevel&
          %Rmatrix(lumpedMassMatrix)%NEQ, NVAR2D, ddensity1,&
          ddensity2, ddensity3, p_DvertexCoords, p_ML, p_u)
      
    case(SYSTEM_BLOCKFORMAT)
      call doRedistributeBlockFormat(rproblemLevel&
          %Rmatrix(lumpedMassMatrix)%NEQ, NVAR2D, ddensity1,&
          ddensity2, ddensity3, p_DvertexCoords, p_ML, p_u)

    case DEFAULT
      call output_line('Invalid system format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'zpinch_redistributeMassEuler')
      call sys_halt()
    end select

  contains

    ! Here, the working routine follow
    
    !**************************************************************
    ! Redistribute mass density stored in interleave format
    
    subroutine doRedistributeInterleaveFormat(neq, nvar, ddensity1,&
        ddensity2, ddensity3, DvertexCoords, ML, u)

      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), dimension(:), intent(in) :: ML
      real(DP), intent(in) :: ddensity1,ddensity2,ddensity3
      integer, intent(in) :: neq, nvar

      real(DP), dimension(nvar,neq), intent(inout) :: u
      
      ! local vairables
      integer :: ieq
      real(DP) :: dradius

      !$omp parallel do private(dradius)
      do ieq = 1, neq
        
        ! Compute distance to origin
        dradius = sqrt(DvertexCoords(1,ieq)**2 +&
                       DvertexCoords(2,ieq)**2)
        
        if (dradius .lt. 1.0_DP) then
          u(1,ieq) = ddensity1
        elseif (dradius .gt. 1.05_DP) then
          u(1,ieq) = ddensity3
        else
          u(1,ieq) = ddensity2
        end if
      end do
      !$omp end parallel do
      
    end subroutine doRedistributeInterleaveFormat

    !**************************************************************
    ! Redistribute mass density stored in block format
    
    subroutine doRedistributeBlockFormat(neq, nvar, ddensity1,&
        ddensity2, ddensity3, DvertexCoords, ML, u)

      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), dimension(:), intent(in) :: ML
      real(DP), intent(in) :: ddensity1,ddensity2,ddensity3
      integer, intent(in) :: neq, nvar

      real(DP), dimension(neq,nvar), intent(inout) :: u
      
      ! local vairables
      integer :: ieq
      real(DP) :: dradius

      !$omp parallel do private(dradius)
      do ieq = 1, neq
        
        ! Compute distance to origin
        dradius = sqrt(DvertexCoords(1,ieq)**2 +&
                       DvertexCoords(2,ieq)**2)
        
        if (dradius .lt. 1.0_DP) then
          u(ieq,1) = ddensity1
        elseif (dradius .gt. 1.05_DP) then
          u(ieq,1) = ddensity3
        else
          u(ieq,1) = ddensity2
        end if
      end do
      !$omp end parallel do
      
    end subroutine doRedistributeBlockFormat
    
  end subroutine zpinch_redistributeMassEuler

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_redistributeMassTransport(rparlist, ssectionName,&
      rproblemLevel, rsolution, rcollection)

!<description>
    ! This subroutine redistributes the mass so that the total mass
    ! is equal to the total mass of Bank and Shadid
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section names in parameter list
    character(LEN=*), intent(in) :: ssectionName

    ! problem level structure
    type(t_problemLevel), intent(in), target :: rproblemLevel
!</input>

!<inputoutput>
    ! solution vector
    type(t_vectorBlock), intent(inout) :: rsolution

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_triangulation), pointer :: p_rtriangulation
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_ML, p_u
    real(DP) :: ddensity, dmass, dradius
    integer :: lumpedMassMatrix
    integer :: ivt

    ! Get global configuration from parameter list
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'lumpedmassmatrix', lumpedMassMatrix)

    ! Compute total density
    ddensity = SYS_PI * (1.05_DP**2-1._DP) * 1e6_DP
    
    ! Set pointer to triangulation
    p_rtriangulation => rproblemLevel%rtriangulation

    ! Set pointer to coordinate vector
    call storage_getbase_double2D(&
        p_rtriangulation%h_DvertexCoords, p_DvertexCoords)

    ! Set pointer to lumped mass matrix
    call lsyssc_getbase_double(&
        rproblemLevel%Rmatrix(lumpedMassMatrix), p_ML)

    ! Set pointer to conservative variables
    call lsysbl_getbase_double(rsolution, p_u)

    ! Compute total mass
    dmass = 0.0_DP
    
    ! Loop over all vertices
    do ivt = 1, p_rtriangulation%NVT
      
      ! Compute distance to origin
      dradius = sqrt(p_DvertexCoords(1,ivt)**2 +&
                     p_DvertexCoords(2,ivt)**2)
      
      if (dradius .ge. 1.0_DP .and.&
          dradius .le. 1.05_DP) then
        dmass = dmass + p_ML(ivt)
      end if
    end do
    
    ddensity = ddensity/dmass

    !$omp parallel do private(dradius)
    do ivt = 1, p_rtriangulation%NVT

      ! Compute distance to origin
      dradius = sqrt(p_DvertexCoords(1,ivt)**2 +&
                     p_DvertexCoords(2,ivt)**2)
    
      if (dradius .ge. 1.0_DP .and.&
          dradius .le. 1.05_DP) then
        p_u(ivt) = ddensity
      else
        p_u(ivt) = 0.0_DP
      end if
    end do
    !$omp end parallel do
    
  end subroutine zpinch_redistributeMassTransport

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_applySourceTerm(rparlist, ssectionName,&
      ssectionNameEuler, ssectionNameTransport, rproblemLevel,&
      rtimestep, rsolutionTransport, rsolutionEuler, rcollection)

!<description>
    ! This subroutine evaluates the Lorentz force term based on the
    ! given solution vectors and applies it to the Euler model
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section names in parameter list
    character(LEN=*), intent(in) :: ssectionName
    character(LEN=*), intent(in) :: ssectionNameEuler
    character(LEN=*), intent(in) :: ssectionNameTransport

    ! solution vector for transport model
    type(t_vectorBlock), intent(in) :: rsolutionTransport
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel
    
    ! time-stepping structure
    type(t_timestep), intent(inout) :: rtimestep
        
    ! solution vector for Euler model
    type(t_vectorBlock), intent(inout) :: rsolutionEuler
    
    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>


    ! local variable
    type(t_fparser), pointer :: p_rfparser
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_DdataTransport, p_DdataEuler
    real(DP), dimension(:), pointer :: p_DconsistentMassMatrix, p_DlumpedMassMatrix
    integer, dimension(:), pointer :: p_Kld, p_Kcol
    character(LEN=SYS_STRLEN) :: slorentzforceName
    real(DP) :: dscale
    integer :: isystemFormat, consistentMassMatrix, lumpedMassMatrix
    integer :: neq, nvar, icomp, icoords
    
    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'slorentzforcename', slorentzforceName)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'icoords', icoords)
    call parlst_getvalue_int(rparlist,&
          ssectionNameEuler, 'lumpedmassmatrix', lumpedMassMatrix)
    call parlst_getvalue_int(rparlist,&
          ssectionNameTransport, 'consistentmassmatrix', consistentMassMatrix)
    call parlst_getvalue_int(rparlist,&
          ssectionNameEuler, 'isystemformat', isystemFormat)

    ! Get lumped and consistent mass matrix
    call lsyssc_getbase_double(rproblemLevel&
        %Rmatrix(lumpedMassMatrix), p_DlumpedMassMatrix)
    call lsyssc_getbase_double(rproblemLevel&
        %Rmatrix(consistentMassMatrix), p_DconsistentMassMatrix)
    call lsyssc_getbase_Kld(rproblemLevel&
        %Rmatrix(consistentMassMatrix), p_Kld)
    call lsyssc_getbase_Kcol(rproblemLevel&
        %Rmatrix(consistentMassMatrix), p_Kcol)
    
    ! Set pointer to global solution vectors
    call lsysbl_getbase_double(rsolutionEuler, p_DdataEuler)
    call lsysbl_getbase_double(rsolutionTransport, p_DdataTransport)
    
    ! Set pointer to the vertex coordinates
    call storage_getbase_double2D(rproblemLevel%rtriangulation&
        %h_DvertexCoords, p_DvertexCoords)
    
    ! Set dimensions
    neq  = rsolutionTransport%NEQ
    nvar = euler_getNVAR(rproblemLevel)

    ! Get function parser from collection structure
    p_rfparser => collct_getvalue_pars(rcollection, 'rfparser')

    ! Get the number of the component used for 
    ! evaluating the Lorentz force term
    icomp = fparser_getFunctionNumber(p_rfparser, slorentzforceName)

    ! Evaluate the function parser
    call fparser_evalFunction(p_rfparser, icomp, (/rtimestep%dTime/), dscale)
    
    ! Multiply scaling parameter by the time step
    dscale = dscale * rtimestep%dStep

    ! What type of system format are we?
    select case(isystemFormat)

    case (SYSTEM_INTERLEAVEFORMAT)

      ! What type of coordinate system are we?
      select case(icoords)
      case(1)
        call calcSourceXYInterleaveFormat(dscale, neq, nvar,&
            p_DvertexCoords, p_Kld, p_Kcol, p_DconsistentMassMatrix,&
            p_DlumpedMassMatrix, p_DdataTransport, p_DdataEuler)

      case(2)
        call calcSourceRZInterleaveFormat(dscale, neq, nvar,&
            p_DvertexCoords, p_Kld, p_Kcol, p_DconsistentMassMatrix,&
            p_DlumpedMassMatrix, p_DdataTransport, p_DdataEuler)
        
      case default
        call output_line('Invalid type of coordinate system!',&
            OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcSourceTerm')
        call sys_halt()
      end select

    case (SYSTEM_BLOCKFORMAT)
      ! What type of coordinate system are we?
      select case(icoords)
      case(1)
        call calcSourceXYBlockFormat(dscale, neq, nvar,&
            p_DvertexCoords, p_Kld, p_Kcol, p_DconsistentMassMatrix,&
            p_DlumpedMassMatrix, p_DdataTransport, p_DdataEuler)
        
      case(2)
        call calcSourceRZBlockFormat(dscale, neq, nvar,&
            p_DvertexCoords, p_Kld, p_Kcol, p_DconsistentMassMatrix,&
            p_DlumpedMassMatrix, p_DdataTransport, p_DdataEuler)
        
      case default
        call output_line('Invalid type of coordinate system!',&
            OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcSourceTerm')
        call sys_halt()
      end select

    case DEFAULT
      call output_line('Invalid system format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcSourceTerm')
      call sys_halt()
    end select
    
  contains
    
    ! Here, the real working routines follow
    
    !**************************************************************
    ! Calculate the source term in x-y coordinates.
    ! The system is stored in interleave format.
    
    subroutine calcSourceXYInterleaveFormat(dscale, neq, nvar,&
        DvertexCoords, Kld, Kcol, MC, ML, DdataTransport, DdataEuler)

      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), dimension(:), intent(in) :: MC,ML,DdataTransport
      real(DP), intent(in) :: dscale
      integer, dimension(:), intent(in) :: Kld, Kcol
      integer, intent(in) :: neq, nvar
      
      real(DP), dimension(nvar,neq), intent(inout) :: DdataEuler
      
      ! local variables
      real(DP), dimension(:,:), allocatable :: DsourceTerm
      real(DP) :: drad, dang, daux, x1, x2, p, rq
      integer :: i,j,ij


      ! Create temporal memory initialized by zero
      allocate(DsourceTerm(2,neq))
      call lalg_clearVector(DsourceTerm)

      ! Loop over all rows
      do i = 1, neq

        ! Loop over all columns
        do ij = Kld(i), Kld(i+1)-1

          ! Get columns number
          j = Kcol(ij)
        
          ! Get coordinates at node j
          x1 = DvertexCoords(1, j)
          x2 = DvertexCoords(2, j)
          
          ! Compute polar coordinates
          drad = sqrt(x1*x1 + x2*x2)
          dang = atan2(x2, x1)
          
          ! Compute unit vector into origin
          if (drad .gt. 1e-4) then
            x1 = -cos(dang)
            x2 = -sin(dang)
          else
            x1 = 0.0; x2 = 0.0
          end if
          
          ! Compute source term
          daux = dscale * MC(ij) * DdataTransport(j) / max(drad, 1.0e-4_DP)
                    
          ! Impose source values
          DsourceTerm(1, i) = DsourceTerm(1, i) + daux * x1
          DsourceTerm(2, i) = DsourceTerm(2, i) + daux * x2

        end do
      end do

      do i = 1, neq

        ! Compute kinetic energy from momentum values without source term
        rq = 0.5 * ( DdataEuler(2, i)*DdataEuler(2, i) +&
                     DdataEuler(3, i)*DdataEuler(3, i) ) / DdataEuler(1, i)

        ! Compute pressure value
        p = DdataEuler(4, i) - rq

        ! Update momentum equations
        DdataEuler(2, i) = DdataEuler(2, i) + DsourceTerm(1, i)/ML(i)
        DdataEuler(3, i) = DdataEuler(3, i) + DsourceTerm(2, i)/ML(i)

        ! Compute kinetic energy from momentum values with source term
        rq = 0.5 * ( DdataEuler(2, i)*DdataEuler(2, i) +&
                     DdataEuler(3, i)*DdataEuler(3, i) ) / DdataEuler(1, i)

        ! Update total energy equation
        DdataEuler(4, i) = p + rq

      end do
      
      deallocate(DsourceTerm)

    end subroutine calcSourceXYInterleaveFormat

    !**************************************************************
    ! Calculate the source term in r-z coordinates.
    ! The system is stored in interleave format.
    
    subroutine calcSourceRZInterleaveFormat(dscale, neq, nvar,&
        DvertexCoords, Kld, Kcol, MC, ML, DdataTransport, DdataEuler)

      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), dimension(:), intent(in) :: MC,ML,DdataTransport
      real(DP), intent(in) :: dscale
      integer, dimension(:), intent(in) :: Kld, Kcol
      integer, intent(in) :: neq, nvar
      
      real(DP), dimension(nvar,neq), intent(inout) :: DdataEuler
      
      ! local variables
      real(DP), dimension(:,:), allocatable :: DsourceTerm
      real(DP) :: drad, dang, daux, x1, x2, p, rq
      integer :: i,j,ij


      ! Create temporal memory initialized by zero
      allocate(DsourceTerm(2,neq))
      call lalg_clearVector(DsourceTerm)

      ! Loop over all rows
      do i = 1, neq

        ! Loop over all columns
        do ij = Kld(i), Kld(i+1)-1

          ! Get columns number
          j = Kcol(ij)
        
          ! Get coordinates at node j
          drad = DvertexCoords(1, j)
          
          ! Compute unit vector into origin
          if (drad .gt. 1e-4) then
            x1 = -1.0; x2 = 0.0
          else
            x1 = 0.0; x2 = 0.0
          end if
          
          ! Compute source term
          daux = dscale * MC(ij) * DdataTransport(j) / max(drad, 1.0e-4_DP)
                    
          ! Impose source values
          DsourceTerm(1, i) = DsourceTerm(1, i) + daux * x1
          DsourceTerm(2, i) = DsourceTerm(2, i) + daux * x2

        end do
      end do

      do i = 1, neq

        ! Compute kinetic energy from momentum values without source term
        rq = 0.5 * ( DdataEuler(2, i)*DdataEuler(2, i) +&
                     DdataEuler(3, i)*DdataEuler(3, i) ) / DdataEuler(1, i)

        ! Compute pressure value
        p = DdataEuler(4, i) - rq

        ! Update momentum equations
        DdataEuler(2, i) = DdataEuler(2, i) + DsourceTerm(1, i)/ML(i)
        DdataEuler(3, i) = DdataEuler(3, i) + DsourceTerm(2, i)/ML(i)

        ! Compute kinetic energy from momentum values with source term
        rq = 0.5 * ( DdataEuler(2, i)*DdataEuler(2, i) +&
                     DdataEuler(3, i)*DdataEuler(3, i) ) / DdataEuler(1, i)

        ! Update total energy equation
        DdataEuler(4, i) = p + rq

      end do
      
      deallocate(DsourceTerm)

    end subroutine calcSourceRZInterleaveFormat
    
    !**************************************************************
    ! Calculate the source term in x-y coordinates.
    ! The system is stored in block format.

    subroutine calcSourceXYBlockFormat(dscale, neq, nvar,&
        DvertexCoords, Kld, Kcol, MC, ML, DdataTransport, DdataEuler)

      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), dimension(:), intent(in) :: MC,ML,DdataTransport
      real(DP), intent(in) :: dscale
      integer, dimension(:), intent(in) :: Kld, Kcol
      integer, intent(in) :: neq, nvar
      
      real(DP), dimension(neq,nvar), intent(inout) :: DdataEuler
      
      ! local variables
      real(DP), dimension(:,:), allocatable :: DsourceTerm
      real(DP) :: drad, dang, daux, x1, x2, p, rq
      integer :: i,j,ij

      ! Create temporal memory initialized by zero
      allocate(DsourceTerm(2,neq))
      call lalg_clearVector(DsourceTerm)

      ! Loop over all rows
      do i = 1, neq

        ! Loop over all columns
        do ij = Kld(i), Kld(i+1)-1

          ! Get columns number
          j = Kcol(ij)
        
          ! Get coordinates at node j
          x1 = DvertexCoords(1, j)
          x2 = DvertexCoords(2, j)
          
          ! Compute polar coordinates
          drad = sqrt(x1*x1 + x2*x2)
          dang = atan2(x2, x1)
          
          ! Compute unit vector into origin
          if (drad .gt. 1e-4) then
            x1 = -cos(dang)
            x2 = -sin(dang)
          else
            x1 = 0.0; x2 = 0.0
          end if
          
          ! Compute source term
          if (DdataTransport(j) > SYS_EPSREAL) then
            daux = dscale * MC(ij) * DdataTransport(j) / max(drad, 1.0e-4_DP)
          else
            daux = 0.0_DP
          end if
                    
          ! Impose source values
          DsourceTerm(1, i) = DsourceTerm(1, i) + daux * x1
          DsourceTerm(2, i) = DsourceTerm(2, i) + daux * x2

        end do
      end do

      do i = 1, neq

        ! Compute kinetic energy from momentum values without source term
        rq = 0.5 * ( DdataEuler(i, 2)*DdataEuler(i, 2) +&
                     DdataEuler(i, 3)*DdataEuler(i, 3) ) /&
                     DdataEuler(i, 1)

        ! Compute pressure value
        p = DdataEuler(i, 4) - rq

        ! Update momentum equations
        DdataEuler(i, 2) = DdataEuler(i, 2) + DsourceTerm(1, i)/ML(i)
        DdataEuler(i, 3) = DdataEuler(i, 3) + DsourceTerm(2, i)/ML(i)

        ! Compute kinetic energy from momentum values with source term
        rq = 0.5 * ( DdataEuler(i, 2)*DdataEuler(i, 2) +&
                     DdataEuler(i, 3)*DdataEuler(i, 3) ) /&
                     DdataEuler(i, 1)

        ! Update total energy equation
        DdataEuler(i, 4) = p + rq

      end do
      
      deallocate(DsourceTerm)

    end subroutine calcSourceXYBlockFormat

    !**************************************************************
    ! Calculate the source term in r-z coordinates.
    ! The system is stored in block format.

    subroutine calcSourceRZBlockFormat(dscale, neq, nvar,&
        DvertexCoords, Kld, Kcol, MC, ML, DdataTransport, DdataEuler)

      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), dimension(:), intent(in) :: MC,ML,DdataTransport
      real(DP), intent(in) :: dscale
      integer, dimension(:), intent(in) :: Kld, Kcol
      integer, intent(in) :: neq, nvar
      
      real(DP), dimension(neq,nvar), intent(inout) :: DdataEuler
      
      ! local variables
      real(DP), dimension(:,:), allocatable :: DsourceTerm
      real(DP) :: drad, dang, daux, x1, x2, p, rq
      integer :: i,j,ij

      ! Create temporal memory initialized by zero
      allocate(DsourceTerm(2,neq))
      call lalg_clearVector(DsourceTerm)

      ! Loop over all rows
      do i = 1, neq

        ! Loop over all columns
        do ij = Kld(i), Kld(i+1)-1

          ! Get columns number
          j = Kcol(ij)
        
          ! Get coordinates at node j
          drad = DvertexCoords(1, j)
          
          ! Compute unit vector into origin
          if (drad .gt. 1e-4) then
            x1 = -1.0; x2 = 0.0
          else
            x1 = 0.0; x2 = 0.0
          end if
          
          ! Compute source term
          if (DdataTransport(j) > SYS_EPSREAL) then
            daux = dscale * MC(ij) * DdataTransport(j) / max(drad, 1.0e-4_DP)
          else
            daux = 0.0_DP
          end if
                    
          ! Impose source values
          DsourceTerm(1, i) = DsourceTerm(1, i) + daux * x1
          DsourceTerm(2, i) = DsourceTerm(2, i) + daux * x2

        end do
      end do

      do i = 1, neq

        ! Compute kinetic energy from momentum values without source term
        rq = 0.5 * ( DdataEuler(i, 2)*DdataEuler(i, 2) +&
                     DdataEuler(i, 3)*DdataEuler(i, 3) ) /&
                     DdataEuler(i, 1)

        ! Compute pressure value
        p = DdataEuler(i, 4) - rq

        ! Update momentum equations
        DdataEuler(i, 2) = DdataEuler(i, 2) + DsourceTerm(1, i)/ML(i)
        DdataEuler(i, 3) = DdataEuler(i, 3) + DsourceTerm(2, i)/ML(i)

        ! Compute kinetic energy from momentum values with source term
        rq = 0.5 * ( DdataEuler(i, 2)*DdataEuler(i, 2) +&
                     DdataEuler(i, 3)*DdataEuler(i, 3) ) /&
                     DdataEuler(i, 1)

        ! Update total energy equation
        DdataEuler(i, 4) = p + rq

      end do
      
      deallocate(DsourceTerm)

    end subroutine calcSourceRZBlockFormat

  end subroutine zpinch_applySourceTerm

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_calcAdaptationIndicator(rsolutionEuler,&
      rsolutionTransport, rindicator)

!<description>
    ! This subroutine computes the element-wise indicator for mesh
    ! refinement and re-coarsening based on the detectation of shocks
    ! and contact discontinuities and the strength of gradient variations
!</description>

!<input>
    ! solution vectors
    type(t_vectorBlock), intent(in) :: rsolutionEuler, rsolutionTransport
!</input>
      
!<inputoutput>
    ! local feature indicator
    type(t_vectorScalar), intent(inout) :: rindicator
!</inputoutput>
!</subroutine>
    
    ! local variables
    type(t_vectorScalar) :: rdensity, rpressure
    type(t_triangulation), pointer :: p_rtriangulation
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Ddensity, p_Dpressure, p_Dtracer, p_Dindicator
    integer, dimension(:,:), pointer :: p_IverticesAtElement, p_IneighboursAtElement
    logical, dimension(:), pointer :: p_BisActiveElement
    real(DP) :: dgradient,dlambda_i, dlambda_j
    integer :: iel,jel,ive,nve,i,j,iprotectLayer
    
    real(DP), parameter :: dEpsS = 0.2_DP
    real(DP), parameter :: dEpsC = 0.1_DP


    ! Set pointer to the underlying triangulation
    p_rtriangulation => rsolutionEuler%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation
    
    ! Extract primitive variables from conservative solution values
    call euler_getVariable(rsolutionEuler, 'density', rdensity)
    call euler_getVariable(rsolutionEuler, 'pressure', rpressure)

    ! Create element-wise indicator
    call lsyssc_createVector(rindicator, p_rtriangulation%NEL, .true.)
    
    ! Set pointers
    call storage_getbase_double2D(p_rtriangulation%h_DvertexCoords,&
        p_DvertexCoords)
    call storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    call storage_getbase_int2D(p_rtriangulation&
        %h_IneighboursAtElement, p_IneighboursAtElement)
    call lsyssc_getbase_double(rdensity, p_Ddensity)
    call lsyssc_getbase_double(rpressure, p_Dpressure)
    call lsyssc_getbase_double(rindicator, p_Dindicator)
    call lsysbl_getbase_double(rsolutionTransport, p_Dtracer)
    
    ! Loop over all elements
    elements: do iel = 1, p_rtriangulation%NEL

      ! Get number of vertices at element
      nve = tria_getNVE(p_IverticesAtElement, iel)

      ! Initialize gradient value
      dgradient = 0.0_DP

      ! Check element for shocks and contact discontinuities
      vertices: do ive = 1, nve
        
        ! Get global vertex numbers
        i = p_IverticesAtElement(ive, iel)
        j = p_IverticesAtElement(mod(ive,nve)+1, iel)

        ! Compute nodal values of tracer
        dlambda_i = p_Dtracer(i)/p_Ddensity(i)
        dlambda_j = p_Dtracer(j)/p_Ddensity(j)
        
        ! Update maximum gradient function
        dgradient = max(dgradient,&
            abs(abs(p_Ddensity(i))  - abs(p_Ddensity(j)))  /&
              max(abs(p_Ddensity(i)),  abs(p_Ddensity(j))),&
            abs(abs(p_Dpressure(i)) - abs(p_Dpressure(j))) /&
              max(abs(p_Dpressure(i)), abs(p_Dpressure(j))),&
            abs(abs(dlambda_i)      - abs(dlambda_j))      /&
              max(1e-2, abs(dlambda_i), abs(dlambda_j))     )
      end do vertices
      
      ! If we end up here, then the maximum gradient function is adopted
      p_Dindicator(iel) = dgradient
    end do elements
    
    ! Release temporal memory
    call lsyssc_releaseVector(rdensity)
    call lsyssc_releaseVector(rpressure)
    
    ! Add protection layers
    allocate(p_BisActiveElement(p_rtriangulation%NEL))
    
    do iprotectLayer = 1, 4

      p_BisActiveElement = .false.

      do iel = 1, p_rtriangulation%NEL

        if (p_BisactiveElement(iel)) cycle
        if (p_Dindicator(iel) .le. 0.8) cycle

        do ive = 1, tria_getNVE(p_IverticesAtElement, iel)
          jel = p_IneighboursAtElement(ive, iel)
          if (jel .eq. 0) cycle
          if (p_BisactiveElement(jel)) then
            p_Dindicator(jel) = max(p_Dindicator(jel), p_Dindicator(iel))
          else
            if (p_Dindicator(jel) .lt. 0.8) then
              p_Dindicator(jel) = max(p_Dindicator(jel), p_Dindicator(iel))
              p_BisactiveElement(jel) = .true.
            end if
          end if
        end do
      end do
    end do

    deallocate(p_BisActiveElement)
    
  end subroutine zpinch_calcAdaptationIndicator
 
  !*****************************************************************************

!<subroutine>

  subroutine zpinch_outputSolution(rparlist, ssectionName,&
      rproblemLevel, rsolution, dtime)

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

    ! array of solution vectors
    type(t_vectorBlock), dimension(:), intent(in), optional :: rsolution

    ! OPTIONAL: simulation time
    real(DP), intent(in), optional :: dtime
!</input>
!</subroutine>

    ! section names
    character(LEN=SYS_STRLEN) :: soutputName
    character(LEN=SYS_STRLEN) :: ssectionNameEuler
    character(LEN=SYS_STRLEN) :: ucdsolution

    ! persistent variable
    integer, save :: ifilenumber = 1
    
    ! local variables
    type(t_ucdExport) :: rexport
    type(t_vectorScalar) :: rvector1, rvector2, rvector3
    real(DP), dimension(:), pointer :: p_Dsolution, p_Ddata1, p_Ddata2, p_Ddata3
    integer :: iformatUCD, isystemFormat, isize, ndim


    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'output', soutputName)
    call parlst_getvalue_string(rparlist,&
        trim(soutputName), 'ucdsolution', ucdsolution)
    call parlst_getvalue_int(rparlist,&
        trim(soutputName), 'iformatucd', iformatUCD)

    call parlst_getvalue_string(rparlist,&
        ssectionName, 'application_euler', ssectionNameEuler)
    call parlst_getvalue_int(rparlist,&
        ssectionNameEuler, 'isystemformat', isystemformat)
    call parlst_getvalue_int(rparlist,&
        ssectionNameEuler, 'isystemformat', isystemformat)

    ! Initialize the UCD exporter
    call flagship_initUCDexport(rproblemLevel, ucdsolution,&
        iformatUCD, rexport, ifilenumber)

    ! Increase filenumber by one
    ifilenumber = ifilenumber+1

    ! Set simulation time
    if (present(dtime)) call ucd_setSimulationTime(rexport, dtime)

    ! Add solution vectors
    if (present(rsolution)) then
      
      ! Add solution for Euler model
      call lsysbl_getbase_double(rsolution(1), p_Dsolution)
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
          call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR1D,&
              'velocity_x', p_Dsolution, p_Ddata1)
          call ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1)
          
          call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR1D,&
              'density', p_Dsolution, p_Ddata1)
          call ucd_addVariableVertexBased (rexport, 'density',&
              UCD_VAR_STANDARD, p_Ddata1)
          
          call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR1D,&
              'energy', p_Dsolution, p_Ddata1)
          call ucd_addVariableVertexBased (rexport, 'energy',&
              UCD_VAR_STANDARD, p_Ddata1)
          
          call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR1D,&
              'pressure', p_Dsolution, p_Ddata1)
          call ucd_addVariableVertexBased (rexport, 'pressure',&
              UCD_VAR_STANDARD, p_Ddata1)
          
        case (NDIM2D)
          call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR2D,&
              'velocity_x', p_Dsolution, p_Ddata1)
          call euler_getVarInterleaveFormat(rvector2%NEQ, NVAR2D,&
              'velocity_y', p_Dsolution, p_Ddata2)
          call ucd_addVarVertBasedVec(rexport, 'velocity',&
              p_Ddata1, p_Ddata2)
          
          call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR2D,&
              'density', p_Dsolution, p_Ddata1)
          call ucd_addVariableVertexBased (rexport, 'density',&
              UCD_VAR_STANDARD, p_Ddata1)
          
          call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR2D,&
              'energy', p_Dsolution, p_Ddata1)
          call ucd_addVariableVertexBased (rexport, 'energy',&
              UCD_VAR_STANDARD, p_Ddata1)
          
          call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR2D,&
              'pressure', p_Dsolution, p_Ddata1)
          call ucd_addVariableVertexBased (rexport, 'pressure',&
              UCD_VAR_STANDARD, p_Ddata1)
          
        case (NDIM3D)
          call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR3D,&
              'velocity_x', p_Dsolution, p_Ddata1)
          call euler_getVarInterleaveFormat(rvector2%NEQ, NVAR3D,&
              'velocity_y', p_Dsolution, p_Ddata2)
          call euler_getVarInterleaveFormat(rvector3%NEQ, NVAR3D, &
              'velocity_z', p_Dsolution, p_Ddata3)
          call ucd_addVarVertBasedVec(rexport, 'velocity',&
              p_Ddata1, p_Ddata2, p_Ddata3)
          
          call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR3D,&
              'density', p_Dsolution, p_Ddata1)
          call ucd_addVariableVertexBased (rexport, 'density',&
              UCD_VAR_STANDARD, p_Ddata1)
          
          call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR3D,&
              'energy', p_Dsolution, p_Ddata1)
          call ucd_addVariableVertexBased (rexport, 'energy',&
              UCD_VAR_STANDARD, p_Ddata1)
          
          call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR3D,&
              'pressure', p_Dsolution, p_Ddata1)
          call ucd_addVariableVertexBased (rexport, 'pressure',&
              UCD_VAR_STANDARD, p_Ddata1)
          
        end select
        
        
      case (SYSTEM_BLOCKFORMAT)
        
        select case(ndim)
        case (NDIM1D)
          call euler_getVarBlockFormat(rvector1%NEQ, NVAR1D,&
              'velocity_x', p_Dsolution, p_Ddata1)
          call ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1)
          
          call euler_getVarBlockFormat(rvector1%NEQ, NVAR1D,&
              'density', p_Dsolution, p_Ddata1)
          call ucd_addVariableVertexBased (rexport, 'density',&
              UCD_VAR_STANDARD, p_Ddata1)
          
          call euler_getVarBlockFormat(rvector1%NEQ, NVAR1D,&
              'energy', p_Dsolution, p_Ddata1)
          call ucd_addVariableVertexBased (rexport, 'energy',&
              UCD_VAR_STANDARD, p_Ddata1)
          
          call euler_getVarBlockFormat(rvector1%NEQ, NVAR1D,&
              'effective_energy', p_Dsolution, p_Ddata1)
          call ucd_addVariableVertexBased (rexport,&
              'effective_energy', UCD_VAR_STANDARD, p_Ddata1)
          
          call euler_getVarBlockFormat(rvector1%NEQ, NVAR1D,&
              'pressure', p_Dsolution, p_Ddata1)
          call ucd_addVariableVertexBased (rexport, 'pressure',&
              UCD_VAR_STANDARD, p_Ddata1)
          
          call euler_getVarBlockFormat(rvector1%NEQ, NVAR1D,&
              'machnumber', p_Dsolution, p_Ddata1)
          call ucd_addVariableVertexBased (rexport, 'machnumber',&
              UCD_VAR_STANDARD, p_Ddata1)
          
          
        case (NDIM2D)
          call euler_getVarBlockFormat(rvector1%NEQ, NVAR2D,&
              'velocity_x', p_Dsolution, p_Ddata1)
          call euler_getVarBlockFormat(rvector2%NEQ, NVAR2D,&
              'velocity_y', p_Dsolution, p_Ddata2)
          call ucd_addVarVertBasedVec(rexport, 'velocity',&
              p_Ddata1, p_Ddata2)
          
          call euler_getVarBlockFormat(rvector1%NEQ, NVAR2D,&
              'density', p_Dsolution, p_Ddata1)
          call ucd_addVariableVertexBased (rexport, 'density',&
              UCD_VAR_STANDARD, p_Ddata1)
          
          call euler_getVarBlockFormat(rvector1%NEQ, NVAR2D,&
              'energy', p_Dsolution, p_Ddata1)
          call ucd_addVariableVertexBased (rexport, 'energy',&
              UCD_VAR_STANDARD, p_Ddata1)
          
          call euler_getVarBlockFormat(rvector1%NEQ, NVAR2D,&
              'effective_energy', p_Dsolution, p_Ddata1)
          call ucd_addVariableVertexBased (rexport,&
              'effective_energy', UCD_VAR_STANDARD, p_Ddata1)
          
          call euler_getVarBlockFormat(rvector1%NEQ, NVAR2D,&
              'pressure', p_Dsolution, p_Ddata1)
          call ucd_addVariableVertexBased (rexport, 'pressure',&
              UCD_VAR_STANDARD, p_Ddata1)
          
          call euler_getVarBlockFormat(rvector1%NEQ, NVAR2D,&
              'machnumber', p_Dsolution, p_Ddata1)
          call ucd_addVariableVertexBased (rexport, 'machnumber',&
              UCD_VAR_STANDARD, p_Ddata1)
          
        case (NDIM3D)
          call euler_getVarBlockFormat(rvector1%NEQ, NVAR3D,&
              'velocity_x', p_Dsolution, p_Ddata1)
          call euler_getVarBlockFormat(rvector2%NEQ, NVAR3D,&
              'velocity_y', p_Dsolution, p_Ddata2)
          call euler_getVarBlockFormat(rvector3%NEQ, NVAR3D,&
              'velocity_z', p_Dsolution, p_Ddata3)
          call ucd_addVarVertBasedVec(rexport, 'velocity',&
              p_Ddata1, p_Ddata2, p_Ddata3)
          
          call euler_getVarBlockFormat(rvector1%NEQ, NVAR3D,&
              'density', p_Dsolution, p_Ddata1)
          call ucd_addVariableVertexBased (rexport, 'density',&
              UCD_VAR_STANDARD, p_Ddata1)
          
          call euler_getVarBlockFormat(rvector1%NEQ, NVAR3D,&
              'energy', p_Dsolution, p_Ddata1)
          call ucd_addVariableVertexBased (rexport, 'energy',&
              UCD_VAR_STANDARD, p_Ddata1)
          
          call euler_getVarBlockFormat(rvector1%NEQ, NVAR3D, 'effectiv' // &
              'e_energy', p_Dsolution, p_Ddata1)
          call ucd_addVariableVertexBased (rexport,&
              'effective_energy', UCD_VAR_STANDARD, p_Ddata1)
          
          call euler_getVarBlockFormat(rvector1%NEQ, NVAR3D,&
              'pressure', p_Dsolution, p_Ddata1)
          call ucd_addVariableVertexBased (rexport, 'pressure',&
              UCD_VAR_STANDARD, p_Ddata1)
          
          call euler_getVarBlockFormat(rvector1%NEQ, NVAR3D,&
              'machnumber', p_Dsolution, p_Ddata1)
          call ucd_addVariableVertexBased (rexport, 'machnumber',&
              UCD_VAR_STANDARD, p_Ddata1)
          
        end select
        
      end select
      
      ! Release temporal memory
      call lsyssc_releaseVector(rvector1)
      call lsyssc_releaseVector(rvector2)
      call lsyssc_releaseVector(rvector3)

      ! Add solution for transport model
      call lsysbl_getbase_double(rsolution(2), p_Ddata1)
      call ucd_addVariableVertexBased (rexport, 'advect',&
          UCD_VAR_STANDARD, p_Ddata1)
    end if
    
    ! Write UCD file
    call ucd_write  (rexport)
    call ucd_release(rexport)

  end subroutine zpinch_outputSolution

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_adaptTriangulation(rparlist, ssectionName, rhadapt,&
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

    ! Get parameters from parameter list
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'isystemformat', isystemFormat)

    ! What type of system format are we?
    select case (isystemFormat)

    case (SYSTEM_INTERLEAVEFORMAT)

      if (rtriangulationSrc%ndim .eq. NDIM2D) then
        call zpinch_hadaptCallbackScalar2D(rcollection,&
            HADAPT_OPR_INITCALLBACK, Ivalue, Ivalue)
        call hadapt_performAdaptation(rhadapt, rindicator,&
            rcollection, zpinch_hadaptCallbackScalar2D)
      else
        call output_line('Invalid spatial dimension!',&
            OU_CLASS_ERROR,OU_MODE_STD,'zpinch_adaptTriangulation')
        call sys_halt()
      end if

    case (SYSTEM_BLOCKFORMAT)

      if (rtriangulationSrc%ndim .eq. NDIM2D) then
        call zpinch_hadaptCallbackBlock2D(rcollection,&
            HADAPT_OPR_INITCALLBACK, Ivalue, Ivalue)
        call hadapt_performAdaptation(rhadapt, rindicator,&
            rcollection, zpinch_hadaptCallbackBlock2D)
      else
        call output_line('Invalid spatial dimension!',&
            OU_CLASS_ERROR,OU_MODE_STD,'zpinch_adaptTriangulation')
        call sys_halt()
      end if
      
    case DEFAULT

      call output_line('Invalid type of system format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'zpinch_adaptTriangulation')
      call sys_halt()
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

  end subroutine zpinch_adaptTriangulation

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_solveTransientPrimal(rparlist, ssectionName,&
      ssectionNameEuler, ssectionNameTransport, rbdrCondEuler,&
      rbdrCondTransport, rproblem, rtimestep, rsolver, rsolution,&
      rcollection)

!<description>
      ! This subroutine solves the transient primal simplified MHD problem.
!</description>

!<input>
    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName
    character(LEN=*), intent(in) :: ssectionNameEuler
    character(LEN=*), intent(in) :: ssectionNameTransport

    ! boundary condition structure
    type(t_boundaryCondition), intent(in) :: rbdrCondEuler
    type(t_boundaryCondition), intent(in) :: rbdrCondTransport
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

    ! primal solution vector
    type(t_vectorBlock), dimension(:), intent(inout), target :: rsolution

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! Pointer to the current multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel
    
    ! Pointer to the discretisation structure
    type(t_blockDiscretisation), pointer :: p_rdiscretisation

    ! Pointer to the solver structures
    type(t_solver), pointer :: p_rsolverEuler, p_rsolverTransport

    ! Pointer to the solution vectors
    type(t_vectorBlock), pointer :: p_rsolutionEuler, p_rsolutionTransport

    ! Vector for the element-wise feature indicator
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
    
    ! Vector for source term
    type(t_vectorBlock), dimension(2) :: rforce

    ! Section names
    character(LEN=SYS_STRLEN) :: sadaptivityName
    character(LEN=SYS_STRLEN) :: soutputName
    character(LEN=SYS_STRLEN) :: ucdimport

    ! local variables
    type(t_ucdExport) :: rimport
    real(dp) :: dstepUCD, dtimeUCD, dstepAdapt, dtimeAdapt, dscale
    integer :: templateMatrix, systemMatrix, isystemFormat
    integer :: discretisationEuler, discretisationTransport
    integer :: isize, ipreadapt, npreadapt, ndimension
    integer, external :: signal_SIGINT

    real(DP) :: dmassEuler0, dmassEuler
    real(DP) :: dmassTransport0, dmassTransport
    
    ! Get timer structures
    p_rtimerPrePostprocess =>&
        collct_getvalue_timer(rcollection, 'rtimerPrePostprocess')
    p_rtimerSolution =>&
        collct_getvalue_timer(rcollection, 'rtimerSolution')
    p_rtimerErrorEstimation =>&
        collct_getvalue_timer(rcollection, 'rtimerErrorEstimation')
    p_rtimerAdaptation =>&
        collct_getvalue_timer(rcollection, 'rtimerAdaptation')
    p_rtimerTriangulation =>&
        collct_getvalue_timer(rcollection, 'rtimerTriangulation')
    p_rtimerAssemblyCoeff =>&
        collct_getvalue_timer(rcollection, 'rtimerAssemblyCoeff')
    
    ! Set pointers to solver structures
    p_rsolverEuler     => solver_getNextSolver(rsolver, 1)
    p_rsolverTransport => solver_getNextSolver(rsolver, 2)
    
    ! Set pointers to solution vectors
    p_rsolutionEuler => rsolution(1)
    p_rsolutionTransport => rsolution(2)

    ! Start time measurement for pre-processing
    call stat_startTimer(p_rtimerPrePostprocess, STAT_TIMERSHORT)

    ! Get global parameters
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'ndimension', ndimension)
    
    ! Set pointer to maximum problem level
    p_rproblemLevel => rproblem%p_rproblemLevelMax
    
    
    !--- compressible Euler model ----------------------------------------------

    ! Get position of discretisation structure
    call parlst_getvalue_int(rparlist,&
        ssectionNameEuler, 'discretisation', discretisationEuler)
    
    ! Set pointer to discretisation structure
    p_rdiscretisation =>&
        p_rproblemLevel%Rdiscretisation(discretisationEuler)

    ! Create the solution vector
    call lsysbl_createVectorBlock(p_rdiscretisation,&
        p_rsolutionEuler, .false., ST_DOUBLE)
    if (p_rdiscretisation%ncomponents .ne.&
        euler_getNVAR(p_rproblemLevel)) then
      p_rsolutionEuler%RvectorBlock(1)%NVAR =&
          euler_getNVAR(p_rproblemLevel)
      isize = p_rsolutionEuler%NEQ*euler_getNVAR(p_rproblemLevel)
      call lsysbl_resizeVectorBlock(p_rsolutionEuler,&
          isize, .false., .false.)
    end if

    ! Initialize the solution vector
    call euler_initSolution(rparlist, ssectionNameEuler,&
        p_rproblemLevel, rtimestep%dinitialTime, p_rsolutionEuler,&
        rcollection)

    ! Redistribute the density
    call zpinch_redistributeMassEuler(rparlist, ssectionNameEuler,&
        p_rproblemLevel, p_rsolutionEuler, rcollection)

    ! Impose boundary conditions
    select case(ndimension)
    case (NDIM1D)
      call bdrf_filterVectorExplicit(rbdrCondEuler, p_rsolutionEuler,&
          rtimestep%dinitialTime, euler_calcBoundaryvalues1d)
    case (NDIM2D)
      call bdrf_filterVectorExplicit(rbdrCondEuler, p_rsolutionEuler,&
          rtimestep%dinitialTime, euler_calcBoundaryvalues2d)
    case (NDIM3D)
      call bdrf_filterVectorExplicit(rbdrCondEuler, p_rsolutionEuler,&
          rtimestep%dinitialTime, euler_calcBoundaryvalues3d)
    end select

    ! Attach the boundary condition
    call solver_setBoundaryCondition(p_rsolverEuler, rbdrCondEuler, .true.)

    ! Set collection to primal problem mode
    call parlst_addvalue(rparlist, ssectionNameEuler, 'mode', 'primal')
    

    !--- transport model -------------------------------------------------------

    ! Get position of discretisation structure
    call parlst_getvalue_int(rparlist, ssectionNameTransport,&
        'discretisation', discretisationTransport)

    ! Set pointer to discretisation structure
    p_rdiscretisation =>&
        p_rproblemLevel%Rdiscretisation(discretisationTransport)

    ! Create the solution vector
    call lsysbl_createVectorBlock(p_rdiscretisation,&
        p_rsolutionTransport, .false., ST_DOUBLE)

    ! Initialize the solution vector and impose boundary conditions
    call transp_initSolution(rparlist, ssectionNameTransport,&
        p_rproblemLevel, rtimestep%dinitialTime,&
        p_rsolutionTransport, rcollection)

    ! Redistribute the tracer
    call zpinch_redistributeMassTransport(rparlist,&
        ssectionNameTransport, p_rproblemLevel, p_rsolutionTransport,&
        rcollection)
    
    ! Impose boundary conditions
    call bdrf_filterVectorExplicit(rbdrCondTransport,&
        p_rsolutionTransport, rtimestep%dinitialTime)
    
    ! Attach the boundary condition
    call solver_setBoundaryCondition(p_rsolverTransport,&
        rbdrCondTransport, .true.)

    ! Set collection to primal problem mode
    call parlst_addvalue(rparlist,&
        ssectionNameTransport, 'mode', 'primal')
       
    
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
        call hadapt_initFromParameterlist(rhadapt,&
            rparlist, sadaptivityName)

        ! Generate a dynamic graph for the sparsity pattern and attach
        ! it to the collection structure which is used to pass this
        ! type to the callback function for h-adaptation. Note that
        ! the template matrix is the same for the compressible Euler
        ! model and the transport model. Therefore, the
        ! template matrix is adopted from the transport Euler model.
        call parlst_getvalue_int(rparlist,&
            ssectionNameEuler, 'templateMatrix', templateMatrix)
        call grph_createGraphFromMatrix(p_rproblemLevel&
            %Rmatrix(templateMatrix), rgraph)
        call collct_setvalue_graph(rcollection,&
            'sparsitypattern', rgraph, .true.)
        
        ! Perform pre-adaptation?
        if (npreadapt > 0) then

          ! Perform number of pre-adaptation steps
          do ipreadapt = 1, npreadapt
            
            ! Compute the error estimator based on the tracer
            call zpinch_calcAdaptationIndicator(p_rsolutionEuler,&
                p_rsolutionTransport, relementError)

            ! Set the names of the template matrix
            rcollection%SquickAccess(1) = 'sparsitypattern'

            ! Attach the solution vectors to the collection structure
            rcollection%p_rvectorQuickAccess1 => rsolution(1)
            rcollection%p_rvectorQuickAccess2 => rsolution(2)

            ! Perform h-adaptation and update the triangulation structure
            call zpinch_adaptTriangulation(rparlist,&
                ssectionnameEuler, rhadapt, p_rproblemLevel&
                %rtriangulation, relementError, rcollection)
            
            ! Release element-wise error distribution
            call lsyssc_releaseVector(rElementError)

            ! Generate standard mesh from raw mesh
            call tria_initStandardMeshFromRaw(p_rproblemLevel&
                %rtriangulation, rproblem%rboundary)
            
            ! Update the template matrix according to the sparsity pattern
            call grph_generateMatrix(rgraph, p_rproblemLevel&
                %Rmatrix(templateMatrix))

            ! Re-initialize all constant coefficient matrices
            call euler_initProblemLevel(rparlist,&
                ssectionNameEuler, p_rproblemLevel, rcollection)
            call transp_initProblemLevel(rparlist,&
                ssectionNameTransport, p_rproblemLevel, rcollection)

            ! Resize the solution vector for the Euler model accordingly
            call parlst_getvalue_int(rparlist,&
                ssectionNameEuler, 'systemMatrix', systemMatrix)
            call lsysbl_resizeVecBlockIndMat(p_rproblemLevel&
                %RmatrixBlock(systemMatrix), p_rsolutionEuler, .false., .true.)
            
            ! Resize the solution vector for the transport model accordingly
            call lsysbl_resizeVectorBlock(p_rsolutionTransport, &
                p_rproblemLevel%Rmatrix(templateMatrix)%NEQ, .false.)

            ! Re-generate the initial solution vector for the Euler model
            call euler_initSolution(rparlist, ssectionnameEuler,&
                p_rproblemLevel, rtimestep%dinitialTime,&
                p_rsolutionEuler, rcollection)

            ! Redistribute the density
            call zpinch_redistributeMassEuler(rparlist,&
                ssectionNameEuler, p_rproblemLevel, p_rsolutionEuler,&
                rcollection)
            
            ! Impose boundary conditions
            select case(ndimension)
            case (NDIM1D)
              call bdrf_filterVectorExplicit(rbdrCondEuler,&
                  p_rsolutionEuler, rtimestep%dinitialTime,&
                  euler_calcBoundaryvalues1d)
              
            case (NDIM2D)
              call bdrf_filterVectorExplicit(rbdrCondEuler,&
                  p_rsolutionEuler, rtimestep%dinitialTime,&
                  euler_calcBoundaryvalues2d)

            case (NDIM3D)
              call bdrf_filterVectorExplicit(rbdrCondEuler,&
                  p_rsolutionEuler, rtimestep%dinitialTime,&
                  euler_calcBoundaryvalues3d)
            end select

            ! Re-generate the initial solution vector for the transport model
            call transp_initSolution(rparlist, ssectionnameTransport,&
                p_rproblemLevel, rtimestep%dinitialTime,&
                p_rsolutionTransport, rcollection)

            ! Redistribute the tracer
            call zpinch_redistributeMassTransport(rparlist,&
                ssectionNameTransport, p_rproblemLevel,&
                p_rsolutionTransport, rcollection)

            ! Impose boundary conditions
            call bdrf_filterVectorExplicit(rbdrCondTransport,&
                p_rsolutionTransport, rtimestep%dinitialTime)
          end do

          ! Prepare internal data arrays of the solver structure
          call parlst_getvalue_int(rparlist,&
              ssectionNameEuler, 'systemMatrix', systemMatrix)
          call parlst_getvalue_int(rparlist,&
              ssectionNameEuler, 'isystemFormat', isystemFormat)
          call flagship_updateSolverMatrix(p_rproblemLevel,&
              p_rsolverEuler, systemMatrix, isystemFormat, UPDMAT_ALL)
          call solver_updateStructure(p_rsolverEuler)
          
          ! Prepare internal data arrays of the solver structure
          call parlst_getvalue_int(rparlist,&
              ssectionNameTransport, 'systemMatrix', systemMatrix)
          call flagship_updateSolverMatrix(p_rproblemLevel,&
              p_rsolverTransport, systemMatrix,&
              SYSTEM_INTERLEAVEFORMAT , UPDMAT_ALL)
          call solver_updateStructure(p_rsolverTransport)

        end if   ! npreadapt > 0
        
      end if   ! dstepAdapt > 0
      
    else
      
      dstepAdapt = 0.0_DP
      
    end if


    ! Get name of import file (if any)
    call parlst_getvalue_string(rparlist,&
        trim(soutputName), 'ucdimport', ucdimport, '')

    ! Do we have to read in a precomdputed solution?
    if (trim(ucdimport) .ne. '') then
      call ucd_readGMV(ucdimport, rimport, p_rproblemLevel%rtriangulation)
      call ucd_getSimulationTime(rimport, rtimestep%dinitialTime)
      call ucd_getSimulationTime(rimport, rtimestep%dTime)
      call euler_setVariables(rimport, p_rsolutionEuler)
      call transp_setVariable(rimport, 'advect', p_rsolutionTransport)
      call ucd_release(rimport)

      ! Set time for solution output
      dtimeUCD = rtimestep%dinitialTime
    end if

    !---------------------------------------------------------------------------
    ! Calculate velocity field (\rho v) for the scalar model problem
    ! based on the initial solution U^0 of the Euler model
    !---------------------------------------------------------------------------

    call zpinch_initVelocityField(rparlist, ssectionNameTransport,&
        p_rproblemLevel, p_rsolutionEuler, rcollection)
    
    ! Stop time measurement for pre-processing
    call stat_stopTimer(p_rtimerPrePostprocess)


    ! CHECKS
    dmassEuler0 = zpinch_checkConservation(rparlist,&
        ssectionNameEuler, p_rproblemLevel, p_rsolutionEuler, 1)
    
    dmassTransport0 = zpinch_checkConservation(rparlist,&
        ssectionNameTransport, p_rproblemLevel, p_rsolutionTransport, 2)
    ! CHECKS

    !---------------------------------------------------------------------------
    ! Infinite time stepping loop
    !---------------------------------------------------------------------------
    
    timeloop: do
      
      ! Check for user interaction
      if (signal_SIGINT(-1) > 0 )&
          call zpinch_outputSolution(rparlist, ssectionName,&
          p_rproblemLevel, rsolution, rtimestep%dTime)

      !-------------------------------------------------------------------------
      ! Compute Euler model + scalar tracer in coupled fashion
      ! for full time step: $U^n \to U^{n+1}$ and $u^n \to u^{n+1}$
      !-------------------------------------------------------------------------

      ! Start time measurement for solution procedure
      call stat_startTimer(p_rtimerSolution, STAT_TIMERSHORT)
      
      ! Compute scaling for explicit part of the Lorentz force
      dscale = (1.0_DP-rtimestep%theta) * rtimestep%dStep

      if (dscale .ne. 0.0_DP) then
        
        ! Calculate explicit part of the Lorentz force term
        call zpinch_initLorentzforceTerm(rparlist, ssectionName,&
            ssectionNameEuler, ssectionNameTransport, p_rproblemLevel,&
            p_rsolutionEuler, p_rsolutionTransport, rtimestep%dTime, &
            dscale, rforce(1), rcollection)
        
      end if
      
      ! Prepare quick access arrays/vectors
      rcollection%SquickAccess(1) = 'null'
      rcollection%SquickAccess(2) = ssectionName
      rcollection%SquickAccess(3) = ssectionNameEuler
      rcollection%SquickAccess(4) = ssectionNameTransport
      rcollection%p_rvectorQuickAccess1 => rsolution(1)
      rcollection%p_rvectorQuickAccess2 => rsolution(2)
      rcollection%p_rvectorQuickAccess3 => rtimestep%RtempVectors(1)
      rcollection%p_rvectorQuickAccess4 => rtimestep%RtempVectors(2)
      
      ! What time-stepping scheme should be used?
      select case(rtimestep%ctimestepType)
        
      case (TSTEP_RK_SCHEME)
        
!!$        ! Adopt explicit Runge-Kutta scheme
!!$        call tstep_performRKStep(p_rproblemLevel, rtimestep,&
!!$            rsolver, rsolution, zpinch_nlsolverCallbackEuler,&
!!$            rcollection, rforce)
        print *, "Not implemented"
        stop
        
      case (TSTEP_THETA_SCHEME)
        
        if (dscale .gt. 0.0_DP) then

          ! Adopt two-level theta-scheme with source term
          call tstep_performThetaStep(p_rproblemLevel, rtimestep,&
              rsolver, rsolution, zpinch_nlsolverCallback,&
              rcollection, rforce)

        else

          ! Adopt two-level theta-scheme without sorce term
          call tstep_performThetaStep(p_rproblemLevel, rtimestep,&
              rsolver, rsolution, zpinch_nlsolverCallback,&
              rcollection)
          
        end if
        
      case DEFAULT
        call output_line('Unsupported time-stepping algorithm!',&
            OU_CLASS_ERROR,OU_MODE_STD,'zpinch_solveTransientPrimal')
        call sys_halt()
      end select
      
      !-------------------------------------------------------------------------
      ! Compute linearized FCT correction for Euler and transport model
      !-------------------------------------------------------------------------
      
      ! Prepare quick access arrays
      rcollection%SquickAccess(1) = ssectionNameEuler
      rcollection%SquickAccess(2) = ssectionNameTransport
      
      ! Apply linearized FCT correction for Euler model
      call zpinch_calcLinearizedFCT(rbdrCondEuler, rbdrCondTransport,&
          p_rproblemLevel, rtimestep, p_rsolutionEuler,&
          p_rsolutionTransport, rcollection)
      
      ! Calculate velocity field (\rho v)
      call zpinch_initVelocityField(rparlist, ssectionNameTransport,&
          p_rproblemLevel, p_rsolutionEuler, rcollection)
      
      ! Stop time measurement for solution procedure
      call stat_stopTimer(p_rtimerSolution)

      ! CHECKS
      dmassEuler = zpinch_checkConservation(rparlist,&
          ssectionNameEuler, p_rproblemLevel, p_rsolutionEuler, 1)
      
      dmassTransport = zpinch_checkConservation(rparlist,&
          ssectionNameTransport, p_rproblemLevel, p_rsolutionTransport, 2)
      ! CHECKS
      
      print *, "################################################################"
      print *, "MASS CONSERVATION"
      print *, (dmassEuler0-dmassEuler)/dmassEuler0
      print *, (dmassTransport0-dmassTransport)/dmassTransport0
      print *, "################################################################"

      ! Reached final time, then exit the infinite time loop?
      if (rtimestep%dTime .ge. rtimestep%dfinalTime) exit timeloop

      !-------------------------------------------------------------------------
      ! Post-process intermediate solution
      !-------------------------------------------------------------------------
      
      if ((dstepUCD .gt. 0.0_DP) .and.&
          (rtimestep%dTime .ge. dtimeUCD)) then
        
        ! Set time for next intermediate solution export
        dtimeUCD = dtimeUCD + dstepUCD
        
        ! Start time measurement for post-processing
        call stat_startTimer(p_rtimerPrepostProcess, STAT_TIMERSHORT)
        
        ! Export the intermediate solution
        call zpinch_outputSolution(rparlist, ssectionName,&
            p_rproblemLevel, rsolution, rtimestep%dTime)

        ! Stop time measurement for post-processing
        call stat_stopTimer(p_rtimerPrepostProcess)

      end if


      !-------------------------------------------------------------------------
      ! Perform adaptation
      !-------------------------------------------------------------------------

      if ((dstepAdapt .gt. 0.0_DP) .and.&
          (rtimestep%dTime .ge. dtimeAdapt)) then
      
        ! Set time for next adaptation step
        dtimeAdapt = dtimeAdapt + dstepAdapt

        !-----------------------------------------------------------------------
        ! Perform error indication
        !-----------------------------------------------------------------------
        
        ! Start time measurement for error estimation
        call stat_startTimer(p_rtimerErrorEstimation, STAT_TIMERSHORT)

!!$        ! Compute the error estimator using recovery techniques
!!$        call euler_estimateRecoveryError(rparlist, ssectionnameEuler,&
!!$            p_rproblemLevel, p_rsolutionEuler, rtimestep&
!!$            %dinitialTime, relementError, derror)

        ! Compute the error indicator based on the tracer
        call zpinch_calcAdaptationIndicator(p_rsolutionEuler,&
            p_rsolutionTransport, relementError)

        ! Stop time measurement for error estimation
        call stat_stopTimer(p_rtimerErrorEstimation)
        
        
        !-------------------------------------------------------------------------
        ! Perform h-adaptation
        !-------------------------------------------------------------------------
        
        ! Start time measurement for mesh adaptation
        call stat_startTimer(p_rtimerAdaptation, STAT_TIMERSHORT)
        
        ! Set the names of the template matrix
        rcollection%SquickAccess(1) = 'sparsitypattern'
        
        ! Attach the solution vector to the collection structure
        rcollection%p_rvectorQuickAccess1 => rsolution(1)
        rcollection%p_rvectorQuickAccess2 => rsolution(2)
        
        ! Perform h-adaptation and update the triangulation structure
        call zpinch_adaptTriangulation(rparlist, ssectionnameEuler,&
            rhadapt, p_rproblemLevel%rtriangulation, relementError,&
            rcollection)
        
        ! Release element-wise error distribution
        call lsyssc_releaseVector(relementError)

        ! Update the template matrix according to the sparsity pattern
        call grph_generateMatrix(rgraph, p_rproblemLevel&
            %Rmatrix(templateMatrix))

        ! Stop time measurement for mesh adaptation
        call stat_stopTimer(p_rtimerAdaptation)
        
        
        !-------------------------------------------------------------------------
        ! Re-generate the discretization and coefficient matrices
        !-------------------------------------------------------------------------
        
        ! Start time measurement for generation of the triangulation
        call stat_startTimer(p_rtimerTriangulation, STAT_TIMERSHORT)
        
        ! Generate standard mesh from raw mesh
        call tria_initStandardMeshFromRaw(p_rproblemLevel&
            %rtriangulation, rproblem%rboundary)
        
        ! Stop time measurement for generation of the triangulation
        call stat_stopTimer(p_rtimerTriangulation)

        ! Start time measurement for generation of constant
        ! coefficient matrices
        call stat_startTimer(p_rtimerAssemblyCoeff, STAT_TIMERSHORT)
        
        ! Re-initialize all constant coefficient matrices
        call euler_initProblemLevel(rparlist, ssectionNameEuler,&
            p_rproblemLevel, rcollection)
        call transp_initProblemLevel(rparlist, ssectionNameTransport,&
            p_rproblemLevel, rcollection)
        
        ! Resize the solution vector for the Euler model accordingly
        call parlst_getvalue_int(rparlist,&
            ssectionNameEuler, 'systemmatrix', systemMatrix)
        call lsysbl_resizeVecBlockIndMat(p_rproblemLevel&
            %RmatrixBlock(systemMatrix), p_rsolutionEuler, .false., .true.)

        ! Resize the solution vector for the transport model accordingly
        call lsysbl_resizeVectorBlock(p_rsolutionTransport, &
            p_rproblemLevel%Rmatrix(templateMatrix)%NEQ, .false.)
        
        ! Prepare internal data arrays of the solver structure
        call parlst_getvalue_int(rparlist,&
            ssectionNameEuler, 'systemmatrix', systemMatrix)
        call parlst_getvalue_int(rparlist,&
            ssectionNameEuler, 'isystemformat', isystemFormat)
        call flagship_updateSolverMatrix(p_rproblemLevel,&
            p_rsolverEuler, systemMatrix, isystemFormat, UPDMAT_ALL)
        call solver_updateStructure(p_rsolverEuler)
        
        ! Prepare internal data arrays of the solver structure
        call parlst_getvalue_int(rparlist,&
            ssectionNameTransport, 'systemmatrix', systemMatrix)
        call flagship_updateSolverMatrix(p_rproblemLevel,&
            p_rsolverTransport, systemMatrix, SYSTEM_INTERLEAVEFORMAT, UPDMAT_ALL)
        call solver_updateStructure(p_rsolverTransport)

        ! Calculate velocity field (\rho v)
        call zpinch_initVelocityField(rparlist, ssectionNameTransport,&
            p_rproblemLevel, p_rsolutionEuler, rcollection)
        
        ! Stop time measurement for generation of constant
        ! coefficient matrices
        call stat_stopTimer(p_rtimerAssemblyCoeff)
        
      end if
      
    end do timeloop

    ! Release Lorentz force vector
    call lsysbl_releaseVector(rforce(1))

    ! Release adaptation structure
    if (trim(adjustl(sadaptivityName)) .ne. '') then
      if ((dstepAdapt > 0.0_DP) .or. (npreadapt > 0)) then
        call hadapt_releaseAdaptation(rhadapt)
        call grph_releaseGraph(rgraph)
      end if
    end if
    
  end subroutine zpinch_solveTransientPrimal

  !*****************************************************************************
  ! AUXILIARY ROUTINES
  !*****************************************************************************

!<subroutine>

  subroutine zpinch_parseCmdlArguments(rparlist)

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

  end subroutine zpinch_parseCmdlArguments

  !*****************************************************************************

  function zpinch_checkConservation(rparlist, ssectionName,&
      rproblemLevel, rvector, isolution) result(dvalue)
  
    type(t_parlist), intent(in) :: rparlist
    character(LEN=*), intent(in) :: ssectionName
    type(t_problemLevel), intent(in) :: rproblemLevel
    type(t_vectorBlock), intent(in) :: rvector
    integer :: isolution

    real(DP) :: dvalue

    ! local variables
    type(t_vectorScalar) :: rvectorScalar
    real(DP), dimension(:), pointer :: p_Dmatrix, p_Dvector
    integer :: lumpedMassMatrix

    ! Get global configuration from parameter list
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'lumpedmassmatrix', lumpedMassMatrix)

    ! Get lumped mass matrix
    call lsyssc_getbase_double(rproblemLevel&
        %Rmatrix(lumpedMassMatrix), p_Dmatrix)
    
    select case(isolution)

    case(1)

      ! Get density distribution from the solution of the Euler model
      call euler_getVariable(rvector, 'density', rvectorScalar)

      ! Get density distribution
      call lsyssc_getbase_double(rvectorScalar, p_Dvector)
      
      ! Compute total mass
      dvalue = sum(p_Dmatrix * p_Dvector)
      
      ! Release temporal vector
      call lsyssc_releaseVector(rvectorScalar)
      
    case(2)
    
      ! Get solution vector
      call lsysbl_getbase_double(rvector, p_Dvector)

      ! Compute total mass
      dvalue = sum(p_Dmatrix * p_Dvector)

    case default
      print *, "Does not exist"
      stop
    end select

  end function zpinch_checkConservation

end module zpinch_application
