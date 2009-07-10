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
!#     -> The application's main routine called from the main problem
!#
!#
!# The following auxiliary routines are available:
!#
!# 1.) zpinch_parseCmdlArguments
!#     -> Parses the list of commandline arguments and overwrites
!#        parameter values from the parameter files
!#
!# 2.) zpinch_initProblem
!#     -> Initializes the global problem structure based on the
!#        parameter settings given by the parameter list
!#
!# 3.) zpinch_outputSolution
!#     -> Outputs the solution vector to file in UCD format
!#
!# 4.) zpinch_solveTransientPrimal
!#     -> Solves the primal formulation of the time-dependent 
!#        simplified MHD equations
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
    type(t_collection) :: rcollectionEuler, rcollectionTransport

    ! Boundary condition structure for the primal problem
    type(t_boundaryCondition) :: rbdrCondEuler, rbdrCondTransport

    ! Problem structure which holds all internal data (vectors/matrices)
    type(t_problem) :: rproblem

    ! Time-stepping structures
    type(t_timestep) :: rtimestepEuler, rtimestepTransport
    
    ! Global solver structure
    type(t_solver) :: rsolverEuler, rsolverTransport

    ! Solution vectors for the primal problem
    type(t_vectorBlock) :: rsolutionEuler, rsolutionTransport

    ! Timer for the total solution process
    type(t_timer) :: rtimerTotal

    ! Parameter file and section names
    character(LEN=SYS_STRLEN) :: sindatfileName
    character(LEN=SYS_STRLEN) :: sbdrcondName
    character(LEN=SYS_STRLEN) :: algorithm
    character(LEN=SYS_STRLEN) :: ssectionNameEuler
    character(LEN=SYS_STRLEN) :: ssectionNameTransport

    ! local variables
    integer :: systemMatrix
    integer :: nlmin, nlmax


    ! Start total time measurement
    call stat_clearTimer(rtimerTotal)
    call stat_startTimer(rtimerTotal)
    
    !---------------------------------------------------------------------------
    ! Pre-processing
    !---------------------------------------------------------------------------
    
    ! Overwrite global configuration from command line arguments. After
    ! this subroutine has been called, the parameter list remains unchanged.
    call zpinch_parseCmdlArguments(rparlist)
    
    call parlst_getvalue_string(rparlist, 'Zpinch', 'application_euler', ssectionNameEuler)
    call parlst_getvalue_string(rparlist, 'Zpinch', 'application_transport', ssectionNameTransport)

    ! Initialize global collection structures
    call collct_init(rcollectionEuler)
    call collct_init(rcollectionTransport)

    ! Initialize the application descriptors
!!!    call euler_initApplication(rparlist, ssectionNameEuler, rappDescrEuler)
!!!    call transp_initApplication(rparlist, ssectionNameTransport, rappDescrTransport)

    ! Start time measurement for pre-processing
!!!    call stat_startTimer(rappDescrEuler%rtimerPrepostProcess, STAT_TIMERSHORT)
!!$    call stat_startTimer(rappDescrTransport%rtimerPrepostProcess, STAT_TIMERSHORT)
    
    ! Initialize the global collections
!!!    call euler_initCollection(rappDescrEuler, rparlist, ssectionNameEuler, rcollectionEuler)
!!!    call transp_initCollection(rparlist, ssectionNameTransport, rcollectionTransport)

    ! Initialize the solver structures
    call euler_initSolvers(rparlist, ssectionNameEuler, rtimestepEuler, rsolverEuler)
    call transp_initSolvers(rparlist, ssectionNameTransport, rtimestepTransport, rsolverTransport)

    ! Initialize the abstract problem structure
    nlmin = min(solver_getMinimumMultigridlevel(rsolverEuler),&
                solver_getMinimumMultigridlevel(rsolverTransport))
    nlmax = max(solver_getMaximumMultigridlevel(rsolverEuler),&
                solver_getMaximumMultigridlevel(rsolverTransport))

    call zpinch_initProblem(rparlist, 'Zpinch', nlmin, nlmax,&
                            rproblem, rcollectionEuler, rcollectionTransport)

    ! Initialize the individual problem levels
    call euler_initAllProblemLevels(rparlist, 'euler', rproblem, rcollectionEuler)
    call transp_initAllProblemLevels(rparlist, 'transport', rproblem, rcollectionTransport)

    ! Prepare internal data arrays of the solver structure for Euler model
    systemMatrix = collct_getvalue_int(rcollectionEuler, 'systemMatrix') 
!!!    isystemFormat = collct_getvalue_int(rcollectionEuler, 'isystemFormat') 
!!!    call flagship_updateSolverMatrix(rproblem%p_rproblemLevelMax, rsolverEuler,&
!!!                                     systemMatrix, isystemFormat, UPDMAT_ALL)
    call solver_updateStructure(rsolverEuler)

    ! Prepare internal data arrays of the solver structure for transport model
    systemMatrix = collct_getvalue_int(rcollectionTransport, 'systemMatrix') 
    call flagship_updateSolverMatrix(rproblem%p_rproblemLevelMax, rsolverTransport,&
                                     systemMatrix, SYSTEM_INTERLEAVEFORMAT, UPDMAT_ALL)
    call solver_updateStructure(rsolverTransport)
    
    ! Stop time measurement for pre-processing
!!!    call stat_stopTimer(rappDescrEuler%rtimerPrePostprocess)
!!!    call stat_stopTimer(rappDescrTransport%rtimerPrePostprocess)


    !---------------------------------------------------------------------------
    ! Solution algorithm
    !---------------------------------------------------------------------------

    call parlst_getvalue_string(rparlist, 'Zpinch', 'algorithm', algorithm)
    
    ! Initialize the boundary condition for the Euler model
    call parlst_getvalue_string(rparlist, ssectionNameEuler, 'sprimalbdrcondname', sbdrcondName)
    call parlst_getvalue_string(rparlist, ssectionNameEuler, 'indatfile', sindatfileName)
!!$    call bdrf_readBoundaryCondition(rbdrCondEuler, sindatfileName,&
!!$                                    '['//trim(sbdrcondName)//']', rappDescrEuler%ndimension)

    ! Initialize the boundary condition for the scalar transport model
    call parlst_getvalue_string(rparlist, ssectionNameTransport, 'sprimalbdrcondname', sbdrcondName)
    call parlst_getvalue_string(rparlist, ssectionNameTransport, 'indatfile', sindatfileName)
!!$    call bdrf_readBoundaryCondition(rbdrCondTransport, sindatfileName,&
!!$                                    '['//trim(sbdrcondName)//']', rappDescrTransport%ndimension)

    
    ! What solution algorithm should be applied?
    select case(trim(algorithm))

    case ('transient_primal')
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Solve the primal formulation for the time-dependent problem
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!$      call zpinch_solveTransientPrimal(rappDescrEuler, rappDescrTransport, rparlist,&
!!$                                       ssectionNameEuler, ssectionNameTransport,&
!!$                                       rbdrCondEuler, rbdrCondTransport, rproblem,&
!!$                                       rtimestepEuler, rsolverEuler,&
!!$                                       rtimestepTransport, rsolverTransport,&
!!$                                       rsolutionEuler, rsolutionTransport,&
!!$                                       rcollectionEuler, rcollectionTransport)

      call zpinch_outputSolution(rparlist, 'Zpinch', rproblem%p_rproblemLevelMax,&
                                 rsolutionEuler, rsolutionTransport, rtimestepEuler%dTime)
      
    case DEFAULT
      call output_line(trim(algorithm)//' is not a valid solution algorithm!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'zpinch_app')
      call sys_halt()
    end select


    !---------------------------------------------------------------------------
    ! Post-processing
    !---------------------------------------------------------------------------
    
    ! Start time measurement for pre-processing
!!!    call stat_startTimer(rappDescrEuler%rtimerPrepostProcess, STAT_TIMERSHORT)
    
    ! Release time-stepping
    call tstep_releaseTimestep(rtimestepEuler)

    ! Release solver
    call solver_releaseSolver(rsolverEuler)
    
    ! Release boundary conditions
    call bdrf_release(rbdrCondEuler)

    ! Release vectors
    call lsysbl_releaseVector(rsolutionEuler)

    ! Release collection
    call collct_done(rcollectionEuler)
    
    ! Stop time measurement for pre-processing
!!!    call stat_stopTimer(rappDescrEuler%rtimerPrePostprocess)

    !---------------------------------------------------------------------------

    ! Start time measurement for pre-processing
!!!    call stat_startTimer(rappDescrTransport%rtimerPrepostProcess, STAT_TIMERSHORT)
    
    ! Release time-stepping
    call tstep_releaseTimestep(rtimestepTransport)

    ! Release solver
    call solver_releaseSolver(rsolverTransport)
    
    ! Release application descriptor
!!!    call transp_doneApplication(rappDescrTransport)

    ! Release boundary conditions
    call bdrf_release(rbdrCondTransport)

    ! Release vectors
    call lsysbl_releaseVector(rsolutionTransport)

    ! Release collection
    call collct_done(rcollectionTransport)
    
    ! Stop time measurement for pre-processing
!!!!    call stat_stopTimer(rappDescrTransport%rtimerPrePostprocess)

    !---------------------------------------------------------------------------

    ! Release function parser
    call fparser_done()
    
    ! Release problem structure
    call problem_releaseProblem(rproblem)
    
    ! Stop time measurement for total time measurement
    call stat_stopTimer(rtimerTotal)

    ! Output statistics
    call output_lbrk()
    call output_separator(OU_SEP_MINUS)
    call output_line('Compressible Euler model')
!!!!    call euler_outputStatistics(rappDescrEuler, rtimerTotal)
    call output_lbrk()
    call output_separator(OU_SEP_MINUS)
    call output_line('Scalar transport model')       
!!!    call transp_outputStatistics(rappDescrTransport, rtimerTotal)

  end subroutine zpinch_app
  
  !*****************************************************************************

!<subroutine>

  subroutine zpinch_initProblem(rparlist, ssectionName, nlmin, nlmax, rproblem,&
                                rcollectionEuler, rcollectionTransport)

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
    ! collection for the Euler model
    type(t_collection), intent(inout) :: rcollectionEuler

    ! collection for the scalar transport model
    type(t_collection), intent(inout) :: rcollectionTransport
!</intputoutput>

!<output>
    ! problem structure
    type(t_problem), intent(out) :: rproblem
!</output>
!</subroutine>

    ! section names
    character(LEN=SYS_STRLEN) :: ssectionNameEuler
    character(LEN=SYS_STRLEN) :: ssectionNameTransport
    character(LEN=SYS_STRLEN) :: sconvectionName
    character(LEN=SYS_STRLEN) :: sinviscidName

    ! abstract problem descriptor
    type(t_problemDescriptor) :: rproblemDescriptor

    ! pointer to the problem level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! local variables
    integer :: convectionAFC
    integer :: inviscidAFC
    integer :: iconvToTria
    

    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist, ssectionName, 'application_euler', ssectionNameEuler)
    call parlst_getvalue_string(rparlist, ssectionName, 'application_transport', ssectionNameTransport)

    call parlst_getvalue_string(rparlist, ssectionName, 'trifile', rproblemDescriptor%trifile)
    call parlst_getvalue_string(rparlist, ssectionName, 'prmfile', rproblemDescriptor%prmfile, '')
    call parlst_getvalue_int(rparlist, ssectionName, 'ndimension', rproblemDescriptor%ndimension)

    call parlst_getvalue_string(rparlist, ssectionNameTransport, 'convection', sconvectionName)
    call parlst_getvalue_string(rparlist, ssectionNameEuler, 'inviscid', sinviscidName)
    
    ! Set additional problem descriptor
    rproblemDescriptor%ndiscretisation = 2
    rproblemDescriptor%nafcstab        = 2   ! for inviscid and convective stabilization
    rproblemDescriptor%nlmin           = nlmin
    rproblemDescriptor%nlmax           = nlmax
    rproblemDescriptor%nmatrixScalar   = rproblemDescriptor%ndimension + 8
    rproblemDescriptor%nmatrixBlock    = 2
    rproblemDescriptor%nvectorScalar   = 0
    rproblemDescriptor%nvectorBlock    = 1   ! external velocity field

    ! Check if quadrilaterals should be converted to triangles
    call parlst_getvalue_int(rparlist, ssectionName, 'iconvtotria', iconvToTria, 0)
    if (iconvToTria .ne. 0)&
        rproblemDescriptor%iproblemSpec = rproblemDescriptor%iproblemSpec &
                                        + PROBDESC_MSPEC_CONVTRIANGLES   

    ! Initialize problem structure
    call problem_initProblem(rproblemDescriptor, rproblem)
    

    ! Initialize the stabilisation structure
    convectionAFC = collct_getvalue_int(rcollectionTransport, 'convectionAFC')
    inviscidAFC   = collct_getvalue_int(rcollectionEuler, 'inviscidAFC')
        
    ! loop over all problem levels
    p_rproblemLevel => rproblem%p_rproblemLevelMax
    do while(associated(p_rproblemLevel))

      if (convectionAFC > 0) then
        call afcstab_initFromParameterlist(rparlist, sconvectionName,&
                                           p_rproblemLevel%Rafcstab(convectionAFC))
      end if
      
      if (inviscidAFC > 0) then
        call afcstab_initFromParameterlist(rparlist, sinviscidName,&
                                           p_rproblemLevel%Rafcstab(inviscidAFC))
      end if

      ! Switch to next coarser level
      p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
    end do

  end subroutine zpinch_initProblem

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_outputSolution(rparlist, ssectionName, rproblemLevel,&
                                   rsolutionEuler, rsolutionTransport, dtime)

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

    ! solution vector for compressible Euler model
    type(t_vectorBlock), intent(in) :: rsolutionEuler

    ! solution vector for scalar transport model
    type(t_vectorBlock), intent(in) :: rsolutionTransport

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
    call parlst_getvalue_string(rparlist, ssectionName, 'output', soutputName)
    call parlst_getvalue_string(rparlist, trim(soutputName), 'ucdsolution', ucdsolution)
    call parlst_getvalue_int(rparlist, trim(soutputName), 'iformatucd', iformatUCD)

    call parlst_getvalue_string(rparlist, ssectionName, 'application_euler', ssectionNameEuler)
    call parlst_getvalue_int(rparlist, ssectionNameEuler, 'isystemformat', isystemformat)
    call parlst_getvalue_int(rparlist, ssectionNameEuler, 'isystemformat', isystemformat)

    ! Initialize the UCD exporter
    call flagship_initUCDexport(rproblemLevel, ucdsolution,&
                                iformatUCD, rexport, ifilenumber)

    ! Increase filenumber by one
    ifilenumber = ifilenumber+1

    ! Set simulation time
    if (present(dtime)) call ucd_setSimulationTime(rexport, dtime)

    ! Set pointers
    call lsysbl_getbase_double(rsolutionEuler, p_Dsolution)
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

        call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR1D, 'effective_energy', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'effective_energy', UCD_VAR_STANDARD, p_Ddata1)
        
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
        
        call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR2D, 'effective_energy', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'effective_energy', UCD_VAR_STANDARD, p_Ddata1)

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
        
        call euler_getVarInterleaveFormat(rvector1%NEQ, NVAR3D, 'effective_energy', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'effective_energy', UCD_VAR_STANDARD, p_Ddata1)

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

        call euler_getVarBlockFormat(rvector1%NEQ, NVAR1D, 'effective_energy', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'effective_energy', UCD_VAR_STANDARD, p_Ddata1)
        
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
        
        call euler_getVarBlockFormat(rvector1%NEQ, NVAR2D, 'effective_energy', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'effective_energy', UCD_VAR_STANDARD, p_Ddata1)

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

        call euler_getVarBlockFormat(rvector1%NEQ, NVAR3D, 'effective_energy', p_Dsolution, p_Ddata1)
        call ucd_addVariableVertexBased (rexport, 'effective_energy', UCD_VAR_STANDARD, p_Ddata1)

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

    ! Add solution for scalar transport model
    call lsysbl_getbase_double(rsolutionTransport, p_Ddata1)
    call ucd_addVariableVertexBased (rexport, 'advect', UCD_VAR_STANDARD, p_Ddata1)
    
    ! Write UCD file
    call ucd_write  (rexport)
    call ucd_release(rexport)

  end subroutine zpinch_outputSolution

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_adaptTriangulation(rhadapt, rtriangulationSrc, rindicator,&
                                    rcollection, rtriangulationDest)

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

      if (rtriangulationSrc%ndim .eq. NDIM2D) then
        call zpinch_hadaptCallbackScalar2D(rcollection, HADAPT_OPR_INITCALLBACK, Ivalue, Ivalue)
        call hadapt_performAdaptation(rhadapt, rindicator, rcollection, zpinch_hadaptCallbackScalar2D)
      else
        call output_line('Invalid spatial dimension!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'zpinch_adaptTriangulation')
        call sys_halt()
      end if

    case (SYSTEM_BLOCKFORMAT)

      if (rtriangulationSrc%ndim .eq. NDIM2D) then
        call zpinch_hadaptCallbackBlock2D(rcollection, HADAPT_OPR_INITCALLBACK, Ivalue, Ivalue)
        call hadapt_performAdaptation(rhadapt, rindicator, rcollection, zpinch_hadaptCallbackBlock2D)
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

    subroutine zpinch_solveTransientPrimal(rparlist, ssectionNameEuler, ssectionNameTransport,&
                                           rbdrCondEuler, rbdrCondTransport, rproblem,&
                                           rtimestepEuler, rsolverEuler,&
                                           rtimestepTransport, rsolverTransport,&
                                           rsolutionEuler, rsolutionTransport,&
                                           rcollectionEuler, rcollectionTransport)

!<description>
      ! This subroutine solves the transient primal simplified MHD problem.
!</description>

!<input>
    ! section name in parameter list
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
    type(t_timestep), intent(inout) :: rtimestepEuler
    type(t_timestep), intent(inout) :: rtimestepTransport

    ! solver struchture
    type(t_solver), intent(inout), target :: rsolverEuler
    type(t_solver), intent(inout), target :: rsolverTransport

    ! primal solution vector
    type(t_vectorBlock), intent(inout), target :: rsolutionEuler
    type(t_vectorBlock), intent(inout), target :: rsolutionTransport

    ! collection structure
    type(t_collection), intent(inout) :: rcollectionEuler
    type(t_collection), intent(inout) :: rcollectionTransport
!</inputoutput>
!</subroutine>

    ! Pointer to the current multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel
    
    ! Pointer to the discretisation structure
    type(t_blockDiscretisation), pointer :: p_rdiscretisation

    ! Vector for the element-wise feature indicator
    type(t_vectorScalar) :: rindicator

    ! Structure for h-adaptation
    type(t_hadapt) :: rhadapt

    ! Structure for the sparsity pattern
    type(t_graph) :: rgraph

    ! Collection structure for h-adaptation
    type(t_collection) :: rcollection

    ! section names
    character(LEN=SYS_STRLEN) :: sadaptivityName
    character(LEN=SYS_STRLEN) :: soutputName

    ! local variables
    real(dp) :: derror, dstepUCD, dtimeUCD, dstepAdapt, dtimeAdapt
    integer :: templateMatrix, systemMatrix, isystemFormat
    integer :: discretisationEuler, discretisationTransport
    integer :: isize, ipreadapt, npreadapt, nerrorvariable
    integer, external :: signal_SIGINT

    
    ! Set pointer to maximum problem level
    p_rproblemLevel => rproblem%p_rproblemLevelMax


    !--- compressible Euler model ----------------------------------------------

    ! Start time measurement for pre-processing
!!$    call stat_startTimer(rappDescrEuler%rtimerPrePostprocess, STAT_TIMERSHORT)

    ! Get position of discretisation structure
    call parlst_getvalue_int(rparlist, ssectionNameEuler, 'discretisation', discretisationEuler)

    ! Set pointer to discretisation structure
    p_rdiscretisation => p_rproblemLevel%Rdiscretisation(discretisationEuler)

    ! Create the solution vector
    call lsysbl_createVectorBlock(p_rdiscretisation, rsolutionEuler, .false., ST_DOUBLE)
    if (p_rdiscretisation%ncomponents .ne. euler_getNVAR(p_rproblemLevel)) then
      rsolutionEuler%RvectorBlock(1)%NVAR = euler_getNVAR(p_rproblemLevel)
      isize = rsolutionEuler%NEQ*euler_getNVAR(p_rproblemLevel)
      call lsysbl_resizeVectorBlock(rsolutionEuler, isize, .false., .false.)
    end if

    ! Initialize the solution vector and impose boundary conditions explicitly
    call euler_initSolution(rparlist, ssectionNameEuler, p_rproblemLevel,&
                            rtimestepEuler%dinitialTime, rsolutionEuler, rcollectionEuler)
!!$    select case(rappDescrEuler%ndimension)
!!$    case (NDIM1D)
!!$      call bdrf_filterVectorExplicit(rbdrCondEuler, p_rproblemLevel%rtriangulation,&
!!$                                     rsolutionEuler, rtimestepEuler%dinitialTime,&
!!$                                     rproblem%rboundary, euler_calcBoundaryvalues1d)
!!$    case (NDIM2D)
!!$      call bdrf_filterVectorExplicit(rbdrCondEuler, p_rproblemLevel%rtriangulation,&
!!$                                     rsolutionEuler, rtimestepEuler%dinitialTime,&
!!$                                     rproblem%rboundary, euler_calcBoundaryvalues2d)
!!$    case (NDIM3D)
!!$      call bdrf_filterVectorExplicit(rbdrCondEuler, p_rproblemLevel%rtriangulation,&
!!$                                     rsolutionEuler, rtimestepEuler%dinitialTime,&
!!$                                     rproblem%rboundary, euler_calcBoundaryvalues3d)
!!$    end select
    
    ! Initialize timer for intermediate UCD exporter
    dtimeUCD = rtimestepEuler%dinitialTime
    call parlst_getvalue_string(rparlist, 'Zpinch', 'output', soutputName)
    call parlst_getvalue_double(rparlist, trim(soutputName), 'dstepUCD', dstepUCD, 0.0_DP)

    ! Attach the boundary condition
    call solver_setBoundaryCondition(rsolverEuler, rbdrCondEuler, .true.)
    call problem_setBoundaryCondition(rproblem, rbdrCondEuler)

    ! Set collection to primal problem mode
    call collct_setvalue_int(rcollectionEuler, 'primaldual', 1, .true.)
    
    ! Stop time measurement for pre-processing
!!!!    call stat_stopTimer(rappDescrEuler%rtimerPrePostprocess)
    

    !--- scalar transport model ------------------------------------------------

    ! Start time measurement for pre-processing
!!!!    call stat_startTimer(rappDescrTransport%rtimerPrePostprocess, STAT_TIMERSHORT)

    ! Get position of discretisation structure
    call parlst_getvalue_int(rparlist, ssectionNameTransport, 'discretisation', discretisationTransport)

    ! Set pointer to discretisation structure
    p_rdiscretisation => p_rproblemLevel%Rdiscretisation(discretisationTransport)

    ! Initialize the solution vector and impose boundary conditions explicitly
    call lsysbl_createVectorBlock(p_rdiscretisation, rsolutionTransport, .false., ST_DOUBLE)
!!$    call transp_initSolution(rappDescrTransport, rparlist, ssectionNameTransport, p_rproblemLevel,&
!!$                             rtimestepTransport%dinitialTime, rsolutionTransport)
    call bdrf_filterVectorExplicit(rbdrCondTransport, p_rproblemLevel%rtriangulation,&
                                   rsolutionTransport, rtimestepTransport%dinitialTime)
    
    ! Attach the boundary condition
    call solver_setBoundaryCondition(rsolverTransport, rbdrCondTransport, .true.)
    call problem_setBoundaryCondition(rproblem, rbdrCondTransport)

    ! Set collection to primal problem mode
    call collct_setvalue_int(rcollectionTransport, 'primaldual', 1, .true.)
    
    ! Stop time measurement for pre-processing
!!!    call stat_stopTimer(rappDescrTransport%rtimerPrePostprocess)

    
    !---------------------------------------------------------------------------
    ! Initialize the h-adaptation structure and perform pre-adaptation
    !---------------------------------------------------------------------------

    call parlst_getvalue_string(rparlist, ssectionNameTransport, 'adaptivity', sadaptivityName, '')
    if (trim(adjustl(sadaptivityName)) .ne. '') then

      call parlst_getvalue_double(rparlist, trim(sadaptivityName), 'dstepAdapt', dstepAdapt)
      call parlst_getvalue_double(rparlist, trim(sadaptivityName), 'dtimeAdapt', dtimeAdapt)
      call parlst_getvalue_int(rparlist, trim(sadaptivityName), 'npreadapt', npreadapt)

      if ((dstepAdapt > 0.0_DP) .or. (npreadapt > 0)) then

        ! Initialize collection structure for the adaptation process
        call collct_init(rcollection)

        ! Initialize adaptation structure from parameter list
        call hadapt_initFromParameterlist(rhadapt, rparlist, sadaptivityName)

        ! Generate a dynamic graph for the sparsity pattern and attach
        ! it to the collection structure which is used to pass this
        ! type to the callback function for h-adaptation. Note that
        ! the template matrix is the same for the compressible Euler
        ! model and the scalar transport model. Therefore, the
        ! template matrix is adopted from the transport Euler model.
        templateMatrix = collct_getvalue_int(rcollectionTransport, 'templateMatrix')
        call grph_createGraphFromMatrix(p_rproblemLevel%Rmatrix(templateMatrix), rgraph)
        call collct_setvalue_graph(rcollection, 'sparsitypattern', rgraph, .true.)
        
        ! Set the system format of the Euler model in the common collection structire
        isystemFormat = collct_getvalue_int(rcollectionEuler, 'isystemFormat')
        call collct_setvalue_int(rcollection, 'isystemFormat', isystemFormat, .true.)

        ! Attach the primal solution vector from the Euler model to the collection
        call collct_setvalue_vec(rcollection, 'solutionvectorEuler', rsolutionEuler, .true.)

        ! Attach the primal solution vector from the scalar transport model to the collection
        call collct_setvalue_vec(rcollection, 'solutionvectorTransport', rsolutionTransport, .true.)

        
        ! Perform pre-adaptation?
        if (npreadapt > 0) then

          ! Set the names of the template matrix and the solution vector
          rcollection%SquickAccess(1) = 'sparsitypattern'
          rcollection%SquickAccess(2) = 'solutionvectorEuler'
          rcollection%SquickAccess(3) = 'solutionvectorTransport'
          
          ! Attach the primal solution vector to the collection structure
          call collct_setvalue_vec(rcollection, 'solutionvectorEuler', rsolutionEuler, .true.)
          call collct_setvalue_vec(rcollection, 'solutionvectorTransport', rsolutionTransport, .true.)
         

          ! Perform number of pre-adaptation steps
          do ipreadapt = 1, npreadapt
            
!!$            ! Compute the error estimator using recovery techniques
!!$            call euler_estimateRecoveryError(rparlist, ssectionnameEuler, p_rproblemLevel,&
!!$                                             rsolutionEuler, rtimestepEuler%dinitialTime,&
!!$                                             rindicator, derror)

            ! Compute the error estimator based on the tracer
            call zpinch_calcAdaptationIndicator(rsolutionEuler, rsolutionTransport, rindicator)

            ! Perform h-adaptation and update the triangulation structure
            call zpinch_adaptTriangulation(rhadapt, p_rproblemLevel%rtriangulation,&
                                           rindicator, rcollection)
            
            ! Release element-wise error distribution
            call lsyssc_releaseVector(rindicator)

            ! Generate standard mesh from raw mesh
            call tria_initStandardMeshFromRaw(p_rproblemLevel%rtriangulation, rproblem%rboundary)
            
            ! Update the template matrix according to the sparsity pattern
            templateMatrix = collct_getvalue_int(rcollectionTransport, 'templateMatrix')
            call grph_generateMatrix(rgraph, p_rproblemLevel%Rmatrix(templateMatrix))

            ! Re-initialize all constant coefficient matrices
!!!!            call euler_initProblemLevel(rappDescrEuler, p_rproblemLevel, rcollectionEuler)
            call transp_initProblemLevel(rparlist, ssectionNameTransport, p_rproblemLevel, rcollectionTransport)

            ! Resize the solution vector for the Euler model accordingly
            systemMatrix = collct_getvalue_int(rcollectionEuler, 'systemMatrix')
            call lsysbl_resizeVecBlockIndMat(p_rproblemLevel%RmatrixBlock(systemMatrix),&
                                             rsolutionEuler, .false., .true.)
            
            ! Resize the solution vector for the scalar transport model accordingly
            call lsysbl_resizeVectorBlock(rsolutionTransport, &
                p_rproblemLevel%Rmatrix(templateMatrix)%NEQ, .false.)


            ! Re-generate the initial solution vectors
            call euler_initSolution(rparlist, ssectionnameEuler, p_rproblemLevel,&
                                    rtimestepEuler%dinitialTime, rsolutionEuler, rcollectionEuler)
!!!!            call transp_initSolution(rappDescrTransport, rparlist, ssectionnameTransport, p_rproblemLevel,&
!!!!                                     rtimestepTransport%dinitialTime, rsolutionTransport)
            
            

            ! Re-generate the initial solution vector and impose boundary conditions explicitly
!!$            select case(rappDescrEuler%ndimension)
!!$            case (NDIM1D)
!!$              call bdrf_filterVectorExplicit(rbdrCondEuler, p_rproblemLevel%rtriangulation,&
!!$                                             rsolutionEuler, rtimestepEuler%dinitialTime,&
!!$                                             rproblem%rboundary, euler_calcBoundaryvalues1d)
!!$              
!!$            case (NDIM2D)
!!$              call bdrf_filterVectorExplicit(rbdrCondEuler, p_rproblemLevel%rtriangulation,&
!!$                                             rsolutionEuler, rtimestepEuler%dinitialTime,&
!!$                                             rproblem%rboundary, euler_calcBoundaryvalues2d)
!!$            case (NDIM3D)
!!$              call bdrf_filterVectorExplicit(rbdrCondEuler, p_rproblemLevel%rtriangulation,&
!!$                                             rsolutionEuler, rtimestepEuler%dinitialTime,&
!!$                                             rproblem%rboundary, euler_calcBoundaryvalues3d)
!!$            end select
            call bdrf_filterVectorExplicit(rbdrCondTransport, p_rproblemLevel%rtriangulation,&
                                             rsolutionTransport, rtimestepTransport%dinitialTime)
          end do

          ! Prepare internal data arrays of the solver structure
          systemMatrix = collct_getvalue_int(rcollectionEuler, 'systemMatrix')
!!$          call flagship_updateSolverMatrix(p_rproblemLevel, rsolverEuler, systemMatrix,&
!!$                                           rappDescrEuler%isystemFormat, UPDMAT_ALL)
          call solver_updateStructure(rsolverEuler)

          ! Prepare internal data arrays of the solver structure
          systemMatrix = collct_getvalue_int(rcollectionTransport, 'systemMatrix')
          call flagship_updateSolverMatrix(p_rproblemLevel, rsolverTransport, systemMatrix,&
                                           SYSTEM_INTERLEAVEFORMAT, UPDMAT_ALL)
          call solver_updateStructure(rsolverTransport)

        end if   ! npreadapt > 0
        
      end if   ! dstepAdapt > 0
      
    else
      
      dstepAdapt = 0.0_DP
      
    end if


    !---------------------------------------------------------------------------
    ! Infinite time stepping loop
    !---------------------------------------------------------------------------
    
    timeloop: do
      
      ! Check for user interaction
      if (signal_SIGINT(-1) > 0 )&
          call zpinch_outputSolution(rparlist, 'Zpinch', p_rproblemLevel,&
                                     rsolutionEuler, rsolutionTransport, rtimestepEuler%dTime)

      !-------------------------------------------------------------------------
      ! Compute Euler model for full time step: U^n -> U^{n+1}
      !-------------------------------------------------------------------------

      ! Start time measurement for solution procedure
!!$      call stat_startTimer(rappDescrEuler%rtimerSolution, STAT_TIMERSHORT)
      
      ! What time-stepping scheme should be used?
      select case(rtimestepEuler%ctimestepType)
        
      case (TSTEP_RK_SCHEME)
        
        ! Adopt explicit Runge-Kutta scheme
        call tstep_performRKStep(p_rproblemLevel, rtimestepEuler, rsolverEuler,&
                                 rsolutionEuler, euler_nlsolverCallback, rcollectionEuler)
        
      case (TSTEP_THETA_SCHEME)
        
        ! Adopt two-level theta-scheme
        call tstep_performThetaStep(p_rproblemLevel, rtimestepEuler, rsolverEuler,&
                                    rsolutionEuler, euler_nlsolverCallback, rcollectionEuler)
        
      case DEFAULT
        call output_line('Unsupported time-stepping algorithm!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'zpinch_solveTransientPrimal')
        call sys_halt()
      end select

      ! Stop time measurement for solution procedure
!!$      call stat_stopTimer(rappDescrEuler%rtimerSolution)

      
      !-------------------------------------------------------------------------
      ! Compute scalar transport model for full time step: u^n -> u^{n+1}
      !-------------------------------------------------------------------------
      
      ! Start time measurement for solution procedure
!!$      call stat_startTimer(rappDescrTransport%rtimerSolution, STAT_TIMERSHORT)

      ! Set velocity field v^{n+1} for scalar model problem
      call zpinch_calcVelocityField(p_rproblemLevel, rsolutionEuler, rcollectionTransport)

      ! What time-stepping scheme should be used?
      select case(rtimestepTransport%ctimestepType)
        
      case (TSTEP_RK_SCHEME)
        
        ! Adopt explicit Runge-Kutta scheme
        call tstep_performRKStep(p_rproblemLevel, rtimestepTransport, rsolverTransport,&
                                 rsolutionTransport, transp_nlsolverCallback, rcollectionTransport)
        
      case (TSTEP_THETA_SCHEME)
        
        ! Adopt two-level theta-scheme
        call tstep_performThetaStep(p_rproblemLevel, rtimestepTransport, rsolverTransport,&
                                    rsolutionTransport, transp_nlsolverCallback, rcollectionTransport)
          
      case DEFAULT
        call output_line('Unsupported time-stepping algorithm!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'zpinch_solveTransientPrimal')
        call sys_halt()
      end select
      
      ! Stop time measurement for solution procedure
!!$      call stat_stopTimer(rappDescrTransport%rtimerSolution)

      
      ! Perform conservative FCT postprocessing
      call zpinch_calcLinearizedFCT(rbdrCondEuler, rbdrCondTransport, p_rproblemLevel,&
                                    rtimestepEuler, rsolutionEuler, rsolutionTransport, rcollectionEuler)
      

      !-------------------------------------------------------------------------
      ! Compute source term for full time step
      !
      ! U^{n+1} - \tilde U^{n+1} = dt * (S^{n+1} - S^n)
      !-------------------------------------------------------------------------

      ! Start time measurement for solution procedure
!!$      call stat_startTimer(rappDescrEuler%rtimerSolution, STAT_TIMERSHORT)

      call zpinch_calcSourceTerm(p_rproblemLevel, rtimestepTransport,&
                                 rsolutionTransport, rsolutionEuler, rcollectionEuler)

!!$      ! Stop time measurement for solution procedure
!!$      call stat_stopTimer(rappDescrEuler%rtimerSolution)
      

      ! Reached final time, then exit the infinite time loop?
      if ((rtimestepEuler%dTime     .ge. rtimestepEuler%dfinalTime) .or.&
          (rtimestepTransport%dTime .ge. rtimestepTransport%dfinalTime)) exit timeloop


      !-------------------------------------------------------------------------
      ! Post-process intermediate solution
      !-------------------------------------------------------------------------
      
      if ((dstepUCD .gt. 0.0_DP) .and. (rtimestepEuler%dTime .ge. dtimeUCD .or.&
                                        rtimestepTransport%dTime .ge. dtimeUCD)) then
        
        ! Set time for next intermediate solution export
        dtimeUCD = dtimeUCD + dstepUCD
        
!!$        ! Start time measurement for post-processing
!!$        call stat_startTimer(rappDescrEuler%rtimerPrepostProcess, STAT_TIMERSHORT)
        
        ! Export the intermediate solution
        call zpinch_outputSolution(rparlist, 'Zpinch', p_rproblemLevel,&
                                   rsolutionEuler, rsolutionTransport, rtimestepEuler%dTime)        

!!$        ! Stop time measurement for post-processing
!!$        call stat_stopTimer(rappDescrEuler%rtimerPrepostProcess)
        
      end if


      !-------------------------------------------------------------------------
      ! Perform adaptation
      !-------------------------------------------------------------------------

      if ((dstepAdapt .gt. 0.0_DP) .and. (rtimestepTransport%dTime .ge. dtimeAdapt)) then
      
        ! Set time for next adaptation step
        dtimeAdapt = dtimeAdapt + dstepAdapt

        !-----------------------------------------------------------------------
        ! Perform error indication
        !-----------------------------------------------------------------------
        
        ! Start time measurement for error estimation
!!$        call stat_startTimer(rappDescrTransport%rtimerErrorEstimation, STAT_TIMERSHORT)

!!! THIS WAS COMMENTED OUT
!!$        ! Compute the error estimator using recovery techniques
!!$        call euler_estimateRecoveryError(rparlist, ssectionnameEuler, p_rproblemLevel,&
!!$                                         rsolutionEuler, rtimestepEuler%dinitialTime,&
!!$                                         rindicator, derror)
!!! THIS WAS COMMENTED OUT

        ! Compute the error indicator based on the tracer
        call zpinch_calcAdaptationIndicator(rsolutionEuler, rsolutionTransport, rindicator)
        
        ! Stop time measurement for error estimation
!!!!        call stat_stopTimer(rappDescrTransport%rtimerErrorEstimation)
        
        !-------------------------------------------------------------------------
        ! Perform h-adaptation
        !-------------------------------------------------------------------------
        
        ! Start time measurement for mesh adaptation
!!!!        call stat_startTimer(rappDescrTransport%rtimerAdaptation, STAT_TIMERSHORT)
        
        ! Set the names of the template matrix and the solution vector
        rcollection%SquickAccess(1) = 'sparsitypattern'
        rcollection%SquickAccess(2) = 'solutionvectorEuler'
        rcollection%SquickAccess(3) = 'solutionvectorTransport'
        
        ! Attach the primal solution vector to the collection structure
        ! Attach the primal solution vector to the collection structure
          call collct_setvalue_vec(rcollection, 'solutionvectorEuler', rsolutionEuler, .true.)
          call collct_setvalue_vec(rcollection, 'solutionvectorTransport', rsolutionTransport, .true.)
        
        ! Perform h-adaptation and update the triangulation structure
        call zpinch_adaptTriangulation(rhadapt, p_rproblemLevel%rtriangulation,&
                                       rindicator, rcollection)
        
        ! Release element-wise error distribution
        call lsyssc_releaseVector(rindicator)

        ! Update the template matrix according to the sparsity pattern
        call grph_generateMatrix(rgraph, p_rproblemLevel%Rmatrix(templateMatrix))
        
        ! Stop time measurement for mesh adaptation
!!!!        call stat_stopTimer(rappDescrTransport%rtimerAdaptation)


        !-------------------------------------------------------------------------
        ! Re-generate the discretization and coefficient matrices
        !-------------------------------------------------------------------------
        
        ! Start time measurement for generation of the triangulation
!!!!        call stat_startTimer(rappDescrTransport%rtimerTriangulation, STAT_TIMERSHORT)
        
        ! Generate standard mesh from raw mesh
        call tria_initStandardMeshFromRaw(p_rproblemLevel%rtriangulation, rproblem%rboundary)
        
        ! Stop time measurement for generation of the triangulation
!!!!        call stat_stopTimer(rappDescrTransport%rtimerTriangulation)
        
        
        ! Start time measurement for generation of constant coefficient matrices
!!!!        call stat_startTimer(rappDescrTransport%rtimerAssemblyCoeff, STAT_TIMERSHORT)
        
        ! Re-initialize all constant coefficient matrices
!!!!        call euler_initProblemLevel(rappDescrEuler, p_rproblemLevel, rcollectionEuler)
        call transp_initProblemLevel(rparlist, ssectionNameTransport, p_rproblemLevel, rcollectionTransport)
        
        ! Resize the solution vector accordingly
        systemMatrix = collct_getvalue_int(rcollectionEuler, 'systemMatrix')
        call lsysbl_resizeVecBlockIndMat(p_rproblemLevel%RmatrixBlock(systemMatrix),&
                                         rsolutionEuler, .false., .true.)

        ! Resize the solution vector for the scalar transport model accordingly
        call lsysbl_resizeVectorBlock(rsolutionTransport, &
            p_rproblemLevel%Rmatrix(templateMatrix)%NEQ, .false.)
        
        ! Prepare internal data arrays of the solver structure
        systemMatrix = collct_getvalue_int(rcollectionEuler, 'systemMatrix')
!!!!        call flagship_updateSolverMatrix(p_rproblemLevel, rsolverEuler, systemMatrix,&
!!!!                                         rappDescrEuler%isystemFormat, UPDMAT_ALL)
        call solver_updateStructure(rsolverEuler)

        systemMatrix = collct_getvalue_int(rcollectionTransport, 'systemMatrix')
        call flagship_updateSolverMatrix(p_rproblemLevel, rsolverTransport, systemMatrix,&
                                         SYSTEM_INTERLEAVEFORMAT, UPDMAT_ALL)
        call solver_updateStructure(rsolverTransport)
        
        ! Stop time measurement for generation of constant coefficient matrices
!!!!        call stat_stopTimer(rappDescrTransport%rtimerAssemblyCoeff)

      end if

    end do timeloop

    ! Release adaptation structure
    if ((dstepAdapt > 0.0_DP) .or. (npreadapt > 0)) then
      call hadapt_releaseAdaptation(rhadapt)
      call grph_releaseGraph(rgraph)
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
  
end module zpinch_application
