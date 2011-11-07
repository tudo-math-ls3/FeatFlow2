!##############################################################################
!# ****************************************************************************
!# <name> zpinch_application </name>
!# ****************************************************************************
!#
!# <purpose>
!# This application solves the time-dependent magnetohydrodynamic equations
!# in the one-, two- or three-dimensional domain $\Omega$.
!#
!# The spatial discretisation is perform by means of the algebraic
!# flux correction (AFC) paradigm by Kuzmin, Moeller and Turek. In
!# particular, high-resolution finite element schemes of TVD- and
!# FCT-type are available. For the temporal discretisation, the
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
!#        parameter list which needs to be initialised and filled in
!#        the main program. It then works black-box, that is, it
!#        determines the solution algorithm to be used and performs
!#        the simulation. The user should only have to modify this
!#        routine if another solution algorithm is implemented.
!#
!# 2.) zpinch_initSolvers
!#     -> Initializes the solve structures from the parameter list.
!#
!# 3.) zpinch_initProblemDescriptor
!#     -> Initializes the abstract problem descriptor based on the
!#        parameter settings given by the parameter list.
!#
!# 4.) zpinch_calcAdaptationIndicator
!#     -> Calculates the element-wise indicator for refinement/-coarsening
!#        based on the detectation of shocks and contact discontinuities
!#
!# 5.) zpinch_outputSolution
!#     -> Outputs the solution vector to file in UCD format
!#
!# 6.) zpinch_adaptTriangulation
!#      -> Performs h-adaptation for the given triangulation
!#
!# 7.) zpinch_solveTransientPrimal
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

#include "hydro.h"

  use afcstabbase
  use boundarycondaux
  use boundaryfilter
  use collection
  use flagship_basic
  use fparser
  use genoutput
  use graph
  use groupfembase
  use hadaptaux
  use hadaptivity
  use hydro_application
  use hydro_basic
  use hydro_basic1d
  use hydro_basic2d
  use hydro_basic3d
  use hydro_callback
  use hydro_callback1d
  use hydro_callback2d
  use hydro_callback3d
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

  subroutine zpinch_app(rparlist, ssectionName)

!<description>
    ! This is the main application for the simplified MHD equations. It
    ! is a so-called driver routine which can be used to start a
    ! standalone MHD simulation.
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

    ! Boundary condition structures for the primal problem and pointers
    type(t_boundaryCondition), dimension(2), target :: RbdrCond
    type(t_boundaryCondition), pointer :: p_rbdrCondHydro, p_rbdrCondTransport

    ! Problem structure which holds all internal data (vectors/matrices)
    type(t_problem) :: rproblem

    ! Time-stepping structures
    type(t_timestep) :: rtimestep

    ! Global solver structure and point
    type(t_solver) :: rsolver
    type(t_solver), pointer :: p_rsolver => null()

    ! Solution vectors and pointers
    type(t_vectorBlock), dimension(2), target :: Rsolution
    type(t_vectorBlock), pointer :: p_rsolutionHydro, p_rsolutionTransport

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
    character(LEN=SYS_STRLEN) :: ssectionNameHydro
    character(LEN=SYS_STRLEN) :: ssectionNameTransport

    ! local variables
    integer :: isystemFormat, systemMatrix, ndimension

    ! Start total time measurement
    call stat_startTimer(rtimerTotal)

    !---------------------------------------------------------------------------
    ! Pre-processing
    !---------------------------------------------------------------------------

    ! Start time measurement
    call stat_startTimer(rtimerPrePostprocess, STAT_TIMERSHORT)

    ! Set pointer
    p_rbdrCondHydro => RbdrCond(1)
    p_rsolutionHydro => Rsolution(1)

    p_rbdrCondTransport => RbdrCond(2)
    p_rsolutionTransport => Rsolution(2)

    ! Retrieve section names of sub-applications
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'subapplication', ssectionNameHydro, isubstring=1)
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'subapplication', ssectionNameTransport, isubstring=2)

    ! Overwrite configuration from command line arguments. After this
    ! subroutine has been called, the parameter list remains unchanged
    ! unless the used updates some parameter values interactively.
    call zpinch_parseCmdlArguments(rparlist)


    ! Initialize global collection structure
    call collct_init(rcollection)

    ! Create separate sections for the z-pinch problem, the
    ! hydrodynamic model and the scalar transport model
    call collct_addsection(rcollection, ssectionName)
    call collct_addsection(rcollection, ssectionNameHydro)
    call collct_addsection(rcollection, ssectionNameTransport)

    ! Define section names of the z-pinch problem, the hydrodynamic
    ! model and the scalar transport model. Below, we will fill each
    ! of the three section with parameter list and the timer
    ! structures. Why is this complicated task necessary? Well, each
    ! submodel assumes that it has its own section in the global
    ! collection structure and the global parameter list.
    call collct_setvalue_string(rcollection, 'ssectionName',&
        ssectionName, .true.)
    call collct_setvalue_string(rcollection, 'ssectionNameHydro',&
        ssectionNameHydro, .true.)
    call collct_setvalue_string(rcollection, 'ssectionNameTransport',&
        ssectionNameTransport, .true.)


    ! Attach the parameter list and the timers to the
    ! collection. Since we do not measure the time individually for
    ! each submodel, the same timers will be attached to the section
    ! that corresponds to the hydrodynamic model and the scalar
    ! transport model so that the timings are cumulated.
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

    ! Attach the parameter list and the timers to the collection
    ! into the section for the hydrodynamic model
    call collct_setvalue_parlst(rcollection,&
        'rparlist', rparlist, .true.,&
        ssectionName=ssectionNameHydro)
    call collct_setvalue_timer(rcollection,&
        'rtimerSolution', rtimerSolution, .true.,&
        ssectionName=ssectionNameHydro)
    call collct_setvalue_timer(rcollection,&
        'rtimerAdaptation', rtimerAdaptation, .true.,&
        ssectionName=ssectionNameHydro)
    call collct_setvalue_timer(rcollection,&
        'rtimerErrorEstimation', rtimerErrorEstimation, .true.,&
        ssectionName=ssectionNameHydro)
    call collct_setvalue_timer(rcollection,&
        'rtimerTriangulation', rtimerTriangulation, .true.,&
        ssectionName=ssectionNameHydro)
    call collct_setvalue_timer(rcollection,&
        'rtimerAssemblyCoeff', rtimerAssemblyCoeff, .true.,&
        ssectionName=ssectionNameHydro)
    call collct_setvalue_timer(rcollection,&
        'rtimerAssemblyMatrix', rtimerAssemblyMatrix, .true.,&
        ssectionName=ssectionNameHydro)
    call collct_setvalue_timer(rcollection,&
        'rtimerAssemblyVector', rtimerAssemblyVector, .true.,&
        ssectionName=ssectionNameHydro)
    call collct_setvalue_timer(rcollection,&
        'rtimerPrePostprocess', rtimerPrePostprocess, .true.,&
        ssectionName=ssectionNameHydro)

    ! Attach the parameter list and the timers to the collection
    ! into the section for the scalar transport model
    call collct_setvalue_parlst(rcollection,&
        'rparlist', rparlist, .true.,&
        ssectionName=ssectionNameTransport)
    call collct_setvalue_timer(rcollection,&
        'rtimerSolution', rtimerSolution, .true.,&
        ssectionName=ssectionNameTransport)
    call collct_setvalue_timer(rcollection,&
        'rtimerAdaptation', rtimerAdaptation, .true.,&
        ssectionName=ssectionNameTransport)
    call collct_setvalue_timer(rcollection,&
        'rtimerErrorEstimation', rtimerErrorEstimation, .true.,&
        ssectionName=ssectionNameTransport)
    call collct_setvalue_timer(rcollection,&
        'rtimerTriangulation', rtimerTriangulation, .true.,&
        ssectionName=ssectionNameTransport)
    call collct_setvalue_timer(rcollection,&
        'rtimerAssemblyCoeff', rtimerAssemblyCoeff, .true.,&
        ssectionName=ssectionNameTransport)
    call collct_setvalue_timer(rcollection,&
        'rtimerAssemblyMatrix', rtimerAssemblyMatrix, .true.,&
        ssectionName=ssectionNameTransport)
    call collct_setvalue_timer(rcollection,&
        'rtimerAssemblyVector', rtimerAssemblyVector, .true.,&
        ssectionName=ssectionNameTransport)
    call collct_setvalue_timer(rcollection,&
        'rtimerPrePostprocess', rtimerPrePostprocess, .true.,&
        ssectionName=ssectionNameTransport)

    ! Create function parser
    call fparser_create(rfparser, 100)

    ! Read in all constants, predefined expressions
    ! and functions from the parameter files
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'indatfile', sindatfileName)
    call fparser_parseFileForKeyword(rfparser,&
        sindatfileName, 'defconst', FPAR_CONSTANT)
    call fparser_parseFileForKeyword(rfparser,&
        sindatfileName, 'defexpr', FPAR_EXPRESSION)
    call fparser_parseFileForKeyword(rfparser,&
        sindatfileName, 'deffunc', FPAR_FUNCTION)

    call parlst_getvalue_string(rparlist,&
        ssectionNameHydro, 'indatfile', sindatfileName)
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

    ! Attach the function parser to all three sections of the collection
    call collct_setvalue_pars(rcollection, 'rfparser', rfparser, .true.,&
        ssectionName=ssectionName)
    call collct_setvalue_pars(rcollection, 'rfparser', rfparser, .true.,&
        ssectionName=ssectionNameHydro)
    call collct_setvalue_pars(rcollection, 'rfparser', rfparser, .true.,&
        ssectionName=ssectionNameTransport)

    ! Initialize the solver structures
    call zpinch_initSolvers(rparlist, ssectionName, rtimestep, rsolver)

    ! Get the spatial dimension from the parameter list
    call parlst_getvalue_int(rparlist, ssectionName,&
        'ndimension', ndimension)
    
    ! Initialize the boundary condition for the hydrodynamic model
    call parlst_getvalue_string(rparlist, ssectionNameHydro,&
        'sprimalbdrcondname', sbdrcondName)
    call parlst_getvalue_string(rparlist, ssectionNameHydro,&
        'indatfile', sindatfileName)
    
    ! The boundary condition for the primal problem is required for
    ! all solution strategies so initialize it from the parameter file
    call bdrc_readBoundaryCondition(p_rbdrCondHydro,&
        sindatfileName, '['//trim(sbdrcondName)//']',&
        ndimension, hydro_parseBoundaryCondition)
    
    ! Initialize the boundary condition for the transport model
    call parlst_getvalue_string(rparlist, ssectionNameTransport,&
        'sprimalbdrcondname', sbdrcondName)
    call parlst_getvalue_string(rparlist, ssectionNameTransport,&
        'indatfile', sindatfileName)
    
    ! The boundary condition for the primal problem is required for
    ! all solution strategies so initialize it from the parameter file
    call bdrc_readBoundaryCondition(p_rbdrCondTransport,&
        sindatfileName, '['//trim(sbdrcondName)//']',&
        ndimension, transp_parseBoundaryCondition)
    
    ! Initialize the abstract problem structure
    call zpinch_initProblemDescriptor(rparlist, ssectionName,&
        solver_getMinimumMultigridlevel(rsolver),&
        solver_getMaximumMultigridlevel(rsolver),&
        rproblemDescriptor)
    call problem_initProblem(rproblemDescriptor, rproblem)

    ! Initialize the individual problem levels
    call hydro_initAllProblemLevels(rparlist,&
        ssectionNameHydro, rproblem, rcollection, p_rbdrCondHydro)
    call transp_initAllProblemLevels(rparlist,&
        ssectionNameTransport, rproblem, rcollection, p_rbdrCondTransport)

    ! Prepare internal data arrays of the solver structure for hydrodynamic model
    call parlst_getvalue_int(rparlist,&
        ssectionNameHydro, 'systemMatrix', systemMatrix)
    call parlst_getvalue_int(rparlist,&
        ssectionNameHydro, 'isystemFormat', isystemFormat)

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
      call parlst_getvalue_string(rparlist, ssectionName,&
          'algorithm', algorithm)
               
      ! What solution algorithm should be applied?
      if (trim(algorithm) .eq. 'transient_primal') then
        
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Solve the primal formulation for the time-dependent problem
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call zpinch_solveTransientPrimal(rparlist, ssectionName,&
            ssectionNameHydro, ssectionNameTransport, RbdrCond,&
            rproblem, rtimestep, rsolver, Rsolution, rcollection)

        call zpinch_outputSolution(rparlist, ssectionName,&
            rproblem%p_rproblemLevelMax, rsolution, rtimestep%dTime)
        
      else
        call output_line(trim(algorithm)//' is not a valid solution algorithm!',&
            OU_CLASS_ERROR,OU_MODE_STD,'zpinch_app')
        call sys_halt()
      end if

    else

      ! Just output the computational mesh and exit
      call zpinch_outputSolution(rparlist, ssectionName,&
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
    call bdrc_release(p_rbdrCondHydro)
    call bdrc_release(p_rbdrCondTransport)

    ! Release vectors
    call lsysbl_releaseVector(p_rsolutionHydro)
    call lsysbl_releaseVector(p_rsolutionTransport)

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

  subroutine zpinch_initProblemDescriptor(rparlist, ssectionName,&
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

    ! section names
    character(LEN=SYS_STRLEN) :: ssectionNameHydro
    character(LEN=SYS_STRLEN) :: ssectionNameTransport

    ! abstract problem descriptor
    type(t_problemDescriptor) :: rproblemDescriptorHydro
    type(t_problemDescriptor) :: rproblemDescriptorTransport

    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'subapplication', ssectionNameHydro, isubstring=1)
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'subapplication', ssectionNameTransport, isubstring=2)

    call hydro_initProblemDescriptor(rparlist, ssectionNameHydro,&
        nlmin, nlmax, rproblemDescriptorHydro)
    call transp_initProblemDescriptor(rparlist, ssectionNameTransport,&
        nlmin, nlmax, rproblemDescriptorTransport)
    
    rproblemDescriptor = problem_combineDescriptors(rproblemDescriptorHydro,&
                                                    rproblemDescriptorTransport)

  end subroutine zpinch_initProblemDescriptor

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_calcAdaptationIndicator(rsolutionHydro,&
      rsolutionTransport, rindicator)

!<description>
    ! This subroutine computes the element-wise indicator for mesh
    ! refinement and re-coarsening based on the detectation of shocks
    ! and contact discontinuities and the strength of gradient variations
!</description>

!<input>
    ! solution vectors
    type(t_vectorBlock), intent(in) :: rsolutionHydro, rsolutionTransport
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
    p_rtriangulation => rsolutionHydro%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation

    ! Extract primitive variables from conservative solution values
    call hydro_getVariable(rsolutionHydro, 'density', rdensity)
    call hydro_getVariable(rsolutionHydro, 'pressure', rpressure)

    ! Create element-wise indicator
    call lsyssc_createVector(rindicator, p_rtriangulation%NEL, .true.)

    ! Set pointers
    call storage_getbase_double2D(p_rtriangulation%h_DvertexCoords,&
        p_DvertexCoords)
    call storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    call storage_getbase_int2D(p_rtriangulation%h_IneighboursAtElement,&
        p_IneighboursAtElement)
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
    character(LEN=SYS_STRLEN) :: ssectionNameHydro
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

    call parlst_getvalue_string(rparlist,&
        ssectionName, 'subapplication', ssectionNameHydro, isubstring=1)
    call parlst_getvalue_int(rparlist,&
        ssectionNameHydro, 'isystemformat', isystemformat)
    
    ! Initialize the UCD exporter
    call flagship_initUCDexport(rproblemLevel, sucdsolution,&
        iformatUCD, rexport, ifilenumber)

    ! Increase filenumber by one
    ifilenumber = ifilenumber+1

    ! Set simulation time
    if (present(dtime)) call ucd_setSimulationTime(rexport, dtime)

    ! Add solution vectors
    if (present(rsolution)) then

      ! Add solution for hydrodynamic model
      call lsysbl_getbase_double(rsolution(1), p_Dsolution)
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

      case default
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
              call hydro_getVarInterleaveFormat1d(rvector1%NEQ, NVAR1D,&
                  'velocity_x', p_Dsolution, p_Ddata1)
              call ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1)

            case (NDIM2D)
              call hydro_getVarInterleaveFormat2d(rvector1%NEQ, NVAR2D,&
                  'velocity_x', p_Dsolution, p_Ddata1)
              call hydro_getVarInterleaveFormat2d(rvector2%NEQ, NVAR2D,&
                  'velocity_y', p_Dsolution, p_Ddata2)
              call ucd_addVarVertBasedVec(rexport, 'velocity',&
                  p_Ddata1, p_Ddata2)

            case (NDIM3D)
              call hydro_getVarInterleaveFormat3d(rvector1%NEQ, NVAR3D,&
                  'velocity_x', p_Dsolution, p_Ddata1)
              call hydro_getVarInterleaveFormat3d(rvector2%NEQ, NVAR3D,&
                  'velocity_y', p_Dsolution, p_Ddata2)
              call hydro_getVarInterleaveFormat3d(rvector3%NEQ, NVAR3D,&
                  'velocity_z', p_Dsolution, p_Ddata3)
              call ucd_addVarVertBasedVec(rexport, 'velocity',&
                  p_Ddata1, p_Ddata2, p_Ddata3)
            end select

          elseif (trim(cvariable) .eq. 'momentum') then
            
            ! Special treatment of momentum vector
            select case(ndim)
            case (NDIM1D)
              call hydro_getVarInterleaveFormat1d(rvector1%NEQ, NVAR1D,&
                  'momentum_x', p_Dsolution, p_Ddata1)
              call ucd_addVarVertBasedVec(rexport, 'momentum', p_Ddata1)

            case (NDIM2D)
              call hydro_getVarInterleaveFormat2d(rvector1%NEQ, NVAR2D,&
                  'momentum_x', p_Dsolution, p_Ddata1)
              call hydro_getVarInterleaveFormat2d(rvector2%NEQ, NVAR2D,&
                  'momentum_y', p_Dsolution, p_Ddata2)
              call ucd_addVarVertBasedVec(rexport, 'momentum',&
                  p_Ddata1, p_Ddata2)

            case (NDIM3D)
              call hydro_getVarInterleaveFormat3d(rvector1%NEQ, NVAR3D,&
                  'momentum_x', p_Dsolution, p_Ddata1)
              call hydro_getVarInterleaveFormat3d(rvector2%NEQ, NVAR3D,&
                  'momentum_y', p_Dsolution, p_Ddata2)
              call hydro_getVarInterleaveFormat3d(rvector3%NEQ, NVAR3D,&
                  'momentum_z', p_Dsolution, p_Ddata3)
              call ucd_addVarVertBasedVec(rexport, 'momentum',&
                  p_Ddata1, p_Ddata2, p_Ddata3)
            end select

          elseif(trim(cvariable) .eq. 'advect') then

            ! Special treatment for tracer quantity
            call lsysbl_getbase_double(rsolution(2), p_Ddata1)
            call ucd_addVariableVertexBased(rexport, 'advect',&
                UCD_VAR_STANDARD, p_Ddata1)

          else

            ! Standard treatment for scalar quantity
            select case(ndim)
            case (NDIM1D)
              call hydro_getVarInterleaveFormat1d(rvector1%NEQ, nvar,&
                  cvariable, p_Dsolution, p_Ddata1)
            case (NDIM2D)
              call hydro_getVarInterleaveFormat2d(rvector1%NEQ, nvar,&
                  cvariable, p_Dsolution, p_Ddata1)
            case (NDIM3D)
              call hydro_getVarInterleaveFormat3d(rvector1%NEQ, nvar,&
                  cvariable, p_Dsolution, p_Ddata1)
            end select
            call ucd_addVariableVertexBased(rexport, cvariable,&
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
              call hydro_getVarBlockFormat1d(rvector1%NEQ, NVAR1D,&
                  'velocity_x', p_Dsolution, p_Ddata1)
              call ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1)

            case (NDIM2D)
              call hydro_getVarBlockFormat2d(rvector1%NEQ, NVAR2D,&
                  'velocity_x', p_Dsolution, p_Ddata1)
              call hydro_getVarBlockFormat2d(rvector2%NEQ, NVAR2D,&
                  'velocity_y', p_Dsolution, p_Ddata2)
              call ucd_addVarVertBasedVec(rexport, 'velocity',&
                  p_Ddata1, p_Ddata2)

            case (NDIM3D)
              call hydro_getVarBlockFormat3d(rvector1%NEQ, NVAR3D,&
                  'velocity_x', p_Dsolution, p_Ddata1)
              call hydro_getVarBlockFormat3d(rvector2%NEQ, NVAR3D,&
                  'velocity_y', p_Dsolution, p_Ddata2)
              call hydro_getVarBlockFormat3d(rvector3%NEQ, NVAR3D,&
                  'velocity_z', p_Dsolution, p_Ddata3)
              call ucd_addVarVertBasedVec(rexport, 'velocity',&
                  p_Ddata1, p_Ddata2, p_Ddata3)
            end select
        
          elseif (trim(cvariable) .eq. 'momentum') then

            ! Special treatment of momentum vector
            select case(ndim)
            case (NDIM1D)
              call hydro_getVarBlockFormat1d(rvector1%NEQ, NVAR1D,&
                  'momentum_x', p_Dsolution, p_Ddata1)
              call ucd_addVarVertBasedVec(rexport, 'momentum', p_Ddata1)

            case (NDIM2D)
              call hydro_getVarBlockFormat2d(rvector1%NEQ, NVAR2D,&
                  'momentum_x', p_Dsolution, p_Ddata1)
              call hydro_getVarBlockFormat2d(rvector2%NEQ, NVAR2D,&
                  'momentum_y', p_Dsolution, p_Ddata2)
              call ucd_addVarVertBasedVec(rexport, 'momentum',&
                  p_Ddata1, p_Ddata2)

            case (NDIM3D)
              call hydro_getVarBlockFormat3d(rvector1%NEQ, NVAR3D,&
                  'momentum_x', p_Dsolution, p_Ddata1)
              call hydro_getVarBlockFormat3d(rvector2%NEQ, NVAR3D,&
                  'momentum_y', p_Dsolution, p_Ddata2)
              call hydro_getVarBlockFormat3d(rvector3%NEQ, NVAR3D,&
                  'momentum_z', p_Dsolution, p_Ddata3)
              call ucd_addVarVertBasedVec(rexport, 'momentum',&
                  p_Ddata1, p_Ddata2, p_Ddata3)
            end select

          elseif(trim(cvariable) .eq. 'advect') then
            
            ! Special treatment for tracer quantity
            call lsysbl_getbase_double(rsolution(2), p_Ddata1)
            call ucd_addVariableVertexBased(rexport, 'advect',&
                UCD_VAR_STANDARD, p_Ddata1)

          else
            
            ! Standard treatment for scalar quantity
            select case(ndim)
            case (NDIM1D)
              call hydro_getVarBlockFormat1d(rvector1%NEQ, nvar,&
                  cvariable, p_Dsolution, p_Ddata1)
            case (NDIM2D)
              call hydro_getVarBlockFormat2d(rvector1%NEQ, nvar,&
                  cvariable, p_Dsolution, p_Ddata1)
            case (NDIM3D)
              call hydro_getVarBlockFormat3d(rvector1%NEQ, nvar,&
                  cvariable, p_Dsolution, p_Ddata1)
            end select
            call ucd_addVariableVertexBased(rexport, cvariable,&
                UCD_VAR_STANDARD, p_Ddata1)
            
          end if
        end do

      case default
        call output_line('Invalid system format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'zpinch_outputSolution')
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

    ! section name in parameter list and collection structure
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

      if (rtriangulationSrc%ndim .eq. NDIM2D) then
        call zpinch_hadaptCallbackScalar2D(&
            HADAPT_OPR_INITCALLBACK, rcollection)
        call hadapt_performAdaptation(rhadapt, rindicator,&
            rcollection, zpinch_hadaptCallbackScalar2D)
      else
        call output_line('Invalid spatial dimension!',&
            OU_CLASS_ERROR,OU_MODE_STD,'zpinch_adaptTriangulation')
        call sys_halt()
      end if

    case (SYSTEM_BLOCKFORMAT)

      if (rtriangulationSrc%ndim .eq. NDIM2D) then
        call zpinch_hadaptCallbackBlock2D(&
            HADAPT_OPR_INITCALLBACK, rcollection)
        call hadapt_performAdaptation(rhadapt, rindicator,&
            rcollection, zpinch_hadaptCallbackBlock2D)
      else
        call output_line('Invalid spatial dimension!',&
            OU_CLASS_ERROR,OU_MODE_STD,'zpinch_adaptTriangulation')
        call sys_halt()
      end if

    case default

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
      ssectionNameHydro, ssectionNameTransport, RbdrCond,&
      rproblem, rtimestep, rsolver, Rsolution, rcollection)

!<description>
      ! This subroutine solves the transient primal simplified MHD problem.
!</description>

!<input>
    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName
    character(LEN=*), intent(in) :: ssectionNameHydro
    character(LEN=*), intent(in) :: ssectionNameTransport

    ! array of boundary condition structures
    type(t_boundaryCondition), dimension(:), intent(in), target :: RbdrCond
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

    ! array of primal solution vectors
    type(t_vectorBlock), dimension(:), intent(inout), target :: Rsolution

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! Pointer to the current multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! Pointer to the discretisation structure
    type(t_blockDiscretisation), pointer :: p_rdiscretisation

    ! Pointer to the solver structures
    type(t_solver), pointer :: p_rsolverHydro, p_rsolverTransport

    ! Pointer to the solution vectors
    type(t_vectorBlock), pointer :: p_rsolutionHydro, p_rsolutionTransport

    ! Pointer to the boundary condition structure
    type(t_boundaryCondition), pointer :: p_rbdrCondHydro, p_rbdrCondTransport

    ! Matrix for the linearised FCT algorithm
    type(t_matrixScalar) :: rmatrix1

    ! Vectors for the linearised FCT algorithm
    type(t_vectorBlock), dimension(2) :: Rvector1, Rvector2, Rvector3

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
    type(t_vectorBlock), dimension(2) :: Rforce

    ! local variables
    type(t_ucdExport) :: rimport
    character(LEN=SYS_STRLEN) :: sadaptivityName
    character(LEN=SYS_STRLEN) :: soutputName
    character(LEN=SYS_STRLEN) :: sucdimport
    
    real(dp) :: dstepUCD, dtimeUCD, dstepAdapt, dtimeAdapt
    real(dp) :: dscale
    integer :: templateMatrix, systemMatrix, isystemFormat
    integer :: discretisationHydro, discretisationTransport
    integer :: ipreadapt, npreadapt, ndimension
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
    
    ! Set pointers to solver structures
    p_rsolverHydro => solver_getNextSolver(rsolver, 1)
    p_rsolverTransport => solver_getNextSolver(rsolver, 2)

    ! Set pointers to solution vectors
    p_rsolutionHydro => Rsolution(1)
    p_rsolutionTransport => Rsolution(2)

    ! Set pointer to boundary condition structures
    p_rbdrCondHydro => RbdrCond(1)
    p_rbdrCondTransport => RbdrCond(2)

    ! Start time measurement for pre-processing
    call stat_startTimer(p_rtimerPrePostprocess, STAT_TIMERSHORT)

    ! Get global parameters
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'ndimension', ndimension)

    ! Set pointer to maximum problem level
    p_rproblemLevel => rproblem%p_rproblemLevelMax


    !--- compressible hydrodynamic model ---------------------------------------

    ! Get position of discretisation structure
    call parlst_getvalue_int(rparlist,&
        ssectionNameHydro, 'discretisation', discretisationHydro)

    ! Set pointer to discretisation structure
    p_rdiscretisation =>&
        p_rproblemLevel%Rdiscretisation(discretisationHydro)

    ! Create the solution vector
    call lsysbl_createVectorBlock(p_rdiscretisation,&
        p_rsolutionHydro, .false., ST_DOUBLE)
    if (p_rdiscretisation%ncomponents .ne. hydro_getNVAR(p_rproblemLevel)) then
      p_rsolutionHydro%RvectorBlock(1)%NVAR = hydro_getNVAR(p_rproblemLevel)
      call lsysbl_resizeVectorBlock(p_rsolutionHydro,&
          p_rsolutionHydro%NEQ, .false., .false.)
    end if

    ! Initialize the solution vector
    call hydro_initSolution(rparlist, ssectionNameHydro,&
        p_rproblemLevel, rtimestep%dinitialTime, p_rsolutionHydro,&
        rcollection)
    
    ! Impose boundary conditions
    select case(ndimension)
    case (NDIM1D)
      call bdrf_filterVectorExplicit(p_rbdrCondHydro, p_rsolutionHydro,&
          rtimestep%dinitialTime, hydro_calcBoundaryvalues1d)
    case (NDIM2D)
      call bdrf_filterVectorExplicit(p_rbdrCondHydro, p_rsolutionHydro,&
          rtimestep%dinitialTime, hydro_calcBoundaryvalues2d)
    case (NDIM3D)
      call bdrf_filterVectorExplicit(p_rbdrCondHydro, p_rsolutionHydro,&
          rtimestep%dinitialTime, hydro_calcBoundaryvalues3d)
    end select

    ! Attach the boundary condition
    call solver_setBoundaryCondition(p_rsolverHydro, p_rbdrCondHydro, .true.)

    ! Set collection to primal problem mode
    call parlst_addvalue(rparlist, ssectionNameHydro, 'mode', 'primal')


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

    ! Impose boundary conditions
    call bdrf_filterVectorExplicit(p_rbdrCondTransport,&
        p_rsolutionTransport, rtimestep%dinitialTime)

    ! Attach the boundary condition
    call solver_setBoundaryCondition(p_rsolverTransport,&
        p_rbdrCondTransport, .true.)

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
        ! the template matrix is the same for the hydrodynamic model
        ! and the transport model. Therefore, the template matrix is
	! adopted from the transport model.
        call parlst_getvalue_int(rparlist,&
            ssectionNameHydro, 'templateMatrix', templateMatrix)
        call grph_createGraphFromMatrix(&
            p_rproblemLevel%Rmatrix(templateMatrix), rgraph)
        call collct_setvalue_graph(rcollection,&
            'sparsitypattern', rgraph, .true.)

        ! Perform pre-adaptation?
        if (npreadapt > 0) then

          ! Perform number of pre-adaptation steps
          do ipreadapt = 1, npreadapt

            ! Compute the error estimator based on the tracer
            call zpinch_calcAdaptationIndicator(p_rsolutionHydro,&
                p_rsolutionTransport, relementError)

            ! Set the names of the template matrix
            rcollection%SquickAccess(1) = 'sparsitypattern'

            ! Attach the solution vectors to the collection structure
            rcollection%p_rvectorQuickAccess1 => rsolution(1)
            rcollection%p_rvectorQuickAccess2 => rsolution(2)

            ! Perform h-adaptation and update the triangulation structure
            call zpinch_adaptTriangulation(rparlist, ssectionnameHydro,&
                rhadapt, p_rproblemLevel%rtriangulation,&
                relementError, rcollection)

            ! Release element-wise error distribution
            call lsyssc_releaseVector(rElementError)

            ! Generate standard mesh from raw mesh
            call tria_initStandardMeshFromRaw(&
                p_rproblemLevel%rtriangulation, rproblem%rboundary)

            ! Update the template matrix according to the sparsity pattern
            call grph_generateMatrix(rgraph,&
                p_rproblemLevel%Rmatrix(templateMatrix))

            ! Re-initialize all constant coefficient matrices
            call hydro_initProblemLevel(rparlist,&
                ssectionNameHydro, p_rproblemLevel, rcollection)
            call transp_initProblemLevel(rparlist,&
                ssectionNameTransport, p_rproblemLevel, rcollection)

            ! Resize the solution vector for the hydrodynamic model accordingly
            call parlst_getvalue_int(rparlist,&
                ssectionNameHydro, 'systemMatrix', systemMatrix)
            call lsysbl_resizeVecBlockIndMat(&
                p_rproblemLevel%RmatrixBlock(systemMatrix),&
                p_rsolutionHydro, .false., .true.)

            ! Resize the solution vector for the transport model accordingly
            call lsysbl_resizeVectorBlock(p_rsolutionTransport, &
                p_rproblemLevel%Rmatrix(templateMatrix)%NEQ, .false.)

            ! Re-generate the initial solution vector for the hydrodynamic model
            call hydro_initSolution(rparlist, ssectionnameHydro,&
                p_rproblemLevel, rtimestep%dinitialTime,&
                p_rsolutionHydro, rcollection)

            ! Impose boundary conditions
            select case(ndimension)
            case (NDIM1D)
              call bdrf_filterVectorExplicit(p_rbdrCondHydro,&
                  p_rsolutionHydro, rtimestep%dinitialTime,&
                  hydro_calcBoundaryvalues1d)

            case (NDIM2D)
              call bdrf_filterVectorExplicit(p_rbdrCondHydro,&
                  p_rsolutionHydro, rtimestep%dinitialTime,&
                  hydro_calcBoundaryvalues2d)

            case (NDIM3D)
              call bdrf_filterVectorExplicit(p_rbdrCondHydro,&
                  p_rsolutionHydro, rtimestep%dinitialTime,&
                  hydro_calcBoundaryvalues3d)
            end select

            ! Re-generate the initial solution vector for the transport model
            call transp_initSolution(rparlist, ssectionnameTransport,&
                p_rproblemLevel, rtimestep%dinitialTime,&
                p_rsolutionTransport, rcollection)

            ! Impose boundary conditions
            call bdrf_filterVectorExplicit(p_rbdrCondTransport,&
                p_rsolutionTransport, rtimestep%dinitialTime)
          end do

          ! Prepare internal data arrays of the solver structure
          call parlst_getvalue_int(rparlist,&
              ssectionNameHydro, 'systemMatrix', systemMatrix)
          call parlst_getvalue_int(rparlist,&
              ssectionNameHydro, 'isystemFormat', isystemFormat)
          call flagship_updateSolverMatrix(p_rproblemLevel,&
              p_rsolverHydro, systemMatrix, isystemFormat, UPDMAT_ALL)
          call solver_updateStructure(p_rsolverHydro)

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
        trim(soutputName), 'sucdimport', sucdimport, '')

    ! Do we have to read in a precomdputed solution?
    if (trim(sucdimport) .ne. '') then
      call ucd_readGMV(sucdimport, rimport, p_rproblemLevel%rtriangulation)
      call ucd_getSimulationTime(rimport, rtimestep%dinitialTime)
      call ucd_getSimulationTime(rimport, rtimestep%dTime)
      call hydro_setVariables(rimport, p_rsolutionHydro)
      call transp_setVariable(rimport, 'advect', p_rsolutionTransport)
      call ucd_release(rimport)

      ! Set time for solution output
      dtimeUCD = rtimestep%dinitialTime
    end if

    !---------------------------------------------------------------------------
    ! Calculate velocity field (\rho v) for the scalar model problem
    ! based on the initial solution U^0 of the hydrodynamic model
    !---------------------------------------------------------------------------

    call zpinch_calcVelocityField(rparlist, ssectionNameTransport,&
        p_rproblemLevel, p_rsolutionHydro, rcollection)

    ! Stop time measurement for pre-processing
    call stat_stopTimer(p_rtimerPrePostprocess)

    !---------------------------------------------------------------------------
    ! Infinite time stepping loop
    !---------------------------------------------------------------------------

    timeloop: do

      ! Check for user interaction
      if (signal_SIGINT(-1) > 0 )&
          call zpinch_outputSolution(rparlist, ssectionName,&
          p_rproblemLevel, rsolution, rtimestep%dTime)

      ! Start time measurement for solution procedure
      call stat_startTimer(p_rtimerSolution, STAT_TIMERSHORT)

      !-----------------------------------------------------------------------
      ! Compute hydrodynamic model + scalar tracer in coupled fashion
      ! for full time step: $U^n \to U^{n+1}$ and $u^n \to u^{n+1}$
      !-----------------------------------------------------------------------

      ! Start time measurement for solution procedure
      call stat_startTimer(p_rtimerSolution, STAT_TIMERSHORT)

      ! Compute scaling for explicit part of the Lorentz force
      dscale = (1.0_DP-rtimestep%theta) * rtimestep%dStep

      if (dscale .ne. 0.0_DP) then
        ! Calculate explicit part of the Lorentz force term
        call zpinch_calcLorentzforceTerm(rparlist, ssectionName,&
            ssectionNameHydro, ssectionNameTransport, p_rproblemLevel,&
            p_rsolutionHydro, p_rsolutionTransport, rtimestep%dTime, &
            dscale, .true., Rforce(1), rcollection)       
      end if
      
      ! Prepare quick access arrays/vectors
      rcollection%p_rvectorQuickAccess1 => rsolution(1)
      rcollection%p_rvectorQuickAccess2 => rsolution(2)
      rcollection%p_rvectorQuickAccess3 => rtimestep%RtempVectors(1)
      rcollection%p_rvectorQuickAccess4 => rtimestep%RtempVectors(2)

      ! What time-stepping scheme should be used?
      select case(rtimestep%ctimestepType)

      case (TSTEP_RK_SCHEME)
        
        ! Adopt explicit Runge-Kutta scheme
        if (dscale .ne. 0.0_DP) then
          
          ! ... with source term
          call tstep_performRKStep(p_rproblemLevel, rtimestep,&
              rsolver, rsolution, zpinch_nlsolverCallback,&
              rcollection, Rforce)
        else
          ! ... without source term
          call tstep_performRKStep(p_rproblemLevel, rtimestep,&
              rsolver, rsolution, zpinch_nlsolverCallback, rcollection)
        end if


      case (TSTEP_THETA_SCHEME)

        ! Adopt two-level theta-scheme
        if (dscale .ne. 0.0_DP) then

          ! ... with source term
          call tstep_performThetaStep(p_rproblemLevel, rtimestep,&
              rsolver, rsolution, zpinch_nlsolverCallback,&
              rcollection, Rforce)
        else
          ! ... without source term
          call tstep_performThetaStep(p_rproblemLevel, rtimestep,&
              rsolver, rsolution, zpinch_nlsolverCallback, rcollection)
        end if


      case default
        call output_line('Unsupported time-stepping algorithm!',&
            OU_CLASS_ERROR,OU_MODE_STD,'zpinch_solveTransientPrimal')
        call sys_halt()
      end select

      !-------------------------------------------------------------------------
      ! Compute linearised FCT correction for hydrodynamic and transport model
      !-------------------------------------------------------------------------

      ! Apply linearised FEM-FCT post-processing
      call zpinch_calcLinearisedFCT(p_rproblemLevel, rtimestep,&
          p_rsolverHydro, p_rsolverTransport, Rsolution, ssectionName,&
          ssectionNameHydro, ssectionNameTransport, rcollection,&
          rmatrix=rmatrix1, Rvector1=Rvector1,&
          Rvector2=Rvector2, Rvector3=Rvector3)

      ! Calculate velocity field (\rho v)
      call zpinch_calcVelocityField(rparlist, ssectionNameTransport,&
          p_rproblemLevel, p_rsolutionHydro, rcollection)
    
      ! Stop time measurement for solution procedure
      call stat_stopTimer(p_rtimerSolution)
      
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
!!$        call hydro_estimateRecoveryError(rparlist,&
!!$            ssectionnameHydro, p_rproblemLevel, p_rsolutionHydro,&
!!$            rtimestep%dinitialTime, relementError, derror)

        ! Compute the error indicator based on the tracer
        call zpinch_calcAdaptationIndicator(p_rsolutionHydro,&
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
        call zpinch_adaptTriangulation(rparlist, ssectionnameHydro,&
            rhadapt, p_rproblemLevel%rtriangulation, relementError,&
            rcollection)

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
        call hydro_initProblemLevel(rparlist, ssectionNameHydro,&
            p_rproblemLevel, rcollection)
        call transp_initProblemLevel(rparlist, ssectionNameTransport,&
            p_rproblemLevel, rcollection)

        ! Resize the solution vector for the hydrodynamic model accordingly
        call parlst_getvalue_int(rparlist,&
            ssectionNameHydro, 'systemmatrix', systemMatrix)
        call lsysbl_resizeVecBlockIndMat(&
            p_rproblemLevel%RmatrixBlock(systemMatrix),&
            p_rsolutionHydro, .false., .true.)

        ! Resize the solution vector for the transport model accordingly
        call lsysbl_resizeVectorBlock(p_rsolutionTransport, &
            p_rproblemLevel%Rmatrix(templateMatrix)%NEQ, .false.)

        ! Prepare internal data arrays of the solver structure
        call parlst_getvalue_int(rparlist,&
            ssectionNameHydro, 'systemmatrix', systemMatrix)
        call parlst_getvalue_int(rparlist,&
            ssectionNameHydro, 'isystemformat', isystemFormat)
        call flagship_updateSolverMatrix(p_rproblemLevel,&
            p_rsolverHydro, systemMatrix, isystemFormat, UPDMAT_ALL)
        call solver_updateStructure(p_rsolverHydro)

        ! Prepare internal data arrays of the solver structure
        call parlst_getvalue_int(rparlist,&
            ssectionNameTransport, 'systemmatrix', systemMatrix)
        call flagship_updateSolverMatrix(p_rproblemLevel,&
            p_rsolverTransport, systemMatrix, SYSTEM_INTERLEAVEFORMAT, UPDMAT_ALL)
        call solver_updateStructure(p_rsolverTransport)

        ! Calculate velocity field (\rho v)
        call zpinch_calcVelocityField(rparlist, ssectionNameTransport,&
            p_rproblemLevel, p_rsolutionHydro, rcollection)

        ! Stop time measurement for generation of constant
        ! coefficient matrices
        call stat_stopTimer(p_rtimerAssemblyCoeff)

      end if

    end do timeloop

    ! Release Lorentz force vector
    call lsysbl_releaseVector(Rforce(1))

    ! Release adaptation structure
    if (trim(adjustl(sadaptivityName)) .ne. '') then
      if ((dstepAdapt > 0.0_DP) .or. (npreadapt > 0)) then
        call hadapt_releaseAdaptation(rhadapt)
        call grph_releaseGraph(rgraph)
      end if
    end if

    ! Release vectors and matrices for the linearised FCT algorithm
    call lsyssc_releaseMatrix(rmatrix1)
    call lsysbl_releaseVector(Rvector1(1))
    call lsysbl_releaseVector(Rvector1(2))
    call lsysbl_releaseVector(Rvector2(1))
    call lsysbl_releaseVector(Rvector2(2))
    call lsysbl_releaseVector(Rvector3(1))
    call lsysbl_releaseVector(Rvector3(2))

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

  end subroutine zpinch_parseCmdlArguments

end module zpinch_application
