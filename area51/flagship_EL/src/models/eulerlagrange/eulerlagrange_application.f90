!##############################################################################
!# ****************************************************************************
!# <name> eulerlagrange_application </name>
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
!# 1.) eulerlagrange_app
!#     -> The main routine of the application called from the main
!#        program. The routine gets all required information from the
!#        parameter list which needs to be initialised and filled in
!#        the main program. It then works black-box, that is, it
!#        determines the solution algorithm to be used and performs
!#        the simulation. The user should only have to modify this
!#        routine if another solution algorithm is implemented.
!#
!# 2.) eulerlagrange_initSolvers
!#     -> Initializes the solve structures from the parameter list.
!#
!# 3.) eulerlagrange_initProblem
!#     -> Initializes the global problem structure based on the
!#        parameter settings given by the parameter list. This routine
!#        is quite universal, that is, it prepares the internal
!#        structure of the global problem and generates a linked
!#        list of problem levels used in the multigrid hierarchy.
!#
!# 4.) eulerlagrange_initProblemLevel
!#     -> Initializes the individual problem level based on the
!#        parameter settings given by the parameter list.
!#        This routine is called repeatedly by the global
!#        initialisation routine eulerlagrange_initAllProblemLevels.
!#
!# 5.) eulerlagrange_initAllProblemLevels
!#     -> Initializes ALL problem levels attached to the global
!#        problem structure based on the parameter settings
!#        given by the parameter list.
!#
!# 6.) eulerlagrange_initSolution
!#     -> Initializes the solution vector based on the parameter
!#        settings given by the parameter list
!#
!# 7.) eulerlagrange_outputSolution
!#     -> Outputs the solution vector to file in UCD format
!#
!# 8.) eulerlagrange_outputStatistics
!#     -> Outputs the application statitics
!#
!# 9.) eulerlagrange_estimateRecoveryError
!#      -> Estimates the solution error using recovery techniques
!#
!# 10.) eulerlagrange_adaptTriangulation
!#      -> Performs h-adaptation for the given triangulation
!#
!# 11.) eulerlagrange_solveTransientPrimal
!#      -> Solves the primal formulation of the time-dependent
!#         compressible Euler equations.
!#
!# 12.) eulerlagrange_init
!#      -> Initialization for the particles
!#
!# 13.) eulerlagrange_step
!#      -> Calculates the particle motion in each timestep
!#
!#
!# The following auxiliary routines are available:
!#
!# 1.) eulerlagrange_parseCmdlArguments
!#     -> Parses the list of commandline arguments and overwrites
!#        parameter values from the parameter files
!# 
!# </purpose>
!##############################################################################

module eulerlagrange_application
    
  use afcstabilisation
  use basicgeometry
  use bilinearformevaluation
  use boundary
  use boundaryfilter
  use collection
  use cubature
  use derivatives
  use element
  use eulerlagrange_basic
  use eulerlagrange_callback
  use eulerlagrange_callback1d
  use eulerlagrange_callback2d
  use eulerlagrange_callback3d
  use flagship_basic
  use fparser
  use fsystem
  use genoutput
  use geometryaux
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
  use thermodynamics
  use timestep
  use timestepaux
  use triasearch
  use ucd

  implicit none

  private
  public :: eulerlagrange_app
  public :: eulerlagrange_adaptTriangulation
  public :: eulerlagrange_adjustParameterlist
  public :: eulerlagrange_estimateRecoveryError
  public :: eulerlagrange_initAllProblemLevels
  public :: eulerlagrange_initProblem
  public :: eulerlagrange_initProblemLevel
  public :: eulerlagrange_initSolution
  public :: eulerlagrange_initSolvers
  public :: eulerlagrange_outputSolution
  public :: eulerlagrange_outputStatistics
  public :: eulerlagrange_solveTransientPrimal
  public :: eulerlagrange_init
  public :: eulerlagrange_step

  contains

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_app(rparlist, ssectionName)

!<description>
    ! This is the main application for the compressible Euler
    ! equations.  It is a so-called driver routine which can be used
    ! to start a standalone Euler simulation.
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
    
    ! Pointer to the multigrid level
    type(t_problemLevel), pointer :: rproblemLevel

    ! Solution vector
    type(t_vectorBlock) :: rsolution

    ! particles
    type(t_Particles) :: rParticles
   
    ! Timer for the total solution process
    type(t_timer) :: rtimerTotal

    ! Timer for the solution process
    type(t_timer) :: rtimerSolution

    ! Timer for the particle phase computation
    type(t_timer) ::rtimerParticlephase

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


    type(t_timer), pointer :: p_rtimerParticlephase



    real(DP), dimension(:), pointer :: p_xpos, p_ypos, p_zpos
    real(DP), dimension(:), pointer :: p_xpos_old, p_ypos_old, p_zpos_old
    real(DP), dimension(:), pointer :: p_xvelo, p_yvelo, p_zvelo
    real(DP), dimension(:), pointer :: p_xvelo_old, p_yvelo_old, p_zvelo_old
    real(DP), dimension(:), pointer :: p_xvelo_gas, p_yvelo_gas, p_zvelo_gas
    real(DP), dimension(:), pointer :: p_xvelo_gas_old, p_yvelo_gas_old, p_zvelo_gas_old
    real(DP), dimension(:), pointer :: p_lambda1, p_lambda2, p_lambda3, p_lambda4
    real(DP), dimension(:), pointer :: p_density
    real(DP), dimension(:), pointer :: p_diam

    integer, dimension(:), pointer :: p_element

    real(DP), dimension(:,:), pointer :: p_midpoints_el

    ! local variables
    integer :: isystemFormat, systemMatrix, ndimension
    
    ! One-way, two-way cpupling or steady state gas solution
    integer :: icouplingpart

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
    call eulerlagrange_parseCmdlArguments(rparlist)
    call eulerlagrange_adjustParameterlist(rparlist, ssectionName)

    ! Initialize global collection structure
    call collct_init(rcollection)

    !  Attach the parameter list and the timers to the collection
    call collct_setvalue_parlst(rcollection,&
        'rparlist', rparlist, .true.)
    call collct_setvalue_timer(rcollection,&
        'rtimerSolution', rtimerSolution, .true.)
    call collct_setvalue_timer(rcollection,&
        'rtimerParticlephase', rtimerParticlephase, .true.)
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
    call eulerlagrange_initSolvers(rparlist, ssectionName, rtimestep, rsolver)

    ! Initialize the abstract problem structure
    call eulerlagrange_initProblem(rparlist, ssectionName,&
        solver_getMinimumMultigridlevel(rsolver),&
        solver_getMaximumMultigridlevel(rsolver), rproblem)

    ! Initialize the individual problem levels
    call eulerlagrange_initAllProblemLevels(rparlist, ssectionName,&
        rproblem, rcollection)

    ! Prepare internal data arrays of the solver structure
    call parlst_getvalue_int(rparlist, ssectionName,&
        'systemMatrix', systemMatrix)
    call parlst_getvalue_int(rparlist, ssectionName,&
        'isystemFormat', isystemFormat)
    call flagship_updateSolverMatrix(rproblem%p_rproblemLevelMax,&
        rsolver, systemMatrix, isystemFormat, UPDMAT_ALL)
    call solver_updateStructure(rsolver)

    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Euler-Lagrange
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    call eulerlagrange_init(rparlist,rproblem%p_rproblemLevelMax,rsolution,&
                            rtimestep,rcollection,rParticles)
    
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
      call bdrf_readBoundaryCondition(rbdrCondPrimal, sindatfileName,&
          '['//trim(sbdrcondName)//']', ndimension)
      
      ! What solution algorithm should be applied?
      if (trim(algorithm) .eq. 'transient_primal') then
        
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Solve the primal formulation for #
        ! the time-dependent problem
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call eulerlagrange_solveTransientPrimal(rparlist, ssectionName,&
            rbdrCondPrimal, rproblem, rtimestep, rsolver,&
            rsolutionPrimal, rcollection,rParticles)

        call eulerlagrange_outputSolution(rparlist, ssectionName, rproblem&
            %p_rproblemLevelMax, rParticles, rsolutionPrimal, dtime=rtimestep&
            %dTime)
 
         ! One way or twoway coupling?
        call parlst_getvalue_int(rparlist, ssectionName, "icouplingpart", icouplingpart)
         
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Euler-Lagrange (for steadystate gas solution)
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if (icouplingpart == -1) then
        
          ! Start time measurement for solution procedure
          !call stat_startTimer(p_rtimerParticlephase)

            !do 
          
              ! Subroutine to compute the particle movement
              call eulerlagrange_step(rparlist,rproblem%p_rproblemLevelMax,rsolution,&
                    rtimestep,rcollection,rParticles)
               
            !end do
                
          ! Stop time measurement for solution procedure
          !call stat_stopTimer(p_rtimerParticlephase)
        
        end if
        
      else
        call output_line(trim(algorithm)//' is not a valid solution algorithm!',&
            OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_app')
        call sys_halt()
      end if

  else

      ! Just output the computational mesh and exit
      call eulerlagrange_outputSolution(rparlist, ssectionName,&
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
    call bdrf_release(rbdrCondPrimal)
    call bdrf_release(rbdrCondDual)

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
    call eulerlagrange_outputStatistics(rtimerTotal, rcollection)

    ! Release collection
    call collct_done(rcollection)

    call storage_free (rParticles%h_xpos)
    call storage_free (rParticles%h_ypos)
    call storage_free (rParticles%h_zpos)
    call storage_free (rParticles%h_xpos_old)
    call storage_free (rParticles%h_ypos_old)
    call storage_free (rParticles%h_zpos_old)
    call storage_free (rParticles%h_xvelo)
    call storage_free (rParticles%h_yvelo)
    call storage_free (rParticles%h_zvelo)
    call storage_free (rParticles%h_xvelo_old)
    call storage_free (rParticles%h_yvelo_old)
    call storage_free (rParticles%h_zvelo_old)
    call storage_free (rParticles%h_xvelo_gas)
    call storage_free (rParticles%h_yvelo_gas)
    call storage_free (rParticles%h_zvelo_gas)
    call storage_free (rParticles%h_xvelo_gas_old)
    call storage_free (rParticles%h_yvelo_gas_old)
    call storage_free (rParticles%h_zvelo_gas_old)
    call storage_free (rParticles%h_lambda1)
    call storage_free (rParticles%h_lambda2)
    call storage_free (rParticles%h_lambda3)
    call storage_free (rParticles%h_lambda4)
    call storage_free (rParticles%h_element)
    call storage_free (rParticles%h_diam)
    call storage_free (rParticles%h_density)
    call storage_free (rParticles%h_mass)
    call storage_free (rParticles%h_temp)
    call storage_free (rParticles%h_midpoints_el)
    call storage_free (rParticles%h_alpha_n)
    call storage_free (rParticles%h_bdy_time)
    call storage_free (rParticles%h_PartVol)
    call storage_free (rParticles%h_PartVolAver)


  end subroutine eulerlagrange_app

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_initSolvers(rparlist, ssectionName,&
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
    
  end subroutine eulerlagrange_initSolvers

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_initProblem(rparlist, ssectionName,&
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

  end subroutine eulerlagrange_initProblem
  
  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_initProblemLevel(rparlist, ssectionName,&
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
    integer :: templateMatrix,systemMatrix,jacobianMatrix,consistentMassMatrix
    integer :: lumpedMassMatrix,coeffMatrix_CX,coeffMatrix_CY, coeffMatrix_CZ
    integer :: inviscidAFC,discretisation,celement,isystemFormat,isystemCoupling
    integer :: imatrixFormat,ivar,jvar,ivariable,nvariable,nvartransformed
    integer :: i,j,nsumcubRefBilForm,nsumcubRefLinForm,nsumcubRefEval

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
              eulerlagrange_getNVAR(rproblemLevel), p_rtriangulation)

        case DEFAULT
          call output_line('Unsupported system format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_initProblemLevel')
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
              OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_initProblemLevel')
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
              OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_initProblemLevel')
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
              OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_initProblemLevel')
          call sys_halt()
        end select
        
      case DEFAULT
        call output_line('Invalid number of spatial dimensions',&
            OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_initProblemLevel')
        call sys_halt()
      end select
    
      ! Duplicate scalar discretisation structure for block matrix format
      if (isystemFormat .eq. SYSTEM_BLOCKFORMAT) then
        do ivar = 2, eulerlagrange_getNVAR(rproblemLevel)
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
          rproblemLevel%Rmatrix(systemMatrix)%NVAR = eulerlagrange_getNVAR(rproblemLevel)
          
          ! What matrix format should be used?
          select case(imatrixFormat)
          case (LSYSSC_MATRIX7)
            rproblemLevel%Rmatrix(systemMatrix)%cmatrixFormat = LSYSSC_MATRIX7INTL
            
          case (LSYSSC_MATRIX9)
            rproblemLevel%Rmatrix(systemMatrix)%cmatrixFormat = LSYSSC_MATRIX9INTL
            
          case DEFAULT
            call output_line('Unsupported matrix format!',&
                OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_initProblemLevel')
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
                             OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_initProblemLevel')
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
            do ivar = 1, eulerlagrange_getNVAR(rproblemLevel)
              call lsyssc_resizeMatrix(&
                  rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                  rproblemLevel%Rmatrix(templateMatrix), .false., .false., .true.)
            end do

          case (SYSTEM_ALLCOUPLED)
            ! Create all NVAR x NVAR blocks
            do ivar = 1, eulerlagrange_getNVAR(rproblemLevel)
              do jvar = 1, eulerlagrange_getNVAR(rproblemLevel)
                call lsyssc_resizeMatrix(&
                    rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar ,jvar),&
                    rproblemLevel%Rmatrix(templateMatrix), .false., .false., .true.)
              end do
            end do

          case DEFAULT
            call output_line('Unsupported block matrix format!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_initProblemLevel')
            call sys_halt()
          end select

        else   ! System matrix has no structure

          ! Create empty NVARxNVAR block matrix directly
          call lsysbl_createEmptyMatrix(&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              eulerlagrange_getNVAR(rproblemLevel),&
              eulerlagrange_getNVAR(rproblemLevel))

          ! Specify matrix as 'group matrix'
          rproblemLevel%RmatrixBlock(systemMatrix)%imatrixSpec = LSYSBS_MSPEC_GROUPMATRIX
          
          ! What kind of global operator should be adopted?
          select case(isystemCoupling)
            
          case (SYSTEM_SEGREGATED)
            ! Create only NVAR diagonal blocks
            do ivar = 1, eulerlagrange_getNVAR(rproblemLevel)
              call lsyssc_duplicateMatrix(&
                  rproblemLevel%Rmatrix(templateMatrix),&
                  rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                  LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
            end do
          
          case (SYSTEM_ALLCOUPLED)
            ! Create all NVAR x NVAR blocks
            do ivar = 1, eulerlagrange_getNVAR(rproblemLevel)
              do jvar = 1, eulerlagrange_getNVAR(rproblemLevel)
                call lsyssc_duplicateMatrix(&
                    rproblemLevel%Rmatrix(templateMatrix),&
                    rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,jvar),&
                    LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
              end do
            end do
            
          case DEFAULT
            call output_line('Unsupported block matrix format!',&
                OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_initProblemLevel')
            call sys_halt()
          end select

        end if
        
        ! Update internal structure of block matrix
        call lsysbl_updateMatStrucInfo(rproblemLevel%RmatrixBlock(systemMatrix))
        
      case DEFAULT
        call output_line('Unsupported system format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_initProblemLevel')
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

    
    ! Create coefficient matrix (phi, dphi/dx) duplicate of the template matrix
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

    
    ! Create coefficient matrix (phi, dphi/dy) duplicate of the template matrix
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

    
    ! Create coefficient matrix (phi, dphi/dz) duplicate of the template matrix
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

    ! Resize stabilisation structures if necessary and remove the
    ! indicator for the subdiagonal edge structure. If they are
    ! needed, then they are re-generated on-the-fly.
    if (inviscidAFC > 0) then
      if (rproblemLevel%Rafcstab(inviscidAFC)%iSpec .eq. AFCSTAB_UNDEFINED) then

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
              eulerlagrange_getNVARtransformed(rproblemLevel, slimitingvariable))
        end do

        ! Initialise stabilisation structure
        call gfsys_initStabilisation(&
            rproblemLevel%RmatrixBlock(systemMatrix),&
            rproblemLevel%Rafcstab(inviscidAFC), nvartransformed)

      else
        ! Resize stabilisation structure
        call afcstab_resizeStabilisation(&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rproblemLevel%Rmatrix(templateMatrix))

        rproblemLevel%Rafcstab(inviscidAFC)%iSpec =&
            iand(rproblemLevel%Rafcstab(inviscidAFC)%iSpec,&
            not(AFCSTAB_HAS_OFFDIAGONALEDGES))
      end if
    end if
    
  end subroutine eulerlagrange_initProblemLevel

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_initAllProblemLevels(rparlist, ssectionName,&
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
      call eulerlagrange_initProblemLevel(rparlist, ssectionName,&
          p_rproblemLevel, rcollection)
      
      ! Switch to next coarser level
      p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
    end do

  end subroutine eulerlagrange_initAllProblemLevels

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_initSolution(rparlist, ssectionName,&
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
          OU_CLASS_ERROR, OU_MODE_STD, 'eulerlagrange_initSolution')
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
            OU_CLASS_ERROR, OU_MODE_STD, 'eulerlagrange_initSolution')
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
            OU_CLASS_ERROR, OU_MODE_STD, 'eulerlagrange_initSolution')
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
      
      !  Initialize number of expressions
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
              rvector%RvectorBlock(iblock), eulerlagrange_coeffVectorAnalytic, rcollection)

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
                rvectorBlock%RvectorBlock(ivar), eulerlagrange_coeffVectorAnalytic,&
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
        rafcstab%iSpec= AFCSTAB_UNDEFINED
        rafcstab%bprelimiting = .false.
        rafcstab%ctypeAFCstabilisation = AFCSTAB_FEMFCT_MASS
        call gfsys_initStabilisation(&
            rproblemLevel%RmatrixBlock(systemMatrix), rafcstab)
        
        ! Compute the raw antidiffusive mass fluxes. Note that we may supply any
        ! callback function for assembling the antidiffusive fluxes since it 
        ! will not be used for assembling antidiffusive mass fluxes !!!
        call gfsys_buildFluxFCT((/p_rconsistentMassMatrix/),&
            rafcstab, rvectorHigh, rvectorHigh, eulerlagrange_calcFluxFCTScalarDiss1d,&
            0.0_DP, 0.0_DP, 1.0_DP, .true., p_rconsistentMassMatrix)
        
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
          call eulerlagrange_calcCorrectionFCT(rproblemLevel, rvector, 1.0_DP,&
              .false., AFCSTAB_FCTALGO_STANDARD-AFCSTAB_FCTALGO_CORRECT,&
              rvector, rcollection, rafcstab, 'ssolutionconstrainvariable')
          
          ! Apply failsafe flux correction
          call afcstab_failsafeLimiting(rafcstab, p_rlumpedMassMatrix,&
              SsolutionFailsafeVariables, 1.0_DP,&
              nsolutionfailsafe, eulerlagrange_getVariable, rvector)
          
          ! Deallocate temporal memory
          deallocate(SsolutionFailsafeVariables)

        else
          
          ! Compute and apply FEM-FCT correction
          call eulerlagrange_calcCorrectionFCT(rproblemLevel, rvector, 1.0_DP,&
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
          OU_CLASS_ERROR, OU_MODE_STD, 'eulerlagrange_initSolution')
      call sys_halt()
    end select

  end subroutine eulerlagrange_initSolution

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_outputSolution(rparlist, ssectionName,&
      rproblemLevel, rParticles, rsolutionPrimal, rsolutionDual, dtime)

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

    ! particles
    type(t_particles), intent(inout), optional :: rParticles

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
    type(t_vectorScalar) :: rvector1, rvector2, rvector3, rvector4, rvector5, rvector6
    real(DP), dimension(:), pointer :: p_Dsolution, p_Ddata1, p_Ddata2, p_Ddata3, p_Ddata4, p_Ddata5, p_Ddata6
    character(len=SYS_NAMELEN) :: cvariable
    integer :: iformatUCD, isystemFormat, isize, ndim, nvar, ivariable, nvariable, ivt, ipartoutput


    ! Local variables
	integer :: i, current
	real(DP) :: c_pi

    ! Current particlenumber
    integer :: iPart
    character(LEN=40) :: sfilename, sfilenamenew

    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'output', soutputName)
    call parlst_getvalue_string(rparlist,&
        trim(soutputName), 'sucdsolution', sucdsolution)
    call parlst_getvalue_int(rparlist,&
        trim(soutputName), 'iformatucd', iformatUCD)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'isystemformat', isystemformat)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'ipartoutput', ipartoutput)

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
      nvar  = eulerlagrange_getNVAR(rproblemLevel)
      isize = size(p_Dsolution)/nvar
      ndim  = rproblemLevel%rtriangulation%ndim
      
      ! Create auxiliary vectors
      select case(ndim)
      case (NDIM1D)
        call lsyssc_createVector(rvector1, isize, .false.)
        call lsyssc_createVector(rvector4, isize, .false.)
        call lsyssc_getbase_double(rvector1, p_Ddata1)
        call lsyssc_getbase_double(rvector4, p_Ddata4)
        
      case (NDIM2D)
        call lsyssc_createVector(rvector1, isize, .false.)
        call lsyssc_createVector(rvector2, isize, .false.)
        call lsyssc_createVector(rvector4, isize, .false.)
        call lsyssc_createVector(rvector5, isize, .false.)
        call lsyssc_getbase_double(rvector1, p_Ddata1)
        call lsyssc_getbase_double(rvector2, p_Ddata2)
        call lsyssc_getbase_double(rvector4, p_Ddata4)
        call lsyssc_getbase_double(rvector5, p_Ddata5)
        
      case (NDIM3D)
        call lsyssc_createVector(rvector1, isize, .false.)
        call lsyssc_createVector(rvector2, isize, .false.)
        call lsyssc_createVector(rvector3, isize, .false.)
        call lsyssc_createVector(rvector4, isize, .false.)
        call lsyssc_createVector(rvector5, isize, .false.)
        call lsyssc_createVector(rvector6, isize, .false.)
        call lsyssc_getbase_double(rvector1, p_Ddata1)
        call lsyssc_getbase_double(rvector2, p_Ddata2)
        call lsyssc_getbase_double(rvector3, p_Ddata3)
        call lsyssc_getbase_double(rvector4, p_Ddata4)
        call lsyssc_getbase_double(rvector5, p_Ddata5)
        call lsyssc_getbase_double(rvector6, p_Ddata6)
        
      case DEFAULT
        call output_line('Invalid number of spatial dimensions',&
                         OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_outputSolution')
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
              call eulerlagrange_getVarInterleaveFormat(rvector1%NEQ, NVAR1D,&
                  'velocity_x', p_Dsolution, p_Ddata1)
              call ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1)

            case (NDIM2D)
              call eulerlagrange_getVarInterleaveFormat(rvector1%NEQ, NVAR2D,&
                  'velocity_x', p_Dsolution, p_Ddata1)
              call eulerlagrange_getVarInterleaveFormat(rvector2%NEQ, NVAR2D,&
                  'velocity_y', p_Dsolution, p_Ddata2)
              call ucd_addVarVertBasedVec(rexport, 'velocity',&
                  p_Ddata1, p_Ddata2)

            case (NDIM3D)
              call eulerlagrange_getVarInterleaveFormat(rvector1%NEQ, NVAR3D,&
                  'velocity_x', p_Dsolution, p_Ddata1)
              call eulerlagrange_getVarInterleaveFormat(rvector2%NEQ, NVAR3D,&
                  'velocity_y', p_Dsolution, p_Ddata2)
              call eulerlagrange_getVarInterleaveFormat(rvector3%NEQ, NVAR3D,&
                  'velocity_z', p_Dsolution, p_Ddata3)
              call ucd_addVarVertBasedVec(rexport, 'velocity',&
                  p_Ddata1, p_Ddata2, p_Ddata3)
            end select

          elseif (trim(cvariable) .eq. 'velo_part') then
            
            ! Special treatment of velocity vector
            select case(ndim)
            case (NDIM1D)
              call eulerlagrange_getVarInterleaveFormat(rvector4%NEQ, NVAR1D,&
                  'velo_part_x', p_Dsolution, p_Ddata4)
              call ucd_addVarVertBasedVec(rexport, 'velo_part', p_Ddata4)

            case (NDIM2D)

              call ucd_addVarVertBasedVec(rexport, 'velo_part',&
                  rParticles%p_PartVelox, rParticles%p_PartVeloy)

            case (NDIM3D)
              call eulerlagrange_getVarInterleaveFormat(rvector4%NEQ, NVAR3D,&
                  'velo_part_x', p_Dsolution, p_Ddata4)
              call eulerlagrange_getVarInterleaveFormat(rvector5%NEQ, NVAR3D,&
                  'velo_part_y', p_Dsolution, p_Ddata5)
              call eulerlagrange_getVarInterleaveFormat(rvector6%NEQ, NVAR3D,&
                  'velo_part_z', p_Dsolution, p_Ddata6)
              call ucd_addVarVertBasedVec(rexport, 'velo_part',&
                  p_Ddata4, p_Ddata5, p_Ddata6)
            end select

          elseif (trim(cvariable) .eq. 'vol_part') then
          
            call ucd_addVariableVertexBased (rexport, cvariable,&
                UCD_VAR_STANDARD, rParticles%p_PartVol)

          elseif (trim(cvariable) .eq. 'vol_partaver') then
          
            do ivt= 1, rproblemLevel%rtriangulation%NVT
                rParticles%p_PartVolaver(ivt)= rParticles%p_PartVolaver(ivt)/rParticles%iPartVolCount
            end do
          
            call ucd_addVariableVertexBased (rexport, cvariable,&
                UCD_VAR_STANDARD, rParticles%p_PartVolaver)

            do ivt= 1, rproblemLevel%rtriangulation%NVT
                rParticles%p_PartVolaver(ivt)= rParticles%p_PartVolaver(ivt)*rParticles%iPartVolCount
            end do

          else

            ! Standard treatment for scalar quantity
            call eulerlagrange_getVarInterleaveFormat(rvector1%NEQ,  nvar,&
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
              call eulerlagrange_getVarBlockFormat(rvector1%NEQ, NVAR1D,&
                  'velocity_x', p_Dsolution, p_Ddata1)
              call ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1)
              
            case (NDIM2D)
              call eulerlagrange_getVarBlockFormat(rvector1%NEQ, NVAR2D,&
                  'velocity_x', p_Dsolution, p_Ddata1)
              call eulerlagrange_getVarBlockFormat(rvector2%NEQ, NVAR2D,&
                  'velocity_y', p_Dsolution, p_Ddata2)
              call ucd_addVarVertBasedVec(rexport, 'velocity',&
                  p_Ddata1, p_Ddata2)
                  
            case (NDIM3D)
              call eulerlagrange_getVarBlockFormat(rvector1%NEQ, NVAR3D,&
                  'velocity_x', p_Dsolution, p_Ddata1)
              call eulerlagrange_getVarBlockFormat(rvector2%NEQ, NVAR3D,&
                  'velocity_y', p_Dsolution, p_Ddata2)
              call eulerlagrange_getVarBlockFormat(rvector3%NEQ, NVAR3D,&
                  'velocity_z', p_Dsolution, p_Ddata3)
              call ucd_addVarVertBasedVec(rexport, 'velocity',&
                  p_Ddata1, p_Ddata2, p_Ddata3)
            end select
   
          elseif (trim(cvariable) .eq. 'velo_part') then

            ! Special treatment of velocity vector
            select case(ndim)
            case (NDIM1D)
              call eulerlagrange_getVarBlockFormat(rvector4%NEQ, NVAR2D,&
                  'velo_part_x', p_Dsolution, p_Ddata4)
              call ucd_addVarVertBasedVec(rexport, 'velo_part', p_Ddata4)   
                   
            case (NDIM2D)

              call ucd_addVarVertBasedVec(rexport, 'velo_part',&
                  rParticles%p_PartVelox, rParticles%p_PartVeloy)
                  
            case (NDIM3D)           
              call eulerlagrange_getVarBlockFormat(rvector4%NEQ, NVAR2D,&
                  'velo_part_x', p_Dsolution, p_Ddata4)
              call eulerlagrange_getVarBlockFormat(rvector5%NEQ, NVAR2D,&
                  'velo_part_y', p_Dsolution, p_Ddata5)
              call eulerlagrange_getVarBlockFormat(rvector6%NEQ, NVAR2D,&
                  'velo_part_z', p_Dsolution, p_Ddata6)
              call ucd_addVarVertBasedVec(rexport, 'velo_part',&
                  p_Ddata4, p_Ddata5, p_Ddata6)

            end select

          elseif (trim(cvariable) .eq. 'vol_part') then
         
              call ucd_addVariableVertexBased (rexport, cvariable,&
                UCD_VAR_STANDARD, rParticles%p_PartVol)
 
           elseif (trim(cvariable) .eq. 'vol_partaver') then
          
            call ucd_addVariableVertexBased (rexport, cvariable,&
                UCD_VAR_STANDARD, rParticles%p_PartVolAver)
               
          else
            
            ! Standard treatment for scalar quantity
            call eulerlagrange_getVarBlockFormat(rvector1%NEQ, nvar,&
                cvariable, p_Dsolution, p_Ddata1)
            call ucd_addVariableVertexBased (rexport, cvariable,&
                UCD_VAR_STANDARD, p_Ddata1)
            
          end if
        end do

      case default
        call output_line('Invalid system format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_outputSolution')
        call sys_halt()
      end select

      ! Release temporal memory
      select case(ndim)
      case (NDIM1D)
          call lsyssc_releaseVector(rvector1)
          call lsyssc_releaseVector(rvector4)
      case (NDIM2D)
          call lsyssc_releaseVector(rvector1)
          call lsyssc_releaseVector(rvector2)
          call lsyssc_releaseVector(rvector4)
          call lsyssc_releaseVector(rvector5)
      case (NDIM3D)           
          call lsyssc_releaseVector(rvector1)
          call lsyssc_releaseVector(rvector2)
          call lsyssc_releaseVector(rvector3)
          call lsyssc_releaseVector(rvector4)
          call lsyssc_releaseVector(rvector5)
          call lsyssc_releaseVector(rvector6)
      end select
      
    end if
    
    
    call ucd_settracers(rexport,reshape((/0.001_dp,0.002_dp,0.003_dp,0.004_dp/),(/2,2/)))

    ! Write UCD file
    call ucd_write  (rexport)
    call ucd_release(rexport)

    write(sfilename,'(i0)') rParticles%iTimestep
    
    sfilenamenew=trim(sucdsolution)//'particles.'//trim(sys_si0(rParticles%iTimestep,5))//'.vtk'
    
    open(20+rParticles%iTimestep,file=trim(sfilenamenew))

    select case(ipartoutput)
      case (0)
      
      case (1)
        ! Store only particle positions
        do iPart = 1, rParticles%nPart
          write(20+rParticles%iTimestep,*) rParticles%p_xpos(iPart), rParticles%p_ypos(iPart)
        end do
          
      case (2) 
        ! Store position, mass, density, diameter, temperature and velocity of the particles          
        write(20+rParticles%iTimestep,*) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
        write(20+rParticles%iTimestep,*) '  <UnstructuredGrid>'
        write(20+rParticles%iTimestep,*) '      <Piece NumberOfPoints="', rParticles%nPart, '" NumberOfCells="0">'
        write(20+rParticles%iTimestep,*) '          <Points>'
        write(20+rParticles%iTimestep,*) '<DataArray name="Position" type="Float32" NumberOfComponents="3" format="ascii">'
        do iPart = 1, rParticles%nPart
          write(20+rParticles%iTimestep,*) rParticles%p_xpos(iPart), rParticles%p_ypos(iPart), '0'
        end do
        write(20+rParticles%iTimestep,*) '</DataArray>'
        write(20+rParticles%iTimestep,*) '</Points>'
        write(20+rParticles%iTimestep,*) '<PointData  Vectors="vector">'
        write(20+rParticles%iTimestep,*) '<DataArray type="Float32" Name="Velocity" NumberOfComponents="3" format="ascii">'
        do iPart = 1, rParticles%nPart
          write(20+rParticles%iTimestep,*) rParticles%p_xvelo(iPart), rParticles%p_yvelo(iPart), '0'
        end do
        write(20+rParticles%iTimestep,*) '</DataArray>'
        write(20+rParticles%iTimestep,*) '<DataArray type="Float32" Name="Diameter" format="ascii">'
        do iPart = 1, rParticles%nPart
          write(20+rParticles%iTimestep,*) rParticles%p_diam(iPart)
        end do
        write(20+rParticles%iTimestep,*) '</DataArray>'
        write(20+rParticles%iTimestep,*) '<DataArray type="Float32" Name="Temperature" format="ascii">'
        do iPart = 1, rParticles%nPart
          write(20+rParticles%iTimestep,*) rParticles%p_temp(iPart)
        end do
        write(20+rParticles%iTimestep,*) '</DataArray>'
        write(20+rParticles%iTimestep,*) '<DataArray type="Float32" Name="Mass" format="ascii">'
        do iPart = 1, rParticles%nPart
          write(20+rParticles%iTimestep,*) rParticles%p_mass(iPart)
        end do
        write(20+rParticles%iTimestep,*) '</DataArray>'
        write(20+rParticles%iTimestep,*) '<DataArray type="Float32" Name="Density" format="ascii">'
        do iPart = 1, rParticles%nPart   
          write(20+rParticles%iTimestep,*) rParticles%p_density(iPart)
        end do
        write(20+rParticles%iTimestep,*) '</DataArray>'
        write(20+rParticles%iTimestep,*) '</PointData>'
        write(20+rParticles%iTimestep,*) '      <Cells>'
        write(20+rParticles%iTimestep,*) '          <DataArray type="Int32" Name="connectivity" format="ascii">'
        write(20+rParticles%iTimestep,*) '          </DataArray>'
        write(20+rParticles%iTimestep,*) '          <DataArray type="Int32" Name="offsets" format="ascii">'
        write(20+rParticles%iTimestep,*) '          </DataArray>'
        write(20+rParticles%iTimestep,*) '       <DataArray type="UInt8" Name="types" format="ascii">'
        write(20+rParticles%iTimestep,*) '       </DataArray>'
        write(20+rParticles%iTimestep,*) '     </Cells>'
        write(20+rParticles%iTimestep,*) '   </Piece>'
        write(20+rParticles%iTimestep,*) '</UnstructuredGrid>'
        write(20+rParticles%iTimestep,*) '</VTKFile>'

    end select


    close(unit=20+rParticles%iTimestep)
    
    rParticles%iTimestep=rParticles%iTimestep+1

  end subroutine eulerlagrange_outputSolution

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_outputStatistics(rtimerTotal, rcollection)

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
    type(t_timer), pointer :: p_rtimerParticlephase 
    real(DP) :: dtotalTime, dfraction

    
    ! Get timer objects from collection
    p_rtimerSolution => collct_getvalue_timer(rcollection, 'rtimerSolution')
    p_rtimerParticlephase => collct_getvalue_timer(rcollection, 'rtimerParticlephase')
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

    call output_line('Time for gas phase            : '//&
                     trim(adjustl(sys_sdE(p_rtimerSolution%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerSolution%delapsedCPU, 5)))//' %')
    call output_line('Time for particle phase       : '//&
                     trim(adjustl(sys_sdE(p_rtimerParticlephase%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerParticlephase%delapsedCPU, 5)))//' %')
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
    call output_line('Time for vector assembly      : '//&
                     trim(adjustl(sys_sdE(p_rtimerAssemblyVector%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerAssemblyVector%delapsedCPU, 5)))//' %')
    call output_line('Time for pre-/post-processing : '//&
                     trim(adjustl(sys_sdE(p_rtimerPrePostprocess%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerPrePostprocess%delapsedCPU, 5)))//' %')
    call output_lbrk()
    call output_line('Time for total simulation     : '//&
                     trim(adjustl(sys_sdE(dtotalTime, 5))))
    call output_lbrk()

  end subroutine eulerlagrange_outputStatistics

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_estimateRecoveryError(rparlist, ssectionName,&
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
      call eulerlagrange_getVariable(rsolution, serrorVariable, rvectorScalar)

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
            OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_estimateRecoveryError')
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

  end subroutine eulerlagrange_estimateRecoveryError

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_adaptTriangulation(rparlist, ssectionName, rhadapt ,&
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
        call eulerlagrange_hadaptCallbackScalar1D(&
            HADAPT_OPR_INITCALLBACK, rcollection)
        call hadapt_performAdaptation(rhadapt, rindicator,&
            rcollection, eulerlagrange_hadaptCallbackScalar1D)
        
      case (NDIM2D)
        call eulerlagrange_hadaptCallbackScalar2D(&
            HADAPT_OPR_INITCALLBACK, rcollection)
        call hadapt_performAdaptation(rhadapt, rindicator,&
            rcollection, eulerlagrange_hadaptCallbackScalar2D)
        
      case (NDIM3D)
        call eulerlagrange_hadaptCallbackScalar3D(&
            HADAPT_OPR_INITCALLBACK, rcollection)
        call hadapt_performAdaptation(rhadapt, rindicator,&
            rcollection, eulerlagrange_hadaptCallbackScalar3D)
      end select

    case (SYSTEM_BLOCKFORMAT)

      ! How many spatial dimensions do we have?
      select case(rtriangulationSrc%ndim)
      case (NDIM1D)
        call eulerlagrange_hadaptCallbackBlock1D(&
            HADAPT_OPR_INITCALLBACK, rcollection)
        call hadapt_performAdaptation(rhadapt, rindicator,&
            rcollection, eulerlagrange_hadaptCallbackBlock1D)
        
      case (NDIM2D)
        call eulerlagrange_hadaptCallbackBlock2D(&
            HADAPT_OPR_INITCALLBACK, rcollection)
        call hadapt_performAdaptation(rhadapt, rindicator,&
            rcollection, eulerlagrange_hadaptCallbackBlock2D)
        
      case (NDIM3D)
        call eulerlagrange_hadaptCallbackBlock3D(&
            HADAPT_OPR_INITCALLBACK, rcollection)
        call hadapt_performAdaptation(rhadapt, rindicator,&
            rcollection, eulerlagrange_hadaptCallbackBlock3D)
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

  end subroutine eulerlagrange_adaptTriangulation

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_solveTransientPrimal(rparlist, ssectionName,&
      rbdrCond, rproblem, rtimestep, rsolver, rsolution, rcollection,rParticles)

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
    
    type(t_Particles), intent(inout) :: rParticles 
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
    type(t_timer), pointer :: p_rtimerParticlephase
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

    integer :: icouplingpart

    ! Get timer structures
    p_rtimerPrePostprocess => collct_getvalue_timer(rcollection, 'rtimerPrePostprocess')
    p_rtimerSolution => collct_getvalue_timer(rcollection, 'rtimerSolution')
    p_rtimerParticlephase => collct_getvalue_timer(rcollection, 'rtimerParticlephase')
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
        eulerlagrange_getNVAR(p_rproblemLevel)) then
      rsolution%RvectorBlock(1)%NVAR = eulerlagrange_getNVAR(p_rproblemLevel)
      isize = rsolution%NEQ*eulerlagrange_getNVAR(p_rproblemLevel)
      call lsysbl_resizeVectorBlock(rsolution, isize, .false., .false.)
    end if
    
    ! Initialize the solution vector and impose boundary conditions
    call eulerlagrange_initSolution(rparlist, ssectionName, p_rproblemLevel,&
        rtimestep%dinitialTime, rsolution, rcollection)

    select case(ndimension)
    case (NDIM1D)
      call bdrf_filterVectorExplicit(rbdrCond, rsolution,&
          rtimestep%dinitialTime, eulerlagrange_calcBoundaryvalues1d)

    case (NDIM2D)
      call bdrf_filterVectorExplicit(rbdrCond, rsolution,&
          rtimestep%dinitialTime, eulerlagrange_calcBoundaryvalues2d)

    case (NDIM3D)
      call bdrf_filterVectorExplicit(rbdrCond, rsolution,&
          rtimestep%dinitialTime, eulerlagrange_calcBoundaryvalues3d)
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
            call eulerlagrange_estimateRecoveryError(rparlist, ssectionname,&
                p_rproblemLevel, rsolution, rtimestep%dinitialTime,&
                relementError, derror)

            ! Set the names of the template matrix
            rcollection%SquickAccess(1) = 'sparsitypattern'

            ! Attach the primal solution vector to the collection structure
            rcollection%p_rvectorQuickAccess1 => rsolution

            ! Perform h-adaptation and update the triangulation structure
            call eulerlagrange_adaptTriangulation(rparlist, ssectionname,&
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
            call eulerlagrange_initProblemLevel(rparlist,&
                ssectionName, p_rproblemLevel, rcollection)

            ! Resize the solution vector accordingly
            call parlst_getvalue_int(rparlist,&
                ssectionName, 'systemMatrix', systemMatrix)
            call lsysbl_resizeVecBlockIndMat(&
                p_rproblemLevel%RmatrixBlock(systemMatrix),&
                rsolution, .false., .true.)

            ! Re-generate the initial solution vector and impose
            ! boundary conditions explicitly
            call eulerlagrange_initSolution(rparlist, ssectionname,&
                p_rproblemLevel, rtimestep%dinitialTime, rsolution,&
                rcollection)

            select case(ndimension)
            case (NDIM1D)
              call bdrf_filterVectorExplicit(rbdrCond, rsolution,&

                  rtimestep%dinitialTime, eulerlagrange_calcBoundaryvalues1d)
            case (NDIM2D)
              call bdrf_filterVectorExplicit(rbdrCond, rsolution,&
                  rtimestep%dinitialTime, eulerlagrange_calcBoundaryvalues2d)

            case (NDIM3D)
              call bdrf_filterVectorExplicit(rbdrCond, rsolution,&
                  rtimestep%dinitialTime, eulerlagrange_calcBoundaryvalues3d)
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
    
    ! Do we have to read in a precomdputed solution?
    if (trim(sucdimport) .ne. '') then
      call ucd_readGMV(sucdimport, rimport, p_rproblemLevel%rtriangulation)
      call ucd_getSimulationTime(rimport, rtimestep%dinitialTime)
      call ucd_getSimulationTime(rimport, rtimestep%dTime)
      call eulerlagrange_setVariables(rimport, rsolution)
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
          call eulerlagrange_outputSolution(rparlist, ssectionName,&
          p_rproblemLevel, rParticles, rsolution, dtime=rtimestep%dTime)
      
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
            rsolution, eulerlagrange_nlsolverCallback, rcollection)

      case (TSTEP_THETA_SCHEME)
        
        ! Adopt two-level theta-scheme
        call tstep_performThetaStep(p_rproblemLevel, rtimestep,&
            rsolver, rsolution, eulerlagrange_nlsolverCallback, rcollection)
       
      case DEFAULT
        call output_line('Unsupported time-stepping algorithm!',&
            OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_solveTransientPrimal')
        call sys_halt()
      end select

      ! Perform linearised FEM-FCT post-processing
      call eulerlagrange_calcLinearisedFCT(rbdrCond, p_rproblemLevel,&
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
        call eulerlagrange_outputSolution(rparlist, ssectionName,&
            p_rproblemLevel, rParticles, rsolution, dtime=rtimestep%dTime)

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
        call eulerlagrange_estimateRecoveryError(rparlist, ssectionname,&
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
        call eulerlagrange_adaptTriangulation(rparlist, ssectionname,&
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
        call eulerlagrange_initProblemLevel(rparlist, ssectionName,&
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
      
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Euler-Lagrange (compute particle phase)
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call parlst_getvalue_int(rparlist, 'Eulerlagrange', "icouplingpart", icouplingpart)

      if (icouplingpart .ne. -1) then
      
          ! Start time measurement for solution procedure
          call stat_startTimer(p_rtimerParticlephase)
          
          ! Subroutine to compute the particle movement
          call eulerlagrange_step(rparlist,rproblem%p_rproblemLevelMax,rsolution,&
                rtimestep,rcollection,rParticles)
                
          ! Stop time measurement for solution procedure
          call stat_stopTimer(p_rtimerParticlephase)
      
      end if


    end do timeloop

    
    ! Release adaptation structure
    if (trim(adjustl(sadaptivityName)) .ne. '') then
      if ((dstepAdapt > 0.0_DP) .or. (npreadapt > 0)) then
        call hadapt_releaseAdaptation(rhadapt)
        call grph_releaseGraph(rgraph)
      end if
    end if
    
  end subroutine eulerlagrange_solveTransientPrimal

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_projectSolution(rsourceVector, rdestVector)
    
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
          rdestVector%RvectorBlock(iblock), eulerlagrange_coeffVectorFE, rcollection)

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
      
      print *, "Density mass", dmass1, dmass2
      
      ! Release matrices
      call lsyssc_releaseMatrix(rmatrix1)
      call lsyssc_releaseMatrix(rmatrix2)

    end do
    
    ! Release the collection structure
    call collct_done(rcollection)

  end subroutine eulerlagrange_projectSolution

  !*****************************************************************************
  ! AUXILIARY ROUTINES
  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_parseCmdlArguments(rparlist)

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

  end subroutine eulerlagrange_parseCmdlArguments

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_adjustParameterlist(rparlist, ssectionName)

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
    
  end subroutine eulerlagrange_adjustParameterlist

  !*****************************************************************************

!<subroutine>

subroutine eulerlagrange_init(rparlist,p_rproblemLevel,rsolution,rtimestep,rcollection,rParticles)


!<description>
    ! This subroutine initialise the particle data. 

!<input>
    ! Parameterlist
    type(t_parlist), intent(inout) :: rparlist
    
    ! Collection structure
    type(t_collection), intent(inout) :: rcollection

    ! Time-stepping structure
    type(t_timestep), intent(inout) :: rtimestep

    ! Primal solution vector
    type(t_vectorBlock), intent(inout), target :: rsolution

    ! Particles
    type(t_Particles), intent(inout) :: rParticles

    ! Local variables
    ! Pointer to the triangulation
    type(t_triangulation), pointer :: p_rtriangulation

    ! Pointer to the multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! pointer to the coordinates of the vertices
    !
    ! A list of all corner(!)-vertices of the elements in the triangulation.
    ! Handle to 
    !       p_RcornerCoordinates = array [1..ndim,1..NVT] of double
    ! with
    !   p_DvertexCoords(1,.) = X-coordinate.
    ! for 1D meshes,
    !   p_DvertexCoords(1,.) = X-coordinate.
    !   p_DvertexCoords(2,.) = Y-coordinate.
    ! for 2D meshes and
    !   p_DvertexCoords(1,.) = X-coordinate.
    !   p_DvertexCoords(2,.) = Y-coordinate.
    !   p_DvertexCoords(3,.) = Z-coordinate.
    ! for 3D meshes.
    ! This is a handle to the old DCORVG-array.
    !
    ! Note that the array may be longer than NVT in general!
    ! (May happen in case of a mesh hierarchy generated by a 2-level
    ! refinement, where the coordinates of the points on the
    ! coarser levels are contained in te coordinates of the
    ! finer levels.)
    ! In such a case, only the first NVT n-tuples in this array are valid!
    real(DP), dimension(:,:), pointer :: p_DvertexCoords

    ! pointer to the vertices adjacent to an element
    !
    ! Handle to h_IverticesAtElement=array [1..NVE,1..NEL] of integer
    ! For each element the node numbers of the corner-vertices
    ! in mathematically positive sense.
    ! On pure triangular meshes, there is NVE=3. On mixed or pure quad
    ! meshes, there is NVE=4. In this case, there is 
    ! IverticesAtElement(4,.)=0 for a triangle in a quad mesh.
    ! This is a handle to the old KVERT array.
    integer, dimension(:,:), pointer :: p_IverticesAtElement
 
    ! Local variables
    integer :: ivt, iPart, iel

    real(DP) :: random1, random2, random3

    ! Midpoints of the elements
    integer(I32) :: h_midpoints
    real(DP), dimension(:,:), pointer :: p_midpoints_el
    integer, dimension(2) :: md_el_length

    ! Startingpostions of the particles
    real(DP) :: partxmin, partxmax, partymin, partymax

    ! Velocity of the particles
    real(DP) :: velopartx, veloparty
    
    ! Variables for particle density, -diameter and -temperature
    real(DP) :: particledensity, particledensitymin, particledensitymax
    real(DP) :: particlediam, particlediammin, particlediammax
    real(DP) :: parttemp, parttempmin, parttempmax

    ! Gravity
    real(DP) :: gravityx, gravityy

    ! Quantity of particles
    integer :: nPart
  
    ! Boundarybehaviour
    integer :: boundbehav

    ! Mode for the startingpositions of the particles
    integer :: istartpos

    ! Mode for the density of the particles
    integer :: idensitypart

    ! Mode for the diameter of the particles
    integer :: idiampart

    ! Mode for the temperature of the particles
    integer :: itemppart
    
    ! Kinematic viscosity of the gas
    real(DP) :: gas_nu
  
    ! Variables for starting position from PGM-file
    integer, dimension(:,:), pointer :: p_Idata
    real(DP) :: x,y,xmin,ymin,xmax,ymax
    integer :: nvt,ix,iy
    type(t_pgm) :: rpgm
    real(DP), dimension(:), pointer :: p_Ddata
    character(LEN=SYS_STRLEN) :: ssolutionname

  
    ! Set pointer to triangulation
    p_rtriangulation => p_rproblemLevel%rtriangulation
  
    ! Get quantity of particles
    call parlst_getvalue_int(rparlist, 'Eulerlagrange', "nPart", nPart)

    ! Number of particles
    rParticles%npart = nPart
  
    ! Storage_new (scall, sname, isize, ctype, ihandle, cinitNewBlock)
    call storage_new ('euler_lagrange', 'Particle:element', rParticles%npart, ST_INT, rParticles%h_element, &
                            ST_NEWBLOCK_NOINIT)

    call storage_new ('euler_lagrange', 'Particle:lambda1', rParticles%npart, ST_DOUBLE, rParticles%h_lambda1, &
                            ST_NEWBLOCK_NOINIT)
    call storage_new ('euler_lagrange', 'Particle:lambda2', rParticles%npart, ST_DOUBLE, rParticles%h_lambda2, &
                            ST_NEWBLOCK_NOINIT)
    call storage_new ('euler_lagrange', 'Particle:lambda3', rParticles%npart, ST_DOUBLE, rParticles%h_lambda3, &
                            ST_NEWBLOCK_NOINIT)
    call storage_new ('euler_lagrange', 'Particle:lambda4', rParticles%npart, ST_DOUBLE, rParticles%h_lambda4, &
                            ST_NEWBLOCK_NOINIT)

    call storage_new ('euler_lagrange', 'Particle:diameter', rParticles%npart, ST_DOUBLE, rParticles%h_diam, &
                            ST_NEWBLOCK_NOINIT)
    call storage_new ('euler_lagrange', 'Particle:temperatur', rParticles%npart, ST_DOUBLE, rParticles%h_temp, &
                            ST_NEWBLOCK_NOINIT)
    call storage_new ('euler_lagrange', 'Particle:mass', rParticles%npart, ST_DOUBLE, rParticles%h_mass, &
                            ST_NEWBLOCK_NOINIT)
    call storage_new ('euler_lagrange', 'Particle:density', rParticles%npart, ST_DOUBLE, rParticles%h_density, &
                            ST_NEWBLOCK_NOINIT)
    call storage_new ('euler_lagrange', 'Particle:alpha_n', rParticles%npart, ST_DOUBLE, rParticles%h_alpha_n, &
                            ST_NEWBLOCK_NOINIT)

    call storage_new ('euler_lagrange', 'Particle:xpos', rParticles%npart, ST_DOUBLE, rParticles%h_xpos, &
                            ST_NEWBLOCK_NOINIT)
    call storage_new ('euler_lagrange', 'Particle:ypos', rParticles%npart, ST_DOUBLE, rParticles%h_ypos, &
                            ST_NEWBLOCK_NOINIT)
    call storage_new ('euler_lagrange', 'Particle:zpos', rParticles%npart, ST_DOUBLE, rParticles%h_zpos, &
                            ST_NEWBLOCK_NOINIT)

    call storage_new ('euler_lagrange', 'Particle:xpos_old', rParticles%npart, ST_DOUBLE, rParticles%h_xpos_old, &
                            ST_NEWBLOCK_NOINIT)
    call storage_new ('euler_lagrange', 'Particle:ypos_old', rParticles%npart, ST_DOUBLE, rParticles%h_ypos_old, &
                            ST_NEWBLOCK_NOINIT)
    call storage_new ('euler_lagrange', 'Particle:zpos_old', rParticles%npart, ST_DOUBLE, rParticles%h_zpos_old, &
                            ST_NEWBLOCK_NOINIT)

    call storage_new ('euler_lagrange', 'Particle:xvelo', rParticles%npart, ST_DOUBLE, rParticles%h_xvelo, &
                            ST_NEWBLOCK_NOINIT)
    call storage_new ('euler_lagrange', 'Particle:yvelo', rParticles%npart, ST_DOUBLE, rParticles%h_yvelo, &
                            ST_NEWBLOCK_NOINIT)
    call storage_new ('euler_lagrange', 'Particle:zvelo', rParticles%npart, ST_DOUBLE, rParticles%h_zvelo, &
                            ST_NEWBLOCK_NOINIT)

    call storage_new ('euler_lagrange', 'Particle:xvelo_old', rParticles%npart, ST_DOUBLE, rParticles%h_xvelo_old, &
                            ST_NEWBLOCK_NOINIT)
    call storage_new ('euler_lagrange', 'Particle:yvelo_old', rParticles%npart, ST_DOUBLE, rParticles%h_yvelo_old, &
                            ST_NEWBLOCK_NOINIT)
    call storage_new ('euler_lagrange', 'Particle:zvelo_old', rParticles%npart, ST_DOUBLE, rParticles%h_zvelo_old, &
                            ST_NEWBLOCK_NOINIT)

    call storage_new ('euler_lagrange', 'Particle:xvelo_gas', rParticles%npart, ST_DOUBLE, rParticles%h_xvelo_gas, &
                            ST_NEWBLOCK_NOINIT)
    call storage_new ('euler_lagrange', 'Particle:yvelo_gas', rParticles%npart, ST_DOUBLE, rParticles%h_yvelo_gas, &
                            ST_NEWBLOCK_NOINIT)
    call storage_new ('euler_lagrange', 'Particle:zvelo_gas', rParticles%npart, ST_DOUBLE, rParticles%h_zvelo_gas, &
                            ST_NEWBLOCK_NOINIT)

    call storage_new ('euler_lagrange', 'Particle:xvelo_gas_old', rParticles%npart, ST_DOUBLE, rParticles%h_xvelo_gas_old, &
                            ST_NEWBLOCK_NOINIT)
    call storage_new ('euler_lagrange', 'Particle:yvelo_gas_old', rParticles%npart, ST_DOUBLE, rParticles%h_yvelo_gas_old, &
                            ST_NEWBLOCK_NOINIT)
    call storage_new ('euler_lagrange', 'Particle:zvelo_gas_old', rParticles%npart, ST_DOUBLE, rParticles%h_zvelo_gas_old, &
                            ST_NEWBLOCK_NOINIT)
    call storage_new ('euler_lagrange', 'Particle:bdy_time', rParticles%npart, ST_DOUBLE, rParticles%h_bdy_time, &
                            ST_NEWBLOCK_NOINIT)

    md_el_length(1)=2
    md_el_length(2)=p_rtriangulation%NEL

    ! Midpoints of the elements
    call storage_new ('euler_lagrange', 'Elements:midpoints', md_el_length, ST_DOUBLE, rParticles%h_midpoints_el, &
                            ST_NEWBLOCK_NOINIT)
    ! Volumepart of the particles
    call storage_new ('euler_lagrange', 'Vertices:particlevolume', p_rtriangulation%NVT, ST_DOUBLE, rParticles%h_PartVol, &
                            ST_NEWBLOCK_NOINIT)
    call storage_new ('euler_lagrange', 'Vertices:particlevolumeaver', p_rtriangulation%NVT, ST_DOUBLE, rParticles%h_PartVolAver, &
                            ST_NEWBLOCK_NOINIT)
    ! Velocity of the particles
    call storage_new ('euler_lagrange', 'Vertices:particlevelocityx', p_rtriangulation%NEL, ST_DOUBLE, rParticles%h_PartVelox, &
                            ST_NEWBLOCK_NOINIT)
    call storage_new ('euler_lagrange', 'Vertices:particlevelocityy', p_rtriangulation%NEL, ST_DOUBLE, rParticles%h_PartVeloy, &
                            ST_NEWBLOCK_NOINIT)
   
    call storage_getbase_double (rParticles%h_xpos, rParticles%p_xpos)
    call storage_getbase_double (rParticles%h_ypos, rParticles%p_ypos)
    call storage_getbase_double (rParticles%h_zpos, rParticles%p_zpos)
    call storage_getbase_double (rParticles%h_xpos_old, rParticles%p_xpos_old)
    call storage_getbase_double (rParticles%h_ypos_old, rParticles%p_ypos_old)
    call storage_getbase_double (rParticles%h_zpos_old, rParticles%p_zpos_old)
    call storage_getbase_double (rParticles%h_xvelo, rParticles%p_xvelo)
    call storage_getbase_double (rParticles%h_yvelo, rParticles%p_yvelo)
    call storage_getbase_double (rParticles%h_zvelo, rParticles%p_zvelo)
    call storage_getbase_double (rParticles%h_xvelo_old, rParticles%p_xvelo_old)
    call storage_getbase_double (rParticles%h_yvelo_old, rParticles%p_yvelo_old)
    call storage_getbase_double (rParticles%h_zvelo_old, rParticles%p_zvelo_old)
    call storage_getbase_double (rParticles%h_xvelo_gas, rParticles%p_xvelo_gas)
    call storage_getbase_double (rParticles%h_yvelo_gas, rParticles%p_yvelo_gas)
    call storage_getbase_double (rParticles%h_zvelo_gas, rParticles%p_zvelo_gas)
    call storage_getbase_double (rParticles%h_xvelo_gas_old, rParticles%p_xvelo_gas_old)
    call storage_getbase_double (rParticles%h_yvelo_gas_old, rParticles%p_yvelo_gas_old)
    call storage_getbase_double (rParticles%h_zvelo_gas_old, rParticles%p_zvelo_gas_old)
    call storage_getbase_int (rParticles%h_element, rParticles%p_element)
    call storage_getbase_double (rParticles%h_lambda1, rParticles%p_lambda1)
    call storage_getbase_double (rParticles%h_lambda2, rParticles%p_lambda2)
    call storage_getbase_double (rParticles%h_lambda3, rParticles%p_lambda3)
    call storage_getbase_double (rParticles%h_lambda4, rParticles%p_lambda4)
    call storage_getbase_double (rParticles%h_diam, rParticles%p_diam)
    call storage_getbase_double (rParticles%h_density, rParticles%p_density)
    call storage_getbase_double (rParticles%h_mass, rParticles%p_mass)
    call storage_getbase_double (rParticles%h_temp, rParticles%p_temp)
    call storage_getbase_double (rParticles%h_alpha_n, rParticles%p_alpha_n)
    call storage_getbase_double (rParticles%h_bdy_time, rParticles%p_bdy_time)

    call storage_getbase_double (rParticles%h_PartVol, rParticles%p_PartVol)
    call storage_getbase_double (rParticles%h_PartVolAver, rParticles%p_PartVolAver)
    call storage_getbase_double (rParticles%h_PartVelox, rParticles%p_PartVelox)
    call storage_getbase_double (rParticles%h_PartVeloy, rParticles%p_PartVeloy)
    call storage_getbase_double2D (rParticles%h_midpoints_el, rParticles%p_midpoints_el)
  
    ! Set pointer to coordinate vector
    call storage_getbase_double2D(&
        p_rtriangulation%h_DvertexCoords, p_DvertexCoords)
 
    ! Set pointer to vertices at element
    call storage_getbase_int2D(&
        p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)

    ! Stores the midpoint for each element
    do iel=1,p_rtriangulation%NEL

      rParticles%p_midpoints_el(1,iel)= &
                                (p_DvertexCoords(1,p_IverticesAtElement(1,iel))+&
                                 p_DvertexCoords(1,p_IverticesAtElement(2,iel))+&
                                 p_DvertexCoords(1,p_IverticesAtElement(3,iel)))/3.0_dp

      rParticles%p_midpoints_el(2,iel)= &
                                (p_DvertexCoords(2,p_IverticesAtElement(1,iel))+&
                                 p_DvertexCoords(2,p_IverticesAtElement(2,iel))+&
                                 p_DvertexCoords(2,p_IverticesAtElement(3,iel)))/3.0_dp

    end do

    ! Get values for the startingpositions of the particles
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "xmin", partxmin)
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "xmax", partxmax)
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "ymin", partymin)
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "ymax", partymax)
 
    ! Get particlevelocity
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "velopartx", velopartx)
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "veloparty", veloparty)

    ! Get particle-density, -temp and -diameter
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "particledensity", particledensity)
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "particlediam", particlediam)
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "parttemp", parttemp)
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "particledensitymin", particledensitymin)
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "particlediammin", particlediammin)
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "parttempmin", parttempmin)
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "particledensitymax", particledensitymax)
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "particlediammax", particlediammax)
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "parttempmax", parttempmax)

    ! Get values for gravity
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "gravityx", gravityx)
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "gravityy", gravityy)

    ! Get value of kinematic viscosity
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "gas_nu", gas_nu)

    ! Get boundarybehaviour
    call parlst_getvalue_int(rparlist, 'Eulerlagrange', "boundbehav", boundbehav)

    ! Get variable for startingposition
    call parlst_getvalue_int(rparlist, 'Eulerlagrange', "startpos", istartpos)

    ! Get variable for density of the particles
    call parlst_getvalue_int(rparlist, 'Eulerlagrange', "idensitypart", idensitypart)

    ! Get variable for diameter of the particles
    call parlst_getvalue_int(rparlist, 'Eulerlagrange', "idiampart", idiampart)

    ! Get variable for temperature of the particles
    call parlst_getvalue_int(rparlist, 'Eulerlagrange', "itemppart", itemppart)

    ! Store particlesnumber, viscosity of the gas and gravity  
    rParticles%npart = nPart
    rParticles%nu_g= gas_nu
    rParticles%gravity(1)= gravityx
    rParticles%gravity(2)= gravityy
    rParticles%iTimestep= 1
    rParticles%maxvalx= maxval(p_DvertexCoords(1,:))
    rParticles%p_PartVol= 0.0_dp
    rParticles%p_PartVolAver= 0.0_dp
    rParticles%iPartVolCount= 1    
    rParticles%p_PartVelox= 0.0_dp
    rParticles%p_PartVeloy= 0.0_dp
    
    random1= 0.0_dp
    random2= 0.0_dp
    random3= 0.0_dp

    ! Set boundaryconditions for the particles
    select case(boundbehav)
    case (0)
        rParticles%tang_val= 1.0_dp
        rParticles%norm_val= -1.0_dp
    case (1)
        rParticles%tang_val= 1.0_dp
        rParticles%norm_val= 0.0_dp
    case default
      call output_line('Invalid boundaryconditions!', &
                       OU_CLASS_ERROR,OU_MODE_STD,'flagship_boundbehav')
      call sys_halt()
    end select

    ! Initialisation for starting position from PGM-file
    if (istartpos == 2) then
        ! Get global configuration from parameter list
        call parlst_getvalue_string(rparlist,&
              'Eulerlagrange', 'filestartpoints', ssolutionName)

        ! Initialize solution from portable graymap image
        call ppsol_readPGM(0, ssolutionName, rpgm)

        ! Set pointer for image data
        call storage_getbase_int2D(rpgm%h_Idata, p_Idata)
        
        ! Determine minimum/maximum values of array
        xmin = huge(DP); xmax = -huge(DP)
        ymin = huge(DP); ymax = -huge(DP)

        do ivt = 1, p_rtriangulation%nvt
            xmin = min(xmin, partxmin)
            xmax = max(xmax, partxmax)
            ymin = min(ymin, partymin)
            ymax = max(ymax, partymax)
        end do
        
        write(*,*) ''
        write(*,*) 'Initialization with graymap!   ...this can take a minute!'

    end if


    ! Initialize data for each particle
    do iPart=1,rParticles%npart
  
        select case(istartpos)
        case(0)
       
          !-------------------------------------------------------------------------
          ! Initialize the starting position by random numbers
          !-------------------------------------------------------------------------

  		  ! Get randomnumber
		  call random_number(random1)
		  call random_number(random2)
		  
          partxmin= minval(p_DvertexCoords(1,:))
          partxmax= minval(p_DvertexCoords(1,:))+&
                    0.2*(maxval(p_DvertexCoords(1,:))-minval(p_DvertexCoords(1,:)))
          partymin= minval(p_DvertexCoords(2,:))
          partymax= maxval(p_DvertexCoords(2,:))

          ! Set startingpositions of the particle
          rParticles%p_xpos(iPart)= partxmin + random1*(partxmax - partxmin)
          rParticles%p_ypos(iPart)= partymin + random2*(partymax - partymin)+(partymax - partymin)*0.001_dp
          rParticles%p_xpos_old(iPart)= partxmin + random1*(partxmax - partxmin)
          rParticles%p_ypos_old(iPart)= partymin + random2*(partymax - partymin)+(partymax - partymin)*0.001_dp


        case(1)
         
          !-------------------------------------------------------------------------
          ! Initialize the starting position by individual starting position
          !-------------------------------------------------------------------------

		  ! Get randomnumber
		  call random_number(random1)
		  call random_number(random2)
		  
          ! Set startingpositions of the particle
          rParticles%p_xpos(iPart)= partxmin + random1*(partxmax - partxmin)
          rParticles%p_ypos(iPart)= partymin + random2*(partymax - partymin)+(partymax - partymin)*0.001_dp
          rParticles%p_xpos_old(iPart)= partxmin + random1*(partxmax - partxmin)
          rParticles%p_ypos_old(iPart)= partymin + random2*(partymax - partymin)+(partymax - partymin)*0.001_dp
        
        case(2)
        
          !-------------------------------------------------------------------------
          ! Initialize the starting position by the data of a graymap image
          !-------------------------------------------------------------------------
  
         if (modulo(10*iPart,rParticles%npart)==0) write(*,*) int(100*iPart/rParticles%npart), '%'
       
         do
            ! Get random numbers
            call random_number(random1)
		    call random_number(random2)
            call random_number(random3)
         
	        partxmin= minval(p_DvertexCoords(1,:))
            partxmax= minval(p_DvertexCoords(1,:))+&
                     (maxval(p_DvertexCoords(1,:))-minval(p_DvertexCoords(1,:)))
            partymin= minval(p_DvertexCoords(2,:))
            partymax= maxval(p_DvertexCoords(2,:))
	    
		    ! Get point in the array
            rParticles%p_xpos(iPart)= partxmin + random1*(partxmax - partxmin)
            rParticles%p_ypos(iPart)= partymin + random2*(partymax - partymin)+(partymax - partymin)*0.001_dp

            ix = 1+(rpgm%width-1)*(rParticles%p_xpos(iPart)-xmin)/(xmax-xmin)
            if (ix .lt. 1 .or. ix .gt. rpgm%width) cycle

            iy = rpgm%height-(rpgm%height-1)*(rParticles%p_ypos(iPart)-ymin)/(ymax-ymin)
            if (iy .lt. 1 .or. iy .gt. rpgm%height) cycle

            ! If there can be particles, exit loop
            if (random3 .le. real(p_Idata(ix,iy),DP)/real(rpgm%maxgray,DP)) exit
          end do
               
          ! Set particle positions
          rParticles%p_xpos_old(iPart)= rParticles%p_xpos(iPart)
          rParticles%p_ypos_old(iPart)= rParticles%p_ypos(iPart)
        
        case(3)
       
          !-------------------------------------------------------------------------
          ! Initialize the starting position uniformly distributed over the domain
          !-------------------------------------------------------------------------

          ! Get random number for element
          call random_number(random1)
          
          ! Set starting element for the particle
          rParticles%p_element(iPart)= int(p_rtriangulation%NEL*random1 + 1)
          
          ! Get random number for barycentric coordinates
          call random_number(random1)
          call random_number(random2)

          ! Set barycentric coordinates          
          rParticles%p_lambda1(iPart)= random1
          rParticles%p_lambda2(iPart)= (1-rParticles%p_lambda1(iPart))*random2
          rParticles%p_lambda3(iPart)= 1-rParticles%p_lambda1(iPart)-rParticles%p_lambda2(iPart)
          
          ! Set startingpositions of the particle
          rParticles%p_xpos(iPart)= &
                 p_DvertexCoords(1,p_IverticesAtElement(1,rParticles%p_element(iPart)))*&
                 rParticles%p_lambda1(iPart)+&
                 p_DvertexCoords(1,p_IverticesAtElement(2,rParticles%p_element(iPart)))*&
                 rParticles%p_lambda2(iPart)+&
                 p_DvertexCoords(1,p_IverticesAtElement(3,rParticles%p_element(iPart)))*&
                 rParticles%p_lambda3(iPart)
          rParticles%p_ypos(iPart)= &
                 p_DvertexCoords(2,p_IverticesAtElement(1,rParticles%p_element(iPart)))*&
                 rParticles%p_lambda1(iPart)+&
                 p_DvertexCoords(2,p_IverticesAtElement(2,rParticles%p_element(iPart)))*&
                 rParticles%p_lambda2(iPart)+&
                 p_DvertexCoords(2,p_IverticesAtElement(3,rParticles%p_element(iPart)))*&
                 rParticles%p_lambda3(iPart)
          rParticles%p_xpos_old(iPart)= rParticles%p_xpos(iPart)
          rParticles%p_ypos_old(iPart)= rParticles%p_ypos(iPart)
          
        case default
            call output_line('Invalid starting position!', &
                       OU_CLASS_ERROR,OU_MODE_STD,'flagship_startpos')
            call sys_halt()
        end select
 
        ! Set diameter of the particles
        select case(idiampart)
        case (0)
            rParticles%p_diam(iPart)= particlediam
            
        case (1)
            ! Get random number
            call random_number(random1)
            
            rParticles%p_diam(iPart)= random1*particlediam
 
        case (2)
            ! Get random number
            call random_number(random1)
            
            rParticles%p_diam(iPart)= particlediammin+random1*(particlediammax-particlediammin)
                      
        case default
          call output_line('Invalid diam type mode!', &
                           OU_CLASS_ERROR,OU_MODE_STD,'flagship_diamtype')
          call sys_halt()
        end select

        ! Set density of the particles
        select case(idensitypart)
        case (0)
            rParticles%p_density(iPart)= particledensity
            
        case (1)
            ! Get random number
            call random_number(random2)
            
            rParticles%p_density(iPart)= random2*particledensity
 
        case (2)
            ! Get random number
            call random_number(random2)
            
            rParticles%p_density(iPart)= particledensitymin+random2*(particledensitymax-particledensitymin)
           
        case default
          call output_line('Invalid density type mode!', &
                           OU_CLASS_ERROR,OU_MODE_STD,'flagship_densitytype')
          call sys_halt()
        end select
 

        ! Set temperature of the particles
        select case(itemppart)
        case (0)
            rParticles%p_temp(iPart)= parttemp
            
        case (1)
            ! Get random number
            call random_number(random2)
            
            rParticles%p_temp(iPart)= random2*parttemp

        case (2)
            ! Get random number
            call random_number(random2)
            
            rParticles%p_temp(iPart)= parttempmin+random2*(parttempmax-parttempmin)
          
        case default
          call output_line('Invalid temp type mode!', &
                           OU_CLASS_ERROR,OU_MODE_STD,'flagship_temptype')
          call sys_halt()
        end select

        ! Set particle mass
        rParticles%p_mass(iPart)= &
               rParticles%p_density(iPart)*(rParticles%p_diam(iPart)**3 * 3.14159265358_dp /6.0_dp)

        ! Set initial values for the particles
        rParticles%p_xvelo(iPart)= velopartx
        rParticles%p_yvelo(iPart)= veloparty
        rParticles%p_xvelo_old(iPart)= velopartx
        rParticles%p_yvelo_old(iPart)= veloparty
        rParticles%p_xvelo_gas(iPart)= 0.0_dp
        rParticles%p_yvelo_gas(iPart)= 0.0_dp
        rParticles%p_xvelo_gas_old(iPart)= 0.0_dp
        rParticles%p_yvelo_gas_old(iPart)= 0.0_dp
        rParticles%p_alpha_n(iPart)= 0.0_dp
        rParticles%p_element(iPart)= 1
        rParticles%p_bdy_time(iPart)= 0.0_dp
                
        ! Find the start element for each particle
        call eulerlagrange_findelement(rparlist,p_rproblemLevel,rParticles,iPart)

        ! Calculate barycentric coordinates
        call eulerlagrange_calcbarycoords(p_rproblemLevel,rParticles,iPart)

        ! Wrong element
        if ((abs(rParticles%p_lambda1(iPart))+abs(rParticles%p_lambda2(iPart))+&
                  abs(rParticles%p_lambda3(iPart))-1) .GE. 0.00001) then
            call eulerlagrange_wrongelement(rparlist,p_rproblemLevel,rParticles,iPart)
        end if

    end do ! Loop over all particles
    
    if (istartpos == 2) then
        ! Release portable graymap image
        call ppsol_releasePGM(rpgm)
    end if
    
    ! Subroutine to calculate the volume part of the particles
    call eulerlagrange_calcvolpart(p_rproblemLevel,rParticles)

    ! Subroutine to calculate the velocity of the particles 
    !call eulerlagrange_calcvelopart(p_rproblemLevel,rParticles)


end subroutine eulerlagrange_init

  !*****************************************************************************

!<subroutine>

subroutine eulerlagrange_step(rparlist,p_rproblemLevel,rsolution,rtimestep,rcollection,rParticles)

!<description>
    ! This subroutine computes 

!<input>
    ! Parameterlist
    type(t_parlist), intent(inout) :: rparlist

    ! Collection structure
    type(t_collection), intent(inout) :: rcollection

    ! Time-stepping structure
    type(t_timestep), intent(inout) :: rtimestep

    ! Particles
    type(t_Particles), intent(inout) :: rParticles

    ! Primal solution vector
    type(t_vectorBlock), intent(inout), target :: rsolution

    ! Pointer to the multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel
    
 
    ! One-way or to-way coppling
    integer :: icouplingpart
    
    ! One way or twoway coupling?
    call parlst_getvalue_int(rparlist, 'Eulerlagrange', "icouplingpart", icouplingpart)

    select case(icouplingpart)
    case(1)
        ! Subroutine to compute the movement of the particles
        call eulerlagrange_moveparticlesoneway(rparlist,p_rproblemLevel,rsolution,rParticles)
    
    case(2)
        ! Subroutine to compute the movement of the particles
        call eulerlagrange_moveparticlestwoway(rparlist,p_rproblemLevel,rsolution,rParticles)

    case default
            call output_line('Invalid coupling mode!', &
                       OU_CLASS_ERROR,OU_MODE_STD,'flagship_coupling')
            call sys_halt()
     
    end select

    ! Subroutine to calculate the volume part of the particles
    call eulerlagrange_calcvolpart(p_rproblemLevel,rParticles)

    ! Subroutine to calculate the velocity of the particles 
    !call eulerlagrange_calcvelopart(p_rproblemLevel,rParticles)

end subroutine eulerlagrange_step

end module eulerlagrange_application
