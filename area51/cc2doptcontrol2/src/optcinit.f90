!##############################################################################
!# ****************************************************************************
!# <name> optcinit </name>
!# ****************************************************************************
!#
!# <purpose>
!# Contains initialisation and preprocessing routines for the optimal
!# control substructures.
!# </purpose>
!##############################################################################

module optcinit

  use fsystem
  use storage
  use basicgeometry
  use boundary
  use triangulation
  use spatialdiscretisation
  use paramlist
  
  use collection
    
  use basicstructures
  use spacediscretisation
  
  use forwardbackwardsimulation

  implicit none

contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine cc_initOptControl (rproblem)
  
!<description>
  ! Basic initialisation of the optimal control structure.
  ! Reads basic parameters from the DAT file to initialise the structure.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  ! The roptcontrol substructure is initialised.
  type(t_problem), intent(INOUT), target :: rproblem
!</inputoutput>

!</subroutine>

    ! Read in the parameters
    !
    ! Alpha/Gamma parameters
    call parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
        'dalphaC',rproblem%roptcontrol%dalphaC,1.0_DP)
        
    call parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
        'dgammaC',rproblem%roptcontrol%dgammaC,0.0_DP)
        
    ! Type of the formulation
    call parlst_getvalue_int (rproblem%rparamList,'OPTIMALCONTROL',&
        'ispaceTimeFormulation',rproblem%roptcontrol%ispaceTimeFormulation,0)
    
    call parlst_getvalue_int (rproblem%rparamList,'OPTIMALCONTROL',&
        'iconvectionExplicit',rproblem%roptcontrol%iconvectionExplicit,0)

    rproblem%roptcontrol%ccontrolConstraints = 0

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_doneOptControl (rproblem)
  
!<description>
  ! Cleans up the structure for the optimal control problem.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  ! The roptcontrol substructure is initialised.
  type(t_problem), intent(INOUT), target :: rproblem
!</inputoutput>

!</subroutine>

    rproblem%roptcontrol%dalphaC = 1.0_DP
    rproblem%roptcontrol%dgammaC = 0.0_DP

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_initTargetFlow (rproblem,roptcontrol,dstartTime,dendTime)
  
!<description>
  ! Reads in and sets up the target flow of the optimal control problem.
!</description>

!<input>
  ! OPTIONAL: Start time in a nonstationary simulation.
  real(DP), intent(IN) :: dstartTime

  ! OPTIONAL: Start time in a nonstationary simulation.
  real(DP), intent(IN) :: dendTime
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem

  ! Structure defining the parameters for the optimal control.
  type(t_problem_optcontrol), intent(inout) :: roptcontrol
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: itypeTargetFlow,iavaillevel
    type(t_collection) :: rcollection
    
    character(len=SYS_STRLEN), dimension(2) :: StargetFlowExpressions
    character(len=SYS_STRLEN) :: smeshTargetFunction,stargetFunction
    integer :: ilevelTargetFlow,ielementTypeTargetFlow,itargetFlowDelta
    integer :: itargetFlowTimesteps
    
    ! Type of the target flow?
    call parlst_getvalue_int (rproblem%rparamList,'OPTIMALCONTROL',&
        'itypeTargetFlow',itypeTargetFlow,0)

    ! Read some more parameters, we may need them.
    call parlst_getvalue_string (rproblem%rparamList,'OPTIMALCONTROL',&
        'stargetFlowExpressionX',StargetFlowExpressions(1),"",bdequote=.true.)

    call parlst_getvalue_string (rproblem%rparamList,'OPTIMALCONTROL',&
        'stargetFlowExpressionY',StargetFlowExpressions(2),"",bdequote=.true.)

    call parlst_getvalue_int (rproblem%rparamList,'OPTIMALCONTROL',&
        'ilevelTargetFlow',ilevelTargetFlow,0)

    call parlst_getvalue_int (rproblem%rparamList,'OPTIMALCONTROL',&
        'ielementTypeTargetFlow',ielementTypeTargetFlow,-1)

    call parlst_getvalue_string (rproblem%rparamList,'OPTIMALCONTROL',&
        'smeshTargetFlow',smeshTargetFunction,"",bdequote=.true.)

    call parlst_getvalue_string (rproblem%rparamList,'OPTIMALCONTROL',&
        'stargetFunction',stargetFunction,"",bdequote=.true.)

    call parlst_getvalue_int (rproblem%rparamList,'OPTIMALCONTROL',&
        'itargetFlowDelta',itargetFlowDelta,1)

    call parlst_getvalue_int (rproblem%rparamList,'OPTIMALCONTROL',&
        'itargetFlowTimesteps',itargetFlowTimesteps,-1)
    
    ! Set it up...
    select case (itypeTargetFlow)
    case (-1)
      ! Zero flow
      call ansol_init(roptcontrol%rtargetFunction)
    
    case (0)
      ! Analytic flow given by a callback routine
      call ansol_init(roptcontrol%rtargetFunction,0)
      
    case (1,3)
      ! Stationary targer flow
      if (ilevelTargetFlow .lt. rproblem%NLMIN) then
        call output_line ('Target flow level invalid!',&
            OU_CLASS_ERROR,OU_MODE_STD,'cc_initTargetFlow')
        call sys_halt()
      end if
      
      ! Initialise the flow. Try to use an existing discretisation,
      ! automatically create a new one if necessary.
      iavaillevel = min(rproblem%NLMAX,ilevelTargetFlow)
      rcollection%IquickAccess(1) = ielementTypeTargetFlow
      rcollection%IquickAccess(2) = NDIM2D+1
      call ansol_init (roptcontrol%rtargetFunction,ilevelTargetFlow,&
          rproblem%RlevelInfo(iavaillevel)%rdiscretisationPrimal,iavaillevel,&
          ielementTypeTargetFlow,&
          fget1LevelDiscretisationNavSt2D,rcollection)
    
      ! Read stationary flow from hard disc
      call ansol_configStationaryFile (roptcontrol%rtargetFunction,&
          stargetFunction,.true.)
      
    case (2,4)
      ! Nonstationary target flow
      if (ilevelTargetFlow .lt. rproblem%NLMIN) then
        call output_line ('Target flow level invalid!',&
            OU_CLASS_ERROR,OU_MODE_STD,'cc_initTargetFlow')
        call sys_halt()
      end if
      
      ! Initialise the flow. Try to use an existing discretisation,
      ! automatically create a new one if necessary.
      iavaillevel = min(rproblem%NLMAX,ilevelTargetFlow)
      rcollection%IquickAccess(1) = ielementTypeTargetFlow
      rcollection%IquickAccess(2) = NDIM2D+1
      call ansol_init (roptcontrol%rtargetFunction,ilevelTargetFlow,&
          rproblem%RlevelInfo(iavaillevel)%rdiscretisationPrimal,iavaillevel,&
          ielementTypeTargetFlow,&
          fget1LevelDiscretisationNavSt2D,rcollection)
    
      ! Read stationary flow from hard disc
      call ansol_configNonstationaryFile (roptcontrol%rtargetFunction, &
          dstartTime,dendTime,itargetFlowTimesteps,&
          '('''//trim(stargetFunction)//'.'',I5.5)',&
          0,itargetFlowDelta,.true.)
          
    case (5)
      ! Analytical expression as target flow.
      call ansol_init (roptcontrol%rtargetFunction)
      call ansol_config (roptcontrol%rtargetFunction,StargetFlowExpressions)

    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_doneTargetFlow (roptcontrol)
  
!<description>
  ! Cleans up the structure for the optimal control problem.
!</description>

!<inputoutput>
  ! Structure defining the parameters for the optimal control.
  type(t_problem_optcontrol), intent(inout) :: roptcontrol
!</inputoutput>

!</subroutine>

    ! Release the target flow.
    call ansol_done(roptcontrol%rtargetFunction)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_initConstraints (rproblem,roptcontrol)
  
!<description>
  ! Initialises constraints in the optimal control.
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</input>

!<inputoutput>
  ! Structure defining the parameters for the optimal control.
  type(t_problem_optcontrol), intent(inout) :: roptcontrol
!</inputoutput>

!</subroutine>
    
    call parlst_getvalue_int (rproblem%rparamList,'OPTIMALCONTROL',&
        'ccontrolConstraints',roptcontrol%ccontrolConstraints,0)

    call parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
        'dumin1',roptcontrol%dumin1,-1.0E10_DP)

    call parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
        'dumax1',roptcontrol%dumax1,1.0E10_DP)

    call parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
        'dumin2',roptcontrol%dumin2,-1.0E10_DP)

    call parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
        'dumax2',roptcontrol%dumax2,1.0E10_DP)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_doneConstraints (roptcontrol)
  
!<description>
  ! Cleans up information about constraints
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  ! The target flow is saved to the roptcontrol substructure.
  type(t_problem_optcontrol), intent(inout) :: roptcontrol
!</inputoutput>

!</subroutine>

    ! Nothing to do here for now.
    ! This may change later if the constraints are given nonstationary...

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine cc_initInitialCondition (rproblem,roptcontrol)
  
!<description>
  ! Initialises the initial solution vector into rvector. Depending on the settings
  ! in the DAT file this is either zero or read from a file.
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(IN) :: rproblem
!</input>

!<inputoutput>
  ! Structure defining the parameters for the optimal control.
  type(t_problem_optcontrol), intent(inout) :: roptcontrol
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ctypeInitCond,iavaillevel
    type(t_collection) :: rcollection
    
    character(len=SYS_STRLEN), dimension(3) :: Sexpressions
    character(len=SYS_STRLEN) :: smeshInitCond,sfileInitialCond
    integer :: ilevelInitCond,ielementTypeInitCond
    
    ! Type of the initial condition?
    call parlst_getvalue_int (rproblem%rparamList,'OPTIMALCONTROL',&
        'ctypeInitCond',ctypeInitCond,0)

    ! Read some more parameters, we may need them.
    call parlst_getvalue_string (rproblem%rparamList,'OPTIMALCONTROL',&
        'ssolutionInitialCondY1',Sexpressions(1),"0",bdequote=.true.)

    call parlst_getvalue_string (rproblem%rparamList,'OPTIMALCONTROL',&
        'ssolutionInitialCondY2',Sexpressions(2),"0",bdequote=.true.)

    call parlst_getvalue_string (rproblem%rparamList,'OPTIMALCONTROL',&
        'ssolutionInitialCondP',Sexpressions(3),"0",bdequote=.true.)

    call parlst_getvalue_int (rproblem%rparamList,'OPTIMALCONTROL',&
        'ilevelInitCond',ilevelInitCond,0)

    call parlst_getvalue_int (rproblem%rparamList,'OPTIMALCONTROL',&
        'ielementTypeInitCond',ielementTypeInitCond,-1)

    call parlst_getvalue_string (rproblem%rparamList,'OPTIMALCONTROL',&
        'smeshInitCond',smeshInitCond,"",bdequote=.true.)

    call parlst_getvalue_string (rproblem%rparamList,'OPTIMALCONTROL',&
        'sfileInitialCond',sfileInitialCond,"",bdequote=.true.)

    ! Set it up...
    select case (ctypeInitCond)
    case (0)
      ! Zero flow
      call ansol_init(roptcontrol%rinitialCondition)
    
    case (1,2)
      ! Initial condition given as file.
      if (ilevelInitCond .lt. rproblem%NLMIN) then
        call output_line ('Level of initial condition invalid!',&
            OU_CLASS_ERROR,OU_MODE_STD,'cc_initInitialSolution')
        call sys_halt()
      end if
      
      ! Initialise the flow. Try to use an existing discretisation,
      ! automatically create a new one if necessary.
      iavaillevel = min(rproblem%NLMAX,ilevelInitCond)
      rcollection%IquickAccess(1) = ielementTypeInitCond
      rcollection%IquickAccess(2) = NDIM2D+1
      call ansol_init (roptcontrol%rinitialCondition,ilevelInitCond,&
          rproblem%RlevelInfo(iavaillevel)%rdiscretisationPrimal,iavaillevel,&
          ielementTypeInitCond,fget1LevelDiscretisationNavSt2D,rcollection)
    
      ! Read stationary flow from hard disc
      call ansol_configStationaryFile (roptcontrol%rinitialCondition,&
          sfileInitialCond,ctypeInitCond .eq. 1)
      
    case (3)
      ! Analytical expression as initial condition
      call ansol_init (roptcontrol%rinitialCondition)
      call ansol_config (roptcontrol%rinitialCondition,Sexpressions)

    end select
    
  end subroutine

  ! ***************************************************************************
  !<subroutine>

    subroutine fgetInitCond (rdiscretisation, rform, &
                  nelements, npointsPerElement, Dpoints, &
                  IdofsTest, rdomainIntSubset, &
                  Dcoefficients, rcollection)
    
    use fsystem
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
  !<description>
    ! Auxiliary routine.
    ! Callback routine that calculates the values of the initial condition.
    ! rcollection%IquickAccess(1) determines the component.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN) :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(IN) :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN) :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN) :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN) :: Dpoints

    ! An array accepting the DOF`s on all elements test in the test space.
    ! DIMENSION(\#local DOF`s in test space,Number of elements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN) :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT), optional :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
      integer :: icomp
      real(dp), dimension(:,:), allocatable :: Dvalues
      
      ! Get the component.
      icomp = rcollection%IquickAccess(1)
      
      ! Evaluate the initial condition.
      allocate(Dvalues(ubound(Dcoefficients,2),ubound(Dcoefficients,3)))
      call ansol_evaluate (rcollection,"INITCOND",icomp,rdiscretisation%p_rtriangulation,&
          npointsPerElement,nelements,Dpoints,rdomainIntSubset%p_Ielements,Dvalues)
      Dcoefficients(1,:,:) = Dvalues(:,:)
      deallocate(Dvalues)
  
    end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_projectInitialCondition (rproblem,roptcontrol,rvector)
  
!<description>
  ! Projects the analytically initial condition in roptcontrol to the discrete
  ! vector rvector using a divergence-free projection based on the current time
  ! step scheme in rproblem.
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT) :: rproblem

  ! Structure defining the parameters for the optimal control.
  type(t_problem_optcontrol), intent(inout) :: roptcontrol
!</input>

!<inputoutput>
  ! Discrete initial condition. The vector must have been initialised
  ! based on the underlying FEM space.
  type(t_vectorBlock), intent(inout), target :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_simSolver) :: rsimsolver
    type(t_ccoptSpaceTimeMatrix) :: rspaceTimeMatrix
    type(t_ccoptSpaceTimeDiscretisation), target :: rspaceTimeDiscr
    type(t_spacetimeVector), target :: rb,rx
    logical :: bsuccess
    integer :: timenlmax

    ! At first clear the vector.
    call lsysbl_clearVector (rvector)
    
    ! Step 1: L2-projection
    ! ---------------------
    ! In the first step, we project the analytically given initial solution
    ! into our FEM space by a standard L2-projection. This will usually lead
    ! to an improper solution which is not divergence-free, what has to be corrected
    ! in a second step.
    !
    ! We do the L2 projection as accurate as possible.
    !
    ! Prepare the evaluation of the initial condition.
    call ansol_prepareEval (roptcontrol%rinitialCondition,&
        rproblem%rcollection,"INITCOND",rproblem%rtimedependence%dtimeInit)
        
    ! Do the projection; use fgetInitCond to return function values.
    !
    ! X-velocity
    rproblem%rcollection%IquickAccess(1) = 1
    call anprj_analytL2projectionByMass (rvector%RvectorBlock(1),&
        rproblem%RlevelInfo(rproblem%nlmax)%rstaticInfo%rmatrixMass,&
        fgetInitCond,rproblem%rcollection)

    ! Y-velocity
    rproblem%rcollection%IquickAccess(1) = 2
    call anprj_analytL2projectionByMass (rvector%RvectorBlock(2),&
        rproblem%RlevelInfo(rproblem%nlmax)%rstaticInfo%rmatrixMass,&
        fgetInitCond,rproblem%rcollection)

    ! Pressure
    rproblem%rcollection%IquickAccess(1) = 3
    call anprj_analytL2projectionByMass (rvector%RvectorBlock(3),&
        rproblem%RlevelInfo(rproblem%nlmax)%rstaticInfo%rmatrixMassPressure,&
        fgetInitCond,rproblem%rcollection)
        
    ! Cleanup
    call ansol_doneEval ("INITCOND",rproblem%rcollection)

    ! Step 2: Create a RHS for the divergence-free projection
    ! -------------------------------------------------------
    ! In this step, we multiply our vector with the operator of the 0th timestep.
    ! This results in a RHS which we later project back to the actual
    ! space.
    
    ! Create a space-time matrix that resembles a stanard forward simulation
    call parlst_getvalue_int (rproblem%rparamList,'TIME-DISCRETISATION',&
                              'TIMENLMAX',timenlmax,1)

    call sptidis_initDiscretisation (rproblem,timenlmax,rproblem%nlmax,rspaceTimeDiscr)
      
    ! Create a sparse space-time vector only holding the 0th rhs and solution.
    call sptivec_initVectorDiscr (rb,rspaceTimeDiscr%rtimeDiscr,&
        rproblem%RlevelInfo(rproblem%nlmax)%rdiscretisationPrimal,1,1)
    call sptivec_initVectorDiscr (rx,rspaceTimeDiscr%rtimeDiscr,&
        rproblem%RlevelInfo(rproblem%nlmax)%rdiscretisationPrimal,1,1)

    ! Put our improper solution vector into rx(0).
    call sptivec_setTimestepData (rx, 1, rvector)

    ! Prepare a space-time matrix
    rspaceTimeMatrix%cmatrixType = 0
    rspaceTimeMatrix%ccontrolConstraints = 0
    rspaceTimeMatrix%p_rspaceTimeDiscr => rspaceTimeDiscr
    rspaceTimeMatrix%p_rsolution => rx
    
    ! Multiply the vector to create the RHS. Note that we created the vector
    ! sparse, containing only the 0th timestep. The multiplication will give
    ! us rhs(0).
    call stlin_spaceTimeMatVec (rproblem, rspaceTimeMatrix, rx, rb, 1.0_DP, 0.0_DP, &
      SPTID_FILTER_NONE)
      
    ! Step 3: Projection into our real space
    ! --------------------------------------
    ! In the 3rd step, we solve a nonlinear system for the 0th timestep
    ! to calculate the actual (divergence-free) initial condition.

    ! Initialise a forward-backward solver for a standard forward simulation
    call fbsim_init (rproblem, rproblem%nlmin, rproblem%nlmax, &
        FBSIM_SOLVER_NLFORWARD, rsimsolver)

    ! Call the solver to calculate the initial condition -- 0th timestep.
    call fbsim_simulate (rsimsolver, rx, rb, 1, 0, rx, bsuccess)

    ! Release the solver
    call fbsim_done (rsimsolver)

    ! Release the discretisation
    call sptidis_doneDiscretisation(rspaceTimeDiscr)

    ! Get back the initial condition
    call sptivec_getTimestepData (rx, 1, rvector)
    
    ! Release rx and rb, that's it.
    call sptivec_releaseVector (rb)
    call sptivec_releaseVector (rx)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_doneInitialCondition (roptcontrol)
  
!<description>
  ! Cleans up the initial condition.
!</description>

!<inputoutput>
  ! Structure defining the parameters for the optimal control.
  type(t_problem_optcontrol), intent(inout) :: roptcontrol
!</inputoutput>

!</subroutine>

    ! Release the initial condition.
    call ansol_done(roptcontrol%rinitialCondition)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_initStartVector (rproblem,roptcontrol,rvector)
  
!<description>
  ! Initialises the space-time vector rvector to be the start vector
  ! for the iteration. Evaluated parameters in the DAT file to create
  ! the vector rvector.
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(in) :: rproblem

  ! Structure defining the parameters for the optimal control.
  type(t_problem_optcontrol), intent(in) :: roptcontrol
!</input>

!<inputoutput>
  ! Space-time vector thaht receives the start vector of the iteration.
  type(t_spacetimeVector), intent(out) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ctypeStartVector

    ! How should the start vector be created?
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
        'ctypeStartVector',ctypeStartVector,0)
    
    select case (ctypeStartVector)
    case (0)
      ! Zero vector -- ok, that's easy.
      call sptivec_clear(rvector)
    end select

  end subroutine

end module
