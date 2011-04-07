!##############################################################################
!# ****************************************************************************
!# <name> user_callback </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains callback functions for the cc2d-mini problem that are
!# used during the matrix/vector assembly for specifying analytical data.
!# There are three callback functions involved, which may be called depending
!# on the situation. All of them correspond to a specific interface for
!# callback functions, defined in 'intf_xxxx.inc' files.
!#
!# 1.) coeff_Stokes
!#     -> Returns the coefficients for the Laplace matrix. This routine is
!#        only used if the problem to calculate has nonconstant coefficients!
!#        Otherwise the routine is dead.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientMatrixSc.inc'
!#
!# 1.) coeff_Pressure
!#     -> Returns the coefficients for the pressure matrix. This routine is
!#        only used if the problem to calculate has nonconstant coefficients!
!#        Otherwise the routine is dead.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientMatrixSc.inc'
!#
!# 2.) coeff_RHS_x
!#     -> Returns analytical values for the right hand side of the X-velocity
!#        equation.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 3.) coeff_RHS_y
!#     -> Returns analytical values for the right hand side of the Y-velocity
!#        equation.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 4.) coeff_TARGET_x
!#     -> Returns analytical values for the desired flow field in X-direction.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 5.) coeff_TARGET_y
!#     -> Returns analytical values for the desired flow field in Y-direction.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 6.) getBoundaryValues
!#     -> Returns analitical values on the boundary of the
!#        problem to solve.
!#     -> Corresponds to the interface defined in the file
!#        'intf_bcassembly.inc'
!#
!# 7.) getBoundaryValuesFBC
!#     -> Returns analytical values on the fictitious boundary components
!#     -> Corresponds to the interface defined in the file
!#        'intf_fbcassembly.inc'
!#
!# 8.) cc_initCollectForAssembly
!#     -> Is called prior to the assembly process.
!#     -> Stores some information from the problem structure to a collection
!#        such that it can be accessed in callback routines
!#
!# 9.) cc_doneCollectForAssembly
!#     -> Is called after the assembly process.
!#     -> Releases information stored in the collection by 
!#        cc_initCollectForAssembly.
!#
!# For nonstationary simulation, it might be neccessary in these routines
!# to access the current simulation time. Before the assembly process, the cc2d
!# framework calls cc_initCollectForAssembly to stores the current point 
!# in time (and probably other necessary information) to the quickaccess-array
!# in the collection which is passed to the callback routines. The callback
!# routines can access this as follows:
!#
!# -> rcollection%IquickAccess(1)   = 0: stationary, 
!#                                    1: nonstationary with explicit time stepping
!# -> rcollection%DquickAccess(1)   = current simulation time
!# -> rcollection%DquickAccess(2)   = minimum simulation time
!# -> rcollection%DquickAccess(3)   = maximum simulation time
!#
!# After the assembly, cc_doneCollectForAssembly is called to clean up.
!# Note: Information stored in the quick-access array are of temporary
!# nature and does not have to be cleaned up.
!#
!# </purpose>
!##############################################################################

module user_callback

  use fsystem
  use storage
  use boundary
  use cubature
  use derivatives
  use matrixfilters
  use vectorfilters
  use discretebc
  use bcassembly
  use mprimitives
  use feevaluation
  use spacetimevectors
  use triasearch
  use bilinearformevaluation
  use linearformevaluation
  use linearsolver
  use fparser
  use element
  use elementpreprocessing
  
  use timeevaluation
  
  use basicstructures
  
  implicit none

contains

! ***************************************************************************

!<subroutine>

  subroutine cc_initCollectForAssembly (rproblem,dtime,rcollection)
  
!<description>
  ! This subroutine is an auxiliary subroutine called by the CC2D framework
  ! and has usually not to be changed by the user.
  !
  ! The subroutine prepares the collection rcollection to be passed to callback
  ! routines for assembling boundary conditions or RHS vectors. It is
  ! called directly prior to the assembly to store problem-specific information
  ! in the quick-access arrays of the collection.
  ! Basically speaking, the routine stores information about whether thw problem
  ! is stationary, nonstationary or about the current simulation time.
!</description>

!<input>
  ! Problem structure with all problem relevant data.
  type(t_problem), intent(INOUT) :: rproblem
  
  ! Current simulation time.
  real(dp), intent(in) :: dtime
!</input>

!<inputoutput>
  ! Collection structure to be initialised.
  type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
  
!</subroutine>
    real(DP) :: dreltime

    ! In a nonstationary simulation, save the simulation time as well as the
    ! minimum and maximum time to the quick-access array of the collection,
    ! so it can be accessed in the callback routines!
    rcollection%Iquickaccess(1) = rproblem%itimedependence
    select case (rproblem%itimedependence)
    case (0)
      ! Stationary simulation
      rcollection%Dquickaccess(1) = 0.0_DP
      rcollection%Dquickaccess(2) = 0.0_DP
      rcollection%Dquickaccess(3) = 0.0_DP
    case (1)
      rcollection%Dquickaccess(1) = dtime
      rcollection%Dquickaccess(2) = rproblem%rtimedependence%dtimeInit
      rcollection%Dquickaccess(3) = rproblem%rtimedependence%dtimeMax
    end select
    
    rcollection%Iquickaccess(2) = rproblem%roptcontrol%itypeTargetFlow
    call collct_setvalue_vec (rcollection, 'TARGETFLOW', &
        rproblem%roptcontrol%rtargetFlow, .true.) 

    call collct_setvalue_pars (rcollection, 'TARGETFLOWPARSER', &
        rproblem%roptcontrol%rparserTargetFlowExpression, .true.) 
    
    ! In case the target vector changes in time, load the current target
    ! from the space-time vector.    
    if ((rproblem%roptcontrol%itypeTargetFlow .eq. 2) .or. &
        (rproblem%roptcontrol%itypeTargetFlow .eq. 4))  then
      if (rproblem%roptcontrol%rtargetFlowNonstat%NEQtime .gt. 1) then
        
        !CALL sptivec_getTimestepData (rproblem%roptcontrol%rtargetFlowNonstat, &
        !    rproblem%rtimedependence%itimeStep, rproblem%roptcontrol%rtargetFlow)
        
        ! Get the timestep based on the time stamp. If necessary, the routine will
        ! interpolate between timesteps. drelTime is the 'relative' time in the range
        ! 0..1 with 0 corresponding to dtimeInit and 1 corresponding to dtimeMax.
        !dreltime = (rproblem%rtimedependence%dtime - &
        !            rproblem%rtimedependence%dtimeInit) / &
        !           (rproblem%rtimedependence%dtimeMax - &
        !            rproblem%rtimedependence%dtimeInit) 
        !CALL sptivec_getTimestepDataByTime (rproblem%roptcontrol%rtargetFlowNonstat, &
        !    dreltime, rproblem%roptcontrol%rtargetFlow)
        
        ! New implementation: Use tmevl_evaluate!
        
        call tmevl_evaluate(rproblem%roptcontrol%rtargetFlowNonstat,&
            dtime,rproblem%roptcontrol%rtargetFlow)
            
      end if
      ! Otherwise, there is no vector.
    end if

    ! Assembly data for the RHS.
    rcollection%Iquickaccess(3) = rproblem%irhs
    call collct_setvalue_pars (rcollection, 'RHSPARSER', &
        rproblem%rrhsParser, .true.) 

  end subroutine
  
! ***************************************************************************
  
!<subroutine>

  subroutine cc_doneCollectForAssembly (rproblem,rcollection)
  
!<description>
  ! This subroutine is an auxiliary subroutine called by the CC2D framework
  ! and has usually not to be changed by the user.
  !
  ! After the assembly process, this subroutine is called to release temporary
  ! information from the collection which was stored there by 
  ! cc_initCollectForAssembly.
!</description>
  
!<input>
  ! Problem structure with all problem relevant data.
  type(t_problem), intent(IN) :: rproblem
!</input>

!<inputoutput>
  ! Collection structure to be cleaned up.
  type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
  
!</subroutine>

    ! Currently, this subroutine is empty as all information stored in
    ! the collection in cc_initCollectForAssembly is put to the quick-access
    ! arrays -- which do not have to be cleaned up. 
    ! This might change in future...

    call collct_deletevalue (rcollection, 'TARGETFLOW') 
    call collct_deletevalue (rcollection, 'TARGETFLOWPARSER') 
    call collct_deletevalue (rcollection, 'RHSPARSER') 

  end subroutine
  
! ***************************************************************************
  !<subroutine>

  subroutine coeff_Stokes (rdiscretisationTrial,rdiscretisationTest,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial,IdofsTest,rdomainIntSubset,&
                  Dcoefficients, rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the matrix assembly. It has to compute
    ! the coefficients in front of the terms of the bilinear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the bilinear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the bilinear form
    ! the corresponding coefficients in front of the terms.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; trial space.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTrial
    
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTest
    
    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(IN)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in trial space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTrial
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    ! Routine is not called anyway in case of constant coefficients!
    Dcoefficients = 1.0_DP

  end subroutine

! ***************************************************************************
  !<subroutine>

  subroutine coeff_Pressure (rdiscretisationTrial,rdiscretisationTest,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial,IdofsTest,rdomainIntSubset, &
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the matrix assembly. It has to compute
    ! the coefficients in front of the terms of the bilinear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the bilinear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the bilinear form
    ! the corresponding coefficients in front of the terms.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; trial space.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTrial
    
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTest
    
    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(IN)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in trial space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTrial
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    ! Routine is not called anyway in case of constant coefficients!
    Dcoefficients = 1.0_DP

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_x (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset, &
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form of the
    ! X-velocity part of the right hand side vector.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(IN)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    real(DP) :: dtime
    integer :: irhs,i,j
    type(t_fparser), pointer :: p_rparser
    real(dp), dimension(:,:), allocatable :: p_Dval
    real(DP), dimension(:), allocatable :: DvaluesAct
    
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      dtime = rcollection%Dquickaccess(1)
    else
      dtime = 0.0_DP
    end if

    ! Type of RHS.
    irhs = rcollection%IquickAccess(3)
    
    ! Depending on the irhs identifier, choose the RHS to evaluate
    select case (irhs)
    case (-1)
      ! Given as analytical expression. Evaluate the expression.
      !
      ! Get the parser object with the RHS expressions from the collection
      p_rparser => collct_getvalue_pars (rcollection, 'RHSPARSER') 
      
      ! Prepare the array with the values for the function.
      ! X-coordinate, Y-coordinate, time.
      allocate(p_Dval(3,npointsPerElement*nelements))
      do i=1,nelements
        do j=1,npointsPerElement
          p_Dval (1,(i-1)*npointsPerElement+j) = Dpoints(1,j,i)
          p_Dval (2,(i-1)*npointsPerElement+j) = Dpoints(2,j,i)
          p_Dval (3,(i-1)*npointsPerElement+j) = dtime
        end do
      end do
      
      ! Evaluate the 1st expression for the X-rhs
      allocate(DvaluesAct(npointsPerElement*nelements))
      call fparser_evalFunction (p_rparser, 1, 2, p_Dval, DvaluesAct)

      ! Reshape the data, that's it.
      do i=0,nelements-1
        do j=1,npointsPerElement
          Dcoefficients(1,j,i+1) = DvaluesAct(i*npointsPerElement+j)
        end do
      end do
      
      deallocate(DvaluesAct)
      deallocate(p_Dval)
      
    case (0)
      ! Analytically given. =0.
      Dcoefficients(:,:,:) = 0.0_DP
      
    ! Other cases may be defined here.
      !Dcoefficients(1,:,:) = Dpoints(1,:,:)
      !Dcoefficients(1,:,:) = -18.0*sin(3.0*SYS_PI*Dpoints(1,:,:))*SYS_PI**2 &
      !                     *sin(3.0*SYS_PI*Dpoints(2,:,:)) &
      !                     + .5*SYS_PI*cos(.5*SYS_PI*(Dpoints(1,:,:)-Dpoints(2,:,:)))
      !Dcoefficients (1,:,:) = -(1./10.0_DP)*(-Dpoints(1,:,:))
      
      ! Without coupling:
      !Dcoefficients (1,:,:) = (1/5.0_DP - dtime/50.0_DP)*(Dpoints(1,:,:))
      !Dcoefficients (1,:,:) = Dpoints(1,:,:)
    end select
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_y (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset, &
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form of the
    ! Y-velocity part of the right hand side vector.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(IN)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    real(DP) :: dtime
    integer :: irhs,i,j
    type(t_fparser), pointer :: p_rparser
    real(dp), dimension(:,:), allocatable :: p_Dval
    real(DP), dimension(:), allocatable :: DvaluesAct
    
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      dtime = rcollection%Dquickaccess(1)
    else
      dtime = 0.0_DP
    end if
    
    ! Type of RHS.
    irhs = rcollection%IquickAccess(3)

    ! Depending on the irhs identifier, choose the RHS to evaluate
    select case (irhs)
    case (-1)
      ! Given as analytical expression. Evaluate the expression.
      !
      ! Get the parser object with the RHS expressions from the collection
      p_rparser => collct_getvalue_pars (rcollection, 'RHSPARSER') 
      
      ! Prepare the array with the values for the function.
      ! X-coordinate, Y-coordinate, time.
      allocate(p_Dval(3,npointsPerElement*nelements))
      do i=1,nelements
        do j=1,npointsPerElement
          p_Dval (1,(i-1)*npointsPerElement+j) = Dpoints(1,j,i)
          p_Dval (2,(i-1)*npointsPerElement+j) = Dpoints(2,j,i)
          p_Dval (3,(i-1)*npointsPerElement+j) = dtime
        end do
      end do
      
      ! Evaluate the 2nd expression for the Y-rhs
      allocate(DvaluesAct(npointsPerElement*nelements))
      call fparser_evalFunction (p_rparser, 2, 2, p_Dval, DvaluesAct)

      ! Reshape the data, that's it.
      do i=0,nelements-1
        do j=1,npointsPerElement
          Dcoefficients(1,j,i+1) = DvaluesAct(i*npointsPerElement+j)
        end do
      end do
      
      deallocate(DvaluesAct)
      deallocate(p_Dval)
    
    case (0)
      Dcoefficients(:,:,:) = 0.0_DP
      
    ! Other cases may be defined here.
      !Dcoefficients(1,:,:) = -Dpoints(2,:,:)
      !Dcoefficients(1,:,:) = -18.0*cos(3.0*SYS_PI*Dpoints(1,:,:))*SYS_PI**2 &
      !                     *cos(3.0*SYS_PI*Dpoints(2,:,:)) &
      !                     - .5*SYS_PI*cos(.5*SYS_PI*(Dpoints(1,:,:)-Dpoints(2,:,:)))
      !Dcoefficients (1,:,:) = -(1./10.0_DP)*(Dpoints(2,:,:))
      
      ! Without coupling:
      !Dcoefficients (1,:,:) = (1/5.0_DP - dtime/50.0_DP)*(-Dpoints(2,:,:))
      !Dcoefficients (1,:,:) = -Dpoints(2,:,:)
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_TARGET_x (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset, &
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form of the
    ! target X-velocity, i.e. the desired flow field z, which enters the dual
    ! equation.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(IN)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    real(DP) :: dtime
  
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      dtime = rcollection%Dquickaccess(1)
    else
      dtime = 0.0_DP
    end if

    Dcoefficients(:,:,:) = 0.0_DP

    ! Call ffunction_TargetX to calculate the analytic function. Store the results
    ! in  Dcoefficients(1,:,:).
    call ffunction_TargetX (DER_FUNC,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dcoefficients(1,:,:),rcollection)
               
    ! Without coupling:
    !Dcoefficients (1,:,:) = - (1/5.0_DP - (10.-dtime)/50.0_DP)*(Dpoints(1,:,:))

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_TARGET_y (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form of the
    ! target X-velocity, i.e. the desired flow field z, which enters the dual
    ! equation.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(IN)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    real(DP) :: dtime

    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      dtime = rcollection%Dquickaccess(1)
    else
      dtime = 0.0_DP
    end if

    Dcoefficients(:,:,:) = 0.0_DP

    ! Call ffunction_TargetX to calculate the analytic function. Store the results
    ! in  Dcoefficients(1,:,:).
    call ffunction_TargetY (DER_FUNC,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dcoefficients(1,:,:),rcollection)

    ! Without coupling:
    !Dcoefficients (1,:,:) = - (1/5.0_DP - (10.-dtime)/50.0_DP)*(-Dpoints(2,:,:))
               
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine ffunction_TargetX (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  use basicgeometry
  use triangulation
  use collection
  use scalarpde
  use domainintegration
  
!<description>
  ! This routine calculates the X-contribution of the target function $z$
  ! in the optimal control problem.
  !
  ! The routine accepts a set of elements and a set of points on these
  ! elements (cubature points) in in real coordinates.
  ! According to the terms in the linear form, the routine has to compute
  ! simultaneously for all these points.
!</description>
  
!<input>
  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
  ! The result must be written to the Dvalue-array below.
  integer, intent(IN)                                         :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(IN)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(IN)                                         :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! DIMENSION(NDIM2D,npointsPerElement,nelements)
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

  ! An array accepting the DOF's on all elements trial in the trial space.
  ! DIMENSION(\#local DOF's in trial space,Number of elements)
  integer, dimension(:,:), intent(IN) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It's usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

  ! Optional: A collection structure to provide additional 
  ! information to the coefficient routine. 
  type(t_collection), intent(INOUT), optional      :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(OUT)                      :: Dvalues
!</output>
  
!</subroutine>

    real(DP) :: dtime,dtimeMax
    integer :: itimedependence,itypeTargetFlow,ieltype,i,j
    type(t_vectorBlock), pointer :: p_rvector
    type(t_spaceTimeVector) :: rspaceTimeVector
    
    integer :: iel,NEL
    real(DP), dimension(:), allocatable :: DvaluesAct
    real(DP), dimension(:,:), allocatable :: DpointsAct
    integer, dimension(:), allocatable :: IelementsAct
    
    type(t_evalElementSet) :: revalElementSet
    integer :: cevaluationTag

    type(t_fparser), pointer :: p_rparser
    real(dp), dimension(:,:), allocatable :: p_Dval
    
    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Ddata
    integer(I32), dimension(1,1) :: ItwistIndex

    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      dtime = rcollection%Dquickaccess(1)
      dtimeMax = rcollection%Dquickaccess(3)
      itimedependence = rcollection%Iquickaccess(1)
      itypeTargetFlow = rcollection%Iquickaccess(2)
    else
      itimedependence = 0
      itypeTargetFlow = 0
      dtime = 0.0_DP
      dtimeMax = 0.0_DP
    end if

!!    Dcoefficients(1,:,:) = -&        
!!    !  mprim_getParabolicProfile (Dpoints(2,:,:),0.41_DP,0.3_DP)
!!      mprim_getParabolicProfile (Dpoints(2,:,:),1.0_DP,0.3_DP)
!    !  -1.0_DP*(Dpoints(1,:,:)/2.2_DP)
!    Dvalues(:,:) = &        
!!           -mprim_signum(Dpoints(2,:,:)-0.5_DP)*mprim_getParabolicProfile(&
!!           MIN(0.5_DP,ABS(Dpoints(2,:,:)-0.5_DP)), &
!!           !*( 0.5_DP - SQRT((Dpoints(1,:,:)-0.5_DP)**2+(Dpoints(2,:,:)-0.5_DP)**2) ) ,&
!!           0.5_DP,1.0_DP) ) 
!    -(Dpoints(2,:,:)-0.5_DP)/SQRT((Dpoints(1,:,:)-0.5_DP)**2+(Dpoints(2,:,:)-0.5_DP)**2) * &
!    mprim_getParabolicProfile( &
!      MIN(0.5_DP,SQRT((Dpoints(1,:,:)-0.5_DP)**2+(Dpoints(2,:,:)-0.5_DP)**2)),&
!      0.5_DP,1.0_DP)
!  
!    IF (ASSOCIATED(p_rcollection)) THEN
!      Dvalues(:,:) = Dvalues(:,:)*mprim_getParabolicProfile(dtime,10.0_DP,1.0_DP)
!    END IF
!    !Dvalues(:,:) = 1.0_DP
    
    if (itimedependence .ne. 0) then  
    
      select case (itypeTargetFlow)
      case (-1)
        ! Zero target flow
        Dvalues(:,:) = 0.0_DP
      case (0)
        ! Analytically given target flow
    
        !Dvalues(:,:) = Dvalues(:,:)*dtime/dtimeMax
        !Dvalues(:,:) = Dvalues(:,:)*dtime
        !Dvalues(:,:) = (-(dtime**2)/100._DP + dtime/5._DP) * Dpoints(1,:,:)
        !Dvalues(:,:) = ((10._DP-dtime)/50._DP - 1._DP/5._DP) * Dpoints(1,:,:)
        Dvalues(:,:) = & ! 1._DP/50._DP * Dpoints(1,:,:) + &
                      (-(dtime**2)/100._DP + dtime/5._DP) * Dpoints(1,:,:)
        !Dvalues(:,:) = ( ((10._DP-dtime)/50._DP - 1._DP/5._DP) + &
        !                 (-(dtime**2)/100._DP + dtime/5._DP)) * Dpoints(1,:,:)
        !Dvalues(:,:) = 0.0_DP
        !IF (dtime .gt. 10._DP) THEN
        !  Dvalues(:,:) = (-(10._DP**2)/100._DP + 10._DP/5._DP) * Dpoints(1,:,:)
        !END IF
        Dvalues(:,:) = Dpoints(1,:,:)
        
        !Dvalues (:,:) = 4.0_DP * Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) * dtime
        !Dvalues (:,:) = 2.0_DP/3.0_DP
        Dvalues (:,:) = 0.3_DP * 4.0_DP * Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:))
        
!      CASE (1:2)
!        ! Target flow is specified by a block vector.
!        !
!        ! Fetch the block vector from the collection
!        p_rvector => collct_getvalue_vec (rcollection,'TARGETFLOW')
!        
!        IF (dof_igetNDofGlob(rdiscretisation) .NE. p_rvector%RvectorBlock(1)%NEQ) THEN
!          CALL output_line ('Target flow vector invalid, NEQ wrong!',&
!              OU_CLASS_ERROR,OU_MODE_STD,'ffunction_TargetX')
!          CALL sys_halt()
!        END IF
!        
!        ! Element type
!        ieltype = &
!            rdiscretisation%RelementDistr(rdomainIntSubset%ielementDistribution)% &
!            itrialElement
!        
!        ! DEBUG!!!
!        CALL lsyssc_getbase_double (p_rvector%RvectorBlock(1),p_Ddata)
!        
!        CALL fevl_evaluate_sim(p_rvector%RvectorBlock(1),rdomainIntSubset%p_Dcoords,&
!            rdomainIntSubset%p_Djac, rdomainIntSubset%p_Ddetj, &
!            ieltype, IdofsTest, npointsPerElement,  nelements, &
!            Dpoints, DER_FUNC, Dvalues, ItwistIndex)
!        
!      CASE (3:4)
      case(1:4)
        ! Target flow is specified by a block vector.
        !
        ! Fetch the block vector from the collection
        p_rvector => collct_getvalue_vec (rcollection,'TARGETFLOW')
        
        ! For every point, find the element of an element nearby the point.
        ! The evaluation routine uses this as hint to speed up the evaluation.
        allocate(IelementsAct(npointsPerElement*nelements))
        allocate(DpointsAct(NDIM2D,npointsPerElement*nelements))
        allocate(DvaluesAct(npointsPerElement*nelements))
        NEL = p_rvector%p_rblockDiscr%p_rtriangulation%NEL
        do i=0,nelements-1
          do j=1,npointsPerElement
            iel = rdomainIntSubset%p_Ielements(i+1)
            
            ! If the FE function is evaluated on a level higher than the
            ! discretisation allows, the element number may be out of bounds.
            ! Reduce it by dividing by 4 which simulates coarsening for a mesh
            ! which is refined by 2-level ordering.
            do while (iel .gt. NEL)
              iel = iel / 4
            end do
            
            IelementsAct(i*npointsPerElement+j) = iel
            
            DpointsAct(1,i*npointsPerElement+j) = Dpoints(1,j,i+1)
            DpointsAct(2,i*npointsPerElement+j) = Dpoints(2,j,i+1)
          end do
        end do
        
        ! Evaluate at the given points - X-coordinate
        call fevl_evaluate (DER_FUNC, DvaluesAct, p_rvector%RvectorBlock(1), &
          DpointsAct, IelementsHint=IelementsAct,cnonmeshPoints=FEVL_NONMESHPTS_ZERO)
          
        do i=0,nelements-1
          do j=1,npointsPerElement
            Dvalues(j,i+1) = DvaluesAct(i*npointsPerElement+j)
          end do
        end do
        
        deallocate(IelementsAct)
        deallocate(DpointsAct)
        deallocate(DvaluesAct)
        
      case (5)
        ! Get the parser object with the RHS expressions from the collection
        p_rparser => collct_getvalue_pars (rcollection, 'TARGETFLOWPARSER') 
        
        ! Prepare the array with the values for the function.
        ! X-coordinate, Y-coordinate, time.
        allocate(p_Dval(3,npointsPerElement*nelements))
        do i=1,nelements
          do j=1,npointsPerElement
            p_Dval (1,(i-1)*npointsPerElement+j) = Dpoints(1,j,i)
            p_Dval (2,(i-1)*npointsPerElement+j) = Dpoints(2,j,i)
            p_Dval (3,(i-1)*npointsPerElement+j) = dtime
          end do
        end do
        
        ! Evaluate the 1st expression for the X-rhs
        allocate(DvaluesAct(npointsPerElement*nelements))
        call fparser_evalFunction (p_rparser, 1, 2, p_Dval, DvaluesAct)

        ! Reshape the data, that's it.
        do i=0,nelements-1
          do j=1,npointsPerElement
            Dvalues(j,i+1) = DvaluesAct(i*npointsPerElement+j)
          end do
        end do
        
        deallocate(DvaluesAct)
        deallocate(p_Dval)
        
      case DEFAULT
        Dvalues(:,:) = Dpoints(1,:,:)
      end select
    end if
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine ffunction_TargetY (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  use basicgeometry
  use triangulation
  use collection
  use scalarpde
  use domainintegration
  
!<description>
  ! This routine calculates the Y-contribution of the target function $z$
  ! in the optimal control problem.
  !
  ! The routine accepts a set of elements and a set of points on these
  ! elements (cubature points) in in real coordinates.
  ! According to the terms in the linear form, the routine has to compute
  ! simultaneously for all these points.
!</description>
  
!<input>
  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
  ! The result must be written to the Dvalue-array below.
  integer, intent(IN)                                         :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(IN)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(IN)                                         :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! DIMENSION(NDIM2D,npointsPerElement,nelements)
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

  ! An array accepting the DOF's on all elements trial in the trial space.
  ! DIMENSION(\#local DOF's in trial space,Number of elements)
  integer, dimension(:,:), intent(IN) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It's usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

  ! Optional: A collection structure to provide additional 
  ! information to the coefficient routine. 
  type(t_collection), intent(INOUT), optional      :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(OUT)                      :: Dvalues
!</output>
  
!</subroutine>

    real(DP) :: dtime,dtimeMax
    integer :: itimedependence,itypeTargetFlow,ieltype,i,j
    type(t_vectorBlock), pointer :: p_rvector

    type(t_fparser), pointer :: p_rparser
    real(dp), dimension(:,:), allocatable :: p_Dval

    integer :: iel,NEL
    real(DP), dimension(:), allocatable :: DvaluesAct
    real(DP), dimension(:,:), allocatable :: DpointsAct
    integer, dimension(:), allocatable :: IelementsAct
    integer(I32), dimension(1,1) :: ItwistIndex

    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      dtime = rcollection%Dquickaccess(1)
      dtimeMax = rcollection%Dquickaccess(3)
      itimedependence = rcollection%Iquickaccess(1)
      itypeTargetFlow = rcollection%Iquickaccess(2)
    else
      itimedependence = 0
      dtime = 0.0_DP
      dtimeMax = 0.0_DP
    end if

!    Dvalues(:,:) = &
!!    !  0.0_DP
!!           mprim_signum(Dpoints(1,:,:)-0.5_DP)*&
!!           mprim_getParabolicProfile(&
!!           MIN(0.5_DP,ABS(Dpoints(1,:,:)-0.5_DP)), &
!!           !*( 0.5_DP - SQRT((Dpoints(1,:,:)-0.5_DP)**2+(Dpoints(2,:,:)-0.5_DP)**2) ),&
!!           0.5_DP,1.0_DP)
!!
!    (Dpoints(1,:,:)-0.5_DP)/SQRT((Dpoints(1,:,:)-0.5_DP)**2+(Dpoints(2,:,:)-0.5_DP)**2) * &
!    mprim_getParabolicProfile( &
!      MIN(0.5_DP,SQRT((Dpoints(1,:,:)-0.5_DP)**2+(Dpoints(2,:,:)-0.5_DP)**2)),&
!      0.5_DP,1.0_DP) 
!      
!    IF (ASSOCIATED(p_rcollection)) THEN  
!      Dvalues(:,:) = Dvalues(:,:)*mprim_getParabolicProfile(dtime,10.0_DP,1.0_DP)
!    END IF
!    !Dvalues(:,:) = 1.0_DP

    if (itimedependence .ne. 0) then  
    
      select case (itypeTargetFlow)
      case (-1)
        ! Zero target flow
        Dvalues(:,:) = 0.0_DP
      case (0)
        ! Analytically given target flow
    
        !Dvalues(:,:) = Dvalues(:,:)*dtime/dtimeMax
        !Dvalues(:,:) = Dvalues(:,:)*dtime
        !Dvalues(:,:) = (-(dtime**2)/100._DP + dtime/5._DP) * (-Dpoints(2,:,:))
        !Dvalues(:,:) = ((10._DP-dtime)/50._DP - 1._DP/5._DP) * (-Dpoints(2,:,:))
        Dvalues(:,:) = & !1._DP/50._DP * (-Dpoints(2,:,:)) + &
                      (-(dtime**2)/100._DP + dtime/5._DP) * (-Dpoints(2,:,:))
        !Dvalues(:,:) = ( ((10._DP-dtime)/50._DP - 1._DP/5._DP) + &
        !                 (-(dtime**2)/100._DP + dtime/5._DP)) * (-Dpoints(2,:,:))
        !Dvalues(:,:) = 0.0_DP
        !IF (dtime .gt. 10._DP) THEN
        !  Dvalues(:,:) = (-(10._DP**2)/100._DP + 10._DP/5._DP) * (-Dpoints(2,:,:))
        !END IF
        Dvalues(:,:) = (-Dpoints(2,:,:))
        
        Dvalues(:,:) = 0.0_DP
        
!      CASE (1:2)
!        ! Target flow is specified by a block vector.
!        !
!        ! Fetch the block vector from the collection
!        p_rvector => collct_getvalue_vec (rcollection,'TARGETFLOW')
!        
!        IF (dof_igetNDofGlob(rdiscretisation) .NE. p_rvector%RvectorBlock(1)%NEQ) THEN
!          CALL output_line ('Target flow vector invalid, NEQ wrong!',&
!              OU_CLASS_ERROR,OU_MODE_STD,'ffunction_TargetY')
!          CALL sys_halt()              
!        END IF
!        
!        ! Element type
!        ieltype = &
!            rdiscretisation%RelementDistr(rdomainIntSubset%ielementDistribution)% &
!            itrialElement
!        
!        ! Evaluate at the given points
!        CALL fevl_evaluate_sim(p_rvector%RvectorBlock(2),rdomainIntSubset%p_Dcoords,&
!            rdomainIntSubset%p_Djac, rdomainIntSubset%p_Ddetj, &
!                  ieltype, IdofsTest, npointsPerElement,  nelements, &
!                  Dpoints, DER_FUNC, Dvalues, ItwistIndex)
!        
!      CASE (3:4)
      case (1:4)
        ! Target flow is specified by a block vector.
        !
        ! Fetch the block vector from the collection
        p_rvector => collct_getvalue_vec (rcollection,'TARGETFLOW')
        
        ! For every point, find the element of an element nearby the point.
        ! The evaluation routine uses this as hint to speed up the evaluation.
        allocate(IelementsAct(npointsPerElement*nelements))
        allocate(DpointsAct(NDIM2D,npointsPerElement*nelements))
        allocate(DvaluesAct(npointsPerElement*nelements))
        NEL = p_rvector%p_rblockDiscr%p_rtriangulation%NEL
        do i=0,nelements-1
          do j=1,npointsPerElement
            iel = rdomainIntSubset%p_Ielements(i+1)
            
            ! If the FE function is evaluated on a level higher than the
            ! discretisation allows, the element number may be out of bounds.
            ! Reduce it by dividing by 4 which simulates coarsening for a mesh
            ! which is refined by 2-level ordering.
            do while (iel .gt. NEL)
              iel = iel / 4
            end do
            
            IelementsAct(i*npointsPerElement+j) = iel
            
            DpointsAct(1,i*npointsPerElement+j) = Dpoints(1,j,i+1)
            DpointsAct(2,i*npointsPerElement+j) = Dpoints(2,j,i+1)
          end do
        end do
        
        ! Evaluate at the given points - Y-coordinate
        call fevl_evaluate (DER_FUNC, DvaluesAct, p_rvector%RvectorBlock(2), &
          DpointsAct, IelementsHint=IelementsAct,cnonmeshPoints=FEVL_NONMESHPTS_ZERO)
          
        do i=0,nelements-1
          do j=1,npointsPerElement
            Dvalues(j,i+1) = DvaluesAct(i*npointsPerElement+j)
          end do
        end do
        
        deallocate(IelementsAct)
        deallocate(DpointsAct)
        deallocate(DvaluesAct)
        
      case (5)
        ! Get the parser object with the RHS expressions from the collection
        p_rparser => collct_getvalue_pars (rcollection, 'TARGETFLOWPARSER') 
        
        ! Prepare the array with the values for the function.
        ! X-coordinate, Y-coordinate, time.
        allocate(p_Dval(3,npointsPerElement*nelements))
        do i=1,nelements
          do j=1,npointsPerElement
            p_Dval (1,(i-1)*npointsPerElement+j) = Dpoints(1,j,i)
            p_Dval (2,(i-1)*npointsPerElement+j) = Dpoints(2,j,i)
            p_Dval (3,(i-1)*npointsPerElement+j) = dtime
          end do
        end do
        
        ! Evaluate the 2nd expression for the Y-rhs
        allocate(DvaluesAct(npointsPerElement*nelements))
        call fparser_evalFunction (p_rparser, 2, 2, p_Dval, DvaluesAct)

        ! Reshape the data, that's it.
        do i=0,nelements-1
          do j=1,npointsPerElement
            Dvalues(j,i+1) = DvaluesAct(i*npointsPerElement+j)
          end do
        end do
        
        deallocate(DvaluesAct)
        deallocate(p_Dval)

      case DEFAULT
        Dvalues(:,:) = -Dpoints(2,:,:)
      end select

    end if

  end subroutine

  ! ***************************************************************************
  ! Values on the real boundary.

!<subroutine>

  subroutine getBoundaryValues (sexpressionName,icomponent,rdiscretisation,&
                                rboundaryRegion,dwhere, dvalue, rcollection)
  
  use collection
  use spatialdiscretisation
  use discretebc
  
!<description>
  ! This routine is called for all segments on the real boundary which are
  ! marked in the DAT file with "evaluate expression of type -2".
  ! sexpressionName is the name of the expression from the DAT file.
  ! The other parameters define the current position on the boundary.
  ! The routine has to return a value which is used as the result of an
  ! expression, e.g. as Diichlet value on the boundary.
  !
  ! The routine is allowed to return SYS_INFINITY_DP. The behaviour of this
  ! value depends on the type of boundary conditions. For Dirichlet
  ! boundary segments, all points where SYS_INFINITY_DP is returned are
  ! treated as Neumann points.
!</description>
  
!<input>
  ! Name of the expression to be evaluated. This name is configured in the
  ! DAT file for the boundary conditions.
  character(LEN=*), intent(IN) :: sexpressionName
  
  ! Solution component that is currently being processed. 
  ! 1 = X-velocity, 2 = y-velocity,...
  integer, intent(IN) :: icomponent
  
  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
  
  ! Boundary region that is currently being processed.
  type(t_boundaryRegion), intent(IN)                          :: rboundaryRegion
  
  ! Current parameter value of the point on the boundary.
  ! 0-1-parametrisation.
  real(DP), intent(IN)                                        :: dwhere
    
  ! Optional: A collection structure to provide additional 
  ! information to the coefficient routine. 
  type(t_collection), intent(IN), optional      :: rcollection
!</input>

!<output>
  ! Return value of the expression. May be SYS_INFINITY_DP.
  real(DP), intent(OUT) :: dvalue
!</output>
  
!</subroutine>

    ! REAL(DP) :: dtime
    ! REAL(DP) :: dx,dy
    !
    ! To get the X/Y-coordinates of the boundary point, use:
    !
    ! CALL boundary_getCoords(rdiscretisation%p_rboundary, &
    !     rboundaryRegion%iboundCompIdx, dwhere, dx, dy)
    !
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    !
    ! dtime = 0.0_DP
    ! IF (PRESENT(rcollection)) dtime = rcollection%Dquickaccess(1)

    dvalue = 0.0_DP

  end subroutine

  ! ***************************************************************************
  ! Values in a fictitious boundary component:

  subroutine getBoundaryValuesFBC_2D (Icomponents,rdiscretisation,&
                                      Revaluation, rcollection)
  
  use collection
  use spatialdiscretisation
  use discretefbc
  
!<description>
  ! This subroutine is called during the discretisation of boundary
  ! conditions on fictitious boundary components. It calculates a special quantity 
  ! on the boundary, which is then used by the discretisation routines to 
  ! generate a discrete 'snapshot' of the (actually analytic) boundary conditions.
  !
  ! The routine must calculate the values on all elements of the element
  ! list Ielements simultaneously. Iwhere is a list with vertex or edge numbers
  ! where information is to be retrieved. Dvalues is filled with function values
  ! while Binside is set to TRUE for every vertex/edge that is inside of the
  ! corresponding fictitious boundary region (identified by rbcRegion).
!</description>
  
!<input>
  ! Component specifier.
  ! For Dirichlet boundary: 
  !   Icomponents(1..SIZE(Icomponents)) defines the number of the solution component,
  !   the value should be calculated for 
  !   (e.g. 1=1st solution component, e.g. X-velocity, 
  !         2=2nd solution component, e.g. Y-velocity,...,
  !         3=3rd solution component, e.g. pressure)
  !   Example: Icomponents(:) = [1,2] -> Compute velues for X- and Y-velocity
  !     (1=x, 2=y component)
  integer, dimension(:), intent(IN)                           :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_blockDiscretisation), intent(IN)                     :: rdiscretisation
  
  ! Optional: A collection structure to provide additional 
  ! information to the coefficient routine. 
  type(t_collection), intent(INOUT), optional                 :: rcollection

!</input>

!<inputoutput>
  ! A t_discreteFBCevaluation structure array that defines what to evaluate, 
  ! where to evaluate and which accepts the return values.
  ! This callback routine must check out the cinfoNeeded-entry in this structure
  ! to find out what to evaluate.
  ! The other entries in this structure describe where to evaluate.
  ! The result of the evaluation must be written into the p_Dvalues array entry
  ! in this structure.
  !
  ! The number of structures in this array depend on what to evaluate:
  !
  ! For Dirichlet boundary:
  !   revaluation contains as many entries as Icomponents; every entry in
  !   Icomponent corresponds to one entry in revaluation
  !   (so Icomponent(1)=1 defines to evaluate the X-velocity while the 
  !    values for the X-velocity are written to revaluation(1)\%p_Dvalues;
  !    Icomponent(2)=2 defines to evaluate the Y-velocity while the values 
  !    for the Y-velocity are written to revaluation(2)\%p_Dvalues, etc).
  !
  type(t_discreteFBCevaluation), dimension(:), intent(INOUT) :: Revaluation
!</inputoutput>
  
!</subroutine>

    ! local variables
    real(DP) :: ddistance, dxcenter, dycenter, dradius, dx, dy
    real(DP), dimension(:,:), pointer :: p_DvertexCoordinates
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    type(t_triangulation), pointer :: p_rtriangulation
    integer :: ipoint,idx
    
    ! Note: the definition of (analytic) fictitious boundary components 
    ! is performed in 'cc_parseFBDconditions'.
    
    ! Get the triangulation array for the point coordinates
    p_rtriangulation => rdiscretisation%RspatialDiscr(1)%p_rtriangulation
    call storage_getbase_double2d (p_rtriangulation%h_DvertexCoords,&
                                    p_DvertexCoordinates)
    call storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
                                p_IverticesAtElement)
    call storage_getbase_int2d (p_rtriangulation%h_IverticesAtEdge,&
                                p_IverticesAtEdge)

    ! Definition of the circle
    dxcenter = 0.6
    dycenter = 0.2
    dradius  = 0.05
    
    ! Loop through the points where to evaluate:
    do idx = 1,Revaluation(1)%nvalues
    
      ! Get the number of the point to process; may also be number of an
      ! edge or element...
      ipoint = Revaluation(1)%p_Iwhere(idx)
      
      ! Get x- and y-coordinate
      call getXYcoord (Revaluation(1)%cinfoNeeded,ipoint,&
                        p_DvertexCoordinates,&
                        p_IverticesAtElement,p_IverticesAtEdge,&
                        p_rtriangulation%NVT,&
                        dx,dy)
      
      ! Get the distance to the center
      ddistance = sqrt( (dx-dxcenter)**2 + (dy-dycenter)**2 )
      
      ! Point inside?
      if (ddistance .le. dradius) then
      
        ! Denote in the p_Iinside array that we prescribe a value here:
        Revaluation(1)%p_Iinside (idx) = 1
        Revaluation(2)%p_Iinside (idx) = 1
        
        ! We prescribe 0.0 as Dirichlet value here - for x- and y-velocity
        Revaluation(1)%p_Dvalues (idx,1) = 0.0_DP
        Revaluation(2)%p_Dvalues (idx,1) = 0.0_DP
      
      end if
      
    end do

!    ! Definition of a 2nd circle
!    dxcenter = 1.1
!    dycenter = 0.31
!    dradius  = 0.05
!    
!    ! Loop through the points where to evaluate:
!    DO idx = 1,Revaluation(1)%nvalues
!    
!      ! Get the number of the point to process; may also be number of an
!      ! edge or element...
!      ipoint = Revaluation(1)%p_Iwhere(idx)
!      
!      ! Get x- and y-coordinate
!      CALL getXYcoord (Revaluation(1)%cinfoNeeded,ipoint,&
!                       p_DvertexCoordinates,&
!                       p_IverticesAtElement,p_IverticesAtEdge,&
!                       p_rtriangulation%NVT,&
!                       dx,dy)
!      
!      ! Get the distance to the center
!      ddistance = SQRT( (dx-dxcenter)**2 + (dy-dycenter)**2 )
!      
!      ! Point inside?
!      IF (ddistance .LE. dradius) THEN
!      
!        ! Denote in the p_Iinside array that we prescribe a value here:
!        Revaluation(1)%p_Iinside (idx) = 1
!        Revaluation(2)%p_Iinside (idx) = 1
!        
!        ! We prescribe 0.0 as Dirichlet value here - vor X- and Y-velocity
!        Revaluation(1)%p_Dvalues (idx,1) = 0.0_DP
!        Revaluation(2)%p_Dvalues (idx,1) = 0.0_DP
!      
!      END IF
!      
!    END DO
    
  contains
  
    ! ---------------------------------------------------------------
    ! Auxiliary routine: Get X/Y coordinate of a point.
    ! This is either the X/Y coordinate of a corner point or
    ! of an edge-midpoint, depending of cinfoNeeded.
    
    subroutine getXYcoord (cinfoNeeded,iwhere,DvertexCoords,&
                           IverticesAtElement,IverticesAtEdge,NVT,&
                           dx,dy)
    
    ! One of the DISCFBC_NEEDxxxx constanr. DISCFBC_NEEDFUNC interprets
    ! iwhere as ivt (vertex number) and returns the coordinates of that
    ! vertex. DISCFBC_NEEDFUNCMID/DISCFBC_NEEDINTMEAN interprets iwhere
    ! as imid (edge number) and returns the coordinates of the edge
    ! midpoint. DISCFBC_NEEDFUNCELMID interprets iwhere as element number
    ! and returns the element midpoint.
    integer, intent(IN) :: cinfoNeeded
    
    ! Identifier. Either vertex number, edge number or element number,
    ! depending on cinfoNeeded.
    integer(I32), intent(IN) :: iwhere
    
    ! Array with coordinates of all corner vertices (DCORVG)
    real(DP), dimension(:,:), intent(IN) :: DvertexCoords
    
    ! Array with numbers of corner coordinates for all elements (KVERT)
    integer, dimension(:,:), intent(IN) :: IverticesAtElement
    
    ! Array with numbers of points adjacent to all edges
    integer, dimension(:,:), intent(IN) :: IverticesAtEdge
    
    ! Number of vertices in the triangulation
    integer, intent(IN) :: NVT
    
    ! Output: X-coordinate
    real(DP), intent(OUT) :: dx
    
    ! Output: Y-coordinate
    real(DP), intent(OUT) :: dy
    
      ! local variables
      real(DP) :: dm1,dm2
      integer :: i,j
      integer :: iv1,iv2
    
      ! Let's see, what do we have...
      select case (cinfoNeeded)
      case (DISCFBC_NEEDFUNC)
        ! That's easy
        dx = DvertexCoords(1,iwhere)
        dy = DvertexCoords(2,iwhere)
      
      case (DISCFBC_NEEDFUNCMID,DISCFBC_NEEDINTMEAN)
        ! Not much harder; get the two points on the edge and calculate
        ! the mean coordinates.
        iv1 = IverticesAtEdge(1,iwhere)
        iv2 = IverticesAtEdge(2,iwhere)
        dx = 0.5_DP*(DvertexCoords(1,iv1)+DvertexCoords(1,iv2))
        dy = 0.5_DP*(DvertexCoords(2,iv1)+DvertexCoords(2,iv2))
      
      case (DISCFBC_NEEDFUNCELMID)
        ! A little bit more to do; we have three or four corners
        dx = 0.0_DP
        dy = 0.0_DP
        do i=1,size(IverticesAtElement,1)
          iv1 = IverticesAtElement(i,iwhere)
          if (iv1 .ne. 0) then
            dx = dx + DvertexCoords(1,iv1)
            dy = dy + DvertexCoords(2,iv1)
            j = i
          end if
        end do
        
        dx = dx / real(j,DP)
        dy = dy / real(j,DP)
      
      case DEFAULT
        print *,'getBoundaryValuesFBC: Insupported coordinate type to return.'
        stop
      end select
      
    end subroutine

  end subroutine

end module
