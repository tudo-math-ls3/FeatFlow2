!##############################################################################
!# ****************************************************************************
!# <name> cc2dmedium_callback </name>
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
!# -> p_rcollection%IquickAccess(1)   = 0: stationary, 
!#                                      1: nonstationary with explicit time stepping
!# -> p_rcollection%DquickAccess(1)   = current simulation time
!# -> p_rcollection%DquickAccess(2)   = minimum simulation time
!# -> p_rcollection%DquickAccess(3)   = maximum simulation time
!#
!# After the assembly, cc_doneCollectForAssembly is called to clean up.
!# Note: Information stored in the quick-access array are of temporary
!# nature and does not have to be cleaned up.
!#
!# </purpose>
!##############################################################################

MODULE cc2dmedium_callback

  USE fsystem
  USE storage
  USE linearsolver
  USE boundary
  USE bilinearformevaluation
  USE linearformevaluation
  USE cubature
  USE matrixfilters
  USE vectorfilters
  USE bcassembly
  USE mprimitives
  USE feevaluation
  USE spacetimevectors
  USE triasearch
  
  USE timeevaluation
  
  USE cc2dmediumm2basic
  USE cc2dmediumm2boundarydef
  
  IMPLICIT NONE

CONTAINS

! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_initCollectForAssembly (rproblem,rcollection)
  
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
  TYPE(t_problem), INTENT(INOUT) :: rproblem
!</input>

!<inputoutput>
  ! Collection structure to be initialised.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
!</inputoutput>
  
!</subroutine>
    REAL(DP) :: dreltime

    ! In a nonstationary simulation, save the simulation time as well as the
    ! minimum and maximum time to the quick-access array of the collection,
    ! so it can be accessed in the callback routines!
    rcollection%Iquickaccess(1) = rproblem%itimedependence
    SELECT CASE (rproblem%itimedependence)
    CASE (0)
      ! Stationary simulation
      rcollection%Dquickaccess(1) = 0.0_DP
      rcollection%Dquickaccess(2) = 0.0_DP
      rcollection%Dquickaccess(3) = 0.0_DP
    CASE (1)
      rcollection%Dquickaccess(1) = rproblem%rtimedependence%dtime
      rcollection%Dquickaccess(2) = rproblem%rtimedependence%dtimeInit
      rcollection%Dquickaccess(3) = rproblem%rtimedependence%dtimeMax
    END SELECT
    
    rcollection%Iquickaccess(2) = rproblem%roptcontrol%itypeTargetFlow
    CALL collct_setvalue_vec (rcollection, 'TARGETFLOW', &
        rproblem%roptcontrol%rtargetFlow, .TRUE.) 
    
    ! In case the target vector changes in time, load the current target
    ! from the space-time vector.    
    IF ((rproblem%roptcontrol%itypeTargetFlow .EQ. 2) .OR. &
        (rproblem%roptcontrol%itypeTargetFlow .EQ. 4))  THEN
      IF (rproblem%roptcontrol%rtargetFlowNonstat%NEQtime .GT. 1) THEN
        
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
        
        CALL tmevl_evaluate(rproblem%roptcontrol%rtargetFlowNonstat,&
            rproblem%rtimedependence%dtime,&
            rproblem%roptcontrol%rtargetFlow)
            
      END IF
      ! Otherwise, there is no vector.
    END IF

  END SUBROUTINE
  
! ***************************************************************************
  
!<subroutine>

  SUBROUTINE cc_doneCollectForAssembly (rproblem,rcollection)
  
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
  TYPE(t_problem), INTENT(IN) :: rproblem
!</input>

!<inputoutput>
  ! Collection structure to be cleaned up.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
!</inputoutput>
  
!</subroutine>

    ! Currently, this subroutine is empty as all information stored in
    ! the collection in cc_initCollectForAssembly is put to the quick-access
    ! arrays -- which do not have to be cleaned up. 
    ! This might change in future...

    CALL collct_deletevalue (rcollection, 'TARGETFLOW') 

  END SUBROUTINE
  
! ***************************************************************************
  !<subroutine>

  SUBROUTINE coeff_Stokes (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial,IdofsTest,rdomainIntSubset, p_rcollection,&
                  Dcoefficients)
    
    USE basicgeometry
    USE triangulation
    USE collection
    USE scalarpde
    USE domainintegration
    
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
    ! analytic boundary boundary description etc.
    TYPE(t_spatialDiscretisation), INTENT(IN)                   :: rdiscretisation
    
    ! The bilinear form which is currently being evaluated:
    TYPE(t_bilinearForm), INTENT(IN)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    INTEGER, INTENT(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    INTEGER, INTENT(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: Dpoints
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in trial space,nelements)
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(IN) :: IdofsTrial
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in test space,nelements)
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(IN) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    TYPE(t_domainIntSubset), INTENT(IN)              :: rdomainIntSubset

    ! A pointer to a collection structure to provide additional 
    ! information to the coefficient routine. May point to NULL() if not defined.
    TYPE(t_collection), POINTER                      :: p_rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    REAL(DP), DIMENSION(:,:,:), INTENT(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    ! Routine is not called anyway in case of constant coefficients!
    Dcoefficients = 1.0_DP

  END SUBROUTINE

! ***************************************************************************
  !<subroutine>

  SUBROUTINE coeff_Pressure (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial,IdofsTest,rdomainIntSubset, p_rcollection,&
                  Dcoefficients)
    
    USE basicgeometry
    USE triangulation
    USE collection
    USE scalarpde
    USE domainintegration
    
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
    ! analytic boundary boundary description etc.
    TYPE(t_spatialDiscretisation), INTENT(IN)                   :: rdiscretisation
    
    ! The bilinear form which is currently being evaluated:
    TYPE(t_bilinearForm), INTENT(IN)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    INTEGER, INTENT(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    INTEGER, INTENT(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: Dpoints
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in trial space,nelements)
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(IN) :: IdofsTrial
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in test space,nelements)
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(IN) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    TYPE(t_domainIntSubset), INTENT(IN)              :: rdomainIntSubset

    ! A pointer to a collection structure to provide additional 
    ! information to the coefficient routine. May point to NULL() if not defined.
    TYPE(t_collection), POINTER                      :: p_rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    REAL(DP), DIMENSION(:,:,:), INTENT(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    ! Routine is not called anyway in case of constant coefficients!
    Dcoefficients = 1.0_DP

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE coeff_RHS_x (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,p_rcollection, &
                  Dcoefficients)
    
    USE basicgeometry
    USE triangulation
    USE collection
    USE scalarpde
    USE domainintegration
    
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
    TYPE(t_spatialDiscretisation), INTENT(IN)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    TYPE(t_linearForm), INTENT(IN)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    INTEGER, INTENT(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    INTEGER, INTENT(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in test space,nelements)
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    TYPE(t_domainIntSubset), INTENT(IN)              :: rdomainIntSubset

    ! A pointer to a collection structure to provide additional 
    ! information to the coefficient routine. May point to NULL() if not defined.
    TYPE(t_collection), POINTER                      :: p_rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    REAL(DP), DIMENSION(:,:,:), INTENT(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    REAL(DP) :: dtime
    
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    IF (ASSOCIATED(p_rcollection)) THEN
      dtime = p_rcollection%Dquickaccess(1)
    ELSE
      dtime = 0.0_DP
    END IF
    
    Dcoefficients(:,:,:) = 0.0_DP
    !Dcoefficients(1,:,:) = Dpoints(1,:,:)
    !Dcoefficients(1,:,:) = -18.0*sin(3.0*SYS_PI*Dpoints(1,:,:))*SYS_PI**2 &
    !                     *sin(3.0*SYS_PI*Dpoints(2,:,:)) &
    !                     + .5*SYS_PI*cos(.5*SYS_PI*(Dpoints(1,:,:)-Dpoints(2,:,:)))
    !Dcoefficients (1,:,:) = -(1./10.0_DP)*(-Dpoints(1,:,:))
    
    ! Without coupling:
    !Dcoefficients (1,:,:) = (1/5.0_DP - dtime/50.0_DP)*(Dpoints(1,:,:))
    !Dcoefficients (1,:,:) = Dpoints(1,:,:)
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE coeff_RHS_y (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,p_rcollection, &
                  Dcoefficients)
    
    USE basicgeometry
    USE triangulation
    USE collection
    USE scalarpde
    USE domainintegration
    
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
    TYPE(t_spatialDiscretisation), INTENT(IN)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    TYPE(t_linearForm), INTENT(IN)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    INTEGER, INTENT(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    INTEGER, INTENT(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in test space,nelements)
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    TYPE(t_domainIntSubset), INTENT(IN)              :: rdomainIntSubset

    ! A pointer to a collection structure to provide additional 
    ! information to the coefficient routine. May point to NULL() if not defined.
    TYPE(t_collection), POINTER                      :: p_rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    REAL(DP), DIMENSION(:,:,:), INTENT(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    REAL(DP) :: dtime
    
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    IF (ASSOCIATED(p_rcollection)) THEN
      dtime = p_rcollection%Dquickaccess(1)
    ELSE
      dtime = 0.0_DP
    END IF
    
    Dcoefficients(:,:,:) = 0.0_DP
    !Dcoefficients(1,:,:) = -Dpoints(2,:,:)
    !Dcoefficients(1,:,:) = -18.0*cos(3.0*SYS_PI*Dpoints(1,:,:))*SYS_PI**2 &
    !                     *cos(3.0*SYS_PI*Dpoints(2,:,:)) &
    !                     - .5*SYS_PI*cos(.5*SYS_PI*(Dpoints(1,:,:)-Dpoints(2,:,:)))
    !Dcoefficients (1,:,:) = -(1./10.0_DP)*(Dpoints(2,:,:))
    
    ! Without coupling:
    !Dcoefficients (1,:,:) = (1/5.0_DP - dtime/50.0_DP)*(-Dpoints(2,:,:))
    !Dcoefficients (1,:,:) = -Dpoints(2,:,:)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE coeff_TARGET_x (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,p_rcollection, &
                  Dcoefficients)
    
    USE basicgeometry
    USE triangulation
    USE collection
    USE scalarpde
    USE domainintegration
    
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
    TYPE(t_spatialDiscretisation), INTENT(IN)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    TYPE(t_linearForm), INTENT(IN)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    INTEGER, INTENT(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    INTEGER, INTENT(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in test space,nelements)
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    TYPE(t_domainIntSubset), INTENT(IN)              :: rdomainIntSubset

    ! A pointer to a collection structure to provide additional 
    ! information to the coefficient routine. May point to NULL() if not defined.
    TYPE(t_collection), POINTER                      :: p_rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    REAL(DP), DIMENSION(:,:,:), INTENT(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    REAL(DP) :: dtime
  
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    IF (ASSOCIATED(p_rcollection)) THEN
      dtime = p_rcollection%Dquickaccess(1)
    ELSE
      dtime = 0.0_DP
    END IF

    Dcoefficients(:,:,:) = 0.0_DP

    ! Call ffunction_TargetX to calculate the analytic function. Store the results
    ! in  Dcoefficients(1,:,:).
    CALL ffunction_TargetX (DER_FUNC,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,p_rcollection, &
                Dcoefficients(1,:,:))
               
    ! Without coupling:
    !Dcoefficients (1,:,:) = - (1/5.0_DP - (10.-dtime)/50.0_DP)*(Dpoints(1,:,:))

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE coeff_TARGET_y (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,p_rcollection, &
                  Dcoefficients)
    
    USE basicgeometry
    USE triangulation
    USE collection
    USE scalarpde
    USE domainintegration
    
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
    TYPE(t_spatialDiscretisation), INTENT(IN)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    TYPE(t_linearForm), INTENT(IN)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    INTEGER, INTENT(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    INTEGER, INTENT(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in test space,nelements)
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    TYPE(t_domainIntSubset), INTENT(IN)              :: rdomainIntSubset

    ! A pointer to a collection structure to provide additional 
    ! information to the coefficient routine. May point to NULL() if not defined.
    TYPE(t_collection), POINTER                      :: p_rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    REAL(DP), DIMENSION(:,:,:), INTENT(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    REAL(DP) :: dtime

    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    IF (ASSOCIATED(p_rcollection)) THEN
      dtime = p_rcollection%Dquickaccess(1)
    ELSE
      dtime = 0.0_DP
    END IF

    Dcoefficients(:,:,:) = 0.0_DP

    ! Call ffunction_TargetX to calculate the analytic function. Store the results
    ! in  Dcoefficients(1,:,:).
    CALL ffunction_TargetY (DER_FUNC,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,p_rcollection, &
                Dcoefficients(1,:,:))

    ! Without coupling:
    !Dcoefficients (1,:,:) = - (1/5.0_DP - (10.-dtime)/50.0_DP)*(-Dpoints(2,:,:))
               
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE ffunction_TargetX (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,p_rcollection, &
                Dvalues)
  
  USE basicgeometry
  USE triangulation
  USE collection
  USE scalarpde
  USE domainintegration
  
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
  INTEGER, INTENT(IN)                                         :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  TYPE(t_spatialDiscretisation), INTENT(IN)                   :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  INTEGER, INTENT(IN)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  INTEGER, INTENT(IN)                                         :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! DIMENSION(NDIM2D,npointsPerElement,nelements)
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: Dpoints

  ! An array accepting the DOF's on all elements trial in the trial space.
  ! DIMENSION(\#local DOF's in trial space,Number of elements)
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(IN) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It's usually used in more complex situations (e.g. nonlinear matrices).
  TYPE(t_domainIntSubset), INTENT(IN)              :: rdomainIntSubset

  ! A pointer to a collection structure to provide additional 
  ! information to the coefficient routine. May point to NULL() if not defined.
  TYPE(t_collection), POINTER                      :: p_rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  REAL(DP), DIMENSION(:,:), INTENT(OUT)                      :: Dvalues
!</output>
  
!</subroutine>

    REAL(DP) :: dtime,dtimeMax
    INTEGER :: itimedependence,itypeTargetFlow,ieltype,i,j
    TYPE(t_vectorBlock), POINTER :: p_rvector
    TYPE(t_spaceTimeVector) :: rspaceTimeVector
    
    INTEGER(PREC_ELEMENTIDX) :: iel,NEL
    REAL(DP), DIMENSION(:), ALLOCATABLE :: DvaluesAct
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: DpointsAct
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), ALLOCATABLE :: IelementsAct
    
    ! DEBUG!!!
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata

    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    IF (ASSOCIATED(p_rcollection)) THEN
      dtime = p_rcollection%Dquickaccess(1)
      dtimeMax = p_rcollection%Dquickaccess(3)
      itimedependence = p_rcollection%Iquickaccess(1)
      itypeTargetFlow = p_rcollection%Iquickaccess(2)
    ELSE
      itimedependence = 0
      itypeTargetFlow = 0
      dtime = 0.0_DP
      dtimeMax = 0.0_DP
    END IF

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

    IF (itimedependence .NE. 0) THEN  
    
      SELECT CASE (itypeTargetFlow)
      CASE (-1)
        ! Zero target flow
        Dvalues(:,:) = 0.0_DP
      CASE (0)
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
      CASE (1:2)
        ! Target flow is specified by a block vector.
        !
        ! Fetch the block vector from the collection
        p_rvector => collct_getvalue_vec (p_rcollection,'TARGETFLOW')
        
        IF (dof_igetNDofGlob(rdiscretisation) .NE. p_rvector%RvectorBlock(1)%NEQ) THEN
          CALL output_line ('Target flow vector invalid, NEQ wrong!',&
              OU_CLASS_ERROR,OU_MODE_STD,'ffunction_TargetX')
          CALL sys_halt()
        END IF
        
        ! Element type
        ieltype = &
            rdiscretisation%RelementDistribution(rdomainIntSubset%ielementDistribution)% &
            itrialElement
        
        ! DEBUG!!!
        CALL lsyssc_getbase_double (p_rvector%RvectorBlock(1),p_Ddata)
        
        ! Evaluate at the given points
        CALL fevl_evaluate_sim(p_rvector%RvectorBlock(1),rdomainIntSubset%p_Dcoords,&
            rdomainIntSubset%p_Djac, rdomainIntSubset%p_Ddetj, &
            ieltype, IdofsTest, npointsPerElement,  nelements, &
            Dpoints, DER_FUNC, Dvalues)
        
      CASE (3:4)
        ! Target flow is specified by a block vector.
        !
        ! Fetch the block vector from the collection
        p_rvector => collct_getvalue_vec (p_rcollection,'TARGETFLOW')
        
        ! For every point, find the element of an element nearby the point.
        ! The evaluation routine uses this as hint to speed up the evaluation.
        ALLOCATE(IelementsAct(npointsPerElement*nelements))
        ALLOCATE(DpointsAct(NDIM2D,npointsPerElement*nelements))
        ALLOCATE(DvaluesAct(npointsPerElement*nelements))
        NEL = p_rvector%p_rblockDiscretisation%p_rtriangulation%NEL
        DO i=0,nelements-1
          DO j=1,npointsPerElement
            iel = rdomainIntSubset%p_Ielements(i+1)
            
            ! If the FE function is evaluated on a level higher than the
            ! discretisation allows, the element number may be out of bounds.
            ! Reduce it by dividing by 4 which simulates coarsening for a mesh
            ! which is refined by 2-level ordering.
            DO WHILE (iel .GT. NEL)
              iel = iel / 4
            END DO
            
            IelementsAct(i*npointsPerElement+j) = iel
            
            DpointsAct(1,i*npointsPerElement+j) = Dpoints(1,j,i+1)
            DpointsAct(2,i*npointsPerElement+j) = Dpoints(2,j,i+1)
          END DO
        END DO
        
        ! Evaluate at the given points - X-coordinate
        CALL fevl_evaluate (DER_FUNC, DvaluesAct, p_rvector%RvectorBlock(1), &
          DpointsAct, IelementsHint=IelementsAct,cnonmeshPoints=FEVL_NONMESHPTS_ZERO)
          
        DO i=0,nelements-1
          DO j=1,npointsPerElement
            Dvalues(j,i+1) = DvaluesAct(i*npointsPerElement+j)
          END DO
        END DO
        
        DEALLOCATE(IelementsAct)
        DEALLOCATE(DpointsAct)
        DEALLOCATE(DvaluesAct)
        
      CASE DEFAULT
        Dvalues(:,:) = Dpoints(1,:,:)
      END SELECT
    END IF
  
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE ffunction_TargetY (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,p_rcollection, &
                Dvalues)
  
  USE basicgeometry
  USE triangulation
  USE collection
  USE scalarpde
  USE domainintegration
  
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
  INTEGER, INTENT(IN)                                         :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  TYPE(t_spatialDiscretisation), INTENT(IN)                   :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  INTEGER, INTENT(IN)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  INTEGER, INTENT(IN)                                         :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! DIMENSION(NDIM2D,npointsPerElement,nelements)
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: Dpoints

  ! An array accepting the DOF's on all elements trial in the trial space.
  ! DIMENSION(\#local DOF's in trial space,Number of elements)
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(IN) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It's usually used in more complex situations (e.g. nonlinear matrices).
  TYPE(t_domainIntSubset), INTENT(IN)              :: rdomainIntSubset

  ! A pointer to a collection structure to provide additional 
  ! information to the coefficient routine. May point to NULL() if not defined.
  TYPE(t_collection), POINTER                      :: p_rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  REAL(DP), DIMENSION(:,:), INTENT(OUT)                      :: Dvalues
!</output>
  
!</subroutine>

    REAL(DP) :: dtime,dtimeMax
    INTEGER :: itimedependence,itypeTargetFlow,ieltype,i,j
    TYPE(t_vectorBlock), POINTER :: p_rvector

    INTEGER(PREC_ELEMENTIDX) :: iel,NEL
    REAL(DP), DIMENSION(:), ALLOCATABLE :: DvaluesAct
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: DpointsAct
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), ALLOCATABLE :: IelementsAct

    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    IF (ASSOCIATED(p_rcollection)) THEN
      dtime = p_rcollection%Dquickaccess(1)
      dtimeMax = p_rcollection%Dquickaccess(3)
      itimedependence = p_rcollection%Iquickaccess(1)
      itypeTargetFlow = p_rcollection%Iquickaccess(2)
    ELSE
      itimedependence = 0
      dtime = 0.0_DP
      dtimeMax = 0.0_DP
    END IF

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

    IF (itimedependence .NE. 0) THEN  
    
      SELECT CASE (itypeTargetFlow)
      CASE (-1)
        ! Zero target flow
        Dvalues(:,:) = 0.0_DP
      CASE (0)
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
      CASE (1:2)
        ! Target flow is specified by a block vector.
        !
        ! Fetch the block vector from the collection
        p_rvector => collct_getvalue_vec (p_rcollection,'TARGETFLOW')
        
        IF (dof_igetNDofGlob(rdiscretisation) .NE. p_rvector%RvectorBlock(1)%NEQ) THEN
          CALL output_line ('Target flow vector invalid, NEQ wrong!',&
              OU_CLASS_ERROR,OU_MODE_STD,'ffunction_TargetY')
          CALL sys_halt()              
        END IF
        
        ! Element type
        ieltype = &
            rdiscretisation%RelementDistribution(rdomainIntSubset%ielementDistribution)% &
            itrialElement
        
        ! Evaluate at the given points
        CALL fevl_evaluate_sim(p_rvector%RvectorBlock(2),rdomainIntSubset%p_Dcoords,&
            rdomainIntSubset%p_Djac, rdomainIntSubset%p_Ddetj, &
                  ieltype, IdofsTest, npointsPerElement,  nelements, &
                  Dpoints, DER_FUNC, Dvalues)
        
      CASE (3:4)
        ! Target flow is specified by a block vector.
        !
        ! Fetch the block vector from the collection
        p_rvector => collct_getvalue_vec (p_rcollection,'TARGETFLOW')
        
        ! For every point, find the element of an element nearby the point.
        ! The evaluation routine uses this as hint to speed up the evaluation.
        ALLOCATE(IelementsAct(npointsPerElement*nelements))
        ALLOCATE(DpointsAct(NDIM2D,npointsPerElement*nelements))
        ALLOCATE(DvaluesAct(npointsPerElement*nelements))
        NEL = p_rvector%p_rblockDiscretisation%p_rtriangulation%NEL
        DO i=0,nelements-1
          DO j=1,npointsPerElement
            iel = rdomainIntSubset%p_Ielements(i+1)
            
            ! If the FE function is evaluated on a level higher than the
            ! discretisation allows, the element number may be out of bounds.
            ! Reduce it by dividing by 4 which simulates coarsening for a mesh
            ! which is refined by 2-level ordering.
            DO WHILE (iel .GT. NEL)
              iel = iel / 4
            END DO
            
            IelementsAct(i*npointsPerElement+j) = iel
            
            DpointsAct(1,i*npointsPerElement+j) = Dpoints(1,j,i+1)
            DpointsAct(2,i*npointsPerElement+j) = Dpoints(2,j,i+1)
          END DO
        END DO
        
        ! Evaluate at the given points - Y-coordinate
        CALL fevl_evaluate (DER_FUNC, DvaluesAct, p_rvector%RvectorBlock(2), &
          DpointsAct, IelementsHint=IelementsAct,cnonmeshPoints=FEVL_NONMESHPTS_ZERO)
          
        DO i=0,nelements-1
          DO j=1,npointsPerElement
            Dvalues(j,i+1) = DvaluesAct(i*npointsPerElement+j)
          END DO
        END DO
        
        DEALLOCATE(IelementsAct)
        DEALLOCATE(DpointsAct)
        DEALLOCATE(DvaluesAct)
        
      CASE DEFAULT
        Dvalues(:,:) = -Dpoints(2,:,:)
      END SELECT

    END IF

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE getBoundaryValues (Icomponents,rdiscretisation,rbcRegion,ielement, &
                                cinfoNeeded,iwhere,dwhere, p_rcollection, Dvalues)
  
  USE collection
  USE spatialdiscretisation
  USE discretebc
  
!<description>
  ! This subroutine is called during the discretisation of boundary
  ! conditions. It calculates a special quantity on the boundary, which is
  ! then used by the discretisation routines to generate a discrete
  ! 'snapshot' of the (actually analytic) boundary conditions.
!</description>
  
!<input>
  ! Component specifier.
  ! For Dirichlet boundary: 
  !   Icomponents(1) defines the number of the boundary component, the value
  !   should be calculated for (e.g. 1=1st solution component, e.g. X-velocitry, 
  !   2=2nd solution component, e.g. Y-velocity,...)
  INTEGER, DIMENSION(:), INTENT(IN)                           :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  TYPE(t_spatialDiscretisation), INTENT(IN)                   :: rdiscretisation
  
  ! Boundary condition region that is currently being processed.
  ! (This e.g. defines the type of boundary conditions that are
  !  currently being calculated, as well as information about the current
  !  boundary segment 'where we are at the moment'.)
  TYPE(t_bcRegion), INTENT(IN)                                :: rbcRegion
  
  
  ! The element number on the boundary which is currently being processed
  INTEGER(I32), INTENT(IN)                                    :: ielement
  
  ! The type of information, the routine should calculate. One of the
  ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
  ! to return one or multiple information value in the result array.
  INTEGER, INTENT(IN)                                         :: cinfoNeeded
  
  ! A reference to a geometric object where information should be computed.
  ! cinfoNeeded=DISCBC_NEEDFUNC : 
  !   iwhere = number of the point in the triangulation or
  !          = 0, if only the parameter value of the point is known; this
  !               can be found in dwhere,
  ! cinfoNeeded=DISCBC_NEEDDERIV : 
  !   iwhere = number of the point in the triangulation or
  !          = 0, if only the parameter value of the point is known; this
  !               can be found in dwhere,
  ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
  !   iwhere = number of the edge where the value integral mean value
  !            should be computed
  INTEGER, INTENT(IN)                                         :: iwhere

  ! A reference to a geometric object where information should be computed.
  ! cinfoNeeded=DISCBC_NEEDFUNC : 
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDDERIV : 
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
  !   dwhere = 0 (not used)
  REAL(DP), INTENT(IN)                                        :: dwhere
    
  ! A pointer to a collection structure to provide additional 
  ! information to the coefficient routine. May point to NULL() if not defined.
  TYPE(t_collection), POINTER                  :: p_rcollection

!</input>

!<output>
  ! This array receives the calculated information. If the caller
  ! only needs one value, the computed quantity is put into Dvalues(1). 
  ! If multiple values are needed, they are collected here (e.g. for 
  ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
  !
  ! The function may return SYS_INFINITY as a value. This indicates the
  ! framework to ignore the node and treat it as 'natural boundary condition'
  ! node (Neumann boundary).
  REAL(DP), DIMENSION(:), INTENT(OUT)                         :: Dvalues
!</output>
  
!</subroutine>

    INTEGER :: icomponent,iexprtyp
    
    REAL(DP) :: dtime
    
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    IF (ASSOCIATED(p_rcollection)) THEN
      dtime = p_rcollection%Dquickaccess(1)
    ELSE
      dtime = 0.0_DP
    END IF
    
    ! Use boundary conditions from DAT files.
    SELECT CASE (cinfoNeeded)
    
    CASE (DISCBC_NEEDFUNC,DISCBC_NEEDFUNCMID,DISCBC_NEEDDERIV, &
          DISCBC_NEEDINTMEAN,DISCBC_NEEDNORMALSTRESS)
      
      ! Dirichlet boundary conditions
    
      ! Get from the current component of the PDE we are discretising:
      icomponent = Icomponents(1)
      
      ! -> 1=X-velocity, 2=Y-velocity.
      
      ! Return zero Dirichlet boundary values for all situations by default.
      Dvalues(1) = 0.0_DP

      ! Now, depending on the problem, calculate the return value.
      
      ! Get the type of the expression to evaluate from the 
      ! integer tag of the BC-region - if there is an expression to evaluate
      ! at all.
      iexprtyp = rbcRegion%ibdrexprtype
            
      ! Now, which boundary condition do we have here?                       
      SELECT CASE (rbcRegion%ctype)
      CASE (BC_DIRICHLET)
        ! Simple Dirichlet BC's. Evaluate the expression iexprtyp.
        Dvalues(1) = evalBoundary (rdiscretisation, rbcRegion%rboundaryRegion, &
                                    iexprtyp, rbcRegion%itag, rbcRegion%dtag, &
                                    dwhere, rbcRegion%stag,dtime,&
                                    p_rcollection)
    
      CASE (BC_PRESSUREDROP)
        ! Normal stress / pressure drop. Evaluate Evaluate the 
        ! expression iexprtyp.
        Dvalues(1) = evalBoundary (rdiscretisation, rbcRegion%rboundaryRegion, &
                                    iexprtyp, rbcRegion%itag, rbcRegion%dtag, &
                                    dwhere, rbcRegion%stag,dtime,&
                                    p_rcollection)
      END SELECT
      
    END SELECT
  
  CONTAINS
  
    ! Auxiliary function: Evaluate a scalar expression on the boundary.
    
    REAL(DP) FUNCTION evalBoundary (rdiscretisation, rboundaryRegion, &
                                    ityp, ivalue, dvalue, dpar, stag, dtime, p_rcollection)
    
    ! Discretisation structure of the underlying discretisation
    TYPE(t_spatialDiscretisation), INTENT(IN) :: rdiscretisation
    
    ! Current boundary region
    TYPE(t_boundaryRegion), INTENT(IN) :: rboundaryRegion
    
    ! Type of expression to evaluate.
    ! One of the BDC_xxxx constants from cc2dmediumm2boundarydef.f90.
    INTEGER, INTENT(IN) :: ityp
    
    ! Integer tag. If ityp=BDC_EXPRESSION, this must specify the number of
    ! the expression in the expression object to evaluate.
    ! Otherwise unused.
    INTEGER, INTENT(IN) :: ivalue
    
    ! Double precision parameter for simple expressions
    REAL(DP), INTENT(IN) :: dvalue
    
    ! Current parameter value of the point on the boundary.
    ! 0-1-parametrisation.
    REAL(DP), INTENT(IN) :: dpar
    
    ! Time in a nonstationary simulation. =0 for stationary simulations.
    REAL(DP), INTENT(IN) :: dtime

    ! String tag that defines more complicated BC's.
    CHARACTER(LEN=*), INTENT(IN) :: stag
    
    ! A compiled expression for evaluation at runtime
    TYPE(t_fparser), POINTER :: p_rparser
    
    ! A pointer to a collection structure to provide additional 
    ! information to the coefficient routine. May point to NULL() if not defined.
    TYPE(t_collection), POINTER                  :: p_rcollection

      ! local variables
      REAL(DP) :: d,dx,dy
      CHARACTER(LEN=PARLST_MLDATA) :: sexpr
      REAL(DP), DIMENSION(SIZE(SEC_EXPRVARIABLES)) :: Rval
      
      SELECT CASE (ityp)
      CASE (BDC_USERDEF)
        ! This is a hardcoded, user-defined identifier.
        ! In stag, the name of the identifier is noted.
        ! Get the identifier itself from the collection.
        CALL collct_getvalue_string (p_rcollection, stag, sexpr, &
                                     0, SEC_SBDEXPRESSIONS)
                                     
        ! Now we can decide on 'sexpr' how to evaluate.
        ! By default, we return 0.0. A user defined calculation can be added here!
        evalBoundary = 0.0_DP
      
      CASE (BDC_VALDOUBLE)
        ! A simple constant, given by dvalue
        evalBoundary = dvalue

      CASE (BDC_EXPRESSION)
        ! A complex expression.
        ! Get the expression object from the collection.
        
        p_rparser => collct_getvalue_pars (p_rcollection, BDC_BDPARSER, &
                                   0, SEC_SBDEXPRESSIONS)
                                   
        ! Set up an array with variables for evaluating the expression.
        ! Give the values in exactly the same order as specified
        ! by SEC_EXPRVARIABLES!
        Rval = 0.0_DP
        
        CALL boundary_getCoords(rdiscretisation%p_rboundary, &
                                rboundaryRegion%iboundCompIdx, &
                                dpar, dx, dy)
        
        ! Get the local parameter value 0 <= d <= 1.
        ! Note that if dpar < rboundaryRegion%dminParam, we have to add the maximum
        ! parameter value on the boundary to dpar as normally 0 <= dpar < max.par.
        ! although 0 <= dminpar <= max.par 
        !      and 0 <= dmaxpar <= max.par!
        d = dpar 
        IF (d .LT. rboundaryRegion%dminParam) &
          d = d + boundary_dgetMaxParVal(rdiscretisation%p_rboundary,&
                                         rboundaryRegion%iboundCompIdx)
        d = d - rboundaryRegion%dminParam
        
        Rval(1) = dx
        Rval(2) = dy
        ! Rval(3) = .
        Rval(4) = d
        Rval(5) = dpar
        Rval(6) = boundary_convertParameter(rdiscretisation%p_rboundary, &
                                            rboundaryRegion%iboundCompIdx, dpar, &
                                            BDR_PAR_01, BDR_PAR_LENGTH) 
        Rval(7) = dtime
        
        ! Evaluate the expression. ivalue is the number of
        ! the expression to evaluate.
        CALL fparser_evalFunction (p_rparser, ivalue, Rval, evalBoundary)
        
      CASE (BDC_VALPARPROFILE)
        ! A parabolic profile. dvalue expresses the
        ! maximum value of the profile. 
        !
        ! Get the local parameter value 0 <= d <= 1.
        ! Note that if dpar < rboundaryRegion%dminParam, we have to add the maximum
        ! parameter value on the boundary to dpar as normally 0 <= dpar < max.par.
        ! although 0 <= dminpar <= max.par 
        !      and 0 <= dmaxpar <= max.par!
        d = dpar 
        IF (d .LT. rboundaryRegion%dminParam) &
          d = d + boundary_dgetMaxParVal(rdiscretisation%p_rboundary,&
                                         rboundaryRegion%iboundCompIdx)
        d = d - rboundaryRegion%dminParam
    
        evalBoundary = mprim_getParabolicProfile (d,1.0_DP,dvalue) 
      END SELECT
    
    END FUNCTION
    
  END SUBROUTINE

  ! ***************************************************************************
  ! Values in a fictitious boundary component:

!<subroutine>

  SUBROUTINE getBoundaryValuesFBC (Icomponents,rdiscretisation,rbcRegion, &
                                   Revaluation, p_rcollection)
  
  USE collection
  USE spatialdiscretisation
  USE discretefbc
  
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
  !   Example: Icomponents(:) = [1,2] -> Compute values for X- and Y-velocity
  !     (1=x, 2=y component)
  INTEGER, DIMENSION(:), INTENT(IN)                           :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  TYPE(t_blockDiscretisation), INTENT(IN)                     :: rdiscretisation
  
  ! Boundary condition region that is currently being processed.
  ! (This e.g. defines the type of boundary conditions that are 
  !  currently being calculated (Dirichlet, Robin,...), as well as information
  !  about the current boundary region 'what is discretised at the moment'.)
  TYPE(t_bcRegion), INTENT(IN)                                :: rbcRegion
  
  ! A pointer to a collection structure to provide additional 
  ! information to the coefficient routine. May point to NULL() if not defined.
  TYPE(t_collection), POINTER                                 :: p_rcollection

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
  ! For Dirichlet boudary:
  !   revaluation contains as many entries as Icomponents; every entry in
  !   Icomponent corresponds to one entry in revaluation
  !   (so Icomponent(1)=1 defines to evaluate the X-velocity while the 
  !    values for the X-velocity are written to revaluation(1)\%p_Dvalues;
  !    Icomponent(2)=2 defines to evaluate the Y-velocity while the values 
  !    for the Y-velocity are written to revaluation(2)\%p_Dvalues, etc).
  !
  TYPE(t_discreteFBCevaluation), DIMENSION(:), INTENT(INOUT) :: Revaluation
!</inputoutput>
  
!</subroutine>

    ! local variables
    REAL(DP) :: ddistance, dxcenter, dycenter, dradius, dx, dy
    REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoordinates
    INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtEdge
    TYPE(t_triangulation), POINTER :: p_rtriangulation
    INTEGER :: ipoint,idx
    
    ! Note: the definition of (analytic) fictitious boundary components 
    ! is performed in 'cc_parseFBDconditions'.
    
    ! Are we evaluating our fictitious boundary component?
    IF (rbcRegion%rfictBoundaryRegion%sname .EQ. 'CIRCLE') THEN
    
      ! Get the triangulation array for the point coordinates
      p_rtriangulation => rdiscretisation%RspatialDiscretisation(1)%p_rtriangulation
      CALL storage_getbase_double2d (p_rtriangulation%h_DvertexCoords,&
                                     p_DvertexCoordinates)
      CALL storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
                                  p_IverticesAtElement)
      CALL storage_getbase_int2d (p_rtriangulation%h_IverticesAtEdge,&
                                  p_IverticesAtEdge)

      ! Definition of the circle
      dxcenter = 0.6
      dycenter = 0.2
      dradius  = 0.05
      
      ! Loop through the points where to evaluate:
      DO idx = 1,Revaluation(1)%nvalues
      
        ! Get the number of the point to process; may also be number of an
        ! edge or element...
        ipoint = Revaluation(1)%p_Iwhere(idx)
        
        ! Get x- and y-coordinate
        CALL getXYcoord (Revaluation(1)%cinfoNeeded,ipoint,&
                         p_DvertexCoordinates,&
                         p_IverticesAtElement,p_IverticesAtEdge,&
                         p_rtriangulation%NVT,&
                         dx,dy)
        
        ! Get the distance to the center
        ddistance = SQRT( (dx-dxcenter)**2 + (dy-dycenter)**2 )
        
        ! Point inside?
        IF (ddistance .LE. dradius) THEN
        
          ! Denote in the p_Iinside array that we prescribe a value here:
          Revaluation(1)%p_Iinside (idx) = 1
          Revaluation(2)%p_Iinside (idx) = 1
          
          ! We prescribe 0.0 as Dirichlet value here - for x- and y-velocity
          Revaluation(1)%p_Dvalues (idx,1) = 0.0_DP
          Revaluation(2)%p_Dvalues (idx,1) = 0.0_DP
        
        END IF
        
      END DO

!      ! Definition of a 2nd circle
!      dxcenter = 1.1
!      dycenter = 0.31
!      dradius  = 0.05
!      
!      ! Loop through the points where to evaluate:
!      DO idx = 1,Revaluation(1)%nvalues
!      
!        ! Get the number of the point to process; may also be number of an
!        ! edge or element...
!        ipoint = Revaluation(1)%p_Iwhere(idx)
!        
!        ! Get x- and y-coordinate
!        CALL getXYcoord (Revaluation(1)%cinfoNeeded,ipoint,&
!                         p_DvertexCoordinates,&
!                         p_IverticesAtElement,p_IverticesAtEdge,&
!                         p_rtriangulation%NVT,&
!                         dx,dy)
!        
!        ! Get the distance to the center
!        ddistance = SQRT( (dx-dxcenter)**2 + (dy-dycenter)**2 )
!        
!        ! Point inside?
!        IF (ddistance .LE. dradius) THEN
!        
!          ! Denote in the p_Iinside array that we prescribe a value here:
!          Revaluation(1)%p_Iinside (idx) = 1
!          Revaluation(2)%p_Iinside (idx) = 1
!          
!          ! We prescribe 0.0 as Dirichlet value here - vor X- and Y-velocity
!          Revaluation(1)%p_Dvalues (idx,1) = 0.0_DP
!          Revaluation(2)%p_Dvalues (idx,1) = 0.0_DP
!        
!        END IF
!        
!      END DO
    
    END IF
    
  CONTAINS
  
    ! ---------------------------------------------------------------
    ! Auxiliary routine: Get X/Y coordinate of a point.
    ! This is either the X/Y coordinate of a corner point or
    ! of an edge-midpoint, depending of cinfoNeeded.
    
    SUBROUTINE getXYcoord (cinfoNeeded,iwhere,DvertexCoords,&
                           IverticesAtElement,IverticesAtEdge,NVT,&
                           dx,dy)
    
    ! One of the DISCFBC_NEEDxxxx constanr. DISCFBC_NEEDFUNC interprets
    ! iwhere as ivt (vertex number) and returns the coordinates of that
    ! vertex. DISCFBC_NEEDFUNCMID/DISCFBC_NEEDINTMEAN interprets iwhere
    ! as imid (edge number) and returns the coordinates of the edge
    ! midpoint. DISCFBC_NEEDFUNCELMID interprets iwhere as element number
    ! and returns the element midpoint.
    INTEGER, INTENT(IN) :: cinfoNeeded
    
    ! Identifier. Either vertex number, edge number or element number,
    ! depending on cinfoNeeded.
    INTEGER(I32), INTENT(IN) :: iwhere
    
    ! Array with coordinates of all corner vertices (DCORVG)
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: DvertexCoords
    
    ! Array with numbers of corner coordinates for all elements (KVERT)
    INTEGER(PREC_POINTIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElement
    
    ! Array with numbers of points adjacent to all edges
    INTEGER(PREC_POINTIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtEdge
    
    ! Number of vertices in the triangulation
    INTEGER(PREC_POINTIDX), INTENT(IN) :: NVT
    
    ! Output: X-coordinate
    REAL(DP), INTENT(OUT) :: dx
    
    ! Output: Y-coordinate
    REAL(DP), INTENT(OUT) :: dy
    
      ! local variables
      REAL(DP) :: dm1,dm2
      INTEGER :: i,j
      INTEGER(PREC_POINTIDX) :: iv1,iv2
    
      ! Let's see, what do we have...
      SELECT CASE (cinfoNeeded)
      CASE (DISCFBC_NEEDFUNC)
        ! That's easy
        dx = DvertexCoords(1,iwhere)
        dy = DvertexCoords(2,iwhere)
      
      CASE (DISCFBC_NEEDFUNCMID,DISCFBC_NEEDINTMEAN)
        ! Not much harder; get the two points on the edge and calculate
        ! the mean coordinates.
        iv1 = IverticesAtEdge(1,iwhere)
        iv2 = IverticesAtEdge(2,iwhere)
        dx = 0.5_DP*(DvertexCoords(1,iv1)+DvertexCoords(1,iv2))
        dy = 0.5_DP*(DvertexCoords(2,iv1)+DvertexCoords(2,iv2))
      
      CASE (DISCFBC_NEEDFUNCELMID)
        ! A little bit more to do; we have three or four corners
        dx = 0.0_DP
        dy = 0.0_DP
        DO i=1,SIZE(IverticesAtElement,1)
          iv1 = IverticesAtElement(i,iwhere)
          IF (iv1 .NE. 0) THEN
            dx = dx + DvertexCoords(1,iv1)
            dy = dy + DvertexCoords(2,iv1)
            j = i
          END IF
        END DO
        
        dx = dx / REAL(j,DP)
        dy = dy / REAL(j,DP)
      
      CASE DEFAULT
        PRINT *,'getBoundaryValuesFBC: Insupported coordinate type to return.'
        STOP
      END SELECT
      
    END SUBROUTINE

  END SUBROUTINE

END MODULE
