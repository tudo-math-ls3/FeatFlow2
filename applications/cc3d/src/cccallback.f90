!##############################################################################
!# ****************************************************************************
!# <name> cccallback </name>
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
!# 4.) coeff_RHS_z
!#     -> Returns analytical values for the right hand side of the Z-velocity
!#        equation.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 5.) ffunction_TargetX
!#     -> Returns analytical values for the desired flow field in X-direction.
!#     -> Is used for error analysis during the postprocessing.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 6.) ffunction_TargetY
!#     -> Returns analytical values for the desired flow field in Y-direction.
!#     -> Is used for error analysis during the postprocessing.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 7.) ffunction_TargetZ
!#     -> Returns analytical values for the desired flow field in Z-direction.
!#     -> Is used for error analysis during the postprocessing.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 8.) ffunction_TargetP
!#     -> Returns analytical values for the desired pressure.
!#     -> Is used for error analysis during the postprocessing.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 11.) cc_initCollectForAssembly
!#     -> Is called prior to the assembly process.
!#     -> Stores some information from the problem structure to a collection
!#        such that it can be accessed in callback routines
!#
!# 12.) cc_doneCollectForAssembly
!#      -> Is called after the assembly process.
!#      -> Releases information stored in the collection by 
!#         cc_initCollectForAssembly.
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

MODULE cccallback

  USE fsystem
  USE storage
  USE linearsolver
  USE bilinearformevaluation
  USE linearformevaluation
  USE cubature
  USE matrixfilters
  USE vectorfilters
  USE bcassembly
  USE mprimitives
  USE meshregion
  
  USE ccbasic
  
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
  TYPE(t_problem), INTENT(IN) :: rproblem
!</input>

!<inputoutput>
  ! Collection structure to be initialised.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
!</inputoutput>
  
!</subroutine>

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

  END SUBROUTINE
  
! ***************************************************************************
  !<subroutine>

  SUBROUTINE coeff_Stokes (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial,IdofsTest,rdomainIntSubset, &
                  Dcoefficients,rcollection)
    
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

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    TYPE(t_collection), INTENT(INOUT), OPTIONAL      :: rcollection
    
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
                  IdofsTrial,IdofsTest,rdomainIntSubset, &
                  Dcoefficients,rcollection)
    
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

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    TYPE(t_collection), INTENT(INOUT), OPTIONAL      :: rcollection
    
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
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
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

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    TYPE(t_collection), INTENT(INOUT), OPTIONAL      :: rcollection
    
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
    IF (PRESENT(rcollection)) THEN
      dtime = rcollection%Dquickaccess(1)
    ELSE
      dtime = 0.0_DP
    END IF
    
    Dcoefficients(:,:,:) = 0.0_DP

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE coeff_RHS_y (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
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

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    TYPE(t_collection), INTENT(INOUT), OPTIONAL      :: rcollection
    
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
    IF (PRESENT(rcollection)) THEN
      dtime = rcollection%Dquickaccess(1)
    ELSE
      dtime = 0.0_DP
    END IF
    
    Dcoefficients(:,:,:) = 0.0_DP

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE coeff_RHS_z (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    USE basicgeometry
    USE triangulation
    USE collection
    USE scalarpde
    USE domainintegration
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form of the
    ! Z-velocity part of the right hand side vector.
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

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    TYPE(t_collection), INTENT(INOUT), OPTIONAL      :: rcollection
    
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
    IF (PRESENT(rcollection)) THEN
      dtime = rcollection%Dquickaccess(1)
    ELSE
      dtime = 0.0_DP
    END IF
    
    Dcoefficients(:,:,:) = 0.0_DP

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE ffunction_TargetX (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  USE basicgeometry
  USE triangulation
  USE collection
  USE scalarpde
  USE domainintegration
  
!<description>
  ! This subroutine is called during the postprocessing. 
  ! It should return values of the analytical solution (if it is known).
  ! These are compared with the calculated solution to calculate the
  ! error in the X-velocity.
  !
  ! If the analytical solution is unknown, this routine doesn't make sense.
  ! In this case, error analysis should be deactivated in the .DAT files!
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
  ! DIMENSION(NDIM3D,npointsPerElement,nelements)
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
  ! information to the coefficient routine. 
  TYPE(t_collection), INTENT(INOUT), OPTIONAL      :: rcollection
  
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
    INTEGER :: itimedependence

    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    IF (PRESENT(rcollection)) THEN
      dtime = rcollection%Dquickaccess(1)
      dtimeMax = rcollection%Dquickaccess(3)
      itimedependence = rcollection%Iquickaccess(1)
    ELSE
      itimedependence = 0
      dtime = 0.0_DP
      dtimeMax = 0.0_DP
    END IF

    Dvalues(:,:) = 0.0_DP
    
    ! Example:
    ! IF (cderivative .EQ. DER_FUNC) THEN
    !   Dvalues(:,:) = (-dtime**2/100.+dtime/5.)*(Dpoints(1,:,:))
    ! END IF

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE ffunction_TargetY (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  USE basicgeometry
  USE triangulation
  USE collection
  USE scalarpde
  USE domainintegration
  
!<description>
  ! This subroutine is called during the postprocessing. 
  ! It should return values of the analytical solution (if it is known).
  ! These are compared with the calculated solution to calculate the
  ! error in the Y-velocity.
  !
  ! If the analytical solution is unknown, this routine doesn't make sense.
  ! In this case, error analysis should be deactivated in the .DAT files!
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
  ! DIMENSION(NDIM3D,npointsPerElement,nelements)
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
  ! information to the coefficient routine. 
  TYPE(t_collection), INTENT(INOUT), OPTIONAL      :: rcollection
  
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
    INTEGER :: itimedependence

    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    IF (PRESENT(rcollection)) THEN
      dtime = rcollection%Dquickaccess(1)
      dtimeMax = rcollection%Dquickaccess(3)
      itimedependence = rcollection%Iquickaccess(1)
    ELSE
      itimedependence = 0
      dtime = 0.0_DP
      dtimeMax = 0.0_DP
    END IF

    Dvalues(:,:) = 0.0_DP
    
    ! Example:
    ! IF (cderivative .EQ. DER_FUNC) THEN
    !   Dvalues(:,:) = (-dtime**2/100.+dtime/5.)*(-Dpoints(2,:,:))
    ! END IF

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE ffunction_TargetZ (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  USE basicgeometry
  USE triangulation
  USE collection
  USE scalarpde
  USE domainintegration
  
!<description>
  ! This subroutine is called during the postprocessing. 
  ! It should return values of the analytical solution (if it is known).
  ! These are compared with the calculated solution to calculate the
  ! error in the Z-velocity.
  !
  ! If the analytical solution is unknown, this routine doesn't make sense.
  ! In this case, error analysis should be deactivated in the .DAT files!
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
  ! DIMENSION(NDIM3D,npointsPerElement,nelements)
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
  ! information to the coefficient routine. 
  TYPE(t_collection), INTENT(INOUT), OPTIONAL      :: rcollection
  
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
    INTEGER :: itimedependence

    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    IF (PRESENT(rcollection)) THEN
      dtime = rcollection%Dquickaccess(1)
      dtimeMax = rcollection%Dquickaccess(3)
      itimedependence = rcollection%Iquickaccess(1)
    ELSE
      itimedependence = 0
      dtime = 0.0_DP
      dtimeMax = 0.0_DP
    END IF

    Dvalues(:,:) = 0.0_DP
    
    ! Example:
    ! IF (cderivative .EQ. DER_FUNC) THEN
    !   Dvalues(:,:) = (-dtime**2/100.+dtime/5.)*(-Dpoints(2,:,:))
    ! END IF

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE ffunction_TargetP (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  USE basicgeometry
  USE triangulation
  USE collection
  USE scalarpde
  USE domainintegration
  
!<description>
  ! This subroutine is called during the postprocessing. 
  ! It should return values of the analytical solution (if it is known).
  ! These are compared with the calculated solution to calculate the
  ! error in the pressure
  !
  ! If the analytical solution is unknown, this routine doesn't make sense.
  ! In this case, error analysis should be deactivated in the .DAT files!
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
  ! DIMENSION(NDIM3D,npointsPerElement,nelements)
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
  ! information to the coefficient routine. 
  TYPE(t_collection), INTENT(INOUT), OPTIONAL      :: rcollection
  
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
    INTEGER :: itimedependence

    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    IF (PRESENT(rcollection)) THEN
      dtime = rcollection%Dquickaccess(1)
      dtimeMax = rcollection%Dquickaccess(3)
      itimedependence = rcollection%Iquickaccess(1)
    ELSE
      itimedependence = 0
      dtime = 0.0_DP
      dtimeMax = 0.0_DP
    END IF

    Dvalues(:,:) = 0.0_DP

    ! Example:
    ! IF (cderivative .EQ. DER_FUNC) THEN
    !   ...
    ! END IF

  END SUBROUTINE


  ! ***************************************************************************
  ! Values on the real boundary, described by a mesh region

!<subroutine>

  SUBROUTINE getBoundaryValues (sexpressionName,icomponent,rdiscretisation,&
                                rboundaryRegion,Dcoords, dvalue, rcollection)
  
  USE collection
  USE spatialdiscretisation
  USE discretebc
  
!<description>
  ! This routine is called for all segments on the real boundary which are
  ! marked in the DAT file with "evaluate expression of type -2".
  ! sexpressionName is the name of the expression from the DAT file.
  ! The other parameters define the current position on the boundary.
  ! The routine has to return a value which is used as the result of an
  ! expression, e.g. as Dirichlet value on the boundary.
  !
  ! The routine is allowed to return SYS_INFINITY. The behaviour of this
  ! value depends on the type of boundary conditions. For Dirichlet
  ! boundary segments, all points where SYS_INFINITY is returned are
  ! treated as Neumann points.
!</description>
  
!<input>
  ! Name of the expression to be evaluated. This name is configured in the
  ! DAT file for the boundary conditions.
  CHARACTER(LEN=*), INTENT(IN)                   :: sexpressionName
  
  ! Solution component that is currently being processed. 
  ! 1 = X-velocity, 2 = y-velocity,...
  INTEGER, INTENT(IN)                            :: icomponent
  
  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation, etc.
  TYPE(t_spatialDiscretisation), INTENT(IN)      :: rdiscretisation
  
  ! Mesh region that is currently being processed.
  TYPE(t_meshRegion), INTENT(IN)                 :: rboundaryRegion
  
  ! The coordinates of the point
  REAL(DP), DIMENSION(:), INTENT(IN)             :: Dcoords
    
  ! Optional: A collection structure to provide additional 
  ! information to the coefficient routine. 
  TYPE(t_collection), INTENT(IN), OPTIONAL       :: rcollection
!</input>

!<output>
  ! Return value of the expression. May be SYS_INFINITY.
  REAL(DP), INTENT(OUT)                          :: dvalue
!</output>
  
!</subroutine>

    ! Return zero here
    dvalue = 0.0_DP

  END SUBROUTINE

END MODULE
