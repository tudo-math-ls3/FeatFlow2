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
!# 2.) coeff_Pressure
!#     -> Returns the coefficients for the pressure matrix. This routine is
!#        only used if the problem to calculate has nonconstant coefficients!
!#        Otherwise the routine is dead.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientMatrixSc.inc'
!#
!# 3.) coeff_RHS_x
!#     -> Returns analytical values for the right hand side of the X-velocity
!#        equation.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 4.) coeff_RHS_y
!#     -> Returns analytical values for the right hand side of the Y-velocity
!#        equation.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 5.) coeff_RHS_p
!#     -> Returns analytical values for the right hand side of the pressure/
!#        divergence equation.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 6.) coeff_AnalyticSolution_X
!#     -> Returns analytical values for the desired flow field in X-direction.
!#     -> Is used for setting up the initial solution.
!#     -> In the basic implementation, this calls ffunction_TargetX.
!#
!# 7.) coeff_AnalyticSolution_Y
!#     -> Returns analytical values for the desired flow field in Y-direction.
!#     -> Is used for setting up the initial solution.
!#     -> In the basic implementation, this calls ffunction_TargetY.
!#
!# 8.) coeff_AnalyticSolution_P
!#     -> Returns analytical values for the desired pressure.
!#     -> Is used for setting up the initial solution.
!#     -> In the basic implementation, this calls ffunction_TargetP.
!#
!# 9.) ffunction_TargetX
!#     -> Returns analytical values for the desired flow field in X-direction.
!#     -> Is used for error analysis during the postprocessing.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 10.) ffunction_TargetY
!#      -> Returns analytical values for the desired flow field in Y-direction.
!#      -> Is used for error analysis during the postprocessing.
!#      -> Corresponds to the interface defined in the file
!#         'intf_coefficientVectorSc.inc'
!#
!# 11.) ffunction_TargetP
!#      -> Returns analytical values for the desired pressure.
!#      -> Is used for error analysis during the postprocessing.
!#      -> Corresponds to the interface defined in the file
!#         'intf_coefficientVectorSc.inc'
!#
!# 12.) getBoundaryValues
!#     -> Returns analytical values on the boundary of the
!#        problem to solve.
!#
!# 13.) getBoundaryValuesFBC
!#     -> Returns analytical values on the fictitious boundary components
!#     -> Corresponds to the interface defined in the file
!#        'intf_fbcassembly.inc'
!#
!# 14.) cc_initCollectForAssembly
!#     -> Is called prior to the assembly process.
!#     -> Stores some information from the problem structure to a collection
!#        such that it can be accessed in callback routines
!#
!# 15.) cc_doneCollectForAssembly
!#      -> Is called after the assembly process.
!#      -> Releases information stored in the collection by
!#         cc_initCollectForAssembly.
!#
!# 16.) getMovingFrameVelocity
!#      -> If the moving frame formulation is activated, this routine
!#         returns the velocity and acceleration of the moving frame.
!#
!# 17.) getNonconstantViscosity
!#      -> If nonconstant viscosity is activated, this routine calculates
!#         the viscosity.
!#
!# 18.) calcAdaptiveTimestep
!#      -> Calculate a user defined timestep length.
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
!#  EXAMPLES
!# ------------------------
!# 1.) Moving frames formulation
!#
!#   The moving frames formulation realises a moving domain, which is equivalent
!#   to the Lagrangian formulation of the Navier-Stokes equation. It can be
!#   used e.g. to track falling particles in an infinite long channel or similar.
!#   The moving frames formulation is realised as follows:
!#
!#   Let us assume, the domain moves with velocity v=(vx,vy) and accelleration
!#   a=(ax,ay). This information must be returned by the function
!#   getMovingFrameVelocity below. To activate the moving frame formulation,
!#   the parameter imovingFrame in discretisation.dat must be set to 1.
!#   Then, the equation is changed as follows:
!#
!#     (Navier-)Stokes (u,p) = f + a
!#     div(u)                = 0
!#
!#     u(Dirichlet-boundary) = u0 + v
!#
!#   The velocity v is substracted from the postprocessing data to give
!#   nice pictures. Both, v and a, can be accessed in the boundary condition
!#   expressions (like L,R,S,x,y,TIME) as:
!#
!#     MFVELX = x-velocity
!#     MFVELY = y-velocity
!#     MFACCX = x-acceleration
!#     MFACCY = y-acceleration
!#
!#   You can for example generate a simple fixed particle in a moving
!#   bench1-channel by adding the following code to getMovingFrameVelocity:
!#
!#     Dvelocity(1) = 0.3_DP*(tanh(dtime))
!#     Dacceleration(1) = 0.3_DP*(1.0_DP-tanh(dtime)**2)
!#
!#   and the following data to the master.dat file.
!#
!#     ###################
!#     [CC-DISCRETISATION]
!#     ###################
!#     imovingFrame = 1
!#
!#     ##############
!#     [BDEXPRESSIONS]
!#     ##############
!#     bdExpressions(2) =
!#       'Dirichlet0'     0    0.0
!#       'mpartx'        -1    '-MFVELX'
!#
!#     ##############
!#     [BDCONDITIONS]
!#     ##############
!#     bdComponent1(4)=
!#        1.0  3  1  'Dirichlet0'  'Dirichlet0'
!#        2.0  0  0
!#        3.0  3  1  'Dirichlet0'  'Dirichlet0'
!#        4.0  0  0
!#     bdComponent2(1)=
!#        4.0  3  1  'mpartx'  'Dirichlet0'
!#
!#     #####################
!#     [TIME-DISCRETISATION]
!#     #####################
!#     itimedependence = 1
!#     dtimeMax = 10.0
!#
!# </purpose>
!##############################################################################

module cccallback

  use fsystem
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
  use discretebc
  use bcassembly
  use mprimitives
  use derivatives
  use feevaluation
  
  use ccbasic
  
  implicit none

contains

! ***************************************************************************

!<subroutine>

  subroutine cc_initCollectForAssembly (rproblem,rcollection)
  
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
  type(t_problem), intent(in) :: rproblem
!</input>

!<inputoutput>
  ! Collection structure to be initialised.
  type(t_collection), intent(inout) :: rcollection
!</inputoutput>
  
!</subroutine>

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
      rcollection%Dquickaccess(1) = rproblem%rtimedependence%dtime
      rcollection%Dquickaccess(2) = rproblem%rtimedependence%dtimeInit
      rcollection%Dquickaccess(3) = rproblem%rtimedependence%dtimeMax
    end select

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
  type(t_problem), intent(in) :: rproblem
!</input>

!<inputoutput>
  ! Collection structure to be cleaned up.
  type(t_collection), intent(inout) :: rcollection
!</inputoutput>
  
!</subroutine>

    ! Currently, this subroutine is empty as all information stored in
    ! the collection in cc_initCollectForAssembly is put to the quick-access
    ! arrays -- which do not have to be cleaned up.
    ! This might change in future...

  end subroutine
  
! ***************************************************************************
  !<subroutine>

  subroutine coeff_Stokes (rdiscretisationTrial,rdiscretisationTest,rform, &
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTrial
        
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints
    
    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTrial
    
    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    ! Routine is not called anyway in case of constant coefficients!
    Dcoefficients = 1.0_DP

  end subroutine
  
! ***************************************************************************
  !<subroutine>

  subroutine coeff_K11 (rdiscretisationTrial,rdiscretisationTest,rform, &
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTrial
        
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints
    
    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTrial
    
    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    ! Routine is not called anyway in case of constant coefficients!
    Dcoefficients = 1.0_DP

  end subroutine
  
! ***************************************************************************
  !<subroutine>

  subroutine coeff_K22 (rdiscretisationTrial,rdiscretisationTest,rform, &
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTrial
        
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints
    
    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTrial
    
    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    ! Routine is not called anyway in case of constant coefficients!
    Dcoefficients = 1.0_DP

  end subroutine
  
! ***************************************************************************
  !<subroutine>

  subroutine coeff_K12 (rdiscretisationTrial,rdiscretisationTest,rform, &
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTrial
        
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints
    
    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTrial
    
    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    ! Routine is not called anyway in case of constant coefficients!
    Dcoefficients = 1.0_DP

  end subroutine
  
! ***************************************************************************
  !<subroutine>

  subroutine coeff_K21 (rdiscretisationTrial,rdiscretisationTest,rform, &
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTrial
        
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints
    
    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTrial
    
    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    ! Routine is not called anyway in case of constant coefficients!
    Dcoefficients = 1.0_DP

  end subroutine



! ***************************************************************************
  !<subroutine>

  subroutine coeff_K15 (rdiscretisationTrial,rdiscretisationTest,rform, &
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTrial
    
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints
    
    ! An array accepting the DOF`s on all elements in the trial space.
    ! DIMENSION(#local DOF`s in trial space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTrial
    
    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rvector
    real(DP) :: dWeight
    real(DP) :: drhoFR


    real(DP), dimension(:,:,:), allocatable :: div_vS
    real(DP), dimension(:,:,:), allocatable :: vS1_x
    real(DP), dimension(:,:,:), allocatable :: vS2_y
    integer :: i, j


    dWeight  = rcollection%DquickAccess(1)
    drhoFR   = rcollection%DquickAccess(2)

    p_rvector => rcollection%p_RvectorQuickAccess1

    allocate (vS1_x(ubound(Dcoefficients,2),ubound(Dcoefficients,3),1))
    allocate (vS2_y(ubound(Dcoefficients,2),ubound(Dcoefficients,3),1))
    allocate (div_vS(ubound(Dcoefficients,2),ubound(Dcoefficients,3),1))

    call fevl_evaluate_sim (DER_DERIV_X, vS1_x(:,:,1), &
        p_rvector%Rvectorblock(3), Dpoints, &
        rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

    call fevl_evaluate_sim (DER_DERIV_Y, vS2_y(:,:,1), &
        p_rvector%Rvectorblock(4), Dpoints, &
        rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

!    div_vS = vS1_x + vS2_y
    do i=1,ubound(Dcoefficients,3)
      do j=1,ubound(Dcoefficients,2)
        div_vS(j,i,1) = vS1_x(j,i,1) + vS2_y(j,i,1)
      end do
    end do

    deallocate (vS1_x)
    deallocate (vS2_y)

! Bilinear coefficients
    do i=1,ubound(Dcoefficients,3)
      do j=1,ubound(Dcoefficients,2)
        Dcoefficients(1,j,i) = dWeight*drhoFR*div_vS(j,i,1)
      end do
    end do

    deallocate (div_vS)

  end subroutine

! ***************************************************************************
  !<subroutine>

  subroutine coeff_M13 (rdiscretisationTrial,rdiscretisationTest,rform, &
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTrial
    
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints
    
    ! An array accepting the DOF`s on all elements in the trial space.
    ! DIMENSION(#local DOF`s in trial space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTrial
    
    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rvector
    real(DP) :: dWeight
    real(DP) :: drhoFR
    real(DP) :: drhoSR
    real(DP) :: dnF0
    real(DP) :: dnS0

    integer :: i, j

    real(DP), dimension(:,:,:), allocatable :: uS1_x
    real(DP), dimension(:,:,:), allocatable :: uS1_y
    real(DP), dimension(:,:,:), allocatable :: uS2_x
    real(DP), dimension(:,:,:), allocatable :: uS2_y
    real(DP), dimension(:,:,:), allocatable :: dnF
    real(DP), dimension(:,:,:), allocatable :: div_uS
    real(DP), dimension(:,:,:), allocatable :: det_grad_uS


    dWeight  = rcollection%DquickAccess(1)
    drhoFR   = rcollection%DquickAccess(2)
    drhoSR   = rcollection%DquickAccess(3)
    dnF0     = rcollection%DquickAccess(4)



    dnS0     = 1.0_DP-dnF0

    p_rvector => rcollection%p_RvectorQuickAccess1

    allocate (uS1_x(ubound(Dcoefficients,2),ubound(Dcoefficients,3),1))
    allocate (uS1_y(ubound(Dcoefficients,2),ubound(Dcoefficients,3),1))
    allocate (uS2_x(ubound(Dcoefficients,2),ubound(Dcoefficients,3),1))
    allocate (uS2_y(ubound(Dcoefficients,2),ubound(Dcoefficients,3),1))
    allocate (div_uS(ubound(Dcoefficients,2),ubound(Dcoefficients,3),1))
    allocate (dnF(ubound(Dcoefficients,2),ubound(Dcoefficients,3),1))
    allocate (det_grad_uS(ubound(Dcoefficients,2),ubound(Dcoefficients,3),1))

! uS1_x
    call fevl_evaluate_sim (DER_DERIV_X, uS1_x(:,:,1), &
        p_rvector%Rvectorblock(1), Dpoints, &
        rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

! uS1_y
    call fevl_evaluate_sim (DER_DERIV_Y, uS1_y(:,:,1), &
        p_rvector%Rvectorblock(1), Dpoints, &
        rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
! uS2_x
    call fevl_evaluate_sim (DER_DERIV_X, uS2_x(:,:,1), &
        p_rvector%Rvectorblock(2), Dpoints, &
        rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
! uS2_y
    call fevl_evaluate_sim (DER_DERIV_Y, uS2_y(:,:,1), &
        p_rvector%Rvectorblock(2), Dpoints, &
        rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)



!  |grad uS | :  is the determinant value of grad uS
    do i=1,ubound(Dcoefficients,3)
      do j=1,ubound(Dcoefficients,2)
        det_grad_uS(j,i,1) = uS1_x(j,i,1)*uS2_y(j,i,1) - uS2_x(j,i,1)*uS1_y(j,i,1)
      end do
    end do

    deallocate (uS1_y)
    deallocate (uS2_x)

!     div_uS = uS1_x + uS2_y
    do i=1,ubound(Dcoefficients,3)
      do j=1,ubound(Dcoefficients,2)
        div_uS(j,i,1) = uS1_x(j,i,1) + uS2_y(j,i,1)
      end do
    end do

    deallocate (uS1_x)
    deallocate (uS2_y)


! nF = 1 - nS0 (1-div uS + |grad uS | )
    do i=1,ubound(Dcoefficients,3)
      do j=1,ubound(Dcoefficients,2)
        dnF(j,i,1) = 1.0_DP - dnS0*(1.0_DP - div_uS(j,i,1) + det_grad_uS(j,i,1))
      end do
    end do

    deallocate (det_grad_uS)
    deallocate (div_uS)




    do i=1,ubound(Dcoefficients,3)
      do j=1,ubound(Dcoefficients,2)
        Dcoefficients(1,j,i) = dWeight*(dnF(j,i,1)*drhoFR + (1.0_DP-dnF(j,i,1))*drhoSR)
      end do
    end do

! print*,dkF

    deallocate (dnF)

  end subroutine

! ***************************************************************************
  !<subroutine>

  subroutine coeff_M55 (rdiscretisationTrial,rdiscretisationTest,rform, &
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTrial
    
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints
    
    ! An array accepting the DOF`s on all elements in the trial space.
    ! DIMENSION(#local DOF`s in trial space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTrial
    
    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rvector
    real(DP) :: dWeight
    real(DP) :: drhoFR
    real(DP) :: dnF0
    real(DP) :: dnS0

    integer :: i, j

    real(DP), dimension(:,:,:), allocatable :: uS1_x
    real(DP), dimension(:,:,:), allocatable :: uS1_y
    real(DP), dimension(:,:,:), allocatable :: uS2_x
    real(DP), dimension(:,:,:), allocatable :: uS2_y
    real(DP), dimension(:,:,:), allocatable :: dnF
    real(DP), dimension(:,:,:), allocatable :: div_uS
    real(DP), dimension(:,:,:), allocatable :: det_grad_uS


    dWeight  = rcollection%DquickAccess(1)
    drhoFR   = rcollection%DquickAccess(2)
    dnF0     = rcollection%DquickAccess(4)



    dnS0     = 1.0_DP-dnF0

    p_rvector => rcollection%p_RvectorQuickAccess1

    allocate (uS1_x(ubound(Dcoefficients,2),ubound(Dcoefficients,3),1))
    allocate (uS1_y(ubound(Dcoefficients,2),ubound(Dcoefficients,3),1))
    allocate (uS2_x(ubound(Dcoefficients,2),ubound(Dcoefficients,3),1))
    allocate (uS2_y(ubound(Dcoefficients,2),ubound(Dcoefficients,3),1))
    allocate (div_uS(ubound(Dcoefficients,2),ubound(Dcoefficients,3),1))
    allocate (dnF(ubound(Dcoefficients,2),ubound(Dcoefficients,3),1))
    allocate (det_grad_uS(ubound(Dcoefficients,2),ubound(Dcoefficients,3),1))

! uS1_x
    call fevl_evaluate_sim (DER_DERIV_X, uS1_x(:,:,1), &
        p_rvector%Rvectorblock(1), Dpoints, &
        rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

! uS1_y
    call fevl_evaluate_sim (DER_DERIV_Y, uS1_y(:,:,1), &
        p_rvector%Rvectorblock(1), Dpoints, &
        rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
! uS2_x
    call fevl_evaluate_sim (DER_DERIV_X, uS2_x(:,:,1), &
        p_rvector%Rvectorblock(2), Dpoints, &
        rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
! uS2_y
    call fevl_evaluate_sim (DER_DERIV_Y, uS2_y(:,:,1), &
        p_rvector%Rvectorblock(2), Dpoints, &
        rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)



!  |grad uS | :  is the determinant value of grad uS
    do i=1,ubound(Dcoefficients,3)
      do j=1,ubound(Dcoefficients,2)
        det_grad_uS(j,i,1) = uS1_x(j,i,1)*uS2_y(j,i,1) - uS2_x(j,i,1)*uS1_y(j,i,1)
      end do
    end do

    deallocate (uS1_y)
    deallocate (uS2_x)

!     div_uS = uS1_x + uS2_y
    do i=1,ubound(Dcoefficients,3)
      do j=1,ubound(Dcoefficients,2)
        div_uS(j,i,1) = uS1_x(j,i,1) + uS2_y(j,i,1)
      end do
    end do

    deallocate (uS1_x)
    deallocate (uS2_y)


! nF = 1 - nS0 (1-div uS + |grad uS | )
    do i=1,ubound(Dcoefficients,3)
      do j=1,ubound(Dcoefficients,2)
        dnF(j,i,1) = 1.0_DP - dnS0*(1.0_DP - div_uS(j,i,1) + det_grad_uS(j,i,1))
      end do
    end do

    deallocate (det_grad_uS)
    deallocate (div_uS)




    do i=1,ubound(Dcoefficients,3)
      do j=1,ubound(Dcoefficients,2)
        Dcoefficients(1,j,i) = dWeight*drhoFR/dnF(j,i,1)
      end do
    end do

! print*,dkF

    deallocate (dnF)

  end subroutine

! ***************************************************************************
  !<subroutine>

  subroutine coeff_K55 (rdiscretisationTrial,rdiscretisationTest,rform, &
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTrial
    
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints
    
    ! An array accepting the DOF`s on all elements in the trial space.
    ! DIMENSION(#local DOF`s in trial space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTrial
    
    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rvector
    real(DP) :: dWeight
    real(DP) :: drhoFR
    real(DP) :: dnF0
    real(DP) :: dnS0
    real(DP) :: dkF0
    real(DP) :: dgammaFR

    integer :: i, j, kappa

    real(DP), dimension(:,:,:), allocatable :: div_vS
    real(DP), dimension(:,:,:), allocatable :: uS1_x
    real(DP), dimension(:,:,:), allocatable :: uS1_y
    real(DP), dimension(:,:,:), allocatable :: uS2_x
    real(DP), dimension(:,:,:), allocatable :: uS2_y
    real(DP), dimension(:,:,:), allocatable :: vS1_x
    real(DP), dimension(:,:,:), allocatable :: vS2_y
    real(DP), dimension(:,:,:), allocatable :: dnF
    real(DP), dimension(:,:,:), allocatable :: dkF
    real(DP), dimension(:,:,:), allocatable :: div_uS
    real(DP), dimension(:,:,:), allocatable :: det_grad_uS


    kappa    = rcollection%IquickAccess(1)

    dWeight  = rcollection%DquickAccess(1)
    drhoFR   = rcollection%DquickAccess(2)
    dnF0     = rcollection%DquickAccess(3)
    dkF0     = rcollection%DquickAccess(4)

    dgammaFR = 10.0_DP*drhoFR
    dnS0     = 1.0_DP-dnF0

    p_rvector => rcollection%p_RvectorQuickAccess1

    allocate (uS1_x(ubound(Dcoefficients,2),ubound(Dcoefficients,3),1))
    allocate (uS1_y(ubound(Dcoefficients,2),ubound(Dcoefficients,3),1))
    allocate (uS2_x(ubound(Dcoefficients,2),ubound(Dcoefficients,3),1))
    allocate (uS2_y(ubound(Dcoefficients,2),ubound(Dcoefficients,3),1))
    allocate (vS1_x(ubound(Dcoefficients,2),ubound(Dcoefficients,3),1))
    allocate (vS2_y(ubound(Dcoefficients,2),ubound(Dcoefficients,3),1))
    allocate (div_vS(ubound(Dcoefficients,2),ubound(Dcoefficients,3),1))
    allocate (div_uS(ubound(Dcoefficients,2),ubound(Dcoefficients,3),1))
    allocate (dnF(ubound(Dcoefficients,2),ubound(Dcoefficients,3),1))
    allocate (dkF(ubound(Dcoefficients,2),ubound(Dcoefficients,3),1))
    allocate (det_grad_uS(ubound(Dcoefficients,2),ubound(Dcoefficients,3),1))

! uS1_x
    call fevl_evaluate_sim (DER_DERIV_X, uS1_x(:,:,1), &
        p_rvector%Rvectorblock(1), Dpoints, &
        rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

! uS1_y
    call fevl_evaluate_sim (DER_DERIV_Y, uS1_y(:,:,1), &
        p_rvector%Rvectorblock(1), Dpoints, &
        rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
! uS2_x
    call fevl_evaluate_sim (DER_DERIV_X, uS2_x(:,:,1), &
        p_rvector%Rvectorblock(2), Dpoints, &
        rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
! uS2_y
    call fevl_evaluate_sim (DER_DERIV_Y, uS2_y(:,:,1), &
        p_rvector%Rvectorblock(2), Dpoints, &
        rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
! vS1_x
    call fevl_evaluate_sim (DER_DERIV_X, vS1_x(:,:,1), &
        p_rvector%Rvectorblock(3), Dpoints, &
        rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
! vS1_y
    call fevl_evaluate_sim (DER_DERIV_Y, vS2_y(:,:,1), &
        p_rvector%Rvectorblock(4), Dpoints, &
        rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

!     div_vS = vS1_x + vS2_y
    do i=1,ubound(Dcoefficients,3)
      do j=1,ubound(Dcoefficients,2)
        div_vS(j,i,1) = vS1_x(j,i,1) + vS2_y(j,i,1)
      end do
    end do

    deallocate (vS1_x)
    deallocate (vS2_y)

!     div_uS = uS1_x + uS2_y
    do i=1,ubound(Dcoefficients,3)
      do j=1,ubound(Dcoefficients,2)
        div_uS(j,i,1) = uS1_x(j,i,1) + uS2_y(j,i,1)
      end do
    end do

!  |grad uS | :  is the determinant value of grad uS
    do i=1,ubound(Dcoefficients,3)
      do j=1,ubound(Dcoefficients,2)
        det_grad_uS(j,i,1) = uS1_x(j,i,1)*uS2_y(j,i,1) - uS2_x(j,i,1)*uS1_y(j,i,1)
      end do
    end do

    deallocate (uS1_x)
    deallocate (uS1_y)
    deallocate (uS2_x)
    deallocate (uS2_y)


! nF = 1 - nS0 (1-div uS + |grad uS | )
    do i=1,ubound(Dcoefficients,3)
      do j=1,ubound(Dcoefficients,2)
        dnF(j,i,1) = 1.0_DP - dnS0*(1.0_DP - div_uS(j,i,1) + det_grad_uS(j,i,1))
      end do
    end do

    deallocate (det_grad_uS)
    deallocate (div_uS)


!  kF
    do i=1,ubound(Dcoefficients,3)
      do j=1,ubound(Dcoefficients,2)
        dkF(j,i,1) = dkF0*(dnF(j,i,1)/dnF0)**kappa
      end do
    end do
!  print*,dkF0

    do i=1,ubound(Dcoefficients,3)
      do j=1,ubound(Dcoefficients,2)
        Dcoefficients(1,j,i) = dWeight*(dgammaFR/dkF(j,i,1) + drhoFR/dnF(j,i,1)*div_vS(j,i,1))
      end do
    end do

! print*,dkF

    deallocate (dnF)
    deallocate (dkF)
    deallocate (div_vS)

  end subroutine

! ! ******************************************************************************
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTrial
        
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints
    
    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTrial
    
    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    ! Routine is not called anyway in case of constant coefficients!
    Dcoefficients = 1.0_DP

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_uSx (rdiscretisation,rform, &
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    ! local variables
    real(DP) :: dtime
    
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      dtime = rcollection%Dquickaccess(1)
    else
      dtime = 0.0_DP
    end if
!        Dcoefficients(:,:,:) = 0.0_DP


         Dcoefficients(1,:,:) = -301.0_DP
!        Dcoefficients(:,:,:) = 0.0_DP    

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_uSy (rdiscretisation,rform, &
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    ! local variables
    real(DP) :: dtime
    
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      dtime = rcollection%Dquickaccess(1)
    else
      dtime = 0.0_DP
    end if

!          Dcoefficients(:,:,:) = 0.0_DP

         Dcoefficients(1,:,:) = -300.0_DP
!        Dcoefficients(:,:,:) = 0.0_DP

  end subroutine


  ! ***************************************************************************

!<subroutine> 

  subroutine coeff_RHS_vSx (rdiscretisation,rform, &
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    ! local variables
    real(DP) :: dtime
    
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      dtime = rcollection%Dquickaccess(1)
    else
      dtime = 0.0_DP
    end if
    
    Dcoefficients(:,:,:) = 0.0_DP

  end subroutine


  ! ***************************************************************************

!<subroutine> 

  subroutine coeff_RHS_vSy (rdiscretisation,rform, &
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    ! local variables
    real(DP) :: dtime
    
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      dtime = rcollection%Dquickaccess(1)
    else
      dtime = 0.0_DP
    end if
    
    Dcoefficients(:,:,:) = 0.0_DP

  end subroutine


  ! ***************************************************************************

!<subroutine> better to be renamed to  RHS_wSx (but not necessary)

  subroutine coeff_RHS_vFx (rdiscretisation,rform, &
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    ! local variables
    real(DP) :: dtime
    
!     ! In a nonstationary simulation, one can get the simulation time
!     ! with the quick-access array of the collection.
!     if (present(rcollection)) then
!       dtime = rcollection%Dquickaccess(1)
!     else
!       dtime = 0.0_DP
!     end if

     Dcoefficients(1,:,:) = -1.0_DP + 16000.0_DP*Dpoints(1,:,:)/ &
     (33.0_DP*Dpoints(1,:,:)-3.0_DP*Dpoints(1,:,:)*Dpoints(2,:,:)+ &
     30.0_DP*Dpoints(2,:,:)+70.0_DP)
!    Dcoefficients(:,:,:) = 0.0_DP

  end subroutine


  ! ***************************************************************************

!<subroutine>   better to be renamed to  RHS_wSy (not necessary)

  subroutine coeff_RHS_vFy (rdiscretisation,rform, &
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    ! local variables
    real(DP) :: dtime
    
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      dtime = rcollection%Dquickaccess(1)
    else
      dtime = 0.0_DP
    end if

     Dcoefficients(1,:,:) = -16000.0_DP*Dpoints(2,:,:)/ &
     (33.0_DP*Dpoints(1,:,:)-3.0_DP*Dpoints(1,:,:)*Dpoints(2,:,:)+ &
     30.0_DP*Dpoints(2,:,:)+70.0_DP)
!    Dcoefficients(:,:,:) = 0.0_DP

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_p (rdiscretisation,rform, &
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
    ! pressure/divergence part of the right hand side vector.
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    ! local variables
    real(DP) :: dtime
    
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      dtime = rcollection%Dquickaccess(1)
    else
      dtime = 0.0_DP
    end if
    
    Dcoefficients(:,:,:) = 0.0_DP

  end subroutine

  ! ***************************************************************************


!  ***************************************************************************************

!/***/  added by Obaid for surface force

! ****************************************************************************************


!<subroutine>
  subroutine RHS_2D_surf (rdiscretisation, rform, nel, nptsPerEl, &
                               Dpoints, ibct, DpointPar, IdofsTest, rdomainIntSubset, &
                               Dcoefficients, rcollection)
    
    use fsystem
    use boundary
    use collection
    use scalarpde
    use domainintegration
    use spatialdiscretisation

!<description>
    ! This subroutine is called during the vector assembly. It computes the function g
    ! in the linear form (usually L(v) = (f,v)_0 + (g,v)_N), i.e. the
    ! surface forces in x- or y-direction (nonzero Neumann boundary conditions).
!</description>

!<input>
    ! discretisation structure
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation
    
    ! linear form to be evaluated
    type(t_linearForm), intent(in) :: rform
    
    ! number of elements where the coefficients are to be computed
    integer, intent(in) :: nel
    
    ! number of points per element where the coefficients are to be computed
    integer, intent(in) :: nptsPerEl
    
    ! array of all points on all elements where coefficients are to be computed
    ! DIMENSION(spatial dimension, nptsPerEl, nel)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints
    
    ! index of the boundary component that contains the points
    integer, intent(in) :: ibct
    
    ! parameter values (LENGTH PARAMETRISATION) of the points on the boundary component
    ! DIMENSION(nptsPerEl, nel)
    real(DP), dimension(:,:), intent(in) :: DpointPar
    
    ! array of test space DOF in all elements
    ! DIMENSION(#local DOF in test space, nel)
    integer, dimension(:,:), intent(in) :: IdofsTest
    
    ! structure providing more detailed information about the current element set
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! optional collection structure for additional information provided by the user
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! array of coefficients of the terms in the linear form for all given points on
    ! all given elements, DIMENSION(itermCount, nptsPerEl, nel)
    ! with itermCount the number of terms in the linear form
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients

    ! local variables
    real(DP) :: dtime
    
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      dtime = rcollection%Dquickaccess(1)
    else
      dtime = 0.0_DP
    end if

!</output>
!</subroutine>

    ! local variables
    !/***/ real(DP), dimension(:,:,:), pointer :: DstressTensor


!     ! in rcollection%IquickAccess(1) the current component is stored
!     icomp = rcollection%IquickAccess(1)
! 
!     ! in rcollection%IquickAccess(2) the current segment number is stored
!     iseg = rcollection%IquickAccess(2)


! ! taken from Marker paper: comparison of monolithic and splitting solution schemes ....
! ! Markert figure: 8
! ! you must go to ccgeneraldiscretisation.f90: line 1856 to specify the segment
! ! and line 2095-2097  to specify the corner parameter values of that segment
! ! you may need also to specify the magnification factors to scal the deformation
! ! in line 1623 & 1624 in ccpostprocessing.f90

    if (dtime .ge. 0.0_DP .and. dtime .le. 0.1_DP) then  		! inf-space
      Dcoefficients(1,:,:) = -1.0E4_DP*(1.0_DP-cos(20.0_DP*sys_pi*dtime))
    else
      Dcoefficients(1,:,:) = 0.0_DP
    end if

! Markert figure 3
! 6*y - (5*x)/4 + 5/8     6.0_DP*Dpoints(2,:,:)-5.0_DP*Dpoints(1,:,:)/4.0_DP+5.0_DP/8.0_DP
! x/4 - 9/8
!      Dcoefficients(1,:,:) = -1.0E3_DP*(1.0_DP-cos(20.0_DP*sys_pi*dtime))  ! Rectangle 2x10

!     if (dtime .ge. 0.0_DP .and. dtime .le. 0.04_DP) then                 ! Rectangle 21x10
!       Dcoefficients(1,:,:) = -1.0E5_DP*sin(25.0_DP*sys_pi*dtime)
!     else
!       Dcoefficients(1,:,:) = 0.0_DP
!     end if

  end subroutine RHS_2D_surf


! ****************************************************************************************


!<subroutine>

  subroutine coeff_AnalyticSolution_uSx (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This routine is called upon program start if ctypeInitialSolution=3.
    ! It returns analytical values for the X-velocity in the
    ! initial solution vector.
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    ! local variables
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
    !
    ! -----
    ! In the basic implementation, we just call ffunction_TargetX to get the
    ! values -- so the target function (which is normally used for
    ! calculating the error to a reference function) is also the definition
    ! for the initial solution.
    ! If necessary, this behaviour can be changed here.
    
    call ffunction_TargetuSx (DER_FUNC,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dcoefficients(1,:,:),rcollection)

  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine coeff_AnalyticSolution_vSx (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This routine is called upon program start if ctypeInitialSolution=3.
    ! It returns analytical values for the X-velocity in the
    ! initial solution vector.
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    ! local variables
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
    !
    ! -----
    ! In the basic implementation, we just call ffunction_TargetX to get the
    ! values -- so the target function (which is normally used for
    ! calculating the error to a reference function) is also the definition
    ! for the initial solution.
    ! If necessary, this behaviour can be changed here.
    
    call ffunction_TargetvSx (DER_FUNC,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dcoefficients(1,:,:),rcollection)

  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine coeff_AnalyticSolution_vFx (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This routine is called upon program start if ctypeInitialSolution=3.
    ! It returns analytical values for the X-velocity in the
    ! initial solution vector.
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    ! local variables
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
    !
    ! -----
    ! In the basic implementation, we just call ffunction_TargetX to get the
    ! values -- so the target function (which is normally used for
    ! calculating the error to a reference function) is also the definition
    ! for the initial solution.
    ! If necessary, this behaviour can be changed here.
    
    call ffunction_TargetvFx (DER_FUNC,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dcoefficients(1,:,:),rcollection)

  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine coeff_AnalyticSolution_uSy (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This routine is called upon program start if ctypeInitialSolution=3.
    ! It returns analytical values for the Y-velocity in the
    ! initial solution vector.
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    ! local variables
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
    !
    ! -----
    ! In the basic implementation, we just call ffunction_TargetX to get the
    ! values -- so the target function (which is normally used for
    ! calculating the error to a reference function) is also the definition
    ! for the initial solution.
    ! If necessary, this behaviour can be changed here.
    
    call ffunction_TargetuSy (DER_FUNC,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dcoefficients(1,:,:),rcollection)

  end subroutine



  ! ***************************************************************************

!<subroutine>

  subroutine coeff_AnalyticSolution_vSy (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This routine is called upon program start if ctypeInitialSolution=3.
    ! It returns analytical values for the Y-velocity in the
    ! initial solution vector.
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    ! local variables
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
    !
    ! -----
    ! In the basic implementation, we just call ffunction_TargetX to get the
    ! values -- so the target function (which is normally used for
    ! calculating the error to a reference function) is also the definition
    ! for the initial solution.
    ! If necessary, this behaviour can be changed here.
    
    call ffunction_TargetvSy (DER_FUNC,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dcoefficients(1,:,:),rcollection)

  end subroutine



  ! ***************************************************************************

!<subroutine>

  subroutine coeff_AnalyticSolution_vFy (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This routine is called upon program start if ctypeInitialSolution=3.
    ! It returns analytical values for the Y-velocity in the
    ! initial solution vector.
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    ! local variables
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
    !
    ! -----
    ! In the basic implementation, we just call ffunction_TargetX to get the
    ! values -- so the target function (which is normally used for
    ! calculating the error to a reference function) is also the definition
    ! for the initial solution.
    ! If necessary, this behaviour can be changed here.
    
    call ffunction_TargetvFy (DER_FUNC,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dcoefficients(1,:,:),rcollection)

  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine coeff_AnalyticSolution_P (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This routine is called upon program start if ctypeInitialSolution=3.
    ! It returns analytical values for the pressure in the
    ! initial solution vector.
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    ! local variables
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
    !
    ! -----
    ! In the basic implementation, we just call ffunction_TargetP to get the
    ! values -- so the target function (which is normally used for
    ! calculating the error to a reference function) is also the definition
    ! for the initial solution.
    ! If necessary, this behaviour can be changed here.
    
    call ffunction_TargetP (DER_FUNC,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dcoefficients(1,:,:),rcollection)
    
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine ffunction_TargetuSx (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  use basicgeometry
  use triangulation
  use collection
  use scalarpde
  use domainintegration
  
!<description>
  ! This subroutine is called during the postprocessing.
  ! It should return values of the analytical solution (if it is known).
  ! These are compared with the calculated solution to calculate the
  ! error in the X-velocity.
  !
  ! If the analytical solution is unknown, this routine does not make sense.
  ! In this case, error analysis should be deactivated in the .DAT files!
!</description>
  
!<input>
  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
  ! The result must be written to the Dvalue-array below.
  integer, intent(in)                                         :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(in)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(in)                                         :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! DIMENSION(NDIM2D,npointsPerElement,nelements)
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(in)  :: Dpoints

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
  integer, dimension(:,:), intent(in) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It is usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

  ! A pointer to a collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional      :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(out)                      :: Dvalues
!</output>
  
!</subroutine>

    ! local variables
    real(DP) :: dtime,dtimeMax
    integer :: itimedependence

    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      dtime = rcollection%Dquickaccess(1)
      dtimeMax = rcollection%Dquickaccess(3)
      itimedependence = rcollection%Iquickaccess(1)
    else
      itimedependence = 0
      dtime = 0.0_DP
      dtimeMax = 0.0_DP
    end if

!      Dvalues(:,:) = 0.0_DP
    
      select case (cderivative)
      case (DER_FUNC);     Dvalues(:,:) =  0.05_DP*Dpoints(1,:,:)**2
      case (DER_DERIV_X);  Dvalues(:,:) =  0.1_DP*Dpoints(1,:,:)
      case (DER_DERIV_Y);  Dvalues(:,:) =  0.0_DP
      case (DER_DERIV_XX); Dvalues(:,:) =  0.1_DP
      case (DER_DERIV_XY); Dvalues(:,:) =  0.0_DP
      case (DER_DERIV_YY); Dvalues(:,:) =  0.0_DP
      end select

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine ffunction_TargetvSx (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  use basicgeometry
  use triangulation
  use collection
  use scalarpde
  use domainintegration
  
!<description>
  ! This subroutine is called during the postprocessing.
  ! It should return values of the analytical solution (if it is known).
  ! These are compared with the calculated solution to calculate the
  ! error in the Y-velocity.
  !
  ! If the analytical solution is unknown, this routine does not make sense.
  ! In this case, error analysis should be deactivated in the .DAT files!
!</description>
  
!<input>
  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
  ! The result must be written to the Dvalue-array below.
  integer, intent(in)                                         :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(in)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(in)                                         :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! DIMENSION(NDIM2D,npointsPerElement,nelements)
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(in)  :: Dpoints

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
  integer, dimension(:,:), intent(in) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It is usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

  ! A pointer to a collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional      :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(out)                      :: Dvalues
!</output>
  
!</subroutine>

    ! local variables
    real(DP) :: dtime,dtimeMax
    integer :: itimedependence

    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      dtime = rcollection%Dquickaccess(1)
      dtimeMax = rcollection%Dquickaccess(3)
      itimedependence = rcollection%Iquickaccess(1)
    else
      itimedependence = 0
      dtime = 0.0_DP
      dtimeMax = 0.0_DP
    end if

     Dvalues(:,:) = 0.0_DP
    
!       select case (cderivative)
!       case (DER_FUNC);     Dvalues(:,:) = Dpoints(2,:,:)
!       case (DER_DERIV_X);  Dvalues(:,:) = 0.0_DP
!       case (DER_DERIV_Y);  Dvalues(:,:) = 1.0_DP
!       case (DER_DERIV_XX); Dvalues(:,:) = 0.0_DP
!       case (DER_DERIV_XY); Dvalues(:,:) = 0.0_DP
!       case (DER_DERIV_YY); Dvalues(:,:) = 0.0_DP
!       end select

  end subroutine


  ! ***************************************************************************
  
!<subroutine>  better to be replaced with ffunction_Targetwx (but not necessary)

  subroutine ffunction_TargetvFx (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  use basicgeometry
  use triangulation
  use collection
  use scalarpde
  use domainintegration
  
!<description>
  ! This subroutine is called during the postprocessing.
  ! It should return values of the analytical solution (if it is known).
  ! These are compared with the calculated solution to calculate the
  ! error in the X-velocity.
  !
  ! If the analytical solution is unknown, this routine does not make sense.
  ! In this case, error analysis should be deactivated in the .DAT files!
!</description>
  
!<input>
  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
  ! The result must be written to the Dvalue-array below.
  integer, intent(in)                                         :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(in)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(in)                                         :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! DIMENSION(NDIM2D,npointsPerElement,nelements)
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(in)  :: Dpoints

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
  integer, dimension(:,:), intent(in) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It is usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

  ! A pointer to a collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional      :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(out)                      :: Dvalues
!</output>
  
!</subroutine>

    ! local variables
    real(DP) :: dtime,dtimeMax
    integer :: itimedependence

    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      dtime = rcollection%Dquickaccess(1)
      dtimeMax = rcollection%Dquickaccess(3)
      itimedependence = rcollection%Iquickaccess(1)
    else
      itimedependence = 0
      dtime = 0.0_DP
      dtimeMax = 0.0_DP
    end if

!     Dvalues(:,:) = 0.0_DP

      select case (cderivative)
      case (DER_FUNC);     Dvalues(:,:) =  Dpoints(1,:,:)
      case (DER_DERIV_X);  Dvalues(:,:) =  1.0_DP
      case (DER_DERIV_Y);  Dvalues(:,:) =  0.0_DP
      case (DER_DERIV_XX); Dvalues(:,:) =  0.0_DP
      case (DER_DERIV_XY); Dvalues(:,:) =  0.0_DP
      case (DER_DERIV_YY); Dvalues(:,:) =  0.0_DP
      end select
    
  end subroutine


  ! ***************************************************************************
  
!<subroutine>

  subroutine ffunction_TargetuSy (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  use basicgeometry
  use triangulation
  use collection
  use scalarpde
  use domainintegration
  
!<description>
  ! This subroutine is called during the postprocessing.
  ! It should return values of the analytical solution (if it is known).
  ! These are compared with the calculated solution to calculate the
  ! error in the X-velocity.
  !
  ! If the analytical solution is unknown, this routine does not make sense.
  ! In this case, error analysis should be deactivated in the .DAT files!
!</description>
  
!<input>
  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
  ! The result must be written to the Dvalue-array below.
  integer, intent(in)                                         :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(in)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(in)                                         :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! DIMENSION(NDIM2D,npointsPerElement,nelements)
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(in)  :: Dpoints

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
  integer, dimension(:,:), intent(in) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It is usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

  ! A pointer to a collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional      :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(out)                      :: Dvalues
!</output>
  
!</subroutine>

    ! local variables
    real(DP) :: dtime,dtimeMax
    integer :: itimedependence

    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      dtime = rcollection%Dquickaccess(1)
      dtimeMax = rcollection%Dquickaccess(3)
      itimedependence = rcollection%Iquickaccess(1)
    else
      itimedependence = 0
      dtime = 0.0_DP
      dtimeMax = 0.0_DP
    end if

!      Dvalues(:,:) = 0.0_DP
    
      select case (cderivative)
      case (DER_FUNC);     Dvalues(:,:) =  0.05_DP*Dpoints(2,:,:)**2 - 0.1_DP*Dpoints(2,:,:)
      case (DER_DERIV_Y);  Dvalues(:,:) =  0.1_DP*Dpoints(2,:,:)-0.1_DP
      case (DER_DERIV_X);  Dvalues(:,:) =  0.0_DP
      case (DER_DERIV_YY); Dvalues(:,:) =  0.1_DP
      case (DER_DERIV_XY); Dvalues(:,:) =  0.0_DP
      case (DER_DERIV_XX); Dvalues(:,:) =  0.0_DP
      end select

  end subroutine


  ! ***************************************************************************
  
!<subroutine>

  subroutine ffunction_TargetvSy (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  use basicgeometry
  use triangulation
  use collection
  use scalarpde
  use domainintegration
  
!<description>
  ! This subroutine is called during the postprocessing.
  ! It should return values of the analytical solution (if it is known).
  ! These are compared with the calculated solution to calculate the
  ! error in the X-velocity.
  !
  ! If the analytical solution is unknown, this routine does not make sense.
  ! In this case, error analysis should be deactivated in the .DAT files!
!</description>
  
!<input>
  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
  ! The result must be written to the Dvalue-array below.
  integer, intent(in)                                         :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(in)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(in)                                         :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! DIMENSION(NDIM2D,npointsPerElement,nelements)
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(in)  :: Dpoints

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
  integer, dimension(:,:), intent(in) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It is usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

  ! A pointer to a collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional      :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(out)                      :: Dvalues
!</output>
  
!</subroutine>

    ! local variables
    real(DP) :: dtime,dtimeMax
    integer :: itimedependence

    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      dtime = rcollection%Dquickaccess(1)
      dtimeMax = rcollection%Dquickaccess(3)
      itimedependence = rcollection%Iquickaccess(1)
    else
      itimedependence = 0
      dtime = 0.0_DP
      dtimeMax = 0.0_DP
    end if

    Dvalues(:,:) = 0.0_DP
    
!       select case (cderivative)
!       case (DER_FUNC);     Dvalues(:,:) = Dpoints(1,:,:)
!       case (DER_DERIV_X);  Dvalues(:,:) = 1.0_DP
!       case (DER_DERIV_Y);  Dvalues(:,:) = 0.0_DP
!       case (DER_DERIV_XX); Dvalues(:,:) = 0.0_DP
!       case (DER_DERIV_XY); Dvalues(:,:) = 0.0_DP
!       case (DER_DERIV_YY); Dvalues(:,:) = 0.0_DP
!       end select

  end subroutine


  ! ***************************************************************************
  
!<subroutine>  better to be replaced with ffunction_Targetwy (but not necessary)

  subroutine ffunction_TargetvFy (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  use basicgeometry
  use triangulation
  use collection
  use scalarpde
  use domainintegration
  
!<description>
  ! This subroutine is called during the postprocessing.
  ! It should return values of the analytical solution (if it is known).
  ! These are compared with the calculated solution to calculate the
  ! error in the X-velocity.
  !
  ! If the analytical solution is unknown, this routine does not make sense.
  ! In this case, error analysis should be deactivated in the .DAT files!
!</description>
  
!<input>
  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
  ! The result must be written to the Dvalue-array below.
  integer, intent(in)                                         :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(in)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(in)                                         :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! DIMENSION(NDIM2D,npointsPerElement,nelements)
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(in)  :: Dpoints

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
  integer, dimension(:,:), intent(in) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It is usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

  ! A pointer to a collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional      :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(out)                      :: Dvalues
!</output>
  
!</subroutine>

    ! local variables
    real(DP) :: dtime,dtimeMax
    integer :: itimedependence

    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      dtime = rcollection%Dquickaccess(1)
      dtimeMax = rcollection%Dquickaccess(3)
      itimedependence = rcollection%Iquickaccess(1)
    else
      itimedependence = 0
      dtime = 0.0_DP
      dtimeMax = 0.0_DP
    end if

    !Dvalues(:,:) = 0.0_DP
    
      select case (cderivative)
      case (DER_FUNC);     Dvalues(:,:) =  -Dpoints(2,:,:)
      case (DER_DERIV_Y);  Dvalues(:,:) =  -1.0_DP
      case (DER_DERIV_X);  Dvalues(:,:) =   0.0_DP
      case (DER_DERIV_YY); Dvalues(:,:) =   0.0_DP
      case (DER_DERIV_XY); Dvalues(:,:) =   0.0_DP
      case (DER_DERIV_XX); Dvalues(:,:) =   0.0_DP
      end select

  end subroutine


  ! ***************************************************************************
  
!<subroutine>

  subroutine ffunction_TargetP (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  use basicgeometry
  use triangulation
  use collection
  use scalarpde
  use domainintegration
  
!<description>
  ! This subroutine is called during the postprocessing.
  ! It should return values of the analytical solution (if it is known).
  ! These are compared with the calculated solution to calculate the
  ! error in the pressure
  !
  ! If the analytical solution is unknown, this routine does not make sense.
  ! In this case, error analysis should be deactivated in the .DAT files!
!</description>
  
!<input>
  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
  ! The result must be written to the Dvalue-array below.
  integer, intent(in)                                         :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(in)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(in)                                         :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! DIMENSION(NDIM2D,npointsPerElement,nelements)
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(in)  :: Dpoints

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
  integer, dimension(:,:), intent(in) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It is usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

  ! A pointer to a collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional      :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(out)                      :: Dvalues
!</output>
  
!</subroutine>

    ! local variables
    real(DP) :: dtime,dtimeMax
    integer :: itimedependence

    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      dtime = rcollection%Dquickaccess(1)
      dtimeMax = rcollection%Dquickaccess(3)
      itimedependence = rcollection%Iquickaccess(1)
    else
      itimedependence = 0
      dtime = 0.0_DP
      dtimeMax = 0.0_DP
    end if

    !Dvalues(:,:) = 0.0_DP
    
      select case (cderivative)
      case (DER_FUNC);     Dvalues(:,:) = 0.5_DP-Dpoints(1,:,:)
      case (DER_DERIV_X);  Dvalues(:,:) = -1.0_DP
      case (DER_DERIV_Y);  Dvalues(:,:) = 0.0_DP
      case (DER_DERIV_XX); Dvalues(:,:) = 0.0_DP
      case (DER_DERIV_XY); Dvalues(:,:) = 0.0_DP
      case (DER_DERIV_YY); Dvalues(:,:) = 0.0_DP
      end select

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
  character(LEN=*), intent(in) :: sexpressionName
  
  ! Solution component that is currently being processed.
  ! 1 = uS1, 2 = uS2, 3 = vS1, 4 = vS2, 5 = vF1, 6 = vF2, 
  integer, intent(in) :: icomponent
  
  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  
  ! Boundary region that is currently being processed.
  type(t_boundaryRegion), intent(in)                          :: rboundaryRegion
  
  ! Current parameter value of the point on the boundary.
  ! 0-1-parametrisation.
  real(DP), intent(in)                                        :: dwhere
    
  ! Optional: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional      :: rcollection
!</input>

!<output>
  ! Return value of the expression. May be SYS_INFINITY_DP.
  real(DP), intent(out) :: dvalue
!</output>
  
!</subroutine>

    ! local variables
     REAL(DP) :: dtime
     REAL(DP) :: dx,dy
    !
    ! To get the X/Y-coordinates of the boundary point, use:
    !
     CALL boundary_getCoords(rdiscretisation%p_rboundary, &
         rboundaryRegion%iboundCompIdx, dwhere, dx, dy)
    !
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    !
     dtime = 0.0_DP
     IF (PRESENT(rcollection)) dtime = rcollection%Dquickaccess(1)

     if (icomponent .EQ. 1) then
       dvalue =  0.05_DP*dx**2
     elseif (icomponent .EQ. 2) then
       dvalue  =  0.05_DP*dy**2-0.1_DP*dy
     elseif (icomponent .EQ. 3) then
      dvalue = 0.0_DP
     elseif (icomponent .EQ. 4) then
       dvalue = 0.0_DP
     elseif (icomponent .EQ. 5) then
       dvalue = dx
     elseif (icomponent .EQ. 6) then
        dvalue = -dy
     else
      !dvalue = 0.5_DP-dx
    end if

  end subroutine

  ! ***************************************************************************
  ! Values in a fictitious boundary component:

!<subroutine>

  subroutine getBoundaryValuesFBC (Icomponents,rdiscretisation,&
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
    integer, dimension(:), intent(in)                           :: Icomponents
  
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_blockDiscretisation), intent(in)                     :: rdiscretisation
    
    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), optional, intent(inout)                 :: rcollection

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
    type(t_discreteFBCevaluation), dimension(:), intent(inout) :: Revaluation
  !</inputoutput>
    
  !</subroutine>

    ! Note: the definition of (analytic) fictitious boundary components
    ! is performed in 'cc_parseFBDconditions'.
    !
    ! By default, fictitious boundary handling is switched off!
    ! To switch the handling on, uncomment tge bcasm_XXXX-call
    ! in cc_parseFBDconditions!!!
    
    ! local variables
    real(DP) :: ddistance, dxcenter, dycenter, dradius, dx, dy
    real(DP), dimension(:,:), pointer :: p_DvertexCoordinates
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    type(t_triangulation), pointer :: p_rtriangulation
    integer :: ipoint,idx
    
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
    integer, intent(in) :: cinfoNeeded
    
    ! Identifier. Either vertex number, edge number or element number,
    ! depending on cinfoNeeded.
    integer, intent(in) :: iwhere
    
    ! Array with coordinates of all corner vertices (DCORVG)
    real(DP), dimension(:,:), intent(in) :: DvertexCoords
    
    ! Array with numbers of corner coordinates for all elements (KVERT)
    integer, dimension(:,:), intent(in) :: IverticesAtElement
    
    ! Array with numbers of points adjacent to all edges
    integer, dimension(:,:), intent(in) :: IverticesAtEdge
    
    ! Number of vertices in the triangulation
    integer, intent(in) :: NVT
    
    ! Output: X-coordinate
    real(DP), intent(out) :: dx
    
    ! Output: Y-coordinate
    real(DP), intent(out) :: dy
    
      ! local variables
      real(DP) :: dm1,dm2
      integer :: i,j
      integer :: iv1,iv2
    
      ! Let us see, what do we have...
      select case (cinfoNeeded)
      case (DISCFBC_NEEDFUNC)
        ! That is easy
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
        call output_line ('Unsupported coordinate type to return.!', &
            OU_CLASS_ERROR,OU_MODE_STD,'getBoundaryValuesFBC')
        call sys_halt()
      end select
      
    end subroutine

  end subroutine
! 
!   ! ***************************************************************************
!!\***\ should be removed
  ! Moving frame formulation.

!<subroutine>
!
  subroutine getMovingFrameVelocity (Dvelocity,Dacceleration,rcollection)

!<description>
  ! This routine is called by the framework if the moving frame formulation
  ! is activated by imovingFrame = 1. It has to return the current velocity
  ! and the current acceleration of the moving frame in the current timestep.
!</description>
    
!<inputoutput>
  ! Optional: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), optional, intent(inout) :: rcollection
!</inputoutput>

!<output>
  ! Current velocity of the moving frame.
  real(DP), dimension(:), intent(out) :: Dvelocity
  
  ! Current acceleration of the moving frame.
  real(DP), dimension(:), intent(out) :: Dacceleration
!</output>
    
!</subroutine>
    
    ! local variables
    real(DP) :: dtime
    
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      dtime = rcollection%Dquickaccess(1)
    else
      dtime = 0.0_DP
    end if

    ! Default implementation: Velocity and acceleration is zero in all dimensinons.
    Dvelocity(:) = 0.0_DP
    
    Dacceleration(:) = 0.0_DP

  end subroutine

! *****************************************************************

!<subroutine>

  subroutine getNonconstantViscosity (cterm,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset, &
                Dcoefficients,rvelocity,rcollection)
  
  use basicgeometry
  use triangulation
  use scalarpde
  use domainintegration
  use spatialdiscretisation
  use collection
  
!<description>
  ! This subroutine is called during the calculation of the SD operator.
  ! It allows to calculate a user defined viscosity coefficient
  ! in case of a nonconstant viscosity.
  !
  ! The routine accepts a set of elements and a set of points on these
  ! elements (cubature points) in in real coordinates.
  ! According to the terms in the linear form, the routine has to compute
  ! simultaneously for all these points.
!</description>
  
!<input>
  ! Term which is to be computed.
  ! =0: Calculate the $\nu$ values in front of the Laplace.
  ! =1: Calculate the $\alpha$ values in front of the Mass matrix.
  integer, intent(in) :: cterm

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(in)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(in)                                         :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! DIMENSION(NDIM2D,npointsPerElement,nelements)
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(in)  :: Dpoints

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
  integer, dimension(:,:), intent(in) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It is usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

  ! Current velocity vector.
  type(t_vectorBlock), intent(in) :: rvelocity

  ! Optional: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional      :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the coefficients
  ! in all the points specified in Dpoints.
  ! cterm specifies what to evaluate.
  real(DP), dimension(:,:), intent(out) :: Dcoefficients
!</output>
  
!</subroutine>

    ! Add a viscolity here. The parameter [CC-DISCRETISATION].cviscoModel
    ! decides on whether this coefficient is evaluated or not.
    !
    ! WARNING: For a nonconstant coefficient, the extended assembly
    ! method must be activated!!! (Parameter [CC-DISCRETISATION].iupwind)
    Dcoefficients(:,:) = 0.0_DP

  end subroutine

! *****************************************************************
  
!<subroutine>
  subroutine calcAdaptiveTimestep (dtstep,dtimeInit,dtime,isolverStatus,rcollection)
!</subroutine>
  
!<description>
  ! Is called by the framework to calculate a new timestep size if
  ! cadaptiveTimeStepping=-1.
!</description>

!<input>
  ! Initial time
  real(dp), intent(in) :: dtimeInit
  
  ! Current simulation time
  real(dp), intent(in) :: dtime

  ! Status of the solver
  !  TADTS_SST_NLFAIL           = failure of the nonlinear solver
  !  TADTS_SST_NLPRECFAIL       = failure of the nonlinear solver and
  !                               preconditioner in the nonlinear solver
  !  TADTS_SST_NLINCOMPLETE     = nonlinear solver did not converge completely
  !  TADTS_SST_NLPREDFAIL       = failure of the nonlinear solver and preconditioner
  !                               in the nonlinear solver during the predictor step
  !  TADTS_SST_NLPREDPRECFAIL   = failure of the nonlinear solver during
  !                               the predictor step
  !  TADTS_SST_NLPREDINCOMPLETE = nonlinear solver in the predictor step did
  !                               not converge completely
  integer(I32), intent(in) :: isolverStatus

  ! Optional: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(in), optional :: rcollection
!</input>

!<inputoutput>
  ! Current timestep size. The routine may change this as needed.
  ! Must be > 0!
  real(dp), intent(inout) :: dtstep
!</inputoutput>
    
    ! Current implementation: No nothing.
    ! Then, a fixed timestep is used.
    
  end subroutine


! ***************************************************************************
  !<subroutine>

  ! ***************************************************************************
end module
