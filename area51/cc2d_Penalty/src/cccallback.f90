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
!# 4.) coeff_AnalyticSolution_X
!#     -> Returns analytical values for the desired flow field in X-direction.
!#     -> Is used for setting up the initial solution.
!#     -> In the basic implementation, this calls ffunction_TargetX. 
!#
!# 5.) coeff_AnalyticSolution_Y
!#     -> Returns analytical values for the desired flow field in Y-direction.
!#     -> Is used for setting up the initial solution.
!#     -> In the basic implementation, this calls ffunction_TargetY. 
!#
!# 6.) coeff_AnalyticSolution_P
!#     -> Returns analytical values for the desired pressure.
!#     -> Is used for setting up the initial solution.
!#     -> In the basic implementation, this calls ffunction_TargetP.
!#
!# 7.) ffunction_TargetX
!#     -> Returns analytical values for the desired flow field in X-direction.
!#     -> Is used for error analysis during the postprocessing.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 8.) ffunction_TargetY
!#     -> Returns analytical values for the desired flow field in Y-direction.
!#     -> Is used for error analysis during the postprocessing.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 9.) ffunction_TargetP
!#     -> Returns analytical values for the desired pressure.
!#     -> Is used for error analysis during the postprocessing.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 10.) getBoundaryValues
!#     -> Returns analytical values on the boundary of the
!#        problem to solve.
!#
!# 11.) getBoundaryValuesFBC
!#     -> Returns analytical values on the fictitious boundary components
!#     -> Corresponds to the interface defined in the file
!#        'intf_fbcassembly.inc'
!#
!# 12.) cc_initCollectForAssembly
!#     -> Is called prior to the assembly process.
!#     -> Stores some information from the problem structure to a collection
!#        such that it can be accessed in callback routines
!#
!# 13.) cc_doneCollectForAssembly
!#      -> Is called after the assembly process.
!#      -> Releases information stored in the collection by 
!#         cc_initCollectForAssembly.
!#
!# 14.) getMovingFrameVelocity
!#      -> If the moving frame formulation is activated, this routine
!#         returns the velocity and acceleration of the moving frame.
!#
!# 15.) getViscosity
!#      -> If nonconstant viscosity is activated, this routine calculates
!#         the viscosity.
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
  use basicgeometry
  use bcassembly
  use boundary
  use bilinearformevaluation
  use ccbasic
  use collection
  use cubature
  use derivatives
  use discretebc
  use domainintegration
  use io
  use linearformevaluation
  use linearsolver
  use matrixfilters
  use mprimitives
  use triangulation
  use scalarpde
  use storage
  use paramlist
  use vectorfilters

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
  type(t_problem), intent(in),target :: rproblem
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
    call collct_setvalue_parlst (rcollection, "PARLST", rproblem%rparamlist, .true.)

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

  subroutine coeff_RHS_x (rdiscretisation,rform, &
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

    integer :: iel,icup,ive,nve,iin
    real(dp) :: dxcenter, dycenter, dradius, ddist
    
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      dtime = rcollection%Dquickaccess(1)
    else
      dtime = 0.0_DP
    end if
    
    
    ! Definition of the circle
    dxcenter = 0.2
    dycenter = 0.2
    dradius  = 0.05

    ! loop over the elements and cubature points
    ! and assign the coefficients 
    do iel=1,nelements
      do icup=1,npointsPerElement
        
        ddist = sqrt( (Dpoints(1,icup,iel) - dxcenter)**2 + (Dpoints(2,icup,iel)-dycenter)**2)
        if(ddist .le. dradius)then
          Dcoefficients(1,icup,iel) = 0.0_dp
        else
          Dcoefficients(1,icup,iel) = 0.0_dp
        end if
      end do
    end do
    
    
   

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_y (rdiscretisation,rform, &
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

  subroutine coeff_AnalyticSolution_X (rdiscretisation,rform, &
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
    
    call ffunction_TargetX (DER_FUNC,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dcoefficients(1,:,:),rcollection)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_AnalyticSolution_Y (rdiscretisation,rform, &
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
    
    call ffunction_TargetY (DER_FUNC,rdiscretisation, &
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
    
    ! Example:
    ! IF (cderivative .EQ. DER_FUNC) THEN
    !   Dvalues(:,:) = (-dtime**2/100.+dtime/5.)*(Dpoints(1,:,:))
    ! END IF

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
    
    ! Example:
    ! IF (cderivative .EQ. DER_FUNC) THEN
    !   Dvalues(:,:) = (-dtime**2/100.+dtime/5.)*(-Dpoints(2,:,:))
    ! END IF

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

    Dvalues(:,:) = 0.0_DP

    ! Example:
    ! IF (cderivative .EQ. DER_FUNC) THEN
    !   ...
    ! END IF

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
  ! The routine is allowed to return SYS_INFINITY. The behaviour of this
  ! value depends on the type of boundary conditions. For Dirichlet
  ! boundary segments, all points where SYS_INFINITY is returned are
  ! treated as Neumann points.
!</description>
  
!<input>
  ! Name of the expression to be evaluated. This name is configured in the
  ! DAT file for the boundary conditions.
  character(LEN=*), intent(in) :: sexpressionName
  
  ! Solution component that is currently being processed. 
  ! 1 = X-velocity, 2 = y-velocity,...
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
  ! Return value of the expression. May be SYS_INFINITY.
  real(DP), intent(out) :: dvalue
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

    dvalue = 0.0_DP

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
    dxcenter = 0.2
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

  ! ***************************************************************************
  ! Moving frame formulation.

!<subroutine>

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
  
!****************************************************************************************
! 
  subroutine cc_lambda(rdiscretisationtrial,rdiscretisationtest,rform, &
                       nelements,npointsperelement,dpoints, &
                       idofstrial,idofstest,rdomainintsubset, &
                       dcoefficients,rcollection)
!
!****************************************************************************************
    
!<description>
! this subroutine is called during the matrix assembly. It has to compute the coefficients 
! in front of the terms of the bilinear form. The routine accepts a set of elements and a set
! of points on these elements (cubature points) in real coordinates. According to the terms 
! in the bilinear form, the routine has to compute simultaneously for all these points and all
! the terms in the bilinear form the corresponding coefficients in front of the terms.
!</description>
    
!<input>
! the discretisation structure that defines the basic shape of the
! triangulation with references to the underlying triangulation,
! analytic boundary boundary description etc.; trial space.
  type(t_spatialdiscretisation), intent(in)                   :: rdiscretisationtrial
! the discretisation structure that defines the basic shape of the
! triangulation with references to the underlying triangulation,
! analytic boundary boundary description etc.; test space.
  type(t_spatialdiscretisation), intent(in)                   :: rdiscretisationtest
! the bilinear form which is currently being evaluated:
  type(t_bilinearform), intent(in)                            :: rform
! number of elements, where the coefficients must be computed.
  integer, intent(in)                                         :: nelements
! number of points per element, where the coefficients must be computed
  integer, intent(in)                                         :: npointsperelement
! this is an array of all points on all the elements where coefficients are needed.
! remark: this usually coincides with rdomainsubset%p_dcubptsreal.
! dimension(dimension,npointsperelement,nelements)
  real(dp), dimension(:,:,:), intent(in)  :: dpoints
! an array accepting the dof's on all elements trial in the trial space.
! dimension(\#local dof's in trial space,nelements)
  integer, dimension(:,:), intent(in) :: idofstrial
! an array accepting the dof's on all elements trial in the trial space.
! dimension(\#local dof's in test space,nelements)
  integer, dimension(:,:), intent(in) :: idofstest
! this is a t_domainintsubset structure specifying more detailed information
! about the element set that is currently being integrated.
! it's usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainintsubset), intent(in)              :: rdomainintsubset
! OPTIONAL: a collection structure to provide additional 
! information to the coefficient routine. 
  type(t_collection), intent(inout), optional      :: rcollection
!</input>
  
!<output>
! a list of all coefficients in front of all terms in the bilinear form - for all given points 
! on all given elements. dimension(itermcount,npointsperelement,nelements) with itermcount the
! number of terms in the bilinear form.
  real(dp), dimension(:,:,:), intent(out)                      :: dcoefficients
!</output>
!</subroutine>

!<local variables>
  integer :: iel,ielreal,icup,ive,nve,iin,ipart,ivert,in,icount,i,ipenaltyact
  real(dp) :: dxcenter, dycenter, dradius, ddist,dlambda,dx,dy,dLocAreea,dElAreea
  real(dp), dimension(:,:), pointer :: p_dvertexcoordinates
  integer, dimension(:,:), pointer :: p_iverticesatelement
  type(t_particlecollection), pointer :: p_rparticlecollection
  type(t_geometryobject), pointer :: p_rgeometryobject
  type(t_parlist), pointer :: p_rparlst


!</local variables>

! Get the parameter list.
  p_rparlst => collct_getvalue_parlst (rcollection, "PARLST")

! Get the triangulation array for the point coordinates
    call storage_getbase_double2d (rdiscretisationtrial%p_rtriangulation%h_dvertexcoords,&
                                   p_dvertexcoordinates)
    call storage_getbase_int2d (rdiscretisationtrial%p_rtriangulation%h_iverticesatelement,&
                                   p_iverticesatelement)  
    
! Pointer to collect all variables about the object/s. We will get data about origin, shape,
! number of objects and so on, depending on the type of object (circle, square, ellipse,etc)
  p_rparticlecollection => collct_getvalue_particles(rcollection,'particles')

! For multiple objects, we need to treat them in the same way, so loop over the number of 
! particles. If we have only 1 particle, the loop has no effect.

  do ipart=1,p_rparticlecollection%nparticles

! Find the geometry of the object   
  p_rgeometryobject => p_rparticlecollection%p_rparticles(ipart)%rgeometryobject    

! Which method is used to implement the penalty matrix? That is choosen in the .dat file and 
! there only can be 2 ways at this moment: "full Lambda" and "fractional Lambda" methods.
! "Full" method is the classical 0/1 for the outside/inside cubature points in one element.
! "Fractional" method calculates an average Lambda as a fraction between element areea and 
! common areea between element and particle, This fraction is always in the interval [0;1]

  call parlst_getvalue_int (p_rparlst,'CC-DISCRETISATION','ipenaltyact',ipenaltyact,0)

! Loop over all elements and calculate the corresponding Lambda value
    do iel=1,nelements
      
    select case(ipenaltyact)
    case (0)  
      ! "Full Lambda" method
        do icup=1,npointsperelement 
          ! get the distance to the center
          call geom_isingeometry (p_rgeometryobject, (/dpoints(1,icup,iel),dpoints(2,icup,iel)/), iin)
          ! check if it is inside      
          if(iin .eq. 1)then 
            dcoefficients(1,icup,iel) = rform%dcoefficients(1) 
          else
            dcoefficients(1,icup,iel) = 0.0_dp
          end if
        end do
  
    case (1)
      !"Fractional Lambda" method
      ! The real element
      ielreal = rdomainintsubset%p_ielements(iel) 
      ! A counter for the inside vertex of the element    
      in = 0     
      ! Get vertices coordinates and check how many are in the object
      do ivert=1,rdiscretisationtrial%p_rtriangulation%nnve
      ! Local node coordinate of an element
        dx = p_dvertexcoordinates(1,p_iverticesatelement(ivert,ielreal))
        dy = p_dvertexcoordinates(2,p_iverticesatelement(ivert,ielreal))

      ! Check if the vertice is inside and count all inside points
        call geom_isingeometry (p_rgeometryobject,(/dx,dy/), iin)
      ! Count how many points are inside the object 
        if (iin .eq. 1) then
           in = in +1
        end if
      end do

      ! For an element with at least one node inside the object, calculate the
      ! intersections nodes between element edges and object. Calculate local areea,
      ! element area and the fractional lambda value (default value:
      ! dlambda=dcoefficients(1)
      if (in .gt. 0) then
        call conectelementobject(ielreal,rdiscretisationtrial%p_rtriangulation,p_rgeometryobject,dElAreea,&
                                 dLocAreea)
      ! calculate the equivalent lambda value
        dlambda = rform%dcoefficients(1)*dLocAreea/dElAreea              
      ! give to all cubature points same fractional lambda
        do icup=1,npointsperelement 
           dcoefficients(1,icup,iel) = dlambda
        end do 
      else !(for in > 0)
        do icup=1,npointsperelement 
           dcoefficients(1,icup,iel) = 0.0_dp
        end do
      end if    
    end select 
 
   end do !(loop over elements)
  end do !(loop over particles)
    
  end subroutine

!**********************************************************************************************************************
!
  subroutine conectElementObject(ielement,p_rtriangulation,p_rgeometryObject,ElAreea,LocAreea)
!
!**********************************************************************************************************************

!<input>
  integer, intent(IN) :: ielement
  type(t_triangulation), intent(IN), pointer    :: p_rtriangulation
  type(t_geometryObject),intent(IN), pointer  :: p_rgeometryObject  
!<output>
  real(DP),intent(OUT) :: ElAreea
  real(DP),intent(OUT) :: LocAreea
  
! local variable
  real(DP), dimension(:,:), pointer :: p_DvertexCoordinates
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  real(DP), dimension(:,:), allocatable :: PolyPoints
  
  integer  :: ivert,SlopeType,iin,icount,i
  real(DP) :: dxmax,dymax,dxmin,dymin,dxs,dys,dxe,dye,dxi1,dxi2,dyi1,dyi2, &
              dxcenter,dycenter,dradius,ddist,dslope,da,db,dc,ddiscr

!**********************************************************************************************************************

! Get the triangulation array for the point coordinates
                                 
  call storage_getbase_double2d (p_rtriangulation%h_DvertexCoords,p_DvertexCoordinates)
  call storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,p_IverticesAtElement)

! Get the parameters of the object (circle case)
  dxcenter = p_rgeometryObject%rcoord2d%dorigin(1)
  dycenter = p_rgeometryObject%rcoord2d%dorigin(2)
  dradius = p_rgeometryObject%rcircle%dradius
! Counter for polygon points
  icount = 0
! Initialise the array that saves the local polygon points
  allocate(PolyPoints(2,16))   
  PolyPoints(:,:)=0.0_dp    

! Calculate total area of the real element
  dxmax=maxval(p_DvertexCoordinates(1,p_IverticesAtElement(:,ielement)))
  dxmin=minval(p_DvertexCoordinates(1,p_IverticesAtElement(:,ielement)))
  dymax=maxval(p_DvertexCoordinates(2,p_IverticesAtElement(:,ielement)))
  dymin=minval(p_DvertexCoordinates(2,p_IverticesAtElement(:,ielement)))
  ElAreea=(dxmax-dxmin)*(dymax-dymin)

  do ivert = 1,p_rtriangulation%NNEE
  ! Pick up the values for start point and end point of an edge of the real element
    dxs = p_DvertexCoordinates(1,p_IverticesAtElement(ivert,ielement))
    dys = p_DvertexCoordinates(2,p_IverticesAtElement(ivert,ielement))
    dxe = p_DvertexCoordinates(1,p_IverticesAtElement(mod(ivert,4)+1,ielement))
    dye = p_DvertexCoordinates(2,p_IverticesAtElement(mod(ivert,4)+1,ielement))

  ! Check wheter the start point of the edge is already in the object. In affirmative case,
  ! save the point into the local polygon array and count the point. At the end we will have
  ! a counter of how many points the polygon has together with the coordinates of these points
    call geom_isInGeometry (p_rgeometryObject,(/dxs,dys/), iin)
          
    if (iin.eq.1) then
        icount=icount+1 ! add one point to the polygon
        PolyPoints(1,icount)=dxs ! x coordinate of the point
        PolyPoints(2,icount)=dys ! y coordinate of the point
    end if
              
  ! Check wheter the investigated edge is cutted by the object between start and end points
  ! We can have 3 cases: vertical edge, horizontal edge and oblic. First 2 cases are treated
  ! in here for a cartezian type mesh.
    SlopeType=0
    if (abs(dxe-dxs) .le. 1.0E-7) then
        SlopeType=1 ! vertical edge case
    else if (abs(dye-dys) .le. 1.0E-7) then
        SlopeType=2 ! horizontal edge case
    else
        SlopeType=3 ! oblic edge case                 
    end if
             
    select case (SlopeType)
    case (1) ! vertical edge
      ! First calculate the distance between the start point and center
      ! point of the object. if this distance is smaller then the radius
      ! then we may calculate conection points. otherwise, there is no 
      ! conection possible between the actual edge and object
      ddist=max(dxs,dxcenter)-min(dxs,dxcenter)
      if (ddist.le.dradius) then
          dyi1=dycenter-sqrt(dradius**2-(dxs-dycenter)**2)
          dyi2=dycenter+sqrt(dradius**2-(dxs-dycenter)**2)
        if ((dyi1.le.dymax).and.(dyi1.ge.dymin)) then
             icount=icount+1
             PolyPoints(1,icount)=dxs        
             PolyPoints(2,icount)=dyi1       
        end if
        if ((dyi2.le.dymax).and.(dyi2.ge.dymin)) then
             icount=icount+1
             PolyPoints(1,icount)=dxs        
             PolyPoints(2,icount)=dyi2        
        end if
      end if
               
    case (2) ! horizontal edge
      ddist=max(dys,dycenter)-min(dys,dycenter)
      if (ddist.le.dradius) then
          dxi1=dxcenter-sqrt(dradius**2-(dys-dycenter)**2)
          dxi2=dxcenter+sqrt(dradius**2-(dys-dycenter)**2)
        if ((dxi1.le.dxmax).and.(dxi1.ge.dxmin)) then
             icount=icount+1
             PolyPoints(1,icount)=dxi1        
             PolyPoints(2,icount)=dys        
        end if
        if ((dxi2.le.dxmax).and.(dxi2.ge.dxmin)) then
             icount=icount+1
             PolyPoints(1,icount)=dxi2        
             PolyPoints(2,icount)=dys        
        end if
      end if
               
    case (3) ! oblic edge
         dslope = (dye-dys)/(dxe-dxs)
         da = 1+dslope**2
         db = -2*dxcenter-2*dslope**2*dxs+2*dslope*(dys-dycenter)
         dc = dxcenter**2+dslope**2*dxs**2-2*dslope*dxs*(dys-dycenter)+(dys-dycenter)**2-dradius**2
         ddiscr = db**2-4*da*dc
         if (ddiscr .ge. 1.0E-7) then
           dxi1 = (-db + sqrt(ddiscr))/(2*da)
           dxi2 = (-db - sqrt(ddiscr))/(2*da)
           if ((dxi1.le.dxmax).and.(dxi1.ge.dxmin)) then
             icount=icount+1
             PolyPoints(1,icount)=dxi1        
             PolyPoints(2,icount)=dslope*(dxi1-dxs)+dys        
           end if
           if ((dxi2.le.dxmax).and.(dxi2.ge.dxmin)) then
             icount=icount+1
             PolyPoints(1,icount)=dxi2        
             PolyPoints(2,icount)=dslope*(dxi2-dxs)+dys        
           end if
         else if ((ddiscr .lt. 1.0E-7).and.(ddiscr .ge. 0.0_DP)) then  
           dxi1 = -db/(2*da)  
           if ((dxi1.le.dxmax).and.(dxi1.ge.dxmin)) then
             icount=icount+1
             PolyPoints(1,icount)=dxi1        
             PolyPoints(2,icount)=dslope*(dxi1-dxs)+dys        
           end if
         else
         end if
        
!      call output_line ('Unsupported edge slope type.')
!      call sys_halt()
             
    end select
  end do    

! Calculate the local area determined by the saved polygon points 
  LocAreea = 0.0_DP
  do i = 1,icount
     LocAreea = LocAreea + PolyPoints(1,i)*PolyPoints(2,mod(i,icount)+1)-PolyPoints(2,i)*PolyPoints(1,mod(i,icount)+1)
  end do
  LocAreea = 0.5*LocAreea
  
  deallocate (PolyPoints)              

  end subroutine

end module
