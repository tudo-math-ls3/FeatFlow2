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
!# 1.) user_coeff_Stokes
!#     -> Returns the coefficients for the Laplace matrix. This routine is
!#        only used if the problem to calculate has nonconstant coefficients!
!#        Otherwise the routine is dead.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientMatrixSc.inc'
!#
!# 1.) user_coeff_Pressure
!#     -> Returns the coefficients for the pressure matrix. This routine is
!#        only used if the problem to calculate has nonconstant coefficients!
!#        Otherwise the routine is dead.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientMatrixSc.inc'
!#
!# 2.) user_coeff_RHS
!#     -> Returns analytical values for the right hand side 
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 3.) user_coeff_Target
!#     -> Returns analytical values for the desired target function
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 4.) user_coeff_Reference
!#     -> Returns analytical values for the desired reference function
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 6.) user_getBoundaryValues
!#     -> Returns analitical values on the boundary of the
!#        problem to solve.
!#     -> Corresponds to the interface defined in the file
!#        'intf_bcassembly.inc'
!#
!# 7.) user_getBoundaryValuesFBC
!#     -> Returns analytical values on the fictitious boundary components
!#     -> Corresponds to the interface defined in the file
!#        'intf_fbcassembly.inc'
!#
!# 8.) user_initCollectForAssembly
!#     -> Is called prior to the assembly process.
!#     -> Stores some information from the problem structure to a collection
!#        such that it can be accessed in callback routines
!#
!# 9.) user_doneCollectForAssembly
!#     -> Is called after the assembly process.
!#     -> Releases information stored in the collection by 
!#        user_initCollectForAssembly.
!#
!# For nonstationary simulation, it might be neccessary in these routines
!# to access the current simulation time. Before the assembly process, the cc2d
!# framework calls user_initCollectForAssembly to stores the current point 
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
!# After the assembly, user_doneCollectForAssembly is called to clean up.
!# Note: Information stored in the quick-access array are of temporary
!# nature and does not have to be cleaned up.
!#
!# </purpose>
!##############################################################################

module user_callback

  use fsystem
  use storage
  use derivatives
  use spatialdiscretisation
  use scalarpde
  use linearformevaluation
  use bilinearformevaluation
  use boundary
  use triangulation
  use collection
  
  use structuresoptc
  use structuresoptflow
  
  use analyticsolution
  use spacetimevectors
  
  use user_flows
  
  implicit none
  
  private
  
  public :: user_initGlobalData
  public :: user_doneGlobalData
  public :: user_initCollectForAssembly
  public :: user_initCollectForVecAssembly
  public :: user_doneCollectForAssembly
  public :: user_coeff_Stokes
  public :: user_coeff_Pressure
  public :: user_coeff_RHS
  public :: user_fct_Target
  public :: user_coeff_Target
  public :: user_coeff_Reference
!  public :: user_coeff_RHSprimal_x
!  public :: user_coeff_RHSprimal_y
!  public :: user_coeff_RHSdual_x
!  public :: user_coeff_RHSdual_y
!  public :: user_coeff_TARGET_x
!  public :: user_coeff_TARGET_y
!  public :: user_ffunction_TargetX
!  public :: user_ffunction_TargetY
  public :: user_getBoundaryValues
  public :: user_getBoundaryValuesFBC_2D

contains

! ***************************************************************************

!<subroutine>

  subroutine user_initGlobalData (rsettings,rrhs,rglobalData)
  
!<description>
  ! Initialises the global data structure rglobalData with information from
  ! the problem structure. This global data contains everything which is
  ! available in user_initCollectForAssembly.
!</description>

!<input>
  ! Global program settings
  type(t_settings_optflow), intent(in), target :: rsettings

  ! Right hand side
  type(t_anSolution), intent(in), target :: rrhs

  ! Global data.
  type(t_globalData), intent(inout) :: rglobalData
!</input>

!</subroutine>

    ! Remember the parameter list and the collection.
    rglobalData%p_rparlist => rsettings%p_rparlist
    rglobalData%p_rtimeCoarse => rsettings%rtimeCoarse
    rglobalData%p_rsettingsOptControl => rsettings%rsettingsOptControl
    rglobalData%p_rrhs => rrhs
    rglobalData%p_rtargetFunction => rsettings%rsettingsOptControl%rtargetFunction
    rglobalData%p_rphysics => rsettings%rphysicsPrimal
    
  end subroutine
  
! ***************************************************************************
  
!<subroutine>

  subroutine user_doneGlobalData (rglobalData)
  
!<description>
  ! Cleans up a global data structure.
!</description>
  
!<input>
  ! Global data.
  type(t_globalData), intent(inout) :: rglobalData
!</input>

!</subroutine>

    nullify(rglobalData%p_rparlist)
    nullify(rglobalData%p_rtimeCoarse)
    nullify(rglobalData%p_rsettingsOptControl)
    nullify(rglobalData%p_rrhs)
    nullify(rglobalData%p_rtargetFunction)

  end subroutine
  
! ***************************************************************************

!<subroutine>

  subroutine user_initCollectForAssembly (rglobalData,dtime,rcollection)
  
!<description>
  ! This subroutine is an auxiliary subroutine called by the framework
  ! and has usually not to be changed by the user.
  !
  ! The subroutine prepares the collection rcollection to be passed to callback
  ! routines for assembling boundary conditions or RHS vectors. It is
  ! called directly prior to the assembly to store problem-specific information
  ! in the quick-access arrays of the collection.
!</description>

!<input>
  ! Current simulation time.
  real(dp), intent(in) :: dtime
!</input>

!<inputoutput>
  ! Global data.
  type(t_globalData), intent(inout) :: rglobalData

  ! Collection structure which is passed to the callback routines.
  type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!</subroutine>

    real(DP) :: dreltime

    ! Save the simulation time as well as the minimum and maximum time 
    ! to the quick-access array of the collection,
    ! so it can be accessed in the callback routines!
    rcollection%IquickAccess(1) = rglobalData%p_rrhs%iid
    rcollection%IquickAccess(2) = rglobalData%p_rsettingsOptControl%rtargetFunction%iid
    rcollection%Dquickaccess(1) = dtime
    rcollection%Dquickaccess(2) = rglobalData%p_rtimeCoarse%dtimeInit
    rcollection%Dquickaccess(3) = rglobalData%p_rtimeCoarse%dtimeMax
    rcollection%Dquickaccess(4) = rglobalData%p_rsettingsOptControl%dalphaC

  end subroutine
  
! ***************************************************************************

!<subroutine>

  subroutine user_initCollectForVecAssembly (rglobalData,iid,icomponent,dtime,rcollection)
  
!<description>
  ! This subroutine is an auxiliary subroutine called by the framework
  ! and has usually not to be changed by the user.
  !
  ! The subroutine prepares the collection rcollection to be passed to callback
  ! routines for assembling boundary conditions or RHS vectors. It is
  ! called directly prior to the assembly to store problem-specific information
  ! in the quick-access arrays of the collection.
  !
  ! The routine is called prior to the assembly of a vector.
  ! ivectype defines the type of vector which is to be assembled
  ! while icomponent defines its component. The routine has to respect the current
  ! physics setting in order to choose the correct function which is to be
  ! assembled!
!</description>

!<input>
  ! Global Id of the flow.
  integer, intent(in) :: iid
  
  ! Component to be assembled.
  integer, intent(in) :: icomponent

  ! Current simulation time.
  real(dp), intent(in) :: dtime
!</input>

!<inputoutput>
  ! Global data.
  type(t_globalData), intent(inout) :: rglobalData

  ! Collection structure which is passed to the callback routines.
  type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!</subroutine>

    real(DP) :: dreltime

    ! Save the simulation time as well as the minimum and maximum time 
    ! to the quick-access array of the collection,
    ! so it can be accessed in the callback routines!
    !
    ! IquickAccess(1) = equation type
    ! IquickAccess(2) = component to assemble
    ! IquickAccess(3) = type of the function
    ! IquickAccess(4) = id of the function
    
    rcollection%IquickAccess(1) = rglobalData%p_rphysics%cequation
    rcollection%IquickAccess(2) = icomponent
    rcollection%IquickAccess(3) = iid
    rcollection%Dquickaccess(1) = dtime
    rcollection%Dquickaccess(2) = rglobalData%p_rtimeCoarse%dtimeInit
    rcollection%Dquickaccess(3) = rglobalData%p_rtimeCoarse%dtimeMax
    rcollection%Dquickaccess(4) = rglobalData%p_rsettingsOptControl%dalphaC

  end subroutine
  
! ***************************************************************************
  
!<subroutine>

  subroutine user_doneCollectForAssembly (rglobalData,rcollection)
  
!<description>
  ! This subroutine is an auxiliary subroutine called by the framework
  ! and has usually not to be changed by the user.
  !
  ! After the assembly process, this subroutine is called to release temporary
  ! information from the collection which was stored there by 
  ! user_initCollectForAssembly.
!</description>
  
!<inputoutput>
  ! Global data.
  type(t_globalData), intent(inout) :: rglobalData

  ! Collection structure which is passed to the callback routines.
  type(t_collection), intent(inout) :: rcollection
!</inputoutput>
  
!</subroutine>

    ! Nothing to be done at the moment.

  end subroutine
  
! ***************************************************************************
  !<subroutine>

  subroutine user_coeff_Stokes (rdiscretisationTrial,rdiscretisationTest,rform, &
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

  subroutine user_coeff_Pressure (rdiscretisationTrial,rdiscretisationTest,rform, &
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

  subroutine user_coeff_RHS (rdiscretisation,rform, &
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
    ! the coefficients in front of the terms of the linear form.
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
  
    real(DP) :: dtime,dalpha
    integer :: ierror,iid,cequation,icomponent
    real(DP), dimension(:,:), allocatable :: DvaluesAct
    
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      cequation = rcollection%Iquickaccess(1)
      icomponent = rcollection%Iquickaccess(2)
      iid = rcollection%Iquickaccess(3)
      dtime = rcollection%Dquickaccess(1)
      dalpha = rcollection%Dquickaccess(4)
    else
      cequation = -1
      iid = 0
      icomponent = 0
      dtime = 0.0_DP
      dalpha = 1.0_DP
    end if
    
    ! Evaluate the RHS if possible.
    allocate(DvaluesAct(npointsPerElement,nelements))
    
    select case (cequation)
    case (0,1)
      ! Stokes, Navier-Stokes, 2D
      
      select case (icomponent)
      case (1)
        ! primal X-velocity
        
        select case (iid)
        case (7)
          DvaluesAct(:,:) = fct_stokesF7_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case default
          ! Analytically given. =0.
          DvaluesAct(:,:) = 0.0_DP
          
        end select

      case (2)
        ! primal Y-velocity

        select case (iid)
        case (7)
          DvaluesAct(:,:) = fct_stokesF7_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case default
          ! Analytically given. =0.
          DvaluesAct(:,:) = 0.0_DP
          
        end select

      case (4)
        ! dual X-velocity
        select case (iid)
        case default
          ! Analytically given. =0.
          DvaluesAct(:,:) = 0.0_DP
          
        end select

      case (5)
        ! dual Y-velocity
        select case (iid)
        case default
          ! Analytically given. =0.
          DvaluesAct(:,:) = 0.0_DP
          
        end select
        
      end select

    end select

    Dcoefficients(1,1:npointsPerElement,1:nelements) = DvaluesAct(:,:)
    deallocate(DvaluesAct)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine user_fct_Target (cderivative, rdiscretisation, &
                                   nelements, npointsPerElement, Dpoints, &
                                   IdofsTest, rdomainIntSubset, &
                                   Dvalues, rcollection)
    
    use fsystem
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
  !<description>
    ! This subroutine is called during the calculation of errors. It computes
    ! an analytically given target function in a couple of points on a couple
    ! of elements. These values are compared to those of a computed FE function
    ! and used to calculate an error.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points.
  !</description>
    
  !<input>
    ! This is a DER_xxxx derivative identifier (from derivative.f90) that
    ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
    ! The result must be written to the Dvalue-array below.
    integer, intent(in) :: cderivative
  
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(inout), optional :: rcollection
    
  !</input>
  
  !<output>
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints, or the appropriate derivative
    ! of the function, respectively, according to cderivative.
    !   DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(out) :: Dvalues
  !</output>
    
  !</subroutine>
  
    real(DP) :: dtime,dalpha
    integer :: ierror,iid,cequation,icomponent
    real(DP), dimension(:,:), allocatable :: DvaluesAct
    
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      cequation = rcollection%Iquickaccess(1)
      icomponent = rcollection%Iquickaccess(2)
      iid = rcollection%Iquickaccess(3)
      dtime = rcollection%Dquickaccess(1)
      dalpha = rcollection%Dquickaccess(4)
    else
      cequation = -1
      iid = 0
      icomponent = 0
      dtime = 0.0_DP
      dalpha = 1.0_DP
    end if
    
    select case (cequation)
    case (0,1)
      ! Stokes, Navier-Stokes, 2D
      
      select case (icomponent)
      case (1)
        ! target X-velocity
        
        select case (iid)
        case (1)
          Dvalues(:,:) = fct_stokesZ1_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (2)
          Dvalues(:,:) = fct_stokesZ2_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (4)
          Dvalues(:,:) = fct_stokesZ4_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (7)
          Dvalues(:,:) = fct_stokesZ7_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)
          
        case default
          ! Analytically given. =0.
          Dvalues(:,:) = 0.0_DP
          
        end select

      case (2)
        ! target Y-velocity

        select case (iid)
        case (1)
          Dvalues(:,:) = fct_stokesZ1_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (2)
          Dvalues(:,:) = fct_stokesZ2_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (4)
          Dvalues(:,:) = fct_stokesZ4_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (7)
          Dvalues(:,:) = fct_stokesZ7_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case default
          ! Analytically given. =0.
          Dvalues(:,:) = 0.0_DP
          
        end select

      end select

    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine user_coeff_Target (rdiscretisation,rform, &
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
    ! target X-velocity, i.e. the desired function z, which enters the dual
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
  
    real(DP) :: dtime,dalpha
    integer :: ierror,iid,cequation,icomponent
    real(DP), dimension(:,:), allocatable :: DvaluesAct
    
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      cequation = rcollection%Iquickaccess(1)
      icomponent = rcollection%Iquickaccess(2)
      iid = rcollection%Iquickaccess(3)
      dtime = rcollection%Dquickaccess(1)
      dalpha = rcollection%Dquickaccess(4)
    else
      cequation = -1
      iid = 0
      icomponent = 0
      dtime = 0.0_DP
      dalpha = 1.0_DP
    end if
    
    ! Evaluate the RHS if possible.
    allocate(DvaluesAct(npointsPerElement,nelements))
    
    select case (cequation)
    case (0,1)
      ! Stokes, Navier-Stokes, 2D
      
      select case (icomponent)
      case (1)
        ! target X-velocity
        
        select case (iid)
        case (1)
          DvaluesAct(:,:) = fct_stokesZ1_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (2)
          DvaluesAct(:,:) = fct_stokesZ2_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (4)
          DvaluesAct(:,:) = fct_stokesZ4_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)
          
        case (7)
          DvaluesAct(:,:) = fct_stokesZ7_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)
          
        case default
          ! Analytically given. =0.
          DvaluesAct(:,:) = 0.0_DP
          
        end select

      case (2)
        ! target Y-velocity

        select case (iid)
        case (1)
          DvaluesAct(:,:) = fct_stokesZ1_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (2)
          DvaluesAct(:,:) = fct_stokesZ2_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (4)
          DvaluesAct(:,:) = fct_stokesZ4_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (7)
          DvaluesAct(:,:) = fct_stokesZ7_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case default
          ! Analytically given. =0.
          DvaluesAct(:,:) = 0.0_DP
          
        end select

      end select

    end select

    Dcoefficients(1,1:npointsPerElement,1:nelements) = DvaluesAct(:,:)
    deallocate(DvaluesAct)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine user_coeff_Reference (cderivative, rdiscretisation, &
                                   nelements, npointsPerElement, Dpoints, &
                                   IdofsTest, rdomainIntSubset, &
                                   Dvalues, rcollection)
    
    use fsystem
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
  !<description>
    ! This subroutine is called during the calculation of errors. It has to compute
    ! the values of an analytically given function in a set of points on a couple
    ! of elements. These values are compared to those of a computed FE function
    ! and used to calculate an error.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points.
  !</description>
    
  !<input>
    ! This is a DER_xxxx derivative identifier (from derivative.f90) that
    ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
    ! The result must be written to the Dvalue-array below.
    integer, intent(in) :: cderivative
  
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(inout), optional :: rcollection
    
  !</input>
  
  !<output>
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints, or the appropriate derivative
    ! of the function, respectively, according to cderivative.
    !   DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(out) :: Dvalues
  !</output>
    
  !</subroutine>
  
    real(DP) :: dtime,dalpha
    integer :: ierror,iid,cequation,icomponent
    real(DP), dimension(:,:), allocatable :: DvaluesAct
    
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      cequation = rcollection%Iquickaccess(1)
      icomponent = rcollection%Iquickaccess(2)
      iid = rcollection%Iquickaccess(3)
      dtime = rcollection%Dquickaccess(1)
      dalpha = rcollection%Dquickaccess(4)
    else
      cequation = -1
      iid = 0
      icomponent = 0
      dtime = 0.0_DP
      dalpha = 1.0_DP
    end if
    
    select case (cequation)
    case (0,1)
      ! Stokes, Navier-Stokes, 2D
      
      select case (icomponent)
      case (1)
        ! primal X-velocity
        
        select case (iid)
        case (1)
          Dvalues(:,:) = fct_stokesY1_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (2)
          Dvalues(:,:) = fct_stokesY2_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (4)
          Dvalues(:,:) = fct_stokesY4_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (7)
          Dvalues(:,:) = fct_stokesY7_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case default
          ! Analytically given. =0.
          Dvalues(:,:) = 0.0_DP
          
        end select

      case (2)
        ! primal Y-velocity

        select case (iid)
        case (1)
          Dvalues(:,:) = fct_stokesY1_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (2)
          Dvalues(:,:) = fct_stokesY2_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (4)
          Dvalues(:,:) = fct_stokesY4_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (7)
          Dvalues(:,:) = fct_stokesY7_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case default
          ! Analytically given. =0.
          Dvalues(:,:) = 0.0_DP
          
        end select

      case (3)
        ! primal pressure

        select case (iid)
        case (1)
          Dvalues(:,:) = fct_stokesP1(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (2)
          Dvalues(:,:) = fct_stokesP2(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (4)
          Dvalues(:,:) = fct_stokesP4(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (7)
          Dvalues(:,:) = fct_stokesP7(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case default
          ! Analytically given. =0.
          Dvalues(:,:) = 0.0_DP
          
        end select

      case (4)
        ! dual X-velocity
        select case (iid)
        case (1)
          Dvalues(:,:) = fct_stokesLambda1_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (2)
          Dvalues(:,:) = fct_stokesLambda2_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (4)
          Dvalues(:,:) = fct_stokesLambda4_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (7)
          Dvalues(:,:) = fct_stokesLambda7_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case default
          ! Analytically given. =0.
          Dvalues(:,:) = 0.0_DP
          
        end select

      case (5)
        ! dual Y-velocity
        select case (iid)
        case (1)
          Dvalues(:,:) = fct_stokesLambda1_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (2)
          Dvalues(:,:) = fct_stokesLambda2_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (4)
          Dvalues(:,:) = fct_stokesLambda4_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (7)
          Dvalues(:,:) = fct_stokesLambda7_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case default
          ! Analytically given. =0.
          Dvalues(:,:) = 0.0_DP
          
        end select

      case (6)
        ! dual pressure
        select case (iid)
        case (1)
          Dvalues(:,:) = fct_stokesXi1(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (2)
          Dvalues(:,:) = fct_stokesXi2(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (4)
          Dvalues(:,:) = fct_stokesXi4(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case (7)
          Dvalues(:,:) = fct_stokesXi7(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

        case default
          ! Analytically given. =0.
          Dvalues(:,:) = 0.0_DP
          
        end select
        
      end select

    end select

  end subroutine

!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine user_coeff_RHSprimal_x (rdiscretisation,rform, &
!                  nelements,npointsPerElement,Dpoints, &
!                  IdofsTest,rdomainIntSubset, &
!                  Dcoefficients,rcollection)
!    
!    use basicgeometry
!    use triangulation
!    use collection
!    use scalarpde
!    use domainintegration
!    
!  !<description>
!    ! This subroutine is called during the vector assembly. It has to compute
!    ! the coefficients in front of the terms of the linear form of the
!    ! X-velocity part of the right hand side vector.
!    !
!    ! The routine accepts a set of elements and a set of points on these
!    ! elements (cubature points) in real coordinates.
!    ! According to the terms in the linear form, the routine has to compute
!    ! simultaneously for all these points and all the terms in the linear form
!    ! the corresponding coefficients in front of the terms.
!  !</description>
!    
!  !<input>
!    ! The discretisation structure that defines the basic shape of the
!    ! triangulation with references to the underlying triangulation,
!    ! analytic boundary boundary description etc.
!    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
!    
!    ! The linear form which is currently to be evaluated:
!    type(t_linearForm), intent(IN)                              :: rform
!    
!    ! Number of elements, where the coefficients must be computed.
!    integer, intent(IN)                                         :: nelements
!    
!    ! Number of points per element, where the coefficients must be computed
!    integer, intent(IN)                                         :: npointsPerElement
!    
!    ! This is an array of all points on all the elements where coefficients
!    ! are needed.
!    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
!    ! DIMENSION(dimension,npointsPerElement,nelements)
!    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints
!
!    ! An array accepting the DOF's on all elements trial in the trial space.
!    ! DIMENSION(\#local DOF's in test space,nelements)
!    integer, dimension(:,:), intent(IN) :: IdofsTest
!
!    ! This is a t_domainIntSubset structure specifying more detailed information
!    ! about the element set that is currently being integrated.
!    ! It's usually used in more complex situations (e.g. nonlinear matrices).
!    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset
!
!    ! Optional: A collection structure to provide additional 
!    ! information to the coefficient routine. 
!    type(t_collection), intent(INOUT), optional      :: rcollection
!    
!  !</input>
!  
!  !<output>
!    ! A list of all coefficients in front of all terms in the linear form -
!    ! for all given points on all given elements.
!    !   DIMENSION(itermCount,npointsPerElement,nelements)
!    ! with itermCount the number of terms in the linear form.
!    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
!  !</output>
!    
!  !</subroutine>
!  
!    real(DP) :: dtime
!    integer :: ierror,iid
!    real(DP), dimension(:,:), allocatable :: DvaluesAct
!    
!    ! In a nonstationary simulation, one can get the simulation time
!    ! with the quick-access array of the collection.
!    if (present(rcollection)) then
!      dtime = rcollection%Dquickaccess(1)
!      iid = rcollection%Iquickaccess(1)
!    else
!      dtime = 0.0_DP
!      iid = 0
!    end if
!
!    ! Evaluate the RHS if possible.
!    allocate(DvaluesAct(npointsPerElement,nelements))
!
!    select case (iid)
!    case (0)
!      ! Analytically given. =0.
!      Dcoefficients(:,:,:) = 0.0_DP
!      
!      ! Other cases may be defined here.
!      !Dcoefficients(1,:,:) = Dpoints(1,:,:)
!      !Dcoefficients(1,:,:) = -18.0*sin(3.0*SYS_PI*Dpoints(1,:,:))*SYS_PI**2 &
!      !                     *sin(3.0*SYS_PI*Dpoints(2,:,:)) &
!      !                     + .5*SYS_PI*cos(.5*SYS_PI*(Dpoints(1,:,:)-Dpoints(2,:,:)))
!      !Dcoefficients (1,:,:) = -(1./10.0_DP)*(-Dpoints(1,:,:))
!      
!      ! Without coupling:
!      !Dcoefficients (1,:,:) = (1/5.0_DP - dtime/50.0_DP)*(Dpoints(1,:,:))
!      !Dcoefficients (1,:,:) = Dpoints(1,:,:)
!    end select
!    
!    deallocate(DvaluesAct)
!
!  end subroutine
!
!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine user_coeff_RHSprimal_y (rdiscretisation,rform, &
!                  nelements,npointsPerElement,Dpoints, &
!                  IdofsTest,rdomainIntSubset, &
!                  Dcoefficients,rcollection)
!    
!    use basicgeometry
!    use triangulation
!    use collection
!    use scalarpde
!    use domainintegration
!    
!  !<description>
!    ! This subroutine is called during the vector assembly. It has to compute
!    ! the coefficients in front of the terms of the linear form of the
!    ! Y-velocity part of the right hand side vector.
!    !
!    ! The routine accepts a set of elements and a set of points on these
!    ! elements (cubature points) in real coordinates.
!    ! According to the terms in the linear form, the routine has to compute
!    ! simultaneously for all these points and all the terms in the linear form
!    ! the corresponding coefficients in front of the terms.
!  !</description>
!    
!  !<input>
!    ! The discretisation structure that defines the basic shape of the
!    ! triangulation with references to the underlying triangulation,
!    ! analytic boundary boundary description etc.
!    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
!    
!    ! The linear form which is currently to be evaluated:
!    type(t_linearForm), intent(IN)                              :: rform
!    
!    ! Number of elements, where the coefficients must be computed.
!    integer, intent(IN)                                         :: nelements
!    
!    ! Number of points per element, where the coefficients must be computed
!    integer, intent(IN)                                         :: npointsPerElement
!    
!    ! This is an array of all points on all the elements where coefficients
!    ! are needed.
!    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
!    ! DIMENSION(dimension,npointsPerElement,nelements)
!    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints
!
!    ! An array accepting the DOF's on all elements trial in the trial space.
!    ! DIMENSION(\#local DOF's in test space,nelements)
!    integer, dimension(:,:), intent(IN) :: IdofsTest
!
!    ! This is a t_domainIntSubset structure specifying more detailed information
!    ! about the element set that is currently being integrated.
!    ! It's usually used in more complex situations (e.g. nonlinear matrices).
!    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset
!
!    ! Optional: A collection structure to provide additional 
!    ! information to the coefficient routine. 
!    type(t_collection), intent(INOUT), optional      :: rcollection
!    
!  !</input>
!  
!  !<output>
!    ! A list of all coefficients in front of all terms in the linear form -
!    ! for all given points on all given elements.
!    !   DIMENSION(itermCount,npointsPerElement,nelements)
!    ! with itermCount the number of terms in the linear form.
!    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
!  !</output>
!    
!  !</subroutine>
!
!    real(DP) :: dtime
!    integer :: ierror,iid
!    real(DP), dimension(:,:), allocatable :: DvaluesAct
!    
!    ! In a nonstationary simulation, one can get the simulation time
!    ! with the quick-access array of the collection.
!    if (present(rcollection)) then
!      dtime = rcollection%Dquickaccess(1)
!      iid = rcollection%Iquickaccess(1)
!    else
!      dtime = 0.0_DP
!      iid = 0
!    end if
!
!    ! Evaluate the RHS if possible.
!    allocate(DvaluesAct(npointsPerElement,nelements))
!    
!    ! Evaluate using iid.
!    select case (iid)
!    case (0)
!      ! Analytically given. =0.
!      Dcoefficients(:,:,:) = 0.0_DP
!      
!      ! Other cases may be defined here.
!      !Dcoefficients(1,:,:) = -Dpoints(2,:,:)
!      !Dcoefficients(1,:,:) = -18.0*cos(3.0*SYS_PI*Dpoints(1,:,:))*SYS_PI**2 &
!      !                     *cos(3.0*SYS_PI*Dpoints(2,:,:)) &
!      !                     - .5*SYS_PI*cos(.5*SYS_PI*(Dpoints(1,:,:)-Dpoints(2,:,:)))
!      !Dcoefficients (1,:,:) = -(1./10.0_DP)*(Dpoints(2,:,:))
!      
!      ! Without coupling:
!      !Dcoefficients (1,:,:) = (1/5.0_DP - dtime/50.0_DP)*(-Dpoints(2,:,:))
!      !Dcoefficients (1,:,:) = -Dpoints(2,:,:)
!    end select
!    
!    deallocate(DvaluesAct)
!
!  end subroutine
!
!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine user_coeff_RHSdual_x (rdiscretisation,rform, &
!                  nelements,npointsPerElement,Dpoints, &
!                  IdofsTest,rdomainIntSubset, &
!                  Dcoefficients,rcollection)
!    
!    use basicgeometry
!    use triangulation
!    use collection
!    use scalarpde
!    use domainintegration
!    
!  !<description>
!    ! This subroutine is called during the vector assembly. It has to compute
!    ! the coefficients in front of the terms of the linear form of the
!    ! X-velocity part of the right hand side vector.
!    !
!    ! The routine accepts a set of elements and a set of points on these
!    ! elements (cubature points) in real coordinates.
!    ! According to the terms in the linear form, the routine has to compute
!    ! simultaneously for all these points and all the terms in the linear form
!    ! the corresponding coefficients in front of the terms.
!  !</description>
!    
!  !<input>
!    ! The discretisation structure that defines the basic shape of the
!    ! triangulation with references to the underlying triangulation,
!    ! analytic boundary boundary description etc.
!    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
!    
!    ! The linear form which is currently to be evaluated:
!    type(t_linearForm), intent(IN)                              :: rform
!    
!    ! Number of elements, where the coefficients must be computed.
!    integer, intent(IN)                                         :: nelements
!    
!    ! Number of points per element, where the coefficients must be computed
!    integer, intent(IN)                                         :: npointsPerElement
!    
!    ! This is an array of all points on all the elements where coefficients
!    ! are needed.
!    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
!    ! DIMENSION(dimension,npointsPerElement,nelements)
!    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints
!
!    ! An array accepting the DOF's on all elements trial in the trial space.
!    ! DIMENSION(\#local DOF's in test space,nelements)
!    integer, dimension(:,:), intent(IN) :: IdofsTest
!
!    ! This is a t_domainIntSubset structure specifying more detailed information
!    ! about the element set that is currently being integrated.
!    ! It's usually used in more complex situations (e.g. nonlinear matrices).
!    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset
!
!    ! Optional: A collection structure to provide additional 
!    ! information to the coefficient routine. 
!    type(t_collection), intent(INOUT), optional      :: rcollection
!    
!  !</input>
!  
!  !<output>
!    ! A list of all coefficients in front of all terms in the linear form -
!    ! for all given points on all given elements.
!    !   DIMENSION(itermCount,npointsPerElement,nelements)
!    ! with itermCount the number of terms in the linear form.
!    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
!  !</output>
!    
!  !</subroutine>
!  
!    ! Return the coefficients for the dual RHS.
!    ! These are normally =0.
!    !
!    ! The coefficients of the target flow will be added somewhere else!
!    Dcoefficients(:,:,:) = 0.0_DP
!      
!  end subroutine
!
!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine user_coeff_RHSdual_y (rdiscretisation,rform, &
!                  nelements,npointsPerElement,Dpoints, &
!                  IdofsTest,rdomainIntSubset, &
!                  Dcoefficients,rcollection)
!    
!    use basicgeometry
!    use triangulation
!    use collection
!    use scalarpde
!    use domainintegration
!    
!  !<description>
!    ! This subroutine is called during the vector assembly. It has to compute
!    ! the coefficients in front of the terms of the linear form of the
!    ! Y-velocity part of the right hand side vector.
!    !
!    ! The routine accepts a set of elements and a set of points on these
!    ! elements (cubature points) in real coordinates.
!    ! According to the terms in the linear form, the routine has to compute
!    ! simultaneously for all these points and all the terms in the linear form
!    ! the corresponding coefficients in front of the terms.
!  !</description>
!    
!  !<input>
!    ! The discretisation structure that defines the basic shape of the
!    ! triangulation with references to the underlying triangulation,
!    ! analytic boundary boundary description etc.
!    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
!    
!    ! The linear form which is currently to be evaluated:
!    type(t_linearForm), intent(IN)                              :: rform
!    
!    ! Number of elements, where the coefficients must be computed.
!    integer, intent(IN)                                         :: nelements
!    
!    ! Number of points per element, where the coefficients must be computed
!    integer, intent(IN)                                         :: npointsPerElement
!    
!    ! This is an array of all points on all the elements where coefficients
!    ! are needed.
!    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
!    ! DIMENSION(dimension,npointsPerElement,nelements)
!    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints
!
!    ! An array accepting the DOF's on all elements trial in the trial space.
!    ! DIMENSION(\#local DOF's in test space,nelements)
!    integer, dimension(:,:), intent(IN) :: IdofsTest
!
!    ! This is a t_domainIntSubset structure specifying more detailed information
!    ! about the element set that is currently being integrated.
!    ! It's usually used in more complex situations (e.g. nonlinear matrices).
!    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset
!
!    ! Optional: A collection structure to provide additional 
!    ! information to the coefficient routine. 
!    type(t_collection), intent(INOUT), optional      :: rcollection
!    
!  !</input>
!  
!  !<output>
!    ! A list of all coefficients in front of all terms in the linear form -
!    ! for all given points on all given elements.
!    !   DIMENSION(itermCount,npointsPerElement,nelements)
!    ! with itermCount the number of terms in the linear form.
!    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
!  !</output>
!    
!  !</subroutine>
!
!    ! Return the coefficients for the dual RHS.
!    ! These are normally =0.
!    !
!    ! The coefficients of the target flow will be added somewhere else!
!    Dcoefficients(:,:,:) = 0.0_DP
!
!  end subroutine
!
!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine user_coeff_TARGET_x (rdiscretisation,rform, &
!                  nelements,npointsPerElement,Dpoints, &
!                  IdofsTest,rdomainIntSubset, &
!                  Dcoefficients,rcollection)
!    
!    use basicgeometry
!    use triangulation
!    use collection
!    use scalarpde
!    use domainintegration
!    
!  !<description>
!    ! This subroutine is called during the vector assembly. It has to compute
!    ! the coefficients in front of the terms of the linear form of the
!    ! target X-velocity, i.e. the desired flow field z, which enters the dual
!    ! equation.
!    !
!    ! The routine accepts a set of elements and a set of points on these
!    ! elements (cubature points) in real coordinates.
!    ! According to the terms in the linear form, the routine has to compute
!    ! simultaneously for all these points and all the terms in the linear form
!    ! the corresponding coefficients in front of the terms.
!  !</description>
!    
!  !<input>
!    ! The discretisation structure that defines the basic shape of the
!    ! triangulation with references to the underlying triangulation,
!    ! analytic boundary boundary description etc.
!    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
!    
!    ! The linear form which is currently to be evaluated:
!    type(t_linearForm), intent(IN)                              :: rform
!    
!    ! Number of elements, where the coefficients must be computed.
!    integer, intent(IN)                                         :: nelements
!    
!    ! Number of points per element, where the coefficients must be computed
!    integer, intent(IN)                                         :: npointsPerElement
!    
!    ! This is an array of all points on all the elements where coefficients
!    ! are needed.
!    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
!    ! DIMENSION(dimension,npointsPerElement,nelements)
!    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints
!
!    ! An array accepting the DOF's on all elements trial in the trial space.
!    ! DIMENSION(\#local DOF's in test space,nelements)
!    integer, dimension(:,:), intent(IN) :: IdofsTest
!
!    ! This is a t_domainIntSubset structure specifying more detailed information
!    ! about the element set that is currently being integrated.
!    ! It's usually used in more complex situations (e.g. nonlinear matrices).
!    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset
!
!    ! Optional: A collection structure to provide additional 
!    ! information to the coefficient routine. 
!    type(t_collection), intent(INOUT), optional      :: rcollection
!    
!  !</input>
!  
!  !<output>
!    ! A list of all coefficients in front of all terms in the linear form -
!    ! for all given points on all given elements.
!    !   DIMENSION(itermCount,npointsPerElement,nelements)
!    ! with itermCount the number of terms in the linear form.
!    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
!  !</output>
!    
!  !</subroutine>
!  
!    real(DP) :: dtime
!    real(DP), dimension(:,:), allocatable :: DvaluesAct
!  
!    ! In a nonstationary simulation, one can get the simulation time
!    ! with the quick-access array of the collection.
!    if (present(rcollection)) then
!      dtime = rcollection%Dquickaccess(1)
!    else
!      dtime = 0.0_DP
!    end if
!
!    Dcoefficients(:,:,:) = 0.0_DP
!
!    ! Call user_ffunction_TargetX to calculate the analytic function. Store the results
!    ! in  Dcoefficients(1,:,:).
!    allocate(DvaluesAct(npointsPerElement,nelements))
!    
!    call user_ffunction_TargetX (DER_FUNC,rdiscretisation, &
!        nelements,npointsPerElement,Dpoints, &
!        IdofsTest,rdomainIntSubset,&
!        DvaluesAct,rcollection)
!    Dcoefficients(1,1:npointsPerElement,1:nelements) = DvaluesAct(:,:)               
!               
!    deallocate(DvaluesAct)
!
!    ! Without coupling:
!    !Dcoefficients (1,:,:) = - (1/5.0_DP - (10.-dtime)/50.0_DP)*(Dpoints(1,:,:))
!
!  end subroutine
!
!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine user_coeff_TARGET_y (rdiscretisation,rform, &
!                  nelements,npointsPerElement,Dpoints, &
!                  IdofsTest,rdomainIntSubset,&
!                  Dcoefficients,rcollection)
!    
!    use basicgeometry
!    use triangulation
!    use collection
!    use scalarpde
!    use domainintegration
!    
!  !<description>
!    ! This subroutine is called during the vector assembly. It has to compute
!    ! the coefficients in front of the terms of the linear form of the
!    ! target X-velocity, i.e. the desired flow field z, which enters the dual
!    ! equation.
!    !
!    ! The routine accepts a set of elements and a set of points on these
!    ! elements (cubature points) in real coordinates.
!    ! According to the terms in the linear form, the routine has to compute
!    ! simultaneously for all these points and all the terms in the linear form
!    ! the corresponding coefficients in front of the terms.
!  !</description>
!    
!  !<input>
!    ! The discretisation structure that defines the basic shape of the
!    ! triangulation with references to the underlying triangulation,
!    ! analytic boundary boundary description etc.
!    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
!    
!    ! The linear form which is currently to be evaluated:
!    type(t_linearForm), intent(IN)                              :: rform
!    
!    ! Number of elements, where the coefficients must be computed.
!    integer, intent(IN)                                         :: nelements
!    
!    ! Number of points per element, where the coefficients must be computed
!    integer, intent(IN)                                         :: npointsPerElement
!    
!    ! This is an array of all points on all the elements where coefficients
!    ! are needed.
!    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
!    ! DIMENSION(dimension,npointsPerElement,nelements)
!    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints
!
!    ! An array accepting the DOF's on all elements trial in the trial space.
!    ! DIMENSION(\#local DOF's in test space,nelements)
!    integer, dimension(:,:), intent(IN) :: IdofsTest
!
!    ! This is a t_domainIntSubset structure specifying more detailed information
!    ! about the element set that is currently being integrated.
!    ! It's usually used in more complex situations (e.g. nonlinear matrices).
!    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset
!
!    ! Optional: A collection structure to provide additional 
!    ! information to the coefficient routine. 
!    type(t_collection), intent(INOUT), optional      :: rcollection
!    
!  !</input>
!  
!  !<output>
!    ! A list of all coefficients in front of all terms in the linear form -
!    ! for all given points on all given elements.
!    !   DIMENSION(itermCount,npointsPerElement,nelements)
!    ! with itermCount the number of terms in the linear form.
!    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
!  !</output>
!    
!  !</subroutine>
!  
!    real(DP) :: dtime
!    real(DP), dimension(:,:), allocatable :: DvaluesAct
!  
!    ! In a nonstationary simulation, one can get the simulation time
!    ! with the quick-access array of the collection.
!    if (present(rcollection)) then
!      dtime = rcollection%Dquickaccess(1)
!    else
!      dtime = 0.0_DP
!    end if
!
!    Dcoefficients(:,:,:) = 0.0_DP
!
!    ! Call user_ffunction_TargetX to calculate the analytic function. Store the results
!    ! in  Dcoefficients(1,:,:).
!    allocate(DvaluesAct(npointsPerElement,nelements))
!    
!    call user_ffunction_TargetY (DER_FUNC,rdiscretisation, &
!        nelements,npointsPerElement,Dpoints, &
!        IdofsTest,rdomainIntSubset,&
!        DvaluesAct,rcollection)
!    Dcoefficients(1,1:npointsPerElement,1:nelements) = DvaluesAct(:,:)               
!               
!    deallocate(DvaluesAct)
!
!    ! Without coupling:
!    !Dcoefficients (1,:,:) = - (1/5.0_DP - (10.-dtime)/50.0_DP)*(-Dpoints(2,:,:))
!               
!  end subroutine
!
!  ! ***************************************************************************
!  
!!<subroutine>
!
!  subroutine user_ffunction_TargetX (cderivative,rdiscretisation, &
!                nelements,npointsPerElement,Dpoints, &
!                IdofsTest,rdomainIntSubset,&
!                Dvalues,rcollection)
!  
!  use basicgeometry
!  use triangulation
!  use collection
!  use scalarpde
!  use domainintegration
!  
!!<description>
!  ! This routine calculates the X-contribution of the target function $z$
!  ! in the optimal control problem.
!  !
!  ! The routine accepts a set of elements and a set of points on these
!  ! elements (cubature points) in in real coordinates.
!  ! According to the terms in the linear form, the routine has to compute
!  ! simultaneously for all these points.
!!</description>
!  
!!<input>
!  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
!  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
!  ! The result must be written to the Dvalue-array below.
!  integer, intent(IN)                                         :: cderivative
!
!  ! The discretisation structure that defines the basic shape of the
!  ! triangulation with references to the underlying triangulation,
!  ! analytic boundary boundary description etc.
!  type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
!  
!  ! Number of elements, where the coefficients must be computed.
!  integer, intent(IN)                                         :: nelements
!  
!  ! Number of points per element, where the coefficients must be computed
!  integer, intent(IN)                                         :: npointsPerElement
!  
!  ! This is an array of all points on all the elements where coefficients
!  ! are needed.
!  ! DIMENSION(NDIM2D,npointsPerElement,nelements)
!  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
!  real(DP), dimension(:,:,:), intent(IN)  :: Dpoints
!
!  ! An array accepting the DOF's on all elements trial in the trial space.
!  ! DIMENSION(\#local DOF's in trial space,Number of elements)
!  integer, dimension(:,:), intent(IN) :: IdofsTest
!
!  ! This is a t_domainIntSubset structure specifying more detailed information
!  ! about the element set that is currently being integrated.
!  ! It's usually used in more complex situations (e.g. nonlinear matrices).
!  type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset
!
!  ! Optional: A collection structure to provide additional 
!  ! information to the coefficient routine. 
!  type(t_collection), intent(INOUT), optional      :: rcollection
!  
!!</input>
!
!!<output>
!  ! This array has to receive the values of the (analytical) function
!  ! in all the points specified in Dpoints, or the appropriate derivative
!  ! of the function, respectively, according to cderivative.
!  !   DIMENSION(npointsPerElement,nelements)
!  real(DP), dimension(:,:), intent(out) :: Dvalues
!!</output>
!  
!!</subroutine>
!
!    real(DP) :: dtime,dtimeMax
!    integer :: ieltype,ierror,iid
!    
!    ! DEBUG!!!
!    ! real(DP), dimension(:), pointer :: p_Ddata
!
!    ! In a nonstationary simulation, one can get the simulation time
!    ! with the quick-access array of the collection.
!    if (present(rcollection)) then
!      dtime = rcollection%Dquickaccess(1)
!      dtimeMax = rcollection%Dquickaccess(3)
!      iid = rcollection%Iquickaccess(2)
!    else
!      dtime = 0.0_DP
!      dtimeMax = 0.0_DP
!    end if
!    
!    ! Evaluate using iid.
!    select case (iid)
!    case (0)    
!
!      ! Analytically given target flow
!  
!      !Dvalues(:,:) = Dvalues(:,:)*dtime/dtimeMax
!      !Dvalues(:,:) = Dvalues(:,:)*dtime
!      !Dvalues(:,:) = (-(dtime**2)/100._DP + dtime/5._DP) * Dpoints(1,:,:)
!      !Dvalues(:,:) = ((10._DP-dtime)/50._DP - 1._DP/5._DP) * Dpoints(1,:,:)
!      Dvalues(:,:) = & ! 1._DP/50._DP * Dpoints(1,:,:) + &
!                    (-(dtime**2)/100._DP + dtime/5._DP) * Dpoints(1,:,:)
!      !Dvalues(:,:) = ( ((10._DP-dtime)/50._DP - 1._DP/5._DP) + &
!      !                 (-(dtime**2)/100._DP + dtime/5._DP)) * Dpoints(1,:,:)
!      !Dvalues(:,:) = 0.0_DP
!      !IF (dtime .gt. 10._DP) THEN
!      !  Dvalues(:,:) = (-(10._DP**2)/100._DP + 10._DP/5._DP) * Dpoints(1,:,:)
!      !END IF
!      Dvalues(:,:) = Dpoints(1,:,:)
!      
!      !Dvalues (:,:) = 4.0_DP * Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) * dtime
!      !Dvalues (:,:) = 2.0_DP/3.0_DP
!      Dvalues (:,:) = 0.3_DP * 4.0_DP * Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:))
!      
!    case default
!      Dvalues(:,:) = 0.0_DP
!    end select
!  
!  end subroutine
!
!  ! ***************************************************************************
!  
!!<subroutine>
!
!  subroutine user_ffunction_TargetY (cderivative,rdiscretisation, &
!                nelements,npointsPerElement,Dpoints, &
!                IdofsTest,rdomainIntSubset,&
!                Dvalues,rcollection)
!  
!  use basicgeometry
!  use triangulation
!  use collection
!  use scalarpde
!  use domainintegration
!  
!!<description>
!  ! This routine calculates the Y-contribution of the target function $z$
!  ! in the optimal control problem.
!  !
!  ! The routine accepts a set of elements and a set of points on these
!  ! elements (cubature points) in in real coordinates.
!  ! According to the terms in the linear form, the routine has to compute
!  ! simultaneously for all these points.
!!</description>
!  
!!<input>
!  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
!  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
!  ! The result must be written to the Dvalue-array below.
!  integer, intent(IN)                                         :: cderivative
!
!  ! The discretisation structure that defines the basic shape of the
!  ! triangulation with references to the underlying triangulation,
!  ! analytic boundary boundary description etc.
!  type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
!  
!  ! Number of elements, where the coefficients must be computed.
!  integer, intent(IN)                                         :: nelements
!  
!  ! Number of points per element, where the coefficients must be computed
!  integer, intent(IN)                                         :: npointsPerElement
!  
!  ! This is an array of all points on all the elements where coefficients
!  ! are needed.
!  ! DIMENSION(NDIM2D,npointsPerElement,nelements)
!  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
!  real(DP), dimension(:,:,:), intent(IN)  :: Dpoints
!
!  ! An array accepting the DOF's on all elements trial in the trial space.
!  ! DIMENSION(\#local DOF's in trial space,Number of elements)
!  integer, dimension(:,:), intent(IN) :: IdofsTest
!
!  ! This is a t_domainIntSubset structure specifying more detailed information
!  ! about the element set that is currently being integrated.
!  ! It's usually used in more complex situations (e.g. nonlinear matrices).
!  type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset
!
!  ! Optional: A collection structure to provide additional 
!  ! information to the coefficient routine. 
!  type(t_collection), intent(INOUT), optional      :: rcollection
!  
!!</input>
!
!!<output>
!  ! This array has to receive the values of the (analytical) function
!  ! in all the points specified in Dpoints, or the appropriate derivative
!  ! of the function, respectively, according to cderivative.
!  !   DIMENSION(npointsPerElement,nelements)
!  real(DP), dimension(:,:), intent(out) :: Dvalues
!!</output>
!  
!!</subroutine>
!
!    real(DP) :: dtime,dtimeMax
!    integer :: ieltype,ierror,iid
!    
!    ! DEBUG!!!
!    ! real(DP), dimension(:), pointer :: p_Ddata
!
!    ! In a nonstationary simulation, one can get the simulation time
!    ! with the quick-access array of the collection.
!    if (present(rcollection)) then
!      dtime = rcollection%Dquickaccess(1)
!      dtimeMax = rcollection%Dquickaccess(3)
!      iid = rcollection%Iquickaccess(2)
!    else
!      dtime = 0.0_DP
!      dtimeMax = 0.0_DP
!    end if
!    
!    ! Evaluate using iid.
!    select case (iid)
!    case (0)    
!      ! Analytically given target flow
!  
!      !Dvalues(:,:) = Dvalues(:,:)*dtime/dtimeMax
!      !Dvalues(:,:) = Dvalues(:,:)*dtime
!      !Dvalues(:,:) = (-(dtime**2)/100._DP + dtime/5._DP) * (-Dpoints(2,:,:))
!      !Dvalues(:,:) = ((10._DP-dtime)/50._DP - 1._DP/5._DP) * (-Dpoints(2,:,:))
!      Dvalues(:,:) = & !1._DP/50._DP * (-Dpoints(2,:,:)) + &
!                    (-(dtime**2)/100._DP + dtime/5._DP) * (-Dpoints(2,:,:))
!      !Dvalues(:,:) = ( ((10._DP-dtime)/50._DP - 1._DP/5._DP) + &
!      !                 (-(dtime**2)/100._DP + dtime/5._DP)) * (-Dpoints(2,:,:))
!      !Dvalues(:,:) = 0.0_DP
!      !IF (dtime .gt. 10._DP) THEN
!      !  Dvalues(:,:) = (-(10._DP**2)/100._DP + 10._DP/5._DP) * (-Dpoints(2,:,:))
!      !END IF
!      Dvalues(:,:) = (-Dpoints(2,:,:))
!      
!      Dvalues(:,:) = 0.0_DP
!      
!    case default
!      Dvalues(:,:) = 0.0_DP
!    end select
!
!  end subroutine

  ! ***************************************************************************
  ! Values on the real boundary.

!<subroutine>

  subroutine user_getBoundaryValues (sexpressionName,icomponent,rdiscretisation,&
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
  ! If this is ="", the analytic function is defined by the tags
  ! in the collection.
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

    real(DP) :: dtime,dalpha
    integer :: ierror,iid,cequation
    real(DP), dimension(:,:), allocatable :: DvaluesAct
    real(DP) :: dx,dy
    
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      cequation = rcollection%Iquickaccess(1)
      iid = rcollection%Iquickaccess(3)
      dtime = rcollection%Dquickaccess(1)
      dalpha = rcollection%Dquickaccess(4)
    else
      cequation = -1
      iid = 0
      dtime = 0.0_DP
      dalpha = 1.0_DP
    end if

    if (sexpressionName .ne. "") then
    
      ! Value defined by expression name.
      ! Currently only =0.
      dvalue = 0.0_DP
      
    else
      
      ! Value defined by reference function, identified by iid.
      call boundary_getCoords(rdiscretisation%p_rboundary, &
          rboundaryRegion%iboundCompIdx, dwhere, dx, dy)

      select case (cequation)
      case (0,1)
        select case (iid)
        case (1)
          select case (icomponent)
          case (1)
            dvalue = fct_stokesY1_x (dx,dy,dtime,dalpha)
          case (2)
            dvalue = fct_stokesY1_y (dx,dy,dtime,dalpha)
          case (4)
            dvalue = fct_stokesLambda1_x (dx,dy,dtime,dalpha)
          case (5)
            dvalue = fct_stokesLambda1_y (dx,dy,dtime,dalpha)
          end select

        case (2)
          select case (icomponent)
          case (1)
            dvalue = fct_stokesY2_x (dx,dy,dtime,dalpha)
          case (2)
            dvalue = fct_stokesY2_y (dx,dy,dtime,dalpha)
          case (4)
            dvalue = fct_stokesLambda2_x (dx,dy,dtime,dalpha)
          case (5)
            dvalue = fct_stokesLambda2_y (dx,dy,dtime,dalpha)
          end select

        case (4)
          select case (icomponent)
          case (1)
            dvalue = fct_stokesY4_x (dx,dy,dtime,dalpha)
          case (2)
            dvalue = fct_stokesY4_y (dx,dy,dtime,dalpha)
          case (4)
            dvalue = fct_stokesLambda4_x (dx,dy,dtime,dalpha)
          case (5)
            dvalue = fct_stokesLambda4_y (dx,dy,dtime,dalpha)
          end select

        case (7)
          select case (icomponent)
          case (1)
            dvalue = fct_stokesY7_x (dx,dy,dtime,dalpha)
          case (2)
            dvalue = fct_stokesY7_y (dx,dy,dtime,dalpha)
          case (4)
            dvalue = fct_stokesLambda7_x (dx,dy,dtime,dalpha)
          case (5)
            dvalue = fct_stokesLambda7_y (dx,dy,dtime,dalpha)
          end select

        end select

      end select

    end if

  end subroutine

  ! ***************************************************************************
  ! Values in a fictitious boundary component:

  subroutine user_getBoundaryValuesFBC_2D (Icomponents,rdiscretisation,&
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
    integer, intent(IN) :: iwhere
    
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
      
      case default
        print *,'user_getBoundaryValuesFBC: Insupported coordinate type to return.'
        stop
      end select
      
    end subroutine

  end subroutine

end module
