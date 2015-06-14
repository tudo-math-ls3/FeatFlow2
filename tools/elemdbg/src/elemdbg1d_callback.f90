!##############################################################################
!# ****************************************************************************
!# <name> elemdbg1d_callback </name>
!# ****************************************************************************
!#
!# <purpose>
!# </purpose>
!##############################################################################

module elemdbg1d_callback

  use fsystem
  use storage
  use boundary
  use cubature
  use matrixfilters
  use vectorfilters
  use spatialdiscretisation
  use discretebc
  use bcassembly
  use derivatives
  use bilinearformevaluation
  use linearformevaluation
  use linearsolver
  
  implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS1D (rdiscretisation,rform, &
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

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
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
    real(DP), dimension(:,:,:), intent(OUT) :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    integer :: itest, isolution
    real(DP) :: dnu, dbeta1
    
    ! Get the info from the collection
    itest = rcollection%IquickAccess(1)
    isolution = rcollection%IquickAccess(2)
    dnu = rcollection%DquickAccess(1)
    dbeta1 = rcollection%DquickAccess(2)

    ! Okay, what type of system do we have here?
    select case(itest)
    case(101)
      ! L2-projection
      select case(isolution)
      case(0)
        Dcoefficients(1,:,:) = sin(SYS_PI*Dpoints(1,:,:))
      case(1)
        Dcoefficients(1,:,:) = Dpoints(1,:,:)
      end select
                   
    case(102)
      ! Poisson-System
      select case(isolution)
      case(0)
        Dcoefficients(1,:,:) = SYS_PI**2 * sin(SYS_PI*Dpoints(1,:,:))
      case(1)
        Dcoefficients(1,:,:) = 0.0_DP
      end select
    
    case(103)
      ! Convection-Diffusion-System
      select case(isolution)
      case(0)
        Dcoefficients(1,:,:) = dnu * SYS_PI**2 * sin(SYS_PI*Dpoints(1,:,:)) &
                             + dbeta1 * SYS_PI * cos(SYS_PI*Dpoints(1,:,:))
      case(1)
        Dcoefficients(1,:,:) = dbeta1
      case(2)
        Dcoefficients(1,:,:) = 0.0_DP
      end select
          
    end select
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceFunction1D (icomponent,cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                rdomainIntSubset,Dvalues,rcollection)
  
  use basicgeometry
  use triangulation
  use collection
  use scalarpde
  use domainintegration
  
!<description>
  ! This subroutine is called during the calculation of errors. It has to compute
  ! the (analytical) values of a function in a couple of points on a couple
  ! of elements. These values are compared to those of a computed FE function
  ! and used to calculate an error.
  !
  ! The routine accepts a set of elements and a set of points on these
  ! elements (cubature points) in in real coordinates.
  ! According to the terms in the linear form, the routine has to compute
  ! simultaneously for all these points.
!</description>
  
!<input>
  ! Specifies which component of the vector field is to be evaluated
  integer, intent(IN) :: icomponent

  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
  ! The result must be written to the Dvalue-array below.
  integer, intent(IN) :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(IN) :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(IN) :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(IN) :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  ! DIMENSION(dimension,npointsPerElement,nelements)
  real(DP), dimension(:,:,:), intent(IN) :: Dpoints

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It's usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(IN) :: rdomainIntSubset

  ! Optional: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(INOUT), optional :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(OUT) :: Dvalues
!</output>
  
!</subroutine>

  real(DP) :: deps, deps2
  
    Dvalues = 0.0_DP

    select case(rcollection%IquickAccess(2))
    case(0)
      select case (cderivative)
      case (DER_FUNC1D)
        Dvalues(:,:) = sin(SYS_PI*Dpoints(1,:,:))
      case (DER_DERIV1D_X)
        Dvalues(:,:) = SYS_PI * cos(SYS_PI*Dpoints(1,:,:))
      end select
    
    case(1)
      select case (cderivative)
      case (DER_FUNC1D)
        Dvalues(:,:) = Dpoints(1,:,:)
      case (DER_DERIV1D_X)
        Dvalues(:,:) = 1.0_DP
      end select
    
    case(2)
      ! Get 1/nu and exp(1/nu)
      deps = 1.0_DP / rcollection%DquickAccess(1)
      deps2 = exp(deps)
      select case (cderivative)
      case (DER_FUNC1D)
        Dvalues(:,:) = (deps2 - exp(deps*Dpoints(1,:,:))) / (deps2 - 1.0_DP)
      case (DER_DERIV1D_X)
        Dvalues(:,:) = -deps*exp(deps*Dpoints(1,:,:)) / (deps2 - 1.0_DP)
      end select

    end select
    
  end subroutine

end module
