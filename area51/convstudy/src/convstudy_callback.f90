!##############################################################################
!# ****************************************************************************
!# <name> convstudy_callback </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains callback functions which are used by the main program.
!#
!# </purpose>
!##############################################################################

module convstudy_callback

  implicit none

  public :: convst_refFunction
  private

contains

  !*****************************************************************************

!<subroutine>

  subroutine convst_refFunction(cderivative, rdiscretisation, nelements,&
      npointsPerElement, Dpoints, IdofsTest, rdomainIntSubset, Dvalues, rcollection)

!<description>
    ! This subroutine is called during the calculation of errors. It has to compute
    ! the (analytical) values of a function in a couple of points on a couple
    ! of elements. These values are compared to those of a computed FE function
    ! and used to calculate an error.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points.
!</description>
   
    use collection
    use domainintegration
    use fsystem
    use spatialdiscretisation
 
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

    call convst_evaluateFE(rdiscretisation%ndimension, npointsPerElement*nelements,&
        cderivative, Dvalues, rcollection%p_rvectorQuickAccess1%RvectorBlock(1), Dpoints)
    
  end subroutine convst_refFunction

 !*****************************************************************************

!<subroutine>

  subroutine convst_evaluateFE(ndim,npoints, iderType, Dvalues, rvectorScalar,&
      Dpoints, Ielements, IelementsHint, cnonmeshPoints)

!<description>
    ! This is a wrapper routine which allows to call subroutine
    ! fevl_evaluate1 for (1..npointsPerElement,1..nelements) data
!</description>

    use feevaluation
    use fsystem
    use linearsystemscalar

!<input>
    ! Number of spatial dimension
    integer, intent(in) :: ndim

    ! Number of points
    integer, intent(in) :: npoints

    ! Type of function value to evaluate. One of the DER_xxxx constants,
    ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
    integer, intent(in) :: iderType
    
    ! The scalar solution vector that is to be evaluated.
    type(t_vectorScalar), intent(in) :: rvectorScalar
    
    ! A list of points where to evaluate.
    ! DIMENSION(1..ndim,1..npoints)
    real(DP), dimension(ndim,npoints), intent(in) :: Dpoints

    ! OPTIONAL: A list of elements containing the points Dpoints.
    ! If this is not specified, the element numbers containing the points
    ! are determined automatically.
    integer, dimension(npoints), intent(in), optional :: Ielements

    ! OPTIONAL: A list of elements that are near the points in Dpoints.
    ! This gives only a hint where to start searching for the actual elements
    ! containing the points. This is ignored if Ielements is specified!
    integer, dimension(npoints), intent(in), optional :: IelementsHint
  
    ! OPTIONAL: A FEVL_NONMESHPTS_xxxx constant that defines what happens
    ! if a point is located outside of the domain. May happen e.g. in
    ! nonconvex domains. FEVL_NONMESHPTS_NONE is the default
    ! parameter if cnonmeshPoints is not specified.
    integer, intent(in), optional :: cnonmeshPoints
!</input>

!<output>
    ! Values of the FE function at the points specified by Dpoints.
    real(DP), dimension(npoints), intent(out) :: Dvalues
!</output>

!</subroutine>

    ! Call FE-evaluation routines
    call fevl_evaluate(iderType, Dvalues, rvectorScalar, Dpoints,&
        Ielements, IelementsHint, cnonmeshPoints)

  end subroutine convst_evaluateFE
end module convstudy_callback
