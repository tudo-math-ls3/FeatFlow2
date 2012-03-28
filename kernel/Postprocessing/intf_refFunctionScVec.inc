  ! This file describes the interface for a coefficient function that is called
  ! during the matrix assembly of scalar vectors.

  interface

  !<subroutine>

    subroutine frefFunction (icomponent, cderivative, rdiscretisation, &
                  nelements, npointsPerElement, Dpoints, rdomainIntSubset, &
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
    integer, intent(in)                          :: icomponent

    ! This is a DER_xxxx derivative identifier (from derivative.f90) that
    ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
    ! The result must be written to the Dvalue-array below.
    integer, intent(in)                          :: cderivative

    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in)    :: rdiscretisation

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                          :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                          :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(in)       :: Dpoints

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)          :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional  :: rcollection

  !</input>

  !<output>
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints, or the appropriate derivative
    ! of the function, respectively, according to cderivative.
    !   DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(out)                  :: Dvalues
  !</output>

  !</subroutine>

    end subroutine

  end interface

  !****************************************************************************

  interface

  !<subroutine>

    subroutine fweightFunction (icomponent,rdiscretisation, nelements,&
                                npointsPerElement, Dpoints, IdofsTest,&
                                rdomainIntSubset,  Dvalues, rcollection)

    use fsystem
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection

  !<description>
    ! This subroutine is called during the calculation of errors. It
    ! has to compute the values of a weighting function in a couple of
    ! points on a couple of elements. These values are multiplied by
    ! the calculated error.
    !
    ! The routine accepts a set of elements and a set of points on
    ! these elements (cubature points) in real coordinates.  According
    ! to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points.
  !</description>

  !<input>
    ! Specifies which component of the vector field is to be evaluated
    integer, intent(in)                          :: icomponent

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
    ! This array has to receive the values of the weights in all the
    ! points specified in Dpoints, or the appropriate derivative of
    ! the function, respectively, according to cderivative.
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(out) :: Dvalues
  !</output>

  !</subroutine>

    end subroutine

  end interface

