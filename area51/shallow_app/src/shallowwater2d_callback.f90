module shallowwater2d_callback

  use fsystem
  use storage
  use triangulation
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use boundary
  use linearformevaluation
  use derivatives
  use cubature

  implicit none

  integer, parameter				:: nvar2d = 3


contains


!<subroutine>

  subroutine shlw_SourceTermCB(rdiscretisation, rform,&
      nelements, npointsPerElement, Dpoints, IdofsTest,&
      rdomainIntSubset, Dcoefficients, rcollection)

    use basicgeometry
    use collection
    use domainintegration
    use fparser
    use scalarpde
    use triangulation
    use feevaluation
    
!<description>
    ! This subroutine is called during the vector assembly. It has to
    ! compute the coefficients in front of the terms of the linear
    ! form. This routine can be used universaly for arbitrary linear
    ! forms for which the coefficients are evaluated analytically
    ! using a function parser which is passed using the collection.
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
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in) :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
  
!<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>
    
!</subroutine>

    ! local variables
    real(DP), dimension(NDIM2D) :: Dpoint
    integer :: itermCount, ipoint, iel, ndim, icomp
    real(DP), dimension(:), allocatable :: Dheigthvalues, Dbottomvalues
    type(t_vectorBlock), pointer :: p_rsolBlock, p_rbottomBlock

    ! Direction to derive the bottom function
    integer :: deriv_direction



    ! Get the pointer to the solution vector from the collection
    p_rsolBlock => Rcollection%p_rvectorQuickAccess1
    
    ! Get the pointer to the bottom vector from the collection
    p_rbottomBlock => Rcollection%p_rvectorQuickAccess2
    
    ! Get the direction to derive the bottom profile
    deriv_direction = Rcollection%IquickAccess(1)
    
    ! Use a higher cubature rule
    p_rsolBlock%RvectorBlock(1)%p_rspatialDiscr%Relementdistr(1)%ccubtypelinform = CUB_G5X5
    

    ! Allocate temp arrays
    allocate(Dheigthvalues(npointsPerElement))
    allocate(Dbottomvalues(npointsPerElement))

    ! Loop over all components of the linear form
    do itermCount = 1, ubound(Dcoefficients,1)

        ! Set number of spatial dimensions
        ndim = size(Dpoints, 1)

        do iel = 1, nelements

          ! Now evaluate the heigth of the water surface
          call fevl_evaluate (DER_FUNC, Dheigthvalues, p_rsolBlock%RvectorBlock(1), Dpoints(:,:,iel))

          ! Now evaluate the derivative of the bottom profile
          call fevl_evaluate (deriv_direction, Dbottomvalues, p_rbottomBlock%RvectorBlock(1), Dpoints(:,:,iel))

          do ipoint = 1, npointsPerElement

            ! Save the coefficients, that are given back by the callback function
            Dcoefficients(itermCount,ipoint,iel) = Dbottomvalues(ipoint) * Dheigthvalues(ipoint)
            

          end do
        end do

    end do ! itermCount

    ! Release the temp array
    deallocate(Dheigthvalues,Dbottomvalues)

  end subroutine


end module shallowwater2d_callback
