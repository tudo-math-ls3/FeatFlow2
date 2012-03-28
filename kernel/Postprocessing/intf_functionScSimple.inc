  ! This file describes the interface for a simple callback function
  ! that calculates values in specific points of a domain.

  interface

  !<subroutine>

    subroutine ffunctionRefSimple (nelements,npointsPerElement,Ielements,&
                                   Dpoints,Dvalues,rcollection,&
                                   DpointsRef,Djac,Ddetj)

    use fsystem
    use collection

  !<description>
    ! This function defines a common interface for a simple callback function.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (e.g. cubature points) in real (and reference) coordinates.
    ! It has to compute function values simultaneously for all these points
    ! and to store them in Dvalues.
  !</description>

  !<input>
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                               :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                               :: npointsPerElement

    ! List of elements in the current triangulation where to compute the values.
    integer, dimension(:), intent(in)                 :: Ielements

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! DIMENSION(#dimensions,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)            :: Dpoints

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional       :: rcollection

    ! Optional: Coordinates of the evaluation points in the reference
    ! coordinate system.
    ! If this information is not available, this parameter is omitted.
    ! DIMENSION(#dimensions,npointsPerElement,nelements)
    real(DP), dimension(:,:,:),intent(in), optional   :: DpointsRef

    ! The Jacobian matrix of the mapping between the reference and each
    ! real element, for all points on all elements in progress.
    ! If this information is not available, this parameter is omitted.
    ! DIMENSION(#dimensions*#dimensions,npointsPerElement,nelements)
    real(DP), dimension(:,:,:),intent(in), optional   :: Djac

    ! The Jacobian determinant of the mapping of each point from the
    ! reference element to each real element in progress.
    ! If this information is not available, this parameter is omitted.
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:),intent(in), optional     :: Ddetj
  !</input>

  !<output>
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints.
    !   DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(out)             :: Dvalues
  !</output>

  !</subroutine>

    end subroutine

  end interface

