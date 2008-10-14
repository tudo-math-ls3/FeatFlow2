  ! This file describes the interface for a reference function that
  ! is used in pperr_scalarBoundary2D to calculate errors on the boundary.
  
  interface
  
  !<subroutine>

    subroutine ffunctionReference (cderivative,rdiscretisation, &
                  DpointsRef, Dpoints, Ielements, Dvalues, rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the calculation of errors with boundary
    ! integrals. It has to compute the (analytical) values of a 
    ! function in a couple of points Dpoints on a couple of elements
    ! Ielements. These values are compared to those of a computed FE 
    ! function to calculate an error.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in in real and reference coordinates.
    ! It has to to simultaneously compute the desired values for all these points.
  !</description>
    
  !<input>
    ! This is a DER_xxxx derivative identifier (from derivative.f90) that
    ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
    ! The result must be written to the Dvalue-array below.
    integer, intent(IN)                                         :: cderivative
  
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed. It specifies the coordinates of the points where
    ! information is needed. These coordinates correspond to the reference
    ! element.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: DpointsRef
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed. It specifies the coordinates of the points where
    ! information is needed. These coordinates are world coordinates,
    ! i.e. on the real element.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints
    
    ! This is a list of elements (corresponding to Dpoints) where information
    ! is needed. To an element iel=Ielements(i), the array Dpoints(:,:,i)
    ! specifies the points where information is needed.
    ! DIMENSION(nelements)
    integer(PREC_ELEMENTIDX), dimension(:), intent(IN) :: Ielements

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT), optional      :: rcollection
  !</input>
  
  !<output>
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints, or the appropriate derivative
    ! of the function, respectively, according to cderivative.
    !   DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(OUT)                      :: Dvalues
  !</output>
    
  !</subroutine>
  
    end subroutine
    
  end interface
  
