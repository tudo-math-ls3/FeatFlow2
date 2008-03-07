  ! This file describes the interface for a function that is used in 
  ! bdint_scalarBoundary2D to calculate integral values on the boundary.
  
  INTERFACE
  
  !<subroutine>

    SUBROUTINE ffunction (DpointsRef, Dpoints, ibct, DpointPar, Ielements, &
        Dvalues,rcollection)
    
    USE basicgeometry
    USE triangulation
    USE collection
    USE scalarpde
    USE domainintegration
    
  !<description>
    ! This subroutine is called during the calculation of boundary
    ! integrals. It has to compute the (analytical) values of a 
    ! function in a couple of points Dpoints on a couple of elements
    ! Ielements. The integral evaluation function uses these values
    ! to compute integrals.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in in real and reference coordinates.
    ! It has to to simultaneously compute the desired values for all these points.
  !</description>
    
  !<input>
    ! This is an array of all points on all the elements where coefficients
    ! are needed. It specifies the coordinates of the points where
    ! information is needed. These coordinates correspond to the reference
    ! element.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: DpointsRef
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed. It specifies the coordinates of the points where
    ! information is needed. These coordinates are world coordinates,
    ! i.e. on the real element.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: Dpoints
    
    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    INTEGER, INTENT(IN) :: ibct
    
    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: DpointPar
    
    ! This is a list of elements (corresponding to Dpoints) where information
    ! is needed. To an element iel=Ielements(i), the array Dpoints(:,:,i)
    ! specifies the points where information is needed.
    ! DIMENSION(nelements)
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), INTENT(IN) :: Ielements

    ! Optional: A collection structure to provide additional 
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
  
    END SUBROUTINE
    
  END INTERFACE
  
