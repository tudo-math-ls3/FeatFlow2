  ! -*- mode: f90; -*-
  
  ! This file describes the interface for a simple callback function
  ! that calculates values in specific points of a domain.
  
  interface
  
  !<subroutine>

    subroutine ffunctionRefSimple (nelements,npointsPerElement,Ielements,&
                                   Dpoints,Dvalues,rcollection)
    
    use fsystem
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
  !<description>
    ! This function defines a common interface for a simple callback function.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in in real coordinates.
    ! It has to compute function values simultaneously for all these points
    ! and to store them in Dvalues
  !</description>
    
  !<input>
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! List of elements in the current triangulation where to compute the values.
    integer, dimension(:), intent(in)                           :: Ielements
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints.
    !   DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(OUT)                      :: Dvalues
  !</output>
    
  !</subroutine>
  
    end subroutine
    
  end interface
  
