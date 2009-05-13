  ! -*- mode: f90; -*-
  
  ! This file describes the interface for a coefficient function that is called
  ! during the matrix assembly in the streamline diffusion method.
  
  interface
  
  !<subroutine>

    subroutine ffunctionCoefficient (cterm,rdiscretisation, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset, &
                  Dcoefficients,rcollection)
    
    use fsystem
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
  !<description>
    ! This subroutine is called during the calculation of the SD operator. It has to 
    ! compute the coefficients in front of the terms.
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
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in trial space,Number of elements)
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
    ! This array has to receive the values of the coefficients
    ! in all the points specified in Dpoints.
    ! cterm specifies what to evaluate.
    real(DP), dimension(:,:), intent(out) :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    end subroutine
    
  end interface
  
