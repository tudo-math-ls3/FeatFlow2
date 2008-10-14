!##############################################################################
!# ****************************************************************************
!# <name> bdryintegralerror </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module calculates the boundary integral error
!#
!#    $$ int_{\Gamma} (u - u_h) * d_n(u)  dx$$
!#
!# of a FE function $u_h$ to a given analytical function $u$.
!# The analytical function $u$ is specified in the callback function
!# bdrierr_referenceFunction.
!#
!# To calculate the integral, the following routine can be used:
!#
!# 1.) bdrierr_calcError 
!#     -> Calculates the boundary integral error to the reference function $u$.
!# </purpose>
!##############################################################################

module bdryintegralerror

  use fsystem
  use boundaryintegral

  private :: ffunctionBoundaryError

contains

  !****************************************************************************

!<subroutine>

  subroutine bdrierr_referenceFunction (dx,dy,dvalue,dderivX,dderivY)

!<description>
  ! USER DEFINED CALLBACK ROUTINE.
  ! This routine calculates in the point (dx,dy) the value of the reference
  ! function $u$.
!</description>

!<input>
  ! Coordinates of the point.
  real(DP), intent(IN) :: dx,dy
!</input>

!<output>
  ! Value of the reference function.
  real(DP), intent(OUT) :: dvalue

  ! X-derivative
  real(DP), intent(OUT) :: dderivX

  ! Y-derivative
  real(DP), intent(OUT) :: dderivY
!</output>

    dvalue = 0.0_DP
    dderivX = 0.0_DP
    dderivY = 0.0_DP

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bdrierr_calcError (rvectorScalar,ccubType,dvalue,rboundaryRegion)

!<description>
  ! Calculates the boundary integral error to the reference function $u$.
  ! To calculate the analytical values of the reference function,
  ! the subroutine bdrierr_referenceFunction is called.
  !
  ! rboundaryRegion allows to specify a region on the boundary where
  ! the integral is computed; if not specified, the integral is
  ! computed over the whole boundary.
!</description>

!<input>
  ! A scalar FE function to compute a boundary integral of.
  type(t_vectorScalar), intent(IN) :: rvectorScalar

  ! A line cubature formula CUB_xxxx_1D to be used for line integration.
  integer, intent(IN) :: ccubType
  
  ! OPTIONAL: A t_boundaryRegion specifying the boundary region where
  ! to calculate. If not specified, the computation is done over
  ! the whole boundary.
  type(t_boundaryRegion), intent(IN), optional :: rboundaryRegion
!</input>

!<output>
  ! The calculated value of the integral.
  real(DP), intent(OUT) :: dvalue
!</output>

    ! local variables
    type(t_spatialDiscretisation), pointer :: p_rdiscr
    type(t_collection) :: rcollection
    
    ! Get the underlying discretisation
    p_rdiscr => rvectorScalar%p_rspatialDiscretisation

    ! Create a collection that allows us to pass the vector to
    ! the callback routine
    call collct_init(rcollection)
    
    ! Save the vector to the collection, so we can restore it in
    ! the callback routine
    call collct_setvalue_vecsca (rcollection, 'vector', rvectorScalar, .true.) 

    call bdint_integral2D (p_rdiscr%p_rboundary,p_rdiscr%p_rtriangulation,&
      ccubType,ffunctionBoundaryError,dvalue,rboundaryRegion,rcollection)
      
    ! Release the collection
    call collct_deletevalue (rcollection,'vector')
    call collct_done(rcollection)

  end subroutine

  !****************************************************************************

  !<subroutine>

  subroutine ffunctionBoundaryError (DpointsRef, Dpoints, ibct, DpointPar, Ielements, &
        p_rcollection, Dvalues)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! AUXILIARY FUNCTION.
    ! This subroutine is called during the calculation of the boundary
    ! integral, initiated by bdrierr_calcError. It computes the (analytical)
    ! values of a function in a couple of points Dpoints on a couple of elements
    ! Ielements. The integral evaluation function uses these values
    ! to compute integrals.
    !
    ! To evaluate the reference function, bdrierr_referenceFunction is
    ! called.
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
    real(DP), dimension(:,:,:), intent(IN)  :: DpointsRef
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed. It specifies the coordinates of the points where
    ! information is needed. These coordinates are world coordinates,
    ! i.e. on the real element.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints
    
    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(IN) :: ibct
    
    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(IN) :: DpointPar
    
    ! This is a list of elements (corresponding to Dpoints) where information
    ! is needed. To an element iel=Ielements(i), the array Dpoints(:,:,i)
    ! specifies the points where information is needed.
    ! DIMENSION(nelements)
    integer(PREC_ELEMENTIDX), dimension(:), intent(IN) :: Ielements

    ! A pointer to a collection structure to provide additional 
    ! information to the coefficient routine. May point to NULL() if not defined.
    type(t_collection), pointer                      :: p_rcollection
  !</input>
  
  !<output>
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints, or the appropriate derivative
    ! of the function, respectively, according to cderivative.
    !   DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(OUT)                      :: Dvalues
  !</output>
    
  !</subroutine>

    ! local variables
    type(t_vectorScalar), pointer :: p_rvector
    type(t_boundary), pointer :: p_rboundary
    real(DP), dimension(:,:), allocatable :: Duh
    real(DP), dimension(:,:,:), allocatable :: Du
    integer(PREC_ELEMENTIDX) :: iel
    integer :: ipoint
    real(DP) :: dminPar,dmaxPar
    real(DP) :: dnx,dny,dt
    
    ! Get the vector with the FE function from the collection
    ! and to the boundary definition
    p_rvector => collct_getvalue_vecsca (p_rcollection, 'vector')
    p_rboundary => p_rvector%p_rspatialDiscretisation%p_rboundary

    ! Allocate memory for the FE-function and for the analytical values.
    ! DU(1,:,:) = function value,
    ! DU(2,:,:) = X-derivative,
    ! DU(3,:,:) = Y-derivative
    allocate(Duh(ubound(Dvalues,1),ubound(Dvalues,2)))
    allocate(Du(3,ubound(Dvalues,1),ubound(Dvalues,2)))
    
    do iel=1,size(Ielements)
      ! Evaluate the FE function in the given points.
      call fevl_evaluate_mult (DER_FUNC, Duh(:,iel), p_rvector, &
          Ielements(iel), DpointsRef(:,:,iel), Dpoints(:,:,iel))
          
      ! Evaluate the analytical function in the given points
      do ipoint = 1,ubound(Dvalues,1)
        call bdrierr_referenceFunction (&
            Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel),&
            Du(1,ipoint,iel),Du(2,ipoint,iel),Du(3,ipoint,iel))
      end do
    end do

    ! Get the minimum and maximum parameter value. The point with the minimal
    ! parameter value is the start point of the interval, the point with the
    ! maximum parameter value the endpoint.
    dminPar = DpointPar(1,1)
    dmaxPar = DpointPar(1,1)
    do iel=1,size(Ielements)
      do ipoint = 1,ubound(Dpoints,1)
        dminPar = min(DpointPar(ipoint,iel),dminPar)
        dmaxPar = max(DpointPar(ipoint,iel),dmaxPar)
      end do
    end do

    ! Calculate the value of the term in the integral in that point.
    do iel=1,size(Ielements)
    
      ! Evaluate the analytical function in the given points
      do ipoint = 1,ubound(Dvalues,1)
       
        ! Parameter value of the point -- in length parametrisation!
        dt = DpointPar(ipoint,iel)
      
        ! Get the normal vector in the point from the boundary.
        ! Note that the parameter value is in length parametrisation!
        ! When we are at the left or right endpoint of the interval, we
        ! calculate the normal vector based on the current edge.
        ! Without that, the behaviour of the routine may lead to some
        ! confusion if the endpoints of the interval coincide with
        ! the endpoints of a boundary edge. In such a case, the routine
        ! would normally compute the normal vector as a mean on the
        ! normal vectors of the edges adjacent to such a point!
        if (dt .eq. dminPar) then
          ! Start point
          call boundary_getNormalVec2D(p_rboundary, ibct, &
              dt, dnx, dny,BDR_NORMAL_RIGHT,BDR_PAR_LENGTH)
        else if (dt .eq. dmaxPar) then
          ! End point
          call boundary_getNormalVec2D(p_rboundary, ibct, &
              dt, dnx, dny,BDR_NORMAL_LEFT,BDR_PAR_LENGTH)
        else
          ! Inner point
          call boundary_getNormalVec2D(p_rboundary, ibct, &
              dt, dnx, dny,cparType=BDR_PAR_LENGTH)
        end if
        
        ! Calculate the function value: (u-u_h) d_n(u)
        Dvalues(ipoint,iel) = &
            ( Du(1,ipoint,iel) - Duh(ipoint,iel)) * &
            ( Du(2,ipoint,iel)*dnx + Du(3,ipoint,iel)*dny )  
            
      end do
      
    end do
    
    ! Release memory
    deallocate(Du,Duh)

  end subroutine
  
end module
