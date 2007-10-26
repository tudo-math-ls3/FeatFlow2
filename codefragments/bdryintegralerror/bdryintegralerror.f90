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

MODULE bdryintegralerror

  USE fsystem
  USE boundaryintegral

  PRIVATE :: ffunctionBoundaryError

CONTAINS

  !****************************************************************************

!<subroutine>

  SUBROUTINE bdrierr_referenceFunction (dx,dy,dvalue,dderivX,dderivY)

!<description>
  ! USER DEFINED CALLBACK ROUTINE.
  ! This routine calculates in the point (dx,dy) the value of the reference
  ! function $u$.
!</description>

!<input>
  ! Coordinates of the point.
  REAL(DP), INTENT(IN) :: dx,dy
!</input>

!<output>
  ! Value of the reference function.
  REAL(DP), INTENT(OUT) :: dvalue

  ! X-derivative
  REAL(DP), INTENT(OUT) :: dderivX

  ! Y-derivative
  REAL(DP), INTENT(OUT) :: dderivY
!</output>

    dvalue = 0.0_DP
    dderivX = 0.0_DP
    dderivY = 0.0_DP

  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE bdrierr_calcError (rvectorScalar,ccubType,dvalue,rboundaryRegion)

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
  TYPE(t_vectorScalar), INTENT(IN) :: rvectorScalar

  ! A line cubature formula CUB_xxxx_1D to be used for line integration.
  INTEGER, INTENT(IN) :: ccubType
  
  ! OPTIONAL: A t_boundaryRegion specifying the boundary region where
  ! to calculate. If not specified, the computation is done over
  ! the whole boundary.
  TYPE(t_boundaryRegion), INTENT(IN), OPTIONAL :: rboundaryRegion
!</input>

!<output>
  ! The calculated value of the integral.
  REAL(DP), INTENT(OUT) :: dvalue
!</output>

    ! local variables
    TYPE(t_spatialDiscretisation), POINTER :: p_rdiscr
    TYPE(t_collection) :: rcollection
    
    ! Get the underlying discretisation
    p_rdiscr => rvectorScalar%p_rspatialDiscretisation

    ! Create a collection that allows us to pass the vector to
    ! the callback routine
    CALL collct_init(rcollection)
    
    ! Save the vector to the collection, so we can restore it in
    ! the callback routine
    CALL collct_setvalue_vecsca (rcollection, 'vector', rvectorScalar, .TRUE.) 

    CALL bdint_integral2D (p_rdiscr%p_rboundary,p_rdiscr%p_rtriangulation,&
      ccubType,ffunctionBoundaryError,dvalue,rboundaryRegion,rcollection)
      
    ! Release the collection
    CALL collct_deletevalue (rcollection,'vector')
    CALL collct_done(rcollection)

  END SUBROUTINE

  !****************************************************************************

  !<subroutine>

  SUBROUTINE ffunctionBoundaryError (DpointsRef, Dpoints, ibct, DpointPar, Ielements, &
        p_rcollection, Dvalues)
    
    USE basicgeometry
    USE triangulation
    USE collection
    USE scalarpde
    USE domainintegration
    
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

    ! A pointer to a collection structure to provide additional 
    ! information to the coefficient routine. May point to NULL() if not defined.
    TYPE(t_collection), POINTER                      :: p_rcollection
  !</input>
  
  !<output>
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints, or the appropriate derivative
    ! of the function, respectively, according to cderivative.
    !   DIMENSION(npointsPerElement,nelements)
    REAL(DP), DIMENSION(:,:), INTENT(OUT)                      :: Dvalues
  !</output>
    
  !</subroutine>

    ! local variables
    TYPE(t_vectorScalar), POINTER :: p_rvector
    TYPE(t_boundary), POINTER :: p_rboundary
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Duh
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Du
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER :: ipoint
    REAL(DP) :: dminPar,dmaxPar
    REAL(DP) :: dnx,dny,dt
    
    ! Get the vector with the FE function from the collection
    ! and to the boundary definition
    p_rvector => collct_getvalue_vecsca (p_rcollection, 'vector')
    p_rboundary => p_rvector%p_rspatialDiscretisation%p_rboundary

    ! Allocate memory for the FE-function and for the analytical values.
    ! DU(1,:,:) = function value,
    ! DU(2,:,:) = X-derivative,
    ! DU(3,:,:) = Y-derivative
    ALLOCATE(Duh(UBOUND(Dvalues,1),UBOUND(Dvalues,2)))
    ALLOCATE(Du(3,UBOUND(Dvalues,1),UBOUND(Dvalues,2)))
    
    DO iel=1,SIZE(Ielements)
      ! Evaluate the FE function in the given points.
      CALL fevl_evaluate_mult (DER_FUNC, Duh(:,iel), p_rvector, &
          Ielements(iel), DpointsRef(:,:,iel), Dpoints(:,:,iel))
          
      ! Evaluate the analytical function in the given points
      DO ipoint = 1,UBOUND(Dvalues,1)
        CALL bdrierr_referenceFunction (&
            Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel),&
            Du(1,ipoint,iel),Du(2,ipoint,iel),Du(3,ipoint,iel))
      END DO
    END DO

    ! Get the minimum and maximum parameter value. The point with the minimal
    ! parameter value is the start point of the interval, the point with the
    ! maximum parameter value the endpoint.
    dminPar = DpointPar(1,1)
    dmaxPar = DpointPar(1,1)
    DO iel=1,SIZE(Ielements)
      DO ipoint = 1,UBOUND(Dpoints,1)
        dminPar = MIN(DpointPar(ipoint,iel),dminPar)
        dmaxPar = MAX(DpointPar(ipoint,iel),dmaxPar)
      END DO
    END DO

    ! Calculate the value of the term in the integral in that point.
    DO iel=1,SIZE(Ielements)
    
      ! Evaluate the analytical function in the given points
      DO ipoint = 1,UBOUND(Dvalues,1)
       
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
        IF (dt .EQ. dminPar) THEN
          ! Start point
          CALL boundary_getNormalVec(p_rboundary, ibct, &
              dt, dnx, dny,BDR_NORMAL_RIGHT,BDR_PAR_LENGTH)
        ELSE IF (dt .EQ. dmaxPar) THEN
          ! End point
          CALL boundary_getNormalVec(p_rboundary, ibct, &
              dt, dnx, dny,BDR_NORMAL_LEFT,BDR_PAR_LENGTH)
        ELSE
          ! Inner point
          CALL boundary_getNormalVec(p_rboundary, ibct, &
              dt, dnx, dny,cparType=BDR_PAR_LENGTH)
        END IF
        
        ! Calculate the function value: (u-u_h) d_n(u)
        Dvalues(ipoint,iel) = &
            ( Du(1,ipoint,iel) - Duh(ipoint,iel)) * &
            ( Du(2,ipoint,iel)*dnx + Du(3,ipoint,iel)*dny )  
            
      END DO
      
    END DO
    
    ! Release memory
    DEALLOCATE(Du,Duh)

  END SUBROUTINE
  
END MODULE
