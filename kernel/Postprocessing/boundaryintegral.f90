!#########################################################################
!# ***********************************************************************
!# <name> boundaryintegral </name>
!# ***********************************************************************
!#
!# <purpose>
!# This module allows to calculate integral values over the boundary of
!# a domain. The following routines can be found here:
!#
!# 1.) bdint_integral2D
!#     -> Calculates the boundary integral
!#            $$ int_{\Gamma} f(x) dx $$
!#        for an arbitrary function $f(x)$.
!#
!# 2.) bdint_scalarBoundaryInt2D
!#     -> Calculates the boundary integral
!#            $$ int_{\Gamma} f_h(x) dx $$
!#        for an arbitrary scalar FE function $f_h(x)$.
!#
!# 3.) bdint_normalDerivativeInt2D
!#     -> Calculates the integral of the normal derivative
!#            $$ int_{\Gamma} \partial_n f_h(x) dx $$
!#        for an arbitrary scalar FE function $f_h(x)$.
!# </purpose>
!#########################################################################

MODULE boundaryintegral

  USE fsystem
  USE storage
  USE boundary
  USE cubature
  USE triangulation
  USE linearsystemscalar
  USE linearsystemblock
  USE spatialdiscretisation
  USE domainintegration
  USE collection
  USE feevaluation

  IMPLICIT NONE

  PRIVATE :: ffunctionFEevaluation2D,ffunctionNormalDeriv2D

CONTAINS

  !****************************************************************************

!<subroutine>

  SUBROUTINE bdint_integral2D (rboundary,rtriangulation,ccubType,&
      ffunction,dvalue,rboundaryRegion,rcollection)

!<description>
  ! This routine calculates a boundary integral of an arbitrary
  ! scalar function in 2D specified by the callback routine ffunction.
  !
  ! rboundaryRegion is a t_boundaryRegion object that allows to
  ! specify the boundary region where the integral should be computed.
  ! If not specified, the integral is computed over the whole boundary.
  !
  ! NOTE: The routine is SLOW but general!
!</description>

!<input>
  ! A boundary object that specifies the analytical boundary of the domain.
  TYPE(t_boundary), INTENT(IN) :: rboundary
  
  ! A triangulation object that specifies the mesh in the domain.
  TYPE(t_triangulation), INTENT(IN) :: rtriangulation

  ! A line cubature formula CUB_xxxx_1D to be used for line integration.
  INTEGER, INTENT(IN) :: ccubType
  
  ! A callback function that provides the analytical reference 
  ! function to which the error should be computed.
  ! If not specified, the reference function is assumed to be zero!
  INCLUDE 'intf_functionScBoundary2D.inc'
  
  ! OPTIONAL: A t_boundaryRegion specifying the boundary region where
  ! to calculate. If not specified, the computation is done over
  ! the whole boundary.
  TYPE(t_boundaryRegion), INTENT(IN), OPTIONAL :: rboundaryRegion
  
  ! OPTIONAL: A collection structure. This structure is given to the
  ! callback function to provide additional information. 
  TYPE(t_collection), INTENT(IN), TARGET, OPTIONAL :: rcollection
!</input>

!<output>
  ! The calculated value of the integral.
  REAL(DP), INTENT(OUT) :: dvalue
!</output>

!</subroutine>

    ! local variables
    TYPE(t_collection), POINTER :: p_rcollection
    TYPE(t_boundaryRegion) :: rboundaryReg
    REAL(DP) :: dlocalValue
    INTEGER :: ibdc
    
    ! Let p_rcollection point to rcollection - or NULL if it's not
    ! given.
    IF (PRESENT(rcollection)) THEN
      p_rcollection => rcollection
    ELSE
      p_rcollection => NULL()
    END IF
    
    IF (rtriangulation%ndim .NE. NDIM2D) THEN
      PRINT *,'pperr_scalarBoundary2D: Only 2D discretisations allowed.'
      CALL sys_halt()
    END IF

    ! If the boundary region is specified, call pperr_scalarBoundary2d_conf
    ! for that boundary region. Otherwise, call pperr_scalarBoundary2d_conf
    ! for all possible boundary regions and sum up the errors.
    IF (PRESENT(rboundaryRegion)) THEN
      CALL bdint_integral2D_conf (rboundary,rtriangulation,ccubType,&
        ffunction,dvalue,rboundaryRegion,p_rcollection)
    ELSE
      dvalue = 0.0_DP
      ! Create a boundary region for each boundary component and call
      ! the calculation routine for that.
      DO ibdc=1,boundary_igetNBoundComp(rboundary)
        CALL boundary_createRegion (rboundary, ibdc, 0, rboundaryReg)
        CALL bdint_integral2D_conf (rboundary,rtriangulation,ccubType,&
          ffunction,dlocalValue,rboundaryReg,p_rcollection)
        dvalue = dvalue + dlocalValue
      END DO
    END IF

  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE bdint_integral2D_conf (rboundary,rtriangulation,ccubType,&
      ffunction,dvalue,rboundaryRegion,rcollection)

!<description>
  ! This routine calculates the error of a given finite element function
  ! in rvector to a given analytical callback function ffunctionReference.
  ! 2D version for double-precision vectors.
!</description>

!<input>
  ! A boundary object that specifies the analytical boundary of the domain.
  TYPE(t_boundary), INTENT(IN) :: rboundary
  
  ! A triangulation object that specifies the mesh in the domain.
  TYPE(t_triangulation), INTENT(IN) :: rtriangulation
  
  ! A line cubature formula CUB_xxxx_1D to be used for line integration.
  INTEGER, INTENT(IN)                      :: ccubType

  ! A t_boundaryRegion specifying the boundary region where
  ! to calculate. 
  TYPE(t_boundaryRegion), INTENT(IN), OPTIONAL :: rboundaryRegion

  ! A callback function that provides the analytical reference 
  ! function to which the error should be computed.
  ! If not specified, the reference function is assumed to be zero!
  INCLUDE 'intf_functionScBoundary2D.inc'

  ! Optional: A collection structure to provide additional 
  ! information to the coefficient routine. 
  TYPE(t_collection), INTENT(INOUT), OPTIONAL      :: rcollection
!</input>

!<output>
  ! The calculated value of the integral.
  REAL(DP), INTENT(OUT) :: dvalue
!</output>

!</subroutine>

    ! local variables
    INTEGER(I32), DIMENSION(:), ALLOCATABLE :: IelementOrientation
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), ALLOCATABLE :: Ielements
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: DedgePosition
    
    INTEGER :: ibdc,ibdcoffset,iedge,ilocaledge,nve
    INTEGER(PREC_ELEMENTIDX) :: NEL,NELbdc,iel
    INTEGER(I32) :: ctrafoType
    
    ! The triangulation structure - to shorten some things...
    INTEGER(I32), DIMENSION(:), POINTER :: p_IboundaryCpIdx
    INTEGER(PREC_EDGEIDX), DIMENSION(:), POINTER :: p_IedgesAtBoundary
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementsAtBoundary
    REAL(DP), DIMENSION(:), POINTER :: p_DedgeParameterValue
    REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoordinates
    REAL(DP), DIMENSION(:), POINTER :: p_DvertexParameterValue
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement

    ! Arrays for cubature points
    REAL(DP), DIMENSION(CUB_MAXCUBP, NDIM3D) :: Dxi1D
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Dxi2D,Dpoints,DpointsRef
    REAL(DP), DIMENSION(CUB_MAXCUBP) :: Domega1D
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Dvalues
    REAL(DP), DIMENSION(NDIM2D,TRIA_MAXNVE) :: Dcoord
    INTEGER :: ncubp,ipoint
    INTEGER(I32) :: icoordSystem
    REAL(DP) :: dlen,dpar1,dpar2

    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: DpointPar
    
    ! Get some pointers and arrays for quicker access
    
    CALL storage_getbase_int (rtriangulation%h_IboundaryCpIdx,&
        p_IboundaryCpIdx)
    CALL storage_getbase_int (rtriangulation%h_IedgesAtBoundary,&
        p_IedgesAtBoundary)
    CALL storage_getbase_int (rtriangulation%h_IelementsAtBoundary,&
        p_IelementsAtBoundary)
    CALL storage_getbase_int2d (rtriangulation%h_IedgesAtElement,&
        p_IedgesAtElement)
    CALL storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    CALL storage_getbase_double2d (rtriangulation%h_DvertexCoords,&
        p_DvertexCoordinates)
    CALL storage_getbase_double (rtriangulation%h_DedgeParameterValue,&
        p_DedgeParameterValue)
    CALL storage_getbase_double (rtriangulation%h_DvertexParameterValue,&
        p_DvertexParameterValue)
        
    ! Boundary component?
    ibdc = rboundaryRegion%iboundCompIdx
    
    ! Number of elements on that boundary component?
    NELbdc = p_IboundaryCpIdx(ibdc+1)-p_IboundaryCpIdx(ibdc)
    
    ! Position of the boundary component?
    ibdcoffset = p_IboundaryCpIdx(ibdc)
        
    ! In a first step, we figure out the elements on the boundary and their
    ! orientation. Allocate arrays that are large enough to hold
    ! even all elements on the boundary if necessary.
    ALLOCATE(Ielements(NELbdc), IelementOrientation(NELbdc))
    
    ! Allocate an array saving the start- and end-parameter values
    ! of the edges on the boundary.
    ALLOCATE(DedgePosition(2,NELbdc))
    
    ! Loop through the edges on the boundary component ibdc.
    ! If the edge is inside, remember the element number and figure out
    ! the orientation of the edge.
    ! NEL counts the total number of elements in the region.
    NEL = 0
    DO iedge = 1,NELbdc
      IF (boundary_isInRegion(rboundaryRegion,ibdc,&
          p_DedgeParameterValue(iedge))) THEN
        NEL = NEL + 1
        
        ! Element number
        Ielements(NEL) = p_IelementsAtBoundary(iedge)
        
        ! Element orientation; i.e. the local number of the boundary edge 
        DO ilocaledge = 1,UBOUND(p_IedgesAtElement,1)
          IF (p_IedgesAtElement(ilocaledge,p_IelementsAtBoundary(iedge)) .EQ. &
              p_IedgesAtBoundary(iedge)) EXIT
        END DO
        IelementOrientation(NEL) = ilocaledge
        
        ! Save the start parameter value of the edge -- in length
        ! parametrisation.
        dpar1 = p_DvertexParameterValue(iedge)
        
        ! Save the end parameter value. Be careful: The last edge
        ! must be treated differently!
        IF (iedge .NE. NELbdc) THEN
          dpar2 = p_DvertexParameterValue(iedge+1)
        ELSE
          dpar2 = boundary_dgetMaxParVal(rboundary,ibdc)
        END IF
        
        DedgePosition(1,NEL) = &
          boundary_convertParameter(rboundary, &
            ibdc, dpar1, rboundaryRegion%cparType, BDR_PAR_LENGTH)
            
        DedgePosition(2,NEL) = &
          boundary_convertParameter(rboundary, &
            ibdc, dpar2, rboundaryRegion%cparType, BDR_PAR_LENGTH)
         
      END IF
    END DO
    
    ! Get the parameter values of the 1D cubature formula
    ! as well as the number of cubature points ncubp
    CALL cub_getCubPoints(ccubType, ncubp, Dxi1D, Domega1D)
    
    ! Map the 1D cubature points to the edges in 2D.
    ALLOCATE(Dxi2D(ncubp,NDIM2D+1,NEL))

    ! The type of the coordinate system may change with every element.
    ! So we may have to switch... ItrialElements in the discretisation
    ! structure informs us about the element type.
    DO iel = 1,NEL
    
      ! How many vertices are on the current element?
      nve = UBOUND(p_IverticesAtElement,1)
      DO WHILE (p_IverticesAtElement(nve,iel) .EQ. 0) 
        nve = nve-1
      END DO
      
      SELECT CASE (nve)
      CASE (3)
        ! Coordinates of cubature points on triangles are given in
        ! barycentric coordinates.
        icoordSystem = TRAFO_CS_BARY2DTRI
      CASE (4)
        ! Coordinates of cubature points on triangles are given in
        ! standard coordinates on the reference element
        icoordSystem = TRAFO_CS_REF2DQUAD
      END SELECT
      
      CALL trafo_mapCubPts1Dto2D(icoordSystem, IelementOrientation(iel), &
          ncubp, Dxi1D, Dxi2D(:,:,iel))
    END DO
    
    ! Transpose the coordinate array such that we get coordinates
    ! we can work with.
    ALLOCATE(DpointsRef(NDIM2D+1,ncubp,NEL))
    DO iel=1,NEL
      DpointsRef(:,:,iel) = TRANSPOSE(Dxi2D(:,:,iel))
    END DO
    
    ! Dxi2D is not needed anymore.
    DEALLOCATE(Dxi2D)
    
    ! Calculate the coordinates of the points on world coordinates
    ALLOCATE(Dpoints(ncubp,NDIM2D+1,NEL))
    
    ! Transformation can be different for all elements
    DO iel = 1,NEL

      ! Get the points forming the element
      DO ipoint = 1,UBOUND(p_IverticesAtElement,1)
        Dcoord(1,ipoint) = &
            p_DvertexCoordinates(1,p_IverticesAtElement(ipoint,iel))
        Dcoord(2,ipoint) = &
            p_DvertexCoordinates(2,p_IverticesAtElement(ipoint,iel))
      END DO

      ! Transform the cubature points
      DO ipoint = 1,ncubp

        ! How many vertices are on the current element?
        nve = UBOUND(p_IverticesAtElement,1)
        DO WHILE (p_IverticesAtElement(nve,iel) .EQ. 0) 
          nve = nve-1
        END DO
        
        SELECT CASE (nve)
        CASE (3)
          ! Coordinates of cubature points on triangles are given in
          ! barycentric coordinates, 2D.
          ctrafoType = TRAFO_ID_LINSIMPLEX + TRAFO_DIM_2D
        CASE (4)
          ! Coordinates of cubature points on triangles are given in
          ! standard coordinates on the reference element, 2D.
          ctrafoType = TRAFO_ID_MLINCUBE + TRAFO_DIM_2D
        END SELECT
      
        CALL trafo_calcRealCoords (ctrafoType,Dcoord,&
            DpointsRef(:,ipoint,iel),Dpoints(:,ipoint,iel))
      END DO
    END DO
    
    ! Calculate the parameter values of the points on the boundary
    ALLOCATE (DpointPar(ncubp,NEL))
    DO iel=1,NEL
      DO ipoint = 1,ncubp
        ! Dxi1D is in [-1,1] while the current edge has parmeter values
        ! [DedgePosition(1),DedgePosition(2)]. So do a linear
        ! transformation to transform Dxi1D into that interval, this 
        ! gives the parameter values in length parametrisation
        CALL mprim_linearRescale(Dxi1D(ipoint,1),-1.0_DP,1.0_DP,&
            DedgePosition(1,iel),DedgePosition(2,iel),DpointPar(ipoint,iel))
      END DO
    END DO
    
    ! So Dxi2 defines the coordinates on the reference element for all
    ! elements. Generally said, we have to evaluate the function in these
    ! points now. For that purpose, we call ffunction!
    
    ALLOCATE (Dvalues(ncubp,NEL))
    Dvalues = 0.0_DP
    
    ! Evaluate the function on the boundary
    CALL ffunction (DpointsRef, Dpoints, rboundaryRegion%iboundCompIdx,&
        DpointPar,Ielements(1:NEL), Dvalues(:,:),rcollection)
          
    ! Now, Dvalues1 contains in Dvalues1(:,:,1) the term
    ! "u(x,y)-u_h(x,y)" -- in every cubature point on every
    ! element. We are finally able to calculate the integral!
    ! That means, run over all the edges and sum up...
    ! (ok, if rvectorScalar is not specified, we have
    !  -u_h(x,y) in Dvalues1(:,:,1), but as we take the square,
    !  it doesn't matter if we have u_h or -u_h there!)
    
    dvalue = 0.0_DP
    DO iel = 1,NEL
    
      ! Get the length of the edge. Let's use the parameter values
      ! on the boundary for that purpose; this is a more general
      ! implementation than using simple lines as it will later 
      ! support isoparametric elements.
      !
      ! The length of the current edge serves as a "determinant"
      ! in the cubature, so we have to divide it by 2 as an edge on 
      ! the unit inverval [-1,1] has length 2.
      dlen = 0.5_DP*(DedgePosition(2,iel)-DedgePosition(1,iel))
    
      DO ipoint = 1,ncubp
        dvalue = dvalue + dlen * Domega1D(ipoint) * (Dvalues(ipoint,iel)**2)
      END DO
    END DO
    
    DEALLOCATE(DpointPar)
    DEALLOCATE(Dvalues)  
      
    ! Release memory

    DEALLOCATE(Dpoints)
      
    DEALLOCATE(DedgePosition)
    DEALLOCATE(Ielements, IelementOrientation)
    
  END SUBROUTINE
  
  !****************************************************************************

!<subroutine>

  SUBROUTINE bdint_scalarBoundaryInt2D (rvectorScalar,ccubType,&
      dvalue,rboundaryRegion)

!<description>
  ! This function calculates the boundary integral
  !    $$ int_{\Gamma} f_h(x) dx $$
  ! for the scalar FE function $f_h$ specified in rvectorScalar.
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

!</subroutine>

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
      ccubType,ffunctionFEevaluation2D,dvalue,rboundaryRegion,rcollection)
      
    ! Release the collection
    CALL collct_deletevalue (rcollection,'vector')
    CALL collct_done(rcollection)
      
  END SUBROUTINE

  !****************************************************************************

  SUBROUTINE ffunctionFEevaluation2D (DpointsRef, Dpoints, ibct, DpointPar, Ielements, &
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

    ! local variables
    TYPE(t_vectorScalar), POINTER :: p_rvector
    INTEGER(PREC_ELEMENTIDX) :: iel

    ! Get the vector with the FE function from the collection
    p_rvector => collct_getvalue_vecsca (rcollection, 'vector')
  
    ! Evaluate the FE function in the given points.
    DO iel=1,SIZE(Ielements)
      CALL fevl_evaluate_mult (DER_FUNC, Dvalues(:,iel), p_rvector, &
          Ielements(iel), DpointsRef(:,:,iel), Dpoints(:,:,iel))
    END DO

  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE bdint_normalDerivativeInt2D (rvectorScalar,ccubType,&
      dvalue,rboundaryRegion)

!<description>
  ! This function calculates the boundary integral
  !    $$ int_{\Gamma} grad(f_h(x))*n dx $$
  ! for the block FE function $f_h$ specified in rvectorScalar.
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

!</subroutine>

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
      ccubType,ffunctionNormalDeriv2D,dvalue,rboundaryRegion,rcollection)
      
    ! Release the collection
    CALL collct_deletevalue (rcollection,'vector')
    CALL collct_done(rcollection)
      
  END SUBROUTINE

  !****************************************************************************

  SUBROUTINE ffunctionNormalDeriv2D (DpointsRef, Dpoints, ibct, DpointPar, Ielements, &
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

    ! local variables
    TYPE(t_vectorScalar), POINTER :: p_rvector
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER :: ipoint
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: DderivX,DderivY
    REAL(DP) :: dnx,dny,dminPar,dmaxPar,dt

    ! Get the vector with the FE function from the collection
    p_rvector => collct_getvalue_vecsca (rcollection, 'vector')
  
    ! Allocate memory for the values of the derivative
    ALLOCATE(DderivX(UBOUND(Dvalues,1),UBOUND(Dvalues,2)))
    ALLOCATE(DderivY(UBOUND(Dvalues,1),UBOUND(Dvalues,2)))
  
    ! Evaluate the derivative of the FE function in the given points.
    DO iel=1,SIZE(Ielements)
      CALL fevl_evaluate_mult (DER_DERIV_X, DderivX(:,iel), p_rvector, &
          Ielements(iel), DpointsRef(:,:,iel), Dpoints(:,:,iel))
      CALL fevl_evaluate_mult (DER_DERIV_Y, DderivY(:,iel), p_rvector, &
          Ielements(iel), DpointsRef(:,:,iel), Dpoints(:,:,iel))
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
    
    ! Multiply the derivative with the normal in each point to get the
    ! normal derivative.
    DO iel=1,SIZE(Ielements)
      DO ipoint = 1,UBOUND(Dpoints,1)
      
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
        IF (DpointPar(ipoint,iel) .EQ. dminPar) THEN
          ! Start point
          CALL boundary_getNormalVec(p_rvector%p_rspatialDiscretisation%p_rboundary, ibct, &
              dt, dnx, dny,BDR_NORMAL_RIGHT,BDR_PAR_LENGTH)
        ELSE IF (DpointPar(ipoint,iel) .EQ. dmaxPar) THEN
          ! End point
          CALL boundary_getNormalVec(p_rvector%p_rspatialDiscretisation%p_rboundary, ibct, &
              dt, dnx, dny,BDR_NORMAL_LEFT,BDR_PAR_LENGTH)
        ELSE
          ! Inner point
          CALL boundary_getNormalVec(p_rvector%p_rspatialDiscretisation%p_rboundary, ibct, &
              dt, dnx, dny,cparType=BDR_PAR_LENGTH)
        END IF
      
        ! Calculate the normal derivative
        Dvalues(ipoint,iel) = DderivX(ipoint,iel)*dnx + DderivY(ipoint,iel)*dny
      
      END DO
    END DO
    
    ! Release memory
    DEALLOCATE(DderivX,DderivY)

  END SUBROUTINE

END MODULE
