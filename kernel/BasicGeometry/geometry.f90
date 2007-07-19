!##############################################################################
!# ****************************************************************************
!# <name> geometry </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module implements some basic 2D/3D geometry objects...
!# TODO: Add a more detailed description of this module.
!#
!#
!# The following routines can be found in this module:
!#
!#  1.) geom_init_circle
!#      -> Creates a 2D circle object.
!#
!#  2.) geom_init_square
!#      -> Creates a 2D square object.
!#
!#  3.) geom_init_ellipse
!#      -> Creates a 2D ellipse object.
!#
!#  4.) geom_init_rectangle
!#      -> Creates a 2D rectangle object.
!#
!#  5.) geom_init_polygon
!#      -> Creates a 2D polygon object.
!#
!#  6.) geom_moveto
!#      -> Moves the origin of a 2D/3D object to a specified point.
!#
!#  7.) geom_rotate2D
!#      -> Overwrites the rotation and scaling factor of a 2D object.
!#
!#  8.) geom_isInGeometry
!#      -> Checks whether a given point is inside a geometry object.
!#
!#  9.) geom_isInGeometryArray
!#      -> Checks whether an array of given points is inside a geometry object.
!#
!# 10.) geom_projectToBoundary
!#      -> Projects a point onto the boundary of a geometry object.
!#
!# 11.) geom_calcSignedDistance
!#      -> Calculates the shortest signed distance between a given point
!#         and a geometry object.
!#
!# 12.) geom_calcSignedDistanceArray
!#      -> Calculates the shortest signed distance between an array of given
!#         points and a geometry object.
!#
!# 13.) geom_polygonise
!#      -> Converts a (non-composed) geometry object to a polygon and
!#         optionally converts the vertices into world coordinates.
!#
!# </purpose>
!##############################################################################
 
MODULE geometry

  USE basicgeometry
  USE storage

  IMPLICIT NONE

!<constants>

!<constantblock description="Geometry type identifiers">

  ! Composed geometry object with sub-objects
  INTEGER, PARAMETER :: GEOM_COMPOSED  = -1

  ! No geometry object
  INTEGER, PARAMETER :: GEOM_NONE      = 0

  ! Circle object
  INTEGER, PARAMETER :: GEOM_CIRCLE    = 1

  ! Square object
  INTEGER, PARAMETER :: GEOM_SQUARE    = 2

  ! Ellipse object
  INTEGER, PARAMETER :: GEOM_ELLIPSE   = 3

  ! Rectangle object
  INTEGER, PARAMETER :: GEOM_RECT      = 4

  ! Polygon object
  INTEGER, PARAMETER :: GEOM_POLYGON   = 10

!</constantblock>

!<constantblock description="Polygon type identifiers">

  ! General polygon
  INTEGER, PARAMETER :: GEOM_POLYGON_GENERAL = 0

  ! Convex polygon
  INTEGER, PARAMETER :: GEOM_POLYGON_CONVEX  = 1

!</constantblock>

!<constantblock description="Composed type identifiers">

  ! Union of objects
  INTEGER, PARAMETER :: GEOM_COMP_TYPE_OR = 1
  
  ! Intersection of objects
  INTEGER, PARAMETER :: GEOM_COMP_TYPE_AND = 2
  
!</constantblock>

! *****************************************************************************
! *****************************************************************************
! *****************************************************************************

!<types>
!<typeblock>

  ! This structure realises the subnode for the composed geometry object.
  
  TYPE t_geometryComposed
  
    ! Type of composed object. One of the GEOM_COMP_TYPE_XXXX constants.
    INTEGER :: ccomposedType
    
    ! Number of sub-objects in the composed object.
    INTEGER :: nsubObjects
    
    ! An array of sub-objects with are composed.
    TYPE(t_geometryObject), DIMENSION(:), POINTER :: p_RsubObjects
    
  END TYPE
  
!</typeblock>
  
! *****************************************************************************

!<typeblock>

  ! This structure realises the subnode for the circle object.
  
  TYPE t_geometryCircle
    
    ! Radius of the circle. Must be positive.
    REAL(DP) :: dradius
    
  END TYPE
  
!</typeblock>
  
! *****************************************************************************

!<typeblock>

  ! This structure realises the subnode for the square object.
  TYPE t_geometrySquare
  
    ! Length of each edge of the square. Must be positive.
    REAL(DP) :: dlength
    
  END TYPE

!</typeblock>

! *****************************************************************************

!<typeblock>

  ! This structure realises the subnode for the ellipse object.
  TYPE t_geometryEllipse
  
    ! X- and Y-radius of the ellipse. Must be positive.
    REAL(DP), DIMENSION(2) :: Dradius
    
  END TYPE

!</typeblock>

! *****************************************************************************

!<typeblock>

  ! This structure realises the subnode for the rectangle object.
  TYPE t_geometryRectangle
  
    ! Lengths of the horizontal and vertical edges of the rectangle. Must be
    ! positive.
    REAL(DP), DIMENSION(2) :: Dlength
    
  END TYPE

!</typeblock>

! *****************************************************************************

!<typeblock>

  ! This structure realises the subnode for the polygon object.
  TYPE t_geometryPolygon

    ! Polygon type. One of the GEOM_POLYGON_XXXX constants defined above.
    INTEGER :: npolyType

    ! Polygon vertices
    REAL(DP), DIMENSION(:,:), POINTER :: p_Dvertices

  END TYPE

!</typeblock>

! *****************************************************************************

!<typeblock>

  ! Geometry object structure for 2D and 3D geometry objects.
  TYPE t_geometryObject

    ! Type identifier of the geometry object.
    ! One of the GEOM_**** constants defined above.
    INTEGER                    :: ctype = GEOM_NONE

    ! Dimension of the coordinate system.
    ! May be NDIM2D or NDIM3D (as defined in basicgeometry.f90).
    INTEGER                    :: ndimension = NDIM2D
    
    ! The 2D coordinate system for this geometry object.
    ! Used for 2D geometry objects - undefined for 3D geometry objects.
    TYPE(t_coordinateSystem2D) :: rcoord2D
    
    ! The 3D coordinate system for this geometry object.
    ! Used for 3D geometry objects - undefined for 2D geometry objects.
    !TYPE(t_coordinateSystem3D) :: rcoord3D
    
    ! A boolean which tells us whether the geometry object is inverted or not.
    LOGICAL                    :: binverted
    
    ! Structure for the composed geometry object
    TYPE(t_geometryComposed)   :: rcomposed
    
    ! -=-=-=-=-=-=-=-=-=-=-=
    ! = 2D object subnodes -
    ! -=-=-=-=-=-=-=-=-=-=-=
    ! Structure for the circle object
    TYPE(t_geometryCircle)     :: rcircle
    
    ! Structure for the square object
    TYPE(t_geometrySquare)     :: rsquare
    
    ! Structure for the ellipse object
    TYPE(t_geometryEllipse)    :: rellipse
    
    ! Structure for the rectangle object
    TYPE(t_geometryRectangle)  :: rrectangle

    ! Structure for the polygon object
    TYPE(t_geometryPolygon)    :: rpolygon
    
    ! -=-=-=-=-=-=-=-=-=-=-=
    ! = 3D object subnodes -
    ! -=-=-=-=-=-=-=-=-=-=-=
    ! Todo...
    
  END TYPE
  
!</typeblock>
!</types>

  INTERFACE geom_init_circle
    MODULE PROCEDURE geom_init_circle_indirect
    MODULE PROCEDURE geom_init_circle_direct
  END INTERFACE
  
  INTERFACE geom_init_square
    MODULE PROCEDURE geom_init_square_indirect
    MODULE PROCEDURE geom_init_square_direct
  END INTERFACE
  
  INTERFACE geom_init_ellipse
    MODULE PROCEDURE geom_init_ellipse_indirect
    MODULE PROCEDURE geom_init_ellipse_direct
  END INTERFACE
  
  INTERFACE geom_init_rectangle
    MODULE PROCEDURE geom_init_rectangle_indirect
    MODULE PROCEDURE geom_init_rectangle_direct
  END INTERFACE

  INTERFACE geom_init_polygon
    MODULE PROCEDURE geom_init_polygon_indirect
    MODULE PROCEDURE geom_init_polygon_direct
  END INTERFACE
  
CONTAINS

  ! ***************************************************************************
  ! *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*
  ! *= Composed Object Routines                                              =*
  ! *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE geom_init_composed(rgeomObject, nsubObjects, ccomposedType, &
                                binverted, ndim)
!<description>
  ! Creates a t_geometryObject representing a composed object, i.e. an union
  ! or an intersection of geometry objects.
!</description>

!<input>
  ! The number of sub-objects. Must be positive.
  INTEGER,                     INTENT(IN)  :: nsubObjects
  
  ! The type of the composed object. One of the GEOM_COMP_TYPE_XXXX constants.
  INTEGER,                     INTENT(IN)  :: ccomposedType
  
  ! OPTIONAL: A boolean telling whether the object is inverted.
  ! Is set to .FALSE. if not given.
  LOGICAL, OPTIONAL,           INTENT(IN)  :: binverted

  ! OPTIONAL: The dimension of the composed object. May be NDIM2D or NDIM3D.
  ! If not given, the dimension is set to NDIM2D.
  INTEGER, OPTIONAL,           INTENT(IN)  :: ndim

!</input>

!<output>
  ! A t_geometryObject structure to be written.
  TYPE(t_geometryObject),      INTENT(OUT) :: rgeomObject

!</output>

!</subroutine>

  ! Check if nsubObjects is valid
  IF (nsubObjects .LE. 0) THEN
    RETURN
  END IF

  ! Set the dimension
  IF (PRESENT(ndim)) THEN
    rgeomObject%ndimension = ndim
  ELSE
    rgeomObject%ndimension = NDIM2D
  END IF
  
  ! Set the composed type
  rgeomObject%ctype = GEOM_COMPOSED
  
  ! Is our object inverted?
  IF (PRESENT(binverted)) THEN
    rgeomObject%binverted = binverted
  ELSE
    rgeomObject%binverted = .FALSE.
  END IF
  
  ! Set the operator
  rgeomObject%rcomposed%ccomposedType = ccomposedType
  
  ! Set the number of sub-objects
  rgeomObject%rcomposed%nsubObjects = nsubObjects
  
  ! Allocate sub-objects
  ALLOCATE(rgeomObject%rcomposed%p_RsubObjects(nsubObjects))
  
  ! That's it
  
END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE geom_composed_addNode(rgeomObject, rsubObject, nindex)
  
!<description>
  ! This routine inserts a geometry object into the composed geometry object
  ! tree.
!</description>

!<input>
  ! The index of the child node of the composed geometry object, where
  ! the geometry object is to be attached.
  INTEGER, INTENT(IN) :: nindex
  
!</input>

!<inputoutput>
  ! The composed geometry object
  TYPE(t_geometryObject), INTENT(INOUT) :: rgeomObject
  
  ! The sub-object
  TYPE(t_geometryObject), INTENT(INOUT) :: rsubObject

!</inputoutput>

!</subroutine>

    ! Make sure the object is composed
    IF (rgeomObject%ctype .NE. GEOM_COMPOSED) THEN
      RETURN
    END IF
    
    ! Make sure the index is in range
    IF ((nindex .LE. 0) .OR. (nindex .GT. rgeomObject%rcomposed%nsubObjects)) &
      THEN
      RETURN
    END IF
    
    ! Insert the sub-node
    rgeomObject%rcomposed%p_RsubObjects(nindex) = rsubObject
    
    ! Reset the sub-node
    rsubObject%ctype = GEOM_NONE
    
    ! That's it

  END SUBROUTINE
  
  ! ***************************************************************************
      
!<subroutine>

  RECURSIVE SUBROUTINE geom_composed_isInGeometry (rgeomObject, Dcoords, &
                                                    iisInObject)

!<description>
  ! This routine checks whether a given point is inside the circle or not.
  !
  ! iisInObject is set to 0 if the point is outside the circle, it is set
  ! to 1 if it is inside the circle and is set to -1 if the point is inside
  ! the circle and the circle is inverted.
!</description>

!<input>
  ! The circle against that the point is to be tested.
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
  
  ! The coordinates of the point that is to be tested.
  REAL(DP), DIMENSION(:), INTENT(IN)  :: Dcoords
  
!</input>

!<output>
  ! An integer for the return value.
  INTEGER(I32),           INTENT(OUT) :: iisInObject
!</output>

!</subroutine>

  ! The relative coordinates
  REAL(DP), DIMENSION(3) :: DrelCoords
  
  INTEGER(I32) :: i, iisInSubObject
  
    ! Transform to local coordinate system
    CALL bgeom_transformBackPoint2D(rgeomObject%rcoord2D, Dcoords, DrelCoords)
  
    ! Go through all sub-objects
    IF (rgeomObject%rcomposed%ccomposedType .EQ. GEOM_COMP_TYPE_AND) THEN
    
      ! Let's assume that the point is inside
      iisInObject = 1
      
      DO i=1, rgeomObject%rcomposed%nsubObjects
      
        CALL geom_isInGeometry(rgeomObject%rcomposed%p_RsubObjects(i), &
                               DrelCoords, iisInSubObject)
      
        IF (iisInSubObject .EQ. 0) THEN
          ! The point is outside the sub-object
          iisInObject = 0
          EXIT
        END IF
        
      END DO
      
    ELSE
    
      ! Let's assume the point is outside
      iisInObject = 0
    
      DO i=1, rgeomObject%rcomposed%nsubObjects
      
        CALL geom_isInGeometry(rgeomObject%rcomposed%p_RsubObjects(i), &
                               DrelCoords, iisInSubObject)
      
        IF (iisInSubObject .EQ. 1) THEN
          ! The point is inside the sub-object
          iisInObject = 1
          EXIT
        END IF
        
      END DO
      
    END IF

    ! Maybe the composed object is inverted?
    IF (rgeomObject%binverted) THEN
      iisInObject = 1 - iisInObject
    END IF
  
    ! That's it

  END SUBROUTINE
  
  ! ***************************************************************************
  ! *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*
  ! *= 2D Circle Routines                                                    =*
  ! *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE geom_init_circle_indirect(rgeomObject, rcoordSys, dradius, &
                                       binverted)

!<description>
  ! Creates a t_geometryObject representing a 2D circle.
!</description>

!<input>
  ! A 2D coordinate system for the circle.
  TYPE(t_coordinateSystem2D),  INTENT(IN)  :: rcoordSys
  
  ! A radius for the circle.
  REAL(DP),                    INTENT(IN)  :: dradius
  
  ! OPTIONAL: A boolean telling us whether the object is inverted.
  ! Is set to .FALSE. if not given.
  LOGICAL, OPTIONAL,           INTENT(IN)  :: binverted

!</input>

!<output>
  ! A t_geometryObject structure to be written.
  TYPE(t_geometryObject),      INTENT(OUT) :: rgeomObject

!</output>

!</subroutine>

    ! The dimension is 2D.
    rgeomObject%ndimension = NDIM2D
    
    ! We want a circle.
    rgeomObject%ctype = GEOM_CIRCLE
    
    ! Store the coordinate system.
    rgeomObject%rcoord2D = rcoordSys
    
    ! Is our object inverted?
    IF (PRESENT(binverted)) THEN
      rgeomObject%binverted = binverted
    ELSE
      rgeomObject%binverted = .FALSE.
    END IF
    
    ! Store the radius of the circle
    rgeomObject%rcircle%dradius = dradius
    
    ! That's it!
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE geom_init_circle_direct(rgeomObject, dradius, Dorigin, &
                                     drotation, dscalingFactor, binverted)

!<description>
  ! Creates a t_geometryObject representing a 2D circle.
!</description>

!<input>
  ! A radius for the circle.
  REAL(DP),                          INTENT(IN)  :: dradius
  
  ! OPTIONAL: The origin of the circle.
  ! Is set to (/ 0.0_DP, 0.0_DP /) if not given.
  REAL(DP), DIMENSION(:), OPTIONAL,  INTENT(IN)  :: Dorigin
  
  ! OPTIONAL: The rotation of the circle.
  ! Is set to 0.0_DP if not given.
  REAL(DP), OPTIONAL,                INTENT(IN)  :: drotation
  
  ! OPTIONAL: The scaling factor of the circle.
  ! Is set to 1.0_DP if not given.
  REAL(DP), OPTIONAL,                INTENT(IN)  :: dscalingFactor
  
  ! OPTIONAL: A boolean telling us whether the object is inverted.
  ! Is set to .FALSE. if not given.
  LOGICAL, OPTIONAL,                 INTENT(IN)  :: binverted

!</input>

!<output>
  ! A t_geometryObject structure to be written.
  TYPE(t_geometryObject),            INTENT(OUT) :: rgeomObject

!</output>

!</subroutine>

    ! The dimension is 2D.
    rgeomObject%ndimension = NDIM2D
    
    ! We want a circle.
    rgeomObject%ctype = GEOM_CIRCLE
    
    ! Now we need to create the coordinate system.
    CALL bgeom_initCoordSys2D (rgeomObject%rcoord2D, Dorigin, drotation, &
                               dscalingFactor)
    
    ! Is our object inverted?
    IF (PRESENT(binverted)) THEN
      rgeomObject%binverted = binverted
    ELSE
      rgeomObject%binverted = .FALSE.
    END IF
    
    ! Store the radius of the circle
    rgeomObject%rcircle%dradius = dradius
    
    ! That's it!
  
  END SUBROUTINE

  ! ***************************************************************************
      
!<subroutine>

  SUBROUTINE geom_circle_isInGeometry (rgeomObject, Dcoords, iisInObject)

!<description>
  ! This routine checks whether a given point is inside the circle or not.
  !
  ! iisInObject is set to 0 if the point is outside the circle, it is set
  ! to 1 if it is inside the circle and is set to -1 if the point is inside
  ! the circle and the circle is inverted.
!</description>

!<input>
  ! The circle against that the point is to be tested.
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
  
  ! The coordinates of the point that is to be tested.
  REAL(DP), DIMENSION(:), INTENT(IN)  :: Dcoords
  
!</input>

!<output>
  ! An integer for the return value.
  INTEGER(I32),           INTENT(OUT) :: iisInObject
!</output>

!</subroutine>

  ! We need one local variable for distance calculation
  REAL(DP) :: ddistance

    ! Checking if a point is inside a circle is quite easy.
    ! We can also improve the performance a bit by not calling the
    ! bgeom_transformBackPoint2D routine for the input point, but take care of
    ! the translation and scaling by hand ( and therefore saving the rotation
    ! of a point around a circle ^_^ ).
    
    ! First we need to calculate the (squared) distance of the coordinate
    ! system's origin and the given point.
    ddistance = ((Dcoords(1) - rgeomObject%rcoord2D%Dorigin(1))**2) &
              + ((Dcoords(2) - rgeomObject%rcoord2D%Dorigin(2))**2)
    
    
    ! Now we check if the squared distance is <= than the scaled radius
    ! squared.
    IF (ddistance .LE. ((rgeomObject%rcoord2D%dscalingFactor * &
                         rgeomObject%rcircle%dradius)**2)) THEN
      ! We are inside the circle
      iisInObject = 1
    ELSE
      ! We are not inside the circle
      iisInObject = 0
    END IF

    ! Maybe the circle is inverted?    
    IF (rgeomObject%binverted) THEN
      iisInObject = 1 - iisInObject
    END IF
        
    ! That's it
    
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>
  
  SUBROUTINE geom_circle_prjToBoundary (rgeomObject, Dcoords, Dproj)
  
!<description>
  ! This routine projects a point onto the boundary of a 2D circle.
!</description>

!<input>
  ! The geometry object to calculate the distance from.
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
  
  ! The coordinates of the point that is to be tested.
  ! The array must hold at least 2 entries for a 2D object, and at least 3
  ! entries for a 3D object.
  REAL(DP), DIMENSION(:), INTENT(IN)  :: Dcoords
  
!</input>

!<output>
  ! The coordinates of the boundary projection.
  ! The array must hold at least 2 entries for a 2D object, and at least 3
  ! entries for a 3D object.
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Dproj
  
!</output>

!</subroutine>
  REAL(DP) :: dlen, drad
  
    ! Calculate scaled radius
    drad = rgeomObject%rcoord2D%dscalingFactor * rgeomObject%rcircle%dradius
  
    ! Set the projection (temporarily) to the vector bewteen the origin of the
    ! circle and the point we're given.
    Dproj(1) = Dcoords(1) - rgeomObject%rcoord2D%Dorigin(1)
    Dproj(2) = Dcoords(2) - rgeomObject%rcoord2D%Dorigin(2)
    
    ! Calculate the length of the vector
    dlen = SQRT(Dproj(1)**2 + Dproj(2)**2)
    
    ! If the length is 0, then the given point is the circle's midpoint.
    IF (dlen == 0.0_DP) THEN
      ! Any point on the circle's boundary is okay
      Dproj(1) = rgeomObject%rcoord2D%Dorigin(1) + drad
      Dproj(2) = rgeomObject%rcoord2D%Dorigin(2)

    ELSE
      ! Normalize the vector and scale it by the radius
      Dproj(1) = (Dproj(1) * drad) / dlen
      Dproj(2) = (Dproj(2) * drad) / dlen

    END IF
    
    ! That's it
  END SUBROUTINE

  ! ***************************************************************************
      
!<subroutine>

  SUBROUTINE geom_circle_calcSignedDistance (rgeomObject, Dcoords, ddistance)

!<description>
  ! This routine calculates the signed distance of a given point and a circle.
!</description>

!<input>
  ! The circle against that the point is to be tested.
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
  
  ! The coordinates of the point that is to be tested.
  REAL(DP), DIMENSION(:), INTENT(IN)  :: Dcoords
  
!</input>

!<output>
  ! The shortest signed distance between the point and the circle's boundary.
  REAL(DP),               INTENT(OUT) :: ddistance
!</output>

!</subroutine>

    ! Once again, we will not call the bgeom_transformBackPoint2D routine,
    ! but we will will calculate the signed distance by hand.
    
    ! Calculate the distance of the coordinate system's origin and the given
    ! point and subtract the scaled radius.
    ddistance = SQRT(((Dcoords(1) - rgeomObject%rcoord2D%Dorigin(1))**2) &
                   + ((Dcoords(2) - rgeomObject%rcoord2D%Dorigin(2))**2)) &
                   - (rgeomObject%rcoord2D%dscalingFactor * &
                   rgeomObject%rcircle%dradius)
    
    ! Now we need to check whether the circle is inverted.
    IF (rgeomObject%binverted) THEN
      ddistance = -ddistance
    END IF
        
    ! That's it
    
  END SUBROUTINE
  
  ! ***************************************************************************
 
!<subroutine>
  
  SUBROUTINE geom_circle_polygonise (rgeomObject, hpolyHandle, &
                                     ndesiredVerticeCount)
  
!<description>
  ! This routine converts a circle to a polygon, so that it can
  ! be printed to an output file via ucd_addPolygon (see ucd.f90 for more
  ! details).
!</description>

!<input>
  ! The geometry object to calculate the distance from.
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
  
  ! The desired number of vertices for the produced polygon.
  ! Is only used for circles and ellipses, and is ignored for all other
  ! geometry objects.
  ! If not given, 32 vertices are generated.
  INTEGER, INTENT(IN) :: ndesiredVerticeCount
  
!</input>

!<output>
  ! Handle to a 2D array holding the vertices of the polygon.
  INTEGER, INTENT(OUT) :: hpolyHandle
  
!</output>

!</subroutine>

  INTEGER :: i
  
  REAL(DP), DIMENSION(:,:), POINTER :: p_Dvertices
  REAL(DP) :: dstep, dangle, dradius
  
  INTEGER(I32), DIMENSION(2) :: Isize
  
    ! Calculate angle step
    dstep = (SYS_PI * 2.0_DP) / REAL(ndesiredVerticeCount, DP)
    
    ! Get radius
    dradius = rgeomObject%rcircle%dradius
  
    ! Allocate desired number of vertices
    Isize = (/ 2, ndesiredVerticeCount /)
    CALL storage_new2D('geom_circle_polygonise', 'hpolyHandle', Isize, &
                       ST_DOUBLE, hpolyHandle, ST_NEWBLOCK_NOINIT)

    ! Get vertice array
    CALL storage_getbase_double2D(hpolyHandle, p_Dvertices)
    
    ! Set vertices
    DO i=1, ndesiredVerticeCount
    
      ! Calculate angle
      dangle = dstep * REAL(i, DP)
      
      ! Calculate vertice position
      p_DVertices(1, i) = dradius * COS(dangle)
      p_DVertices(2, i) = dradius * SIN(dangle)
    
    END DO

  END SUBROUTINE
  
  ! ***************************************************************************
  ! *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*
  ! *= 2D Square Routines                                                    =*
  ! *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*
  ! ***************************************************************************
    
!<subroutine>

  SUBROUTINE geom_init_square_indirect(rgeomObject, rcoordSys, dlength, &
                                       binverted)

!<description>
  ! Creates a t_geometryObject representing a 2D square.
!</description>

!<input>
  ! A 2D coordinate system for the square.
  TYPE(t_coordinateSystem2D),  INTENT(IN)  :: rcoordSys
  
  ! The edge length for the square.
  REAL(DP),                    INTENT(IN)  :: dlength
  
  ! OPTIONAL: A boolean telling us whether the object is inverted.
  ! Is set to .FALSE. if not given.
  LOGICAL, OPTIONAL,           INTENT(IN)  :: binverted

!</input>

!<output>
  ! A t_geometryObject structure to be written.
  TYPE(t_geometryObject),      INTENT(OUT) :: rgeomObject

!</output>

!</subroutine>

    ! The dimension is 2D.
    rgeomObject%ndimension = NDIM2D
    
    ! We want a square.
    rgeomObject%ctype = GEOM_SQUARE
    
    ! Store the coordinate system.
    rgeomObject%rcoord2D = rcoordSys
    
    ! Is our object inverted?
    IF (PRESENT(binverted)) THEN
      rgeomObject%binverted = binverted
    ELSE
      rgeomObject%binverted = .FALSE.
    END IF
    
    ! Store the edge length of the square.
    rgeomObject%rsquare%dlength = dlength
    
    ! That's it!
  
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>
  
  SUBROUTINE geom_init_square_direct(rgeomObject, dlength, Dorigin,drotation, &
                                     dscalingFactor, binverted)

!<description>
  ! Creates a t_geometryObject representing a 2D square.
!</description>

!<input>
  ! The edge length for the square.
  REAL(DP),                          INTENT(IN)  :: dlength
  
  ! OPTIONAL: The origin of the square.
  ! Is set to (/ 0.0_DP, 0.0_DP /) if not given.
  REAL(DP), DIMENSION(:), OPTIONAL,  INTENT(IN)  :: Dorigin
  
  ! OPTIONAL: The rotation of the square.
  ! Is set to 0.0_DP if not given.
  REAL(DP), OPTIONAL,                INTENT(IN)  :: drotation
  
  ! OPTIONAL: The scaling factor of the square.
  ! Is set to 1.0_DP if not given.
  REAL(DP), OPTIONAL,                INTENT(IN)  :: dscalingFactor
  
  ! OPTIONAL: A boolean telling us whether the object is inverted.
  ! Is set to .FALSE. if not given.
  LOGICAL, OPTIONAL,                 INTENT(IN)  :: binverted

!</input>

!<output>
  ! A t_geometryObject structure to be written.
  TYPE(t_geometryObject),            INTENT(OUT) :: rgeomObject

!</output>

!</subroutine>

    ! The dimension is 2D.
    rgeomObject%ndimension = NDIM2D
    
    ! We want a square.
    rgeomObject%ctype = GEOM_SQUARE
    
    ! Now we need to create the coordinate system.
    CALL bgeom_initCoordSys2D (rgeomObject%rcoord2D, Dorigin, drotation, &
                               dscalingFactor)
    
    ! Is our object inverted?
    IF (PRESENT(binverted)) THEN
      rgeomObject%binverted = binverted
    ELSE
      rgeomObject%binverted = .FALSE.
    END IF
    
    ! Store the edge length of the square.
    rgeomObject%rsquare%dlength = dlength
    
    ! That's it!
  
  END SUBROUTINE

  ! ***************************************************************************
      
!<subroutine>

  SUBROUTINE geom_square_isInGeometry (rgeomObject, Dcoords, iisInObject)

!<description>
  ! This routine checks whether a given point is inside the square or not.
  !
  ! iisInObject is set to 0 if the point is outside the square, it is set
  ! to 1 if it is inside the square and is set to -1 if the point is inside
  ! the square and the square is inverted.
!</description>

!<input>
  ! The square against that the point is to be tested.
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
    ! The coordinates of the point that is to be tested.
  REAL(DP), DIMENSION(:), INTENT(IN)  :: Dcoords
  
!</input>

!<output>
  ! An integer for the return value.
  INTEGER(I32),           INTENT(OUT) :: iisInObject
!</output>

!</subroutine>

  ! We need one local variable for distance calculation
  REAL(DP) :: ddistance
  
  ! And an array for the transformed X- and Y-coordinates of our point
  REAL(DP), DIMENSION(2) :: DrelCoords

    ! First transfrom the point's coordinates into the square's local
    ! coordinate system.
    CALL bgeom_transformBackPoint2D(rgeomObject%rcoord2D, Dcoords, DrelCoords)
    
    ! Get the distance
    ddistance = MAX(ABS(DrelCoords(1)), ABS(DrelCoords(2)))
    
    ! Check against half of the edge length
    IF (ddistance .LE. (0.5_DP * rgeomObject%rsquare%dlength)) THEN
      ! We are inside the square
      iisInObject = 1
    ELSE
      ! We are outside the square
      iisInObject = 0
    END IF
    
    ! Maybe the square is inverted
    IF (rgeomObject%binverted) THEN
      iisInObject = 1 - iisInObject
    END IF

    ! That's it
    
  END SUBROUTINE

  ! ***************************************************************************
      
!<subroutine>
  
  SUBROUTINE geom_square_prjToBoundary (rgeomObject, Dcoords, Dproj)

!<description>
  ! This routine calculates the projection of a point onto a 2D square's
  ! boundary.
  !
  ! For a detailed description of the method used here, take a look at
  ! the large comment block in the routine geom_square_calcSignedDistance
!</description>

!<input>
  ! The geometry object to calculate the distance from.
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
  
  ! The coordinates of the point that is to be tested.
  ! The array must hold at least 2 entries for a 2D object, and at least 3
  ! entries for a 3D object.
  REAL(DP), DIMENSION(:), INTENT(IN)  :: Dcoords
  
!</input>

!<output>
  ! The coordinates of the boundary projection.
  ! The array must hold at least 2 entries for a 2D object, and at least 3
  ! entries for a 3D object.
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Dproj
  
!</output>

!</subroutine>

  ! Square's half edge length
  REAL(DP) :: dlen
  
  ! Temporary coords in reference coordinate system
  REAL(DP), DIMENSION(2) :: DcoordsRef
  
  ! mirroring factors
  REAL(DP), DIMENSION(2) :: Dmirror = (/ 1.0_DP, 1.0_DP /)
  
    ! Get square's half egde length
    dlen = (rgeomObject%rsquare%dlength * 0.5_DP)
  
    ! First we need to transform the point's coordinates into the square's
    ! local coordinate system.
    CALL bgeom_transformBackPoint2D(rgeomObject%rcoord2D, Dcoords, DcoordsRef)
    
    ! Now mirror the point, so that the resulting point is inside the positive
    ! quadrant.
    IF (DcoordsRef(1) .LT. 0.0) THEN
      ! X-coord of point is negative, so multiply with -1
      Dmirror(1) = -1.0_DP
      DcoordsRef(1) = -DcoordsRef(1)
    END IF
    
    IF (DcoordsRef(2) .LT. 0.0) THEN
      ! Y-coord of point is negative, so multiply with -1
      Dmirror(2) = -1.0_DP
      DcoordsRef(2) = -DcoordsRef(2)
    END IF
    
    ! If both coordinates are greater than half of the square's egde length,
    ! then the projection is the square's corner.
    
    IF ((DcoordsRef(1) .GE. dlen) .AND. (DcoordsRef(2) .GE. dlen)) THEN
    
      ! Save square's corner
      DcoordsRef(1) = dlen
      DcoordsRef(2) = dlen
    
    ELSE
      ! In all other cases, the projection is on an edge of the square.
      ! Now find out which coordinate is greater.
      IF (DcoordsRef(1) .GE. DcoordsRef(2)) THEN
        ! The X-coordinate is greater than the Y-coordinate.
        ! Now we need to set the X-coordinate to half of the square's edge
        ! length and then we have our projection.
        DcoordsRef(1) = dlen
      
      ELSE
        ! Otherwise correct Y-coordinate
        DcoordsRef(2) = dlen
      
      END IF
    
    END IF
    
    ! The projection itself is calculated, now we need to mirror the projection
    ! into the quadrant where our original point was in.
    DcoordsRef(1) = DcoordsRef(1) * Dmirror(1)
    DcoordsRef(2) = DcoordsRef(2) * Dmirror(2)
    
    ! And transform the projection back into world coordinates
    CALL bgeom_transformPoint2D(rgeomObject%rcoord2D, DcoordsRef, Dproj)
    
    ! That's it

  END SUBROUTINE
  
  ! ***************************************************************************
      
!<subroutine>

  SUBROUTINE geom_square_calcSignedDistance(rgeomObject, Dcoords, ddistance)
  
!<description>
  ! This routine calculates the shortest signed distance between a point and
  ! a square's boundary.
!</description>

!<input>
  ! The square against that the point is to be tested.
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
  
  ! The coordinates of the point that is to be tested.
  REAL(DP), DIMENSION(:), INTENT(IN)  :: Dcoords

!</input>

!<output>
  ! The shortest signed distance between the point and the square's boundary.
  REAL(DP),               INTENT(OUT) :: ddistance
!</output>

!</subroutine>

  ! We need one local array for the transformed X- and Y-coordinates of our
  ! point.
  REAL(DP), DIMENSION(2) :: DrelCoords
  
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Now we are going to use some nasty tricks to calculate the shortest
    ! distance, namely: calculating the distance without calculating a
    ! projection of the point onto the square's boundary.
    !
    ! These are the steps we are going to make:
    ! 1. Transform the point into the square's local coordinate system.
    ! 2. "Rotating" the whole coordinate system by a multiple of 90°, such
    !    that the resulting point has non-negative coordinates.
    !   (In fact, we will not "rotate" the point, but simply get the
    !    absolute values of its coordinates, since the square is symmetrical
    !    in respect to his own local coordinate system axes and therefore we
    !    do not need to rotate the whole coordinate system.)
    ! 3. Subtracting the square's vertice with positive coordinates from
    !    the "rotated" point.
    ! 4. Checking whether the resulting point has positive coordinates.
    !    If yes, the projection of our original point onto the square's
    !    boundary is a vertice - so the distance is simply the length of
    !    our "rotated" vertice.
    !    If no, then the projection of our original point lies on an edge
    !    of our square.
    !      If both "rotated" coordinates are negative then we are inside our
    !      square, and the shortest distance is simply the minimum of the
    !      absolute values of the coordinates multiplied with -1, or in other
    !      words: the maximum of the "rotated" coordinates.
    !      Otherwise the point is outside of our square, and one of the two
    !      coordinates is negative and one is non-negative.
    !      The shortest distance is equal to the non-negative coordinate in
    !      this case, or in other words: the maximum of the "rotated"
    !      coordinates.
    ! 5. Scale the calculated distance by the scaling factor of the square's
    !    local coordinate system scaling factor.
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
    ! We are now going to transform our point into the square's local
    ! coordinate system.
    CALL bgeom_transformBackPoint2D(rgeomObject%rcoord2D, Dcoords, DrelCoords)
  
    ! Now get the absolute values of the relative coords, and subtract one
    ! half of the square's edge length
    DrelCoords(1) = ABS(DrelCoords(1)) - (rgeomObject%rsquare%dlength * 0.5_DP)
    DrelCoords(2) = ABS(DrelCoords(2)) - (rgeomObject%rsquare%dlength * 0.5_DP)
  
    ! Check whether the resulting coordinates are both positive
    IF ((DrelCoords(1) > 0.0_DP) .AND. (DrelCoords(2) > 0.0_DP)) THEN
      ddistance = SQRT(DrelCoords(1)**2 + DrelCoords(2)**2)
    ELSE
      ddistance = MAX(DrelCoords(1), DrelCoords(2))
    END IF
    
    ! Is the square inverted?
    IF (rgeomObject%binverted) THEN
      ddistance = -ddistance
    END IF
    
    ! Make sure to scale the distance by the square's coordinate system's
    ! scaling factor
    ddistance = ddistance * rgeomObject%rcoord2D%dscalingFactor
    
    ! That's it
    
  END SUBROUTINE
  
  ! ***************************************************************************
 
!<subroutine>
  
  SUBROUTINE geom_square_polygonise (rgeomObject, hpolyHandle)
  
!<description>
  ! This routine converts a square to a polygon, so that it can
  ! be printed to an output file via ucd_addPolygon (see ucd.f90 for more
  ! details).
!</description>

!<input>
  ! The geometry object to calculate the distance from.
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
  
!</input>

!<output>
  ! Handle to a 2D array holding the vertices of the polygon.
  INTEGER, INTENT(OUT) :: hpolyHandle
  
!</output>

!</subroutine>

  REAL(DP), DIMENSION(:,:), POINTER :: p_Dvertices
  REAL(DP) :: dedge

  INTEGER(I32), DIMENSION(2), PARAMETER :: Isize = (/ 2, 4 /)  
  
    ! Get edge length
    dedge = rgeomObject%rsquare%dlength * 0.5_DP
  
    ! Allocate desired number of vertices
    CALL storage_new2D('geom_square_polygonise', 'hpolyHandle', Isize, &
                       ST_DOUBLE, hpolyHandle, ST_NEWBLOCK_NOINIT)

    ! Get vertice array
    CALL storage_getbase_double2D(hpolyHandle, p_Dvertices)
    
    ! Set coords
    p_Dvertices(1,1) = -dedge
    p_Dvertices(2,1) = dedge
    p_Dvertices(1,2) = -dedge
    p_Dvertices(2,2) = -dedge
    p_Dvertices(1,3) = dedge
    p_Dvertices(2,3) = -dedge
    p_Dvertices(1,4) = dedge
    p_Dvertices(2,4) = dedge
    
  END SUBROUTINE
  
  ! ***************************************************************************
  ! *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*
  ! *= 2D Ellipse Routines                                                   =*
  ! *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE geom_init_ellipse_indirect(rgeomObject, rcoordSys, Dradius, &
                                        binverted)

!<description>
  ! Creates a t_geometryObject representing a 2D ellipse.
!</description>

!<input>
  ! A 2D coordinate system for the ellipse.
  TYPE(t_coordinateSystem2D),  INTENT(IN)  :: rcoordSys
  
  ! An array holding the X- and Y-radiuses for the ellipse.
  REAL(DP), DIMENSION(:),      INTENT(IN)  :: Dradius
  
  ! OPTIONAL: A boolean telling us whether the object is inverted.
  ! Is set to .FALSE. if not given.
  LOGICAL, OPTIONAL,           INTENT(IN)  :: binverted

!</input>

!<output>
  ! A t_geometryObject structure to be written.
  TYPE(t_geometryObject),      INTENT(OUT) :: rgeomObject

!</output>

!</subroutine>

    ! The dimension is 2D.
    rgeomObject%ndimension = NDIM2D
    
    ! We want an ellipse.
    rgeomObject%ctype = GEOM_ELLIPSE
    
    ! Store the coordinate system.
    rgeomObject%rcoord2D = rcoordSys
    
    ! Is our object inverted?
    IF (PRESENT(binverted)) THEN
      rgeomObject%binverted = binverted
    ELSE
      rgeomObject%binverted = .FALSE.
    END IF
    
    ! Store the X- and Y-radius of the ellipse
    rgeomObject%rellipse%Dradius = Dradius(1:2)
    
    ! That's it!
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE geom_init_ellipse_direct(rgeomObject, Dradius, Dorigin, &
                                      drotation, dscalingFactor, binverted)

!<description>
  ! Creates a t_geometryObject representing a 2D ellipse.
!</description>

!<input>
  ! A array holding the X- and Y-radiuses for the ellipse.
  REAL(DP), DIMENSION(:),            INTENT(IN)  :: Dradius
  
  ! OPTIONAL: The origin of the ellipse.
  ! Is set to (/ 0.0_DP, 0.0_DP /) if not given.
  REAL(DP), DIMENSION(:), OPTIONAL,  INTENT(IN)  :: Dorigin
  
  ! OPTIONAL: The rotation of the ellipse.
  ! Is set to 0.0_DP if not given.
  REAL(DP), OPTIONAL,                INTENT(IN)  :: drotation
  
  ! OPTIONAL: The scaling factor of the ellipse.
  ! Is set to 1.0_DP if not given.
  REAL(DP), OPTIONAL,                INTENT(IN)  :: dscalingFactor
  
  ! OPTIONAL: A boolean telling us whether the object is inverted.
  ! Is set to .FALSE. if not given.
  LOGICAL, OPTIONAL,                 INTENT(IN)  :: binverted

!</input>

!<output>
  ! A t_geometryObject structure to be written.
  TYPE(t_geometryObject),            INTENT(OUT) :: rgeomObject

!</output>

!</subroutine>

    ! The dimension is 2D.
    rgeomObject%ndimension = NDIM2D
    
    ! We want a ellipse.
    rgeomObject%ctype = GEOM_ELLIPSE
    
    ! Now we need to create the coordinate system.
    CALL bgeom_initCoordSys2D (rgeomObject%rcoord2D, Dorigin, drotation, &
                               dscalingFactor)
    
    ! Is our object inverted?
    IF (PRESENT(binverted)) THEN
      rgeomObject%binverted = binverted
    ELSE
      rgeomObject%binverted = .FALSE.
    END IF
    
    ! Store the X- and Y-radius of the ellipse
    rgeomObject%rellipse%Dradius = Dradius(1:2)
    
    ! That's it!
  

  END SUBROUTINE

  ! ***************************************************************************
      
!<subroutine>

  SUBROUTINE geom_ellipse_isInGeometry (rgeomObject, Dcoords, iisInObject)

!<description>
  ! This routine checks whether a given point is inside the ellipse or not.
  !
  ! iisInObject is set to 0 if the point is outside the ellipse, it is set
  ! to 1 if it is inside the ellipse and is set to -1 if the point is inside
  ! the ellipse and the ellipse is inverted.
!</description>

!<input>
  ! The ellipse against that the point is to be tested.
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
  
  ! The coordinates of the point that is to be tested.
  REAL(DP), DIMENSION(:), INTENT(IN)  :: Dcoords
  
!</input>

!<output>
  ! An integer for the return value.
  INTEGER(I32),           INTENT(OUT) :: iisInObject
!</output>

!</subroutine>

  ! An array for the transformed X- and Y-coordinates of our point
  REAL(DP), DIMENSION(2) :: DrelCoords
  
    ! First transfrom the point's coordinates into the ellipse's local
    ! coordinate system
    CALL bgeom_transformBackPoint2D(rgeomObject%rcoord2D, Dcoords, DrelCoords)

    ! And scale the coordinates by the inverses of our radiuses.
    ! By doing this, we "transform" our ellipse into a circle with radius 1
    DrelCoords(1) = DrelCoords(1) / rgeomObject%rellipse%Dradius(1)
    DrelCoords(2) = DrelCoords(2) / rgeomObject%rellipse%Dradius(2)
    
    ! Now check the length of our resulting point
    IF ((DrelCoords(1)**2 + DrelCoords(2)**2) .LE. 1.0_DP) THEN
      ! We're inside the ellipse
      iisInObject = 1
      
    ELSE
      ! We are outside the ellipse
      iisInObject = 0
      
    END IF

    ! Is the ellipse inverted?
    IF (rgeomObject%binverted) THEN
      iisInObject = 1 - iisInObject
    END IF
    
    ! That's it
    
  END SUBROUTINE

  ! ***************************************************************************
      
!<subroutine>
  
  SUBROUTINE geom_ellipse_prjToBoundary (rgeomObject, Dcoords, Dproj)

!<description>
  ! This routine calculates the projection of a point onto an ellipse's
  ! boundary.
  !
  ! This is probably the most complicated routine in this module ^_^
  ! The algorithm for the projection was based on the following paper:
  !
  ! -> David Eberly - Distance from a Point to an Ellipse in 2D
  !
  ! The paper can be found on the following page:
  ! -> http://www.geometrictools.com
  !
!</description>

!<input>
  ! The geometry object to calculate the distance from.
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
  
  ! The coordinates of the point that is to be tested.
  ! The array must hold at least 2 entries for a 2D object, and at least 3
  ! entries for a 3D object.
  REAL(DP), DIMENSION(:), INTENT(IN)  :: Dcoords
  
!</input>

!<output>
  ! The coordinates of the boundary projection.
  ! The array must hold at least 2 entries for a 2D object, and at least 3
  ! entries for a 3D object.
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Dproj
  
!</output>

!</subroutine>

  ! Two variables for the ellipses radiuses
  REAL(DP) :: dradX, dradY
  
  ! Some other temporary variables needed for the Newton iteration
  REAL(DP) :: dT, dF, dFDer, dXDivA, dYDivB, dradXSqr, dradYSqr, dprX, dprY
  REAL(DP) :: dratio, dXDivASqr, dYDivBSqr
  LOGICAL :: btranspose = .FALSE.
  INTEGER :: i

  ! Temporary coords in reference coordinate system
  REAL(DP), DIMENSION(2) :: DcoordsRef
  
  ! mirroring factors
  REAL(DP), DIMENSION(2) :: Dmirror = (/ 1.0_DP, 1.0_DP /)

    ! Get the ellipse's X- and Y-radius
    dradX = rgeomObject%rellipse%Dradius(1)
    dradY = rgeomObject%rellipse%Dradius(2)
  
    ! First we need to transform the point's coordinates into the ellipse's
    ! local coordinate system.
    CALL bgeom_transformBackPoint2D(rgeomObject%rcoord2D, Dcoords, DcoordsRef)
    
    ! First we are going to check whether the ellipse degenerates to a line
    ! or even a point.
    IF(dradY .LE. 0.0_DP) THEN
      ! Y-radius is 0
      ! Maybe the ellipse is even a point?
      IF (dradX .LE. 0.0_DP) THEN
        ! Special case: point
        ! The projection is the ellipse's midpoint
        DcoordsRef = Dcoords - rgeomObject%rcoord2D%Dorigin
        
      ELSE
        ! The ellipse is a line on the X-axis.
        ! So we can set the Y-coordinate of the projection to 0.
        DcoordsRef(2) = 0.0_DP
        
        ! Where is the point's X-coordinate in respect to the line?
        IF (DcoordsRef(1) .LE. -dradX) THEN
          ! Left of the line
          DcoordsRef(1) = -dradX
          
        ELSE IF (DcoordsRef(1) .GE. dradX) THEN
          ! Right of the line
          DcoordsRef(1) = dradX;
          
        END IF
        
        ! If the point was on the line, then the X-coordinate does not need
        ! to be modified.
              
      END IF
    
      ! The projection is calculated - transform it to world coordinates.
      CALL bgeom_transformPoint2D(rgeomObject%rcoord2D, DcoordsRef, Dproj)
      
      RETURN
      
    ELSE IF(dradX .LE. 0.0_DP) THEN
      ! In this case the ellipse is a line on the Y-axis
      ! So we can set the X-coordinate of the projection to 0.
      DcoordsRef(1) = 0.0_DP
        
      ! Where is the point's Y-coordinate in respect to the line?
      IF (DcoordsRef(2) .LE. -dradY) THEN
        ! Below of the line
        DcoordsRef(2) = -dradY
          
      ELSE IF (DcoordsRef(2) .GE. dradY) THEN
        ! Above of the line
        DcoordsRef(2) = dradY;
          
      END IF
        
      ! If the point was on the line, then the Y-coordinate does not need
      ! to be modified.
      ! The projection is calculated - transform it to world coordinates.
      CALL bgeom_transformPoint2D(rgeomObject%rcoord2D, DcoordsRef, Dproj)
      
      RETURN
      
    END IF
    
    ! One more special case: The ellipse is a circle.
    IF (dradX == dradY) THEN
      ! Get vector between circle's midpoint and point
      Dproj = Dcoords - rgeomObject%rcoord2D%Dorigin
      
      ! Get length of vector
      dT = SQRT(Dproj(1)**2 + Dproj(2)**2)
      
      ! Is the point equal to the circle's midpoint?
      IF (dT == 0.0_DP) THEN
        ! Then the projection is the circle's midpoint
        Dproj = Dcoords
        
        RETURN
        
      END IF
      
      ! Scale vector
      Dproj = (Dproj * rgeomObject%rcoord2D%dscalingFactor * dradX) / dT
      
      RETURN
    
    END IF
    
    ! Now here comes the more interesting part - the ellipse is not a circle,
    ! line or point.

    ! The first thing we have to do is to move our point into the positive
    ! quadrant.
    IF (DcoordsRef(1) .LT. 0.0_DP) THEN
      ! X-coord of point is negative, so multiply with -1
      Dmirror(1) = -1.0_DP
      DcoordsRef(1) = -DcoordsRef(1)
    ELSE 
      Dmirror(1) = 1.0_DP
    END IF
    
    IF (DcoordsRef(2) .LT. 0.0_DP) THEN
      ! Y-coord of point is negative, so multiply with -1
      Dmirror(2) = -1.0_DP
      DcoordsRef(2) = -DcoordsRef(2)
    ELSE
      Dmirror(2) = 1.0_DP
    END IF
    
    ! Now maybe we need to transpose the ellipse?
    IF (dradX .LT. dradY) THEN
      btranspose = .TRUE.

      ! Change ellipse radius and point coordinates
      dT = dradX
      dradX = dradY
      dradY = dT
      dT = DcoordsRef(1)
      DcoordsRef(1) = DcoordsRef(2)
      DcoordsRef(2) = dT

    END IF
    
    ! Now where is the point we want to project?
    IF (DcoordsRef(1) .NE. 0.0_DP) THEN
    
      IF (DcoordsRef(2) .NE. 0.0_DP) THEN
      
        ! The point to be projection is not on one of the axes
        ! This is the part where we are going to start that Newton-iteration.
        
        ! initial guess
        dT = dradY * (DcoordsRef(2) - dradY)
        
        ! Some precalculations
        dradXSqr = dradX**2
        dradYSqr = dradY**2
        dprX = dradX * DcoordsRef(1)
        dprY = dradY * DcoordsRef(2)
        
        ! iterate
        DO i = 1, 50
          
          dXDivA = dprX / (dT + dradXSqr)
          dYDivB = dprY / (dT + dradYSqr)
          dXDivASqr = dXDivA**2
          dYDivBSqr = dYDivB**2
          
          dF = dXDivASqr + dYDivBSqr - 1.0_DP
          
          IF (dF .LE. 1D-10) THEN
            EXIT
          END IF
          
          dFDer = 2.0_DP * ((dXDivASqr / (dT + dradXSqr)) + &
                            (dYDivBSqr / (dT + dradYSqr)))
          
          dratio = dF / dFDer
          IF (dRatio .LE. 1D-10) THEN
            EXIT
          END IF

          dT = dT + dRatio
        END DO

        DcoordsRef(1) = dXDivA * dradX
        DcoordsRef(2) = dYDivB * dradY
        
      ELSE
        
        dradYSqr = dradY**2
        IF (DcoordsRef(1) .LT. dradX - (dradYSqr / dradX)) THEN
        
          dradXSqr = dradX**2
          DcoordsRef(1) = dradXSqr * DcoordsRef(1) / (dradXSqr - dradYSqr)
          dXDivA = DcoordsRef(1) / dradX
          DcoordsRef(2) = dradY * SQRT(ABS(1.0 - dXDivA**2))
          
        ELSE
          DcoordsRef(1) = dradX
          DcoordsRef(2) = 0.0_DP
          
        END IF
        
      END IF
      
    ELSE
      DcoordsRef(1) = 0.0_DP
      DcoordsRef(2) = dradY
    END IF
    
    ! Do we need to transpose the result?
    IF (btranspose) THEN
      dT = DcoordsRef(1)
      DcoordsRef(1) = DcoordsRef(2)
      DcoordsRef(2) = dT
    END IF
    
    ! Multiplty with mirror
    DcoordsRef(1) = DcoordsRef(1) * Dmirror(1)
    DcoordsRef(2) = DcoordsRef(2) * Dmirror(2)
    
    ! And transform the projection back into world coordinates
    CALL bgeom_transformPoint2D(rgeomObject%rcoord2D, DcoordsRef, Dproj)
    
    ! That's it
      
  END SUBROUTINE
  
  ! ***************************************************************************
      
!<subroutine>

  SUBROUTINE geom_ellipse_calcSignedDistance (rgeomObject, Dcoords, ddistance)

!<description>
  ! This routine calculates the signed distance of a given point and a ellipse.
!</description>

!<input>
  ! The ellipse against that the point is to be tested.
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
  
  ! The coordinates of the point that is to be tested.
  REAL(DP), DIMENSION(:), INTENT(IN)  :: Dcoords
  
!</input>

!<output>
  ! The shortest signed distance between the point and the ellipse's boundary.
  REAL(DP),               INTENT(OUT) :: ddistance
!</output>

!</subroutine>

  ! We need one local variable for temporary distance calculation
  INTEGER :: iInside
  
  ! And an array for the transformed X- and Y-coordinates of our point
  REAL(DP), DIMENSION(2) :: DrelCoords
  
  ! And one array for the projection of our point onto the ellipses boundary
  REAL(DP), DIMENSION(2) :: Dprojection

    
    ! Project the point onto the ellipse's boundary.
    CALL geom_ellipse_prjToBoundary(rgeomObject, DCoords, Dprojection)
    
    ! Now subtract the projection from our point
    DrelCoords = Dcoords - Dprojection
    
    ! Caluclate the new length
    ddistance = SQRT(DrelCoords(1)**2 + DrelCoords(2)**2)

    ! Is the point inside or outside our ellipse?
    CALL geom_ellipse_isInGeometry(rgeomObject, Dcoords, iInside)
    
    ! Now we still need to know whether we are inside the ellipse or not.
    ! This information is implicitly given by dreldist
    IF (iInside .EQ. 1) THEN
      ddistance = -ddistance
    END IF

    ! That's it
    
  END SUBROUTINE
  
  ! ***************************************************************************
 
!<subroutine>
  
  SUBROUTINE geom_ellipse_polygonise (rgeomObject, hpolyHandle, &
                                      ndesiredVerticeCount)
  
!<description>
  ! This routine converts an ellipse to a polygon, so that it can
  ! be printed to an output file via ucd_addPolygon (see ucd.f90 for more
  ! details).
!</description>

!<input>
  ! The geometry object to calculate the distance from.
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
  
  ! The desired number of vertices for the produced polygon.
  ! Is only used for circles and ellipses, and is ignored for all other
  ! geometry objects.
  INTEGER, INTENT(IN) :: ndesiredVerticeCount
  
!</input>

!<output>
  ! Handle to a 2D array holding the vertices of the polygon.
  INTEGER, INTENT(OUT) :: hpolyHandle
  
!</output>

!</subroutine>

  INTEGER :: i
  
  REAL(DP), DIMENSION(:,:), POINTER :: p_Dvertices
  REAL(DP) :: dstep, dangle, dradiusX, dradiusY
  
  INTEGER(I32), DIMENSION(2) :: Isize
  
    ! Calculate angle step
    dstep = (SYS_PI * 2.0_DP) / REAL(ndesiredVerticeCount, DP)
    
    ! Get radius
    dradiusX = rgeomObject%rellipse%Dradius(1)
    dradiusY = rgeomObject%rellipse%Dradius(2)
  
    ! Allocate desired number of vertices
    Isize = (/ 2, ndesiredVerticeCount /)
    CALL storage_new2D('geom_ellipse_polygonise', 'hpolyHandle', Isize, &
                       ST_DOUBLE, hpolyHandle, ST_NEWBLOCK_NOINIT)

    ! Get vertice array
    CALL storage_getbase_double2D(hpolyHandle, p_Dvertices)
    
    ! Set vertices
    DO i=1, ndesiredVerticeCount
    
      ! Calculate angle
      dangle = dstep * REAL(i, DP)
      
      ! Calculate vertice position
      p_DVertices(1, i) = dradiusX * COS(dangle)
      p_DVertices(2, i) = dradiusY * SIN(dangle)
    
    END DO

  END SUBROUTINE

  ! ***************************************************************************
  ! *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*
  ! *= 2D Rectangle Routines                                                 =*
  ! *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*
  ! ***************************************************************************
    
!<subroutine>

  SUBROUTINE geom_init_rectangle_indirect(rgeomObject, rcoordSys, Dlength, &
                                          binverted)

!<description>
  ! Creates a t_geometryObject representing a 2D rectangle.
!</description>

!<input>
  ! A 2D coordinate system for the rectangle.
  TYPE(t_coordinateSystem2D),  INTENT(IN)  :: rcoordSys
  
  ! The edge length for the rectangle.
  REAL(DP), DIMENSION(:),      INTENT(IN)  :: Dlength
  
  ! OPTIONAL: A boolean telling us whether the object is inverted.
  ! Is set to .FALSE. if not given.
  LOGICAL, OPTIONAL,           INTENT(IN)  :: binverted

!</input>

!<output>
  ! A t_geometryObject structure to be written.
  TYPE(t_geometryObject),      INTENT(OUT) :: rgeomObject

!</output>

!</subroutine>

    ! The dimension is 2D.
    rgeomObject%ndimension = NDIM2D
    
    ! We want a rectangle.
    rgeomObject%ctype = geom_rect
    
    ! Store the coordinate system.
    rgeomObject%rcoord2D = rcoordSys
    
    ! Is our object inverted?
    IF (PRESENT(binverted)) THEN
      rgeomObject%binverted = binverted
    ELSE
      rgeomObject%binverted = .FALSE.
    END IF
    
    ! Store the edge lengths of the rectangle.
    rgeomObject%rrectangle%Dlength = Dlength(1:2)
    
    ! That's it!
  
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE geom_init_rectangle_direct(rgeomObject, Dlength, Dorigin, &
                                        drotation, dscalingFactor, binverted)

!<description>
  ! Creates a t_geometryObject representing a 2D rectangle.
!</description>

!<input>
  ! The edge lengths for the rectangle.
  REAL(DP), DIMENSION(:),            INTENT(IN)  :: Dlength
  
  ! OPTIONAL: The origin of the rectangle.
  ! Is set to (/ 0.0_DP, 0.0_DP /) if not given.
  REAL(DP), DIMENSION(:), OPTIONAL,  INTENT(IN)  :: Dorigin
  
  ! OPTIONAL: The rotation of the rectangle.
  ! Is set to 0.0_DP if not given.
  REAL(DP), OPTIONAL,                INTENT(IN)  :: drotation
  
  ! OPTIONAL: The scaling factor of the rectangle.
  ! Is set to 1.0_DP if not given.
  REAL(DP), OPTIONAL,                INTENT(IN)  :: dscalingFactor
  
  ! OPTIONAL: A boolean telling us whether the object is inverted.
  ! Is set to .FALSE. if not given.
  LOGICAL, OPTIONAL,                 INTENT(IN)  :: binverted

!</input>

!<output>
  ! A t_geometryObject structure to be written.
  TYPE(t_geometryObject),            INTENT(OUT) :: rgeomObject

!</output>

!</subroutine>

    ! The dimension is 2D.
    rgeomObject%ndimension = NDIM2D
    
    ! We want a rectangle.
    rgeomObject%ctype = geom_rect
    
    ! Now we need to create the coordinate system.
    CALL bgeom_initCoordSys2D (rgeomObject%rcoord2D, Dorigin, drotation, &
                               dscalingFactor)
    
    ! Is our object inverted?
    IF (PRESENT(binverted)) THEN
      rgeomObject%binverted = binverted
    ELSE
      rgeomObject%binverted = .FALSE.
    END IF
    
    ! Store the edge length of the rectangle.
    rgeomObject%rrectangle%Dlength = Dlength(1:2)
    
    ! That's it!
  
  END SUBROUTINE
  
  ! ***************************************************************************
      
!<subroutine>

  SUBROUTINE geom_rect_isInGeometry (rgeomObject, Dcoords, iisInObject)

!<description>
  ! This routine checks whether a given point is inside the rectangle or not.
  !
  ! iisInObject is set to 0 if the point is outside the rectangle, it is set
  ! to 1 if it is inside the rectangle and is set to -1 if the point is inside
  ! the rectangle and the rectangle is inverted.
!</description>

!<input>
  ! The rectangle against that the point is to be tested.
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
  
  ! The coordinates of the point that is to be tested.
  REAL(DP), DIMENSION(:), INTENT(IN)  :: Dcoords
  
!</input>

!<output>
  ! An integer for the return value.
  INTEGER(I32),           INTENT(OUT) :: iisInObject
!</output>

!</subroutine>

  ! And an array for the transformed X- and Y-coordinates of our point
  REAL(DP), DIMENSION(2) :: DrelCoords

    ! First transfrom the point's coordinates into the rectangle's local
    ! coordinate system
    CALL bgeom_transformBackPoint2D(rgeomObject%rcoord2D, Dcoords, DrelCoords)
    
    ! Get the absolute values of the coords
    DrelCoords(1) = ABS(DrelCoords(1))
    DrelCoords(2) = ABS(DrelCoords(2))
    
    ! Check against half of the edge lengths
    IF ((DrelCoords(1) .LE. (0.5_DP * rgeomObject%rrectangle%Dlength(1))) .AND. &
        (DrelCoords(2) .LE. (0.5_DP * rgeomObject%rrectangle%Dlength(2)))) THEN
      
      ! We are inside the rectangle
      iisInObject = 1
      
    ELSE
      ! We are outside the rectangle
      iisInObject = 0
      
    END IF

    ! Is the rectangle inverted?
    IF (rgeomObject%binverted) THEN
      iisInObject = 1 - iisInObject
    END IF

    ! That's it
    
  END SUBROUTINE

  ! ***************************************************************************
      
!<subroutine>
  
  SUBROUTINE geom_rect_prjToBoundary (rgeomObject, Dcoords, Dproj)

!<description>
  ! This routine calculates the projection of a point onto a 2D rectangle's
  ! boundary.
  !
  ! This routine is nearly identical to the geom_square_prjToBoundary
  ! routine.
!</description>

!<input>
  ! The geometry object to calculate the distance from.
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
  
  ! The coordinates of the point that is to be tested.
  ! The array must hold at least 2 entries for a 2D object, and at least 3
  ! entries for a 3D object.
  REAL(DP), DIMENSION(:), INTENT(IN)  :: Dcoords
  
!</input>

!<output>
  ! The coordinates of the boundary projection.
  ! The array must hold at least 2 entries for a 2D object, and at least 3
  ! entries for a 3D object.
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Dproj
  
!</output>

!</subroutine>

  ! Rectangle's half edge lengths
  REAL(DP), DIMENSION(2) :: Dlen
  
  ! Temporary coords in reference coordinate system
  REAL(DP), DIMENSION(2) :: DcoordsRef
  
  ! mirroring factors
  REAL(DP), DIMENSION(2) :: Dmirror = (/ 1.0_DP, 1.0_DP /)
  
    ! Get rectangle's half egde lengths
    Dlen = rgeomObject%rrectangle%Dlength * 0.5_DP
  
    ! First we need to transform the point's coordinates into the square's
    ! local coordinate system.
    CALL bgeom_transformBackPoint2D(rgeomObject%rcoord2D, Dcoords, DcoordsRef)
    
    ! Now mirror the point, so that the resulting point is inside the positive
    ! quadrant.
    IF (DcoordsRef(1) .LT. 0.0) THEN
      ! X-coord of point is negative, so multiply with -1
      Dmirror(1) = -1.0_DP
      DcoordsRef(1) = -DcoordsRef(1)
    END IF
    
    IF (DcoordsRef(2) .LT. 0.0) THEN
      ! Y-coord of point is negative, so multiply with -1
      Dmirror(2) = -1.0_DP
      DcoordsRef(2) = -DcoordsRef(2)
    END IF
    
    ! If both coordinates are greater than half of the correspoding
    ! rectangle's egde length, then the projection is the rectangle's corner.
    
    IF ((DcoordsRef(1) .GE. Dlen(1)) .AND. (DcoordsRef(2) .GE. Dlen(2))) THEN
    
      ! Save square's corner
      DcoordsRef = Dlen
    
    ELSE
      ! In all other cases, the projection is on an edge of the square.
      ! Now find out which coordinate is greater.
      IF (DcoordsRef(1) .GE. DcoordsRef(2)) THEN
        ! The X-coordinate is greater than the Y-coordinate.
        ! Now we need to set the X-coordinate to half of the rectangle's edge
        ! length and then we have our projection.
        DcoordsRef(1) = Dlen(1)
      
      ELSE
        ! Otherwise correct Y-coordinate
        DcoordsRef(2) = Dlen(2)
      
      END IF
    
    END IF
    
    ! The projection itself is calculated, now we need to mirror the projection
    ! into the quadrant where our original point was in.
    DcoordsRef(1) = DcoordsRef(1) * Dmirror(1)
    DcoordsRef(2) = DcoordsRef(2) * Dmirror(2)
    
    ! And transform the projection back into world coordinates
    CALL bgeom_transformPoint2D(rgeomObject%rcoord2D, DcoordsRef, Dproj)
    
    ! That's it

  END SUBROUTINE

  ! ***************************************************************************
      
!<subroutine>

  SUBROUTINE geom_rect_calcSignedDistance(rgeomObject, Dcoords, ddistance)
  
!<description>
  ! This routine calculates the shortest signed distance between a point and
  ! a rectangle's boundary.
!</description>

!<input>
  ! The rectangle against that the point is to be tested.
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
  
  ! The coordinates of the point that is to be tested.
  REAL(DP), DIMENSION(:), INTENT(IN)  :: Dcoords

!</input>

!<output>
  ! The shortest signed distance between the point and the rectangle's boundary.
  REAL(DP),               INTENT(OUT) :: ddistance
!</output>

!</subroutine>

  ! We need one local array for the transformed X- and Y-coordinates of our
  ! point.
  REAL(DP), DIMENSION(2) :: DrelCoords
  
    ! We are going to use the same trick as in the routine for the square.
  
    ! We are now going to transform our point into the rectangle's local
    ! coordinate system.
    CALL bgeom_transformBackPoint2D(rgeomObject%rcoord2D, Dcoords, DrelCoords)
  
    ! Now get the absolute values of the relative coords, and subtract one
    ! half of the rectangle's edge length
    DrelCoords(1) = ABS(DrelCoords(1)) - &
                    (rgeomObject%rrectangle%Dlength(1) * 0.5_DP)
    DrelCoords(2) = ABS(DrelCoords(2)) - &
                    (rgeomObject%rrectangle%Dlength(2) * 0.5_DP)
  
    ! Check whether the resulting coordinates are both positive
    IF ((DrelCoords(1) > 0.0_DP) .AND. (DrelCoords(2) > 0.0_DP)) THEN
      ddistance = SQRT(DrelCoords(1)**2 + DrelCoords(2)**2)
    ELSE
      ddistance = MAX(DrelCoords(1), DrelCoords(2))
    END IF
    
    ! Is the rectangle inverted?
    IF (rgeomObject%binverted) THEN
      ddistance = -ddistance
    END IF
    
    ! Make sure to scale the distance by the rectangle's coordinate system's
    ! scaling factor
    ddistance = ddistance * rgeomObject%rcoord2D%dscalingFactor
    
    ! That's it
    
  END SUBROUTINE
  
  ! ***************************************************************************
 
!<subroutine>
  
  SUBROUTINE geom_rect_polygonise (rgeomObject, hpolyHandle)
  
!<description>
  ! This routine converts a rectangle to a polygon, so that it can
  ! be printed to an output file via ucd_addPolygon (see ucd.f90 for more
  ! details).
!</description>

!<input>
  ! The geometry object to calculate the distance from.
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
  
!</input>

!<output>
  ! Handle to a 2D array holding the vertices of the polygon.
  INTEGER, INTENT(OUT) :: hpolyHandle
  
!</output>

!</subroutine>

  REAL(DP), DIMENSION(:,:), POINTER :: p_Dvertices
  REAL(DP) :: dedgeX, dedgeY
  
  INTEGER(I32), DIMENSION(2), PARAMETER :: Isize = (/ 2, 4 /)
  
    ! Get edge lengths
    dedgeX = rgeomObject%rrectangle%Dlength(1) * 0.5_DP
    dedgeY = rgeomObject%rrectangle%Dlength(2) * 0.5_DP
  
    ! Allocate desired number of vertices
    CALL storage_new2D('geom_rect_polygonise', 'hpolyHandle', Isize, &
                       ST_DOUBLE, hpolyHandle, ST_NEWBLOCK_NOINIT)

    ! Get vertice array
    CALL storage_getbase_double2D(hpolyHandle, p_Dvertices)
    
    ! Set coords
    p_Dvertices(1,1) = -dedgeX
    p_Dvertices(2,1) = dedgeY
    p_Dvertices(1,2) = -dedgeX
    p_Dvertices(2,2) = -dedgeY
    p_Dvertices(1,3) = dedgeX
    p_Dvertices(2,3) = -dedgeY
    p_Dvertices(1,4) = dedgeX
    p_Dvertices(2,4) = dedgeY
    
  END SUBROUTINE
  
  ! ***************************************************************************
  ! *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*
  ! *= 2D Polygon Routines                                                   =*
  ! *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*
  ! ***************************************************************************
    
!<subroutine>

  SUBROUTINE geom_init_polygon_indirect(rgeomObject, rcoordSys, Dvertices, &
                                        npolyType, binverted)

!<description>
  ! Creates a t_geometryObject representing a 2D polygon.
!</description>

!<input>
  ! A 2D coordinate system for the polygon.
  TYPE(t_coordinateSystem2D),  INTENT(IN)  :: rcoordSys
  
  ! An 2D array for the vertices of the polygon
  REAL(DP), DIMENSION(:,:), TARGET, INTENT(IN)     :: Dvertices
  
  ! OPTIONAL: One of the GEOM_POLYGON_XXXX constants.
  ! Is set to GEOM_POLYGON_GENERAL if not given.
  INTEGER, OPTIONAL, INTENT(IN)            :: npolyType
  
  ! OPTIONAL: A boolean telling us whether the object is inverted.
  ! Is set to .FALSE. if not given.
  LOGICAL, OPTIONAL,           INTENT(IN)  :: binverted

!</input>  

!<output>
  ! A t_geometryObject structure to be written.
  TYPE(t_geometryObject),      INTENT(OUT) :: rgeomObject

!</output>

!</subroutine>

    ! The dimension is 2D.
    rgeomObject%ndimension = NDIM2D
    
    ! We want a polygon.
    rgeomObject%ctype = GEOM_POLYGON
    
    ! Store the coordinate system.
    rgeomObject%rcoord2D = rcoordSys
    
    ! Is our object inverted?
    IF (PRESENT(binverted)) THEN
      rgeomObject%binverted = binverted
    ELSE
      rgeomObject%binverted = .FALSE.
    END IF
    
    ! Store the vertices of the polygon
    rgeomObject%rpolygon%p_Dvertices => Dvertices
    
    ! Store the polygon type
    IF (PRESENT(npolyType)) THEN
      rgeomObject%rpolygon%npolyType = npolyType
    ELSE
      rgeomObject%rpolygon%npolyType = GEOM_POLYGON_GENERAL
    END IF
    
    ! That's it!
    
  END SUBROUTINE

  ! ***************************************************************************
      
!<subroutine>

  SUBROUTINE geom_init_polygon_direct(rgeomObject, Dvertices, npolyType, &
                                      Dorigin, drotation, dscalingFactor, &
                                      binverted)

!<description>
  ! Creates a t_geometryObject representing a 2D polygon.
!</description>

!<input>
  ! An 2D array for the vertices of the polygon
  REAL(DP), DIMENSION(:,:), TARGET, INTENT(IN)     :: Dvertices
  
  ! OPTIONAL: One of the GEOM_POLYGON_XXXX constants.
  ! Is set to GEOM_POLYGON_GENERAL if not given.
  INTEGER, OPTIONAL, INTENT(IN)            :: npolyType
  
  ! OPTIONAL: The origin of the polygon.
  ! Is set to (/ 0.0_DP, 0.0_DP /) if not given.
  REAL(DP), DIMENSION(:), OPTIONAL,  INTENT(IN)  :: Dorigin
  
  ! OPTIONAL: The rotation of the polygon.
  ! Is set to 0.0_DP if not given.
  REAL(DP), OPTIONAL,                INTENT(IN)  :: drotation
  
  ! OPTIONAL: The scaling factor of the polygon.
  ! Is set to 1.0_DP if not given.
  REAL(DP), OPTIONAL,                INTENT(IN)  :: dscalingFactor
  
  ! OPTIONAL: A boolean telling us whether the object is inverted.
  ! Is set to .FALSE. if not given.
  LOGICAL, OPTIONAL,           INTENT(IN)  :: binverted

!</input>  

!<output>
  ! A t_geometryObject structure to be written.
  TYPE(t_geometryObject),      INTENT(OUT) :: rgeomObject

!</output>

!</subroutine>

    ! The dimension is 2D.
    rgeomObject%ndimension = NDIM2D
    
    ! We want a polygon.
    rgeomObject%ctype = GEOM_POLYGON
    
    ! Now we need to create the coordinate system.
    CALL bgeom_initCoordSys2D (rgeomObject%rcoord2D, Dorigin, drotation, &
                               dscalingFactor)
    
    ! Is our object inverted?
    IF (PRESENT(binverted)) THEN
      rgeomObject%binverted = binverted
    ELSE
      rgeomObject%binverted = .FALSE.
    END IF
    
    ! Store the vertices of the polygon
    rgeomObject%rpolygon%p_Dvertices => Dvertices
    
    ! Store the polygon type
    IF (PRESENT(npolyType)) THEN
      rgeomObject%rpolygon%npolyType = npolyType
    ELSE
      rgeomObject%rpolygon%npolyType = GEOM_POLYGON_GENERAL
    END IF
    
    ! That's it!
    
  END SUBROUTINE

  ! ***************************************************************************
      
!<subroutine>

  SUBROUTINE geom_polygon_isInGeometry (rgeomObject, Dcoords, iisInObject)

!<description>
  ! This routine checks whether a given point is inside the polygon or not.
  !
  ! This routine calls the geom_polygon_iIG_convex routine if the polygon
  ! type is set to GEOM_POLYGON_CONVEX, otherwise it calls the 
  ! geom_polygon_projector method.
  !
  ! iisInObject is set to 0 if the point is outside the polygon, it is set
  ! to 1 if it is inside the polygon and is set to -1 if the point is inside
  ! the polygon and the polygon is inverted.
!</description>

!<input>
  ! The polygon against that the point is to be tested.
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
  
  ! The coordinates of the point that is to be tested.
  REAL(DP), DIMENSION(:), INTENT(IN)  :: Dcoords
  
!</input>

!<output>
  ! An integer for the return value.
  INTEGER(I32),           INTENT(OUT) :: iisInObject
!</output>

!</subroutine>
  
  ! The output of the projector routine
  LOGICAL :: bisInObject
  
  ! 2 Dummys for the projector routine
  REAL(DP) :: ddummy1
  REAL(DP), DIMENSION(2) :: Ddummy2

    SELECT CASE (rgeomObject%rpolygon%npolyType)

    CASE (GEOM_POLYGON_CONVEX)
      ! Special Case: Convex Polygon
      CALL geom_polygon_iIG_convex(rgeomObject, Dcoords, iisInObject)

    CASE DEFAULT
      ! Call Projector
      CALL geom_polygon_projector(rgeomObject, Dcoords, Ddummy2, ddummy1, &
                                  bisInObject)
      ! Are we inside the polygon?
      IF (bisInObject) THEN
        iisInObject = 1
      ELSE
        iisInObject = 0
      END IF
      
      ! Maybe the polygon is inverted?
      IF (rgeomObject%binverted) THEN
        iisInObject = 1 - iisInObject
      END IF
      
    END SELECT
    
    ! That's it
    
  END SUBROUTINE

  ! ***************************************************************************
      
!<subroutine>

  SUBROUTINE geom_polygon_iIG_convex (rgeomObject, Dcoords, iisInObject)

!<description>
  ! This routine checks whether a given point is inside the convex (!) polygon
  ! or not.
  !
  ! iisInObject is set to 0 if the point is outside the polygon, it is set
  ! to 1 if it is inside the polygon and is set to -1 if the point is inside
  ! the polygon and the polygon is inverted.
!</description>

!<input>
  ! The polygon against that the point is to be tested.
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
  
  ! The coordinates of the point that is to be tested.
  REAL(DP), DIMENSION(:), INTENT(IN)  :: Dcoords
  
!</input>

!<output>
  ! An integer for the return value.
  INTEGER(I32),           INTENT(OUT) :: iisInObject
!</output>

!</subroutine>

  ! A pointer to the vertex array
  REAL(DP), DIMENSION(:,:), POINTER :: p_Dvertices

  ! The lower and upper bounds of our vertice array
  INTEGER:: lb, ub, i
  
  ! The point in the reference coordinate system and 2 temporary vectors:
  ! an edge and a ray
  REAL(DP), DIMENSION(2) :: DcoordsRef, Dedge, Dray
  
    ! Let's assume that the point is outside the polygon
    iisInObject = 0
  
    ! First we're going to transform the point into the local reference
    ! coordinate system
    CALL bgeom_transformBackPoint2D(rgeomObject%rcoord2D, Dcoords, DcoordsRef)
  
    ! Get our vertice array
    p_Dvertices => rgeomObject%rpolygon%p_Dvertices
  
    ! Get the bounds of our vector array
    lb = LBOUND(p_Dvertices, 2)
    ub = UBOUND(p_Dvertices, 2)
    
    ! Once again, we will use a fancy trick to find out whether a point is
    ! inside a convex polygon or not.
    ! The trick is simple:
    ! Since we know that the edges of our polygon are given in counter-
    ! clockwise direction, we will simply loop through all edges of the polygon
    ! and calculate the scalar product of the edge and a ray vector.
    ! The ray vector is the clockwise normal vector of the vector from the
    ! edge's start-vertice to the point that is to be checked (transformed into
    ! the polygon's local coordinate system, of course).
    ! If all the scalar products of all edges and their corresponding ray
    ! vectors are non-negative, then the point is inside the polygon.
    
    
    ! Loop through the first n-1 edges of our polygon
    DO i=lb, ub-1
    
      ! Calculate edge vector
      Dedge = p_Dvertices(1:2,i+1) - p_Dvertices(1:2,i)
      
      ! Calculate ray vector
      Dray = DcoordsRef - p_Dvertices(1:2, i)
      
      ! Calculate scalar product of ray vector's normal and the edge
      IF (((Dray(2) * Dedge(1)) - (Dray(1) * Dedge(2))) .LT. 0.0_DP) THEN
        
        ! The point is outside
        RETURN
        
      END IF
    
    END DO
        ! Check last edge
    Dedge = p_Dvertices(1:2,lb) - p_Dvertices(1:2,ub)
      
    ! Calculate ray vector
    Dray = DcoordsRef - p_Dvertices(1:2, ub)
    
    ! Calculate scalar product of ray vector's normal and the edge
    IF (((Dray(2) * Dedge(1)) - (Dray(1) * Dedge(2))) .GE. 0.0_DP) THEN
        
      ! All scalar products are non-negative - so the point is inside the poly
      iisInObject = 1
      
    ELSE
      iisInObject = 0;

    END IF
    
    IF (rgeomObject%binverted) THEN
      iisInObject = 1 - iisInObject
    END IF
  
    ! That's it
    
  END SUBROUTINE

  ! ***************************************************************************
      
!<subroutine>

  SUBROUTINE geom_polygon_projector (rgeomObject, Dcoords, Dproj, ddistance, &
                                     bisInside)

!<description>
  ! This routine calculates the projection of a point onto a 2D polygon's
  ! boundary.
!</description>

!<input>
  ! The geometry object to calculate the distance from.
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
  
  ! The coordinates of the point that is to be tested.
  ! The array must hold at least 2 entries for a 2D object, and at least 3
  ! entries for a 3D object.
  REAL(DP), DIMENSION(:), INTENT(IN)  :: Dcoords
  
!</input>

!<output>
  ! The coordinates of the boundary projection.
  ! The array must hold at least 2 entries for a 2D object, and at least 3
  ! entries for a 3D object.
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Dproj
  
  ! OPTIONAL: The distance between the given point and the projection.
  ! Note: This distance is absolute, not signed!
  REAL(DP), OPTIONAL, INTENT(OUT) :: ddistance
  
  ! OPTIONAL: A boolean deciding whether the point is inside or outside the
  ! polygon.
  LOGICAL, OPTIONAL, INTENT(OUT) :: bisInside
  
!</output>

!</subroutine>

  ! A pointer to the vertex array
  REAL(DP), DIMENSION(:,:), POINTER :: p_Dvertices

  ! The lower and upper bounds of our vertice array
  INTEGER :: lb, ub, i, iminVert, icurVert, iprev, inext
  LOGICAL :: bminVert, bcurVert
  
  ! The point in the reference coordinate system and 5 temporary vectors:
  ! two edges, two rays, and 2 vectors for the projection
  REAL(DP), DIMENSION(2) :: DcoordsRef, Dedge, Dray1, Dray2
  REAL(DP), DIMENSION(2) :: DprojMin, DprojCur
  
  ! The shortest squared distance from the point to the polygon's vertices
  REAL(DP) :: ddistMin, ddistCur
  REAL(DP) :: dalpha, dbeta, dgamma, ddelta
  LOGICAL :: binsideMin, binsideCur
  
    ! First we're going to transform the point into the local reference
    ! coordinate system
    CALL bgeom_transformBackPoint2D(rgeomObject%rcoord2D, Dcoords, DcoordsRef)
  
    ! Get our vertice array
    p_Dvertices => rgeomObject%rpolygon%p_Dvertices
  
    ! Get the bounds of our vector array
    lb = LBOUND(p_Dvertices, 2)
    ub = UBOUND(p_Dvertices, 2)
    
    ! Calculate the last edge
    Dedge = p_Dvertices(1:2, lb) - p_Dvertices(1:2, ub)
    
    ! Calculate both rays
    Dray1 = DcoordsRef - p_Dvertices(1:2, ub)
    Dray2 = DcoordsRef - p_Dvertices(1:2, lb)
    
    ! Calculate (squared) vector lengths
    dalpha = Dray1(1)**2 + Dray1(2)**2
    dbeta = Dray2(1)**2 + Dray2(2)**2
    dgamma = Dedge(1)**2 + Dedge(2)**2
    
    ! Calculate interpolation factor
    ddelta = (dalpha - dbeta + dgamma) / (2.0_DP * dgamma)
    
    ! clip interpolation factor to [0, 1]
    bminVert = .FALSE.
    IF (ddelta .LE. 0.0_DP) THEN
      ddelta = 0.0_DP
      bminVert = .TRUE.
      iminVert = ub
    ELSE IF (ddelta .GE. 1.0_DP) THEN
      ddelta = 1.0_DP
      bminVert = .TRUE.
      iminVert = lb
    END IF
    
    ! Assume that this is the minimal projection
    DprojMin = p_Dvertices(1:2, ub) + (ddelta * Dedge)
    
    ! abuse ray1 for the vector between the point and its projection
    Dray1 = DcoordsRef - DprojMin
    
    ! calculate distance
    ddistMin = Dray1(1)**2 + Dray1(2)**2
    
    ! Decide whether the point is inside or outside the polygon in
    ! respect to this edge
    binsideMin = (((Dray1(2) * Dedge(1)) - (Dray1(1) * Dedge(2))) &
                  .GE. 0.0_DP)
    
    ! Now loop though all other edges
    DO i = lb, ub-1
      
      ! Calculate i-th edge
      Dedge = p_Dvertices(1:2, i+1) - p_Dvertices(1:2, i)
      
      ! Calculate both rays
      Dray1 = Dray2
      Dray2 = DcoordsRef - p_Dvertices(1:2, i+1)
      
      ! Calculate (squared) vector lengths
      dalpha = Dray1(1)**2 + Dray1(2)**2
      dbeta = Dray2(1)**2 + Dray2(2)**2
      dgamma = Dedge(1)**2 + Dedge(2)**2
        
      ! Calculate interpolation factor
      ddelta = (dalpha - dbeta + dgamma) / (2.0_DP * dgamma)
        
      ! clip interpolation factor to [0, 1]
      bcurVert = .FALSE.
      icurVert = 0
      IF (ddelta .LE. 0.0_DP) THEN
        ddelta = 0.0_DP
        bcurVert = .TRUE.
        icurVert = i
      ELSE IF (ddelta .GE. 1.0_DP) THEN
        ddelta = 1.0_DP
        bcurVert = .TRUE.
        icurVert = i+1
      END IF
        
      ! Calculate current projection
      DprojCur = p_Dvertices(1:2, i) + (ddelta * Dedge)
        
      ! abuse ray1 for the vector between the point and its projection
      Dray1 = DcoordsRef - DprojCur
        
      ! calculate distance
      ddistCur = Dray1(1)**2 + Dray1(2)**2

      ! Decide whether the point is inside or outside the polygon in
      ! respect to this edge
      binsideCur = (((Dray1(2) * Dedge(1)) - (Dray1(1) * Dedge(2))) &
                    .GE. 0.0_DP)
      
      ! Maybe the distance is minimal?
      IF (ddistCur .LT. ddistMin) THEN
      
        ! Save the current projection as the minimal one
        DprojMin = DprojCur
        ddistMin = ddistCur
        binsideMin = binsideCur
        bminVert = bcurVert
        iminVert = icurVert
      
      END IF
      
    END DO

    ! Transform projection to world coordinates
    CALL bgeom_transformPoint2D(rgeomObject%rcoord2D, DprojMin, Dproj)
    
    IF (PRESENT(ddistance)) THEN
      ddistance = SQRT(ddistMin) * rgeomObject%rcoord2D%dscalingFactor
    END IF
    
    IF (PRESENT(bisInside)) THEN
    
      ! First we need to check if the projection is a vertice
      IF (bminVert) THEN
      
        ! We now need to find out whether the projection vertice is convex or
        ! (strictly) concav.
        IF (iminVert .LT. ub) THEN
          inext = iminVert + 1
        ELSE
          inext = lb
        END IF
        IF (iminVert .GT. lb) THEN
          iprev = iminVert - 1
        ELSE
          iprev = ub
        END IF 

        Dray1 = p_Dvertices(1:2, iminVert) - p_Dvertices(1:2, iprev)
        Dray2 = p_Dvertices(1:2, inext) - p_Dvertices(1:2, iminVert)
        
        ! Calculate the scalar product of both edges to find out whether the
        ! vertice is convex.
        dalpha = (Dray2(2) * Dray1(1)) - (Dray2(1) * Dray1(2))
        
        ! Calculate the scalar product of the ray and the 2 edges
        Dedge = DcoordsRef - p_Dvertices(1:2, iminVert)
        dbeta  = (Dedge(2) * Dray1(1)) - (Dedge(1) * Dray1(2))
        dgamma = (Dedge(2) * Dray2(1)) - (Dedge(1) * Dray2(2))
        
        IF (dalpha .GE. 0.0_DP) THEN
          ! The vertice is convex.
          ! The point is inside the polygon if and only if both beta and gamma
          ! are greater or equal to 0.
          bisInside = (dbeta .GE. 0.0_DP) .AND. (dgamma .GE. 0.0_DP)
        
        ELSE
          ! The vertice is stricly concav.
          ! The point is outside the polygon if and only if both beta and gamma
          ! are strictly less than 0.
          bisInside = (dbeta .GE. 0.0_DP) .OR. (dgamma .GE. 0.0_DP)
        
        END IF
      
      ELSE
        ! The projection is not a vertice.
        bisInside = binsideMin
      END IF
    END IF
  
    ! That's it
    
  END SUBROUTINE

  ! ***************************************************************************
      
!<subroutine>

  SUBROUTINE geom_polygon_calcSignedDistance (rgeomObject, Dcoords, ddistance)

!<description>
  ! This routine calculates the signed distance of a given point and a polygon.
  !
  ! This routine is not a wrapper - it calls geom_polygon_prjToBoundary to get
  ! a projection onto the boundary of the polygon and calculates the distance
  ! to the projection.
!</description>

!<input>
  ! The polygon against that the point is to be tested.
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
  
  ! The coordinates of the point that is to be tested.
  REAL(DP), DIMENSION(:), INTENT(IN)  :: Dcoords
  
!</input>

!<output>
  ! The shortest signed distance between the point and the circle's boundary.
  REAL(DP),               INTENT(OUT) :: ddistance
!</output>

!</subroutine>


  ! The projection
  REAL(DP), DIMENSION(2) :: Dproj
  
  ! Inside the polygon?
  LOGICAL :: bisInside
  
    ! Calculate projection
    CALL geom_polygon_projector(rgeomObject, Dcoords, Dproj, ddistance, bisInside)
    
    IF (bisInside .AND. (.NOT. rgeomObject%binverted)) THEN
      ddistance = -ddistance
    END IF
  
    ! That's it
    
  END SUBROUTINE
  
  ! ***************************************************************************
 
!<subroutine>
  
  SUBROUTINE geom_polygon_polygonise (rgeomObject, hpolyHandle)
  
!<description>
  ! This routine converts a polygon to a polygon, so that it can
  ! be printed to an output file via ucd_addPolygon (see ucd.f90 for more
  ! details).
  ! This routine simply copies the vertices stored in the geometry object
  ! into a new array of vertices without changing them.
!</description>

!<input>
  ! The geometry object to calculate the distance from.
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
  
!</input>

!<output>
  ! Handle to a 2D array holding the vertices of the polygon.
  INTEGER, INTENT(OUT) :: hpolyHandle
  
!</output>

!</subroutine>

  INTEGER :: i, inumVerts
  
  REAL(DP), DIMENSION(:,:), POINTER :: p_Dvertices
  
  INTEGER(I32), DIMENSION(2) :: Isize
  
    ! Get number of vertices
    inumVerts = UBOUND(rgeomObject%rpolygon%p_Dvertices, 2)
  
    ! Allocate desired number of vertices
    Isize = (/ 2, inumVerts /)
    CALL storage_new2D('geom_polygon_polygonise', 'hpolyHandle', Isize, &
                       ST_DOUBLE, hpolyHandle, ST_NEWBLOCK_NOINIT)

    ! Get vertice array
    CALL storage_getbase_double2D(hpolyHandle, p_Dvertices)
    
    ! Copy all vertices
    DO i=1, inumVerts

      p_Dvertices(1, i) = rgeomObject%rpolygon%p_Dvertices(1, i)
      p_Dvertices(2, i) = rgeomObject%rpolygon%p_Dvertices(2, i)
    
    END DO
    
    
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  PURE SUBROUTINE geom_moveto(rgeomObject, DnewOrigin)
  
!<description>
  ! This subroutine moves a given geometry object (2D or 3D) to a new origin
  ! by simply overwriting the coordinate system's old origin.
  !
  ! If the geometry object is a 2D object, this routine expects DnewOrigin to
  ! be an array with at least 2 entries, if it is a 3D object,  DnewOrigin has
  ! to have at least 3 entries.
!</description>

!<input>
  ! The new origin of the geometric object.
  REAL(DP), DIMENSION(:), INTENT(IN)    :: DnewOrigin
  
!</input>

!<inputoutput>
  ! The geometry object that is to be moved.
  TYPE(t_geometryObject), INTENT(INOUT) :: rgeomObject
  
!</inputoutput>

!</subroutine>

    ! We need to check whether we use a 2D or 3D coordinate system first.
    ! Then we simply overwrite the old origin.
    SELECT CASE (rgeomObject%ndimension)
    CASE (NDIM2D)
      rgeomObject%rcoord2D%Dorigin = DnewOrigin(1 : 2)
    !CASE (NDIM3D)
      !rgeomObject%rcoord3D%Dorigin = DnewOrigin(1 : 3)
    END SELECT
    
    ! That's it
    
  END SUBROUTINE

  ! ***************************************************************************
      
!<subroutine>

  ELEMENTAL SUBROUTINE geom_rotate2D(rgeomObject, dnewRotation, dscalingFactor)
  
!<description>
  ! This routine overwrites the rotation and scaling factor a 2D geometry
  ! object.
!</description>

!<input>
  ! A new rotation angle, range: 0..2*PI
  REAL(DP),               INTENT(IN)  :: dnewRotation
  
  ! OPTIONAL: A new scaling factor
  REAL(DP), OPTIONAL,     INTENT(IN)  :: dscalingFactor
  
!</input>

!<inputoutput>
  ! The geometry object that is to be rotated and/or scaled.
  TYPE(t_geometryObject), INTENT(INOUT) :: rgeomObject
  
!</inputoutput>

!</subroutine>

    ! *************************************************************************
    ! * We do not check the geometry object's dimension - we assume that the  *
    ! * caller calls this routine only for 2D objects.                        *
    ! *************************************************************************
    
    ! Overwrite rotation angle
    rgeomObject%rcoord2D%drotation = dnewRotation
    
    ! Recalculate SIN and COS values of angle
    rgeomObject%rcoord2D%dsin_rotation = SIN(dnewRotation)
    rgeomObject%rcoord2D%dcos_rotation = COS(dnewRotation)
    
    ! Overwrite scaling factor, if given.
    IF (PRESENT(dscalingFactor)) THEN
      rgeomObject%rcoord2D%dscalingFactor = dscalingFactor
    END IF
    
    ! That's it!
    
  END SUBROUTINE

  ! ***************************************************************************
      
!<subroutine>

  RECURSIVE SUBROUTINE geom_isInGeometry (rgeomObject, Dcoords, iisInObject)

!<description>
  ! This routine is a wrapper for the geom_****_isInGeometry routines.
  !
  ! The return value iisInObject holds the number of objects where the given
  ! point is inside the object's geometry.
  ! The returned integer is 0, if the point is outside of the object's
  ! geometry. It is 1 if the point is inside of a (simple) object's geometry,
  ! and may be > 1 if the object is a composed geometry object.
  ! For composed objects, iisInObject holds the number of subobjects, where
  ! the point is inside the subobject's geometry.
!</description>

!<input>
  ! The geometry object against that the point is to be tested.
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
  
  ! The coordinates of the point that is to be tested.
  ! The array must hold at least 2 entries for a 2D object, and at least 3
  ! entries for a 3D object.
  REAL(DP), DIMENSION(:), INTENT(IN)  :: Dcoords
  
!</input>

!<output>
  ! An integer storing the number of objects where the point is inside
  ! the object's geometry.
  INTEGER(I32),           INTENT(OUT) :: iisInObject
!</output>

!</subroutine>

    ! Check what type of object we have and call the object's corresponding
    ! subroutine.
    SELECT CASE(rgeomObject%ctype)
    CASE (GEOM_COMPOSED)
      CALL geom_composed_isInGeometry(rgeomObject, Dcoords, iisInObject)
    CASE (GEOM_CIRCLE)
      CALL geom_circle_isInGeometry(rgeomObject, Dcoords, iisInObject)
    CASE (GEOM_SQUARE)
      CALL geom_square_isInGeometry(rgeomObject, Dcoords, iisInObject)
    CASE (GEOM_ELLIPSE)
      CALL geom_ellipse_isInGeometry(rgeomObject, Dcoords, iisInObject)
    CASE (GEOM_RECT)
      CALL geom_rect_isInGeometry(rgeomObject, Dcoords, iisInObject)
    CASE (GEOM_POLYGON)
      CALL geom_polygon_isInGeometry(rgeomObject, Dcoords, iisInObject)
    CASE DEFAULT
      iisInObject = 0
    END SELECT
    
    ! That's it
    
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE geom_isInGeometryArray (rgeomObject, Dcoords, IisInObject)

!<description>
  ! This routine check whether an array of given points is inside a given
  ! geometry object or not.
  !
!</description>

!<input>
  ! The geometry object against that the points are to be tested.
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
  
  ! An array holding the coordinates of the points that are to be tested.
  REAL(DP), DIMENSION(:,:), INTENT(IN)  :: Dcoords
  
!</input>

!<output>
  ! An array of integers storing the number of objects where the point is inside
  ! the object's geometry.
  ! The lower and upper bounds of the array are assumed to be the same as the ones
  ! for the coordinate array.
  INTEGER(I32), DIMENSION(:), INTENT(OUT) :: IisInObject
!</output>

!</subroutine>

  INTEGER :: i, lb,ub

    ! Until now, this routine is a simple DO-loop
    lb = LBOUND(Dcoords, 2)
    ub = UBOUND(Dcoords, 2)

    DO i = lb, ub

      ! Call the geom_isInGeometry routine
      CALL geom_isInGeometry(rgeomObject, Dcoords(:,i), IisInObject(i))

    END DO
    
    ! That's it
    
  END SUBROUTINE
  
  ! ***************************************************************************
      
  !<subroutine>
  
  SUBROUTINE geom_projectToBoundary (rgeomObject, Dcoords, Dproj)

!<description>
  ! This routine is a wrapper for the geom_****_prjToBoundary routines.
!</description>

!<input>
  ! The geometry object to calculate the distance from.
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
  
  ! The coordinates of the point that is to be tested.
  ! The array must hold at least 2 entries for a 2D object, and at least 3
  ! entries for a 3D object.
  REAL(DP), DIMENSION(:), INTENT(IN)  :: Dcoords
  
!</input>

!<output>
  ! The coordinates of the boundary projection.
  ! The array must hold at least 2 entries for a 2D object, and at least 3
  ! entries for a 3D object.
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Dproj
  
!</output>

!</subroutine>

   ! Check what type of object we have and call the object's corresponding
   ! subroutine.
   SELECT CASE(rgeomObject%ctype)
   !CASE (GEOM_COMPOSED)
     !CALL geom_composed_prjToBoundary(rgeomObject, Dcoords, Dproj)
   CASE (GEOM_CIRCLE)
     CALL geom_circle_prjToBoundary(rgeomObject, Dcoords, Dproj)
   CASE (GEOM_SQUARE)
     CALL geom_square_prjToBoundary(rgeomObject, Dcoords, Dproj)
   CASE (GEOM_ELLIPSE)
     CALL geom_ellipse_prjToBoundary(rgeomObject, Dcoords, Dproj)
   CASE (GEOM_RECT)
     CALL geom_rect_prjToBoundary(rgeomObject, Dcoords, Dproj)
   CASE (GEOM_POLYGON)
     CALL geom_polygon_projector(rgeomObject, Dcoords, Dproj)
   END SELECT
   
   ! That's it!
   
 END SUBROUTINE
 
  ! ***************************************************************************
      
!<subroutine>
  
  SUBROUTINE geom_calcSignedDistance (rgeomObject, Dcoords, ddistance)

!<description>
  ! This routine is a wrapper for the geom_****_calcSignedDistance routines.
  !
  ! The return value ddistance holds the absolute shortest distance from a
  ! given point to the geometry object.
  ! If ddistance is < 0, then the given point is inside the object, and
  ! ddistance is > 0, if the given point if outside the object.
!</description>

!<input>
  ! The geometry object to calculate the distance from.
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
  
  ! The coordinates of the point that is to be tested.
  ! The array must hold at least 2 entries for a 2D object, and at least 3
  ! entries for a 3D object.
  REAL(DP), DIMENSION(:), INTENT(IN)  :: Dcoords
  
!</input>

!<output>
  ! The calculated signed distance
  REAL(DP),               INTENT(OUT) :: ddistance
  
!</output>

!</subroutine>

   ! Check what type of object we have and call the object's corresponding
   ! subroutine.
   SELECT CASE(rgeomObject%ctype)
   !CASE (GEOM_COMPOSED)
     !CALL geom_composed_calcSignedDistance(rgeomObject, Dcoords, ddistance)
   CASE (GEOM_CIRCLE)
     CALL geom_circle_calcSignedDistance(rgeomObject, Dcoords, ddistance)
   CASE (GEOM_SQUARE)
     CALL geom_square_calcSignedDistance(rgeomObject, Dcoords, ddistance)
   CASE (GEOM_ELLIPSE)
     CALL geom_ellipse_calcSignedDistance(rgeomObject, Dcoords, ddistance)
   CASE (GEOM_RECT)
     CALL geom_rect_calcSignedDistance(rgeomObject, Dcoords, ddistance)
   CASE (GEOM_POLYGON)
     CALL geom_polygon_calcSignedDistance(rgeomObject, Dcoords, ddistance)
   CASE DEFAULT
     ddistance = 0.0_DP
   END SELECT
   
   ! That's it!
   
 END SUBROUTINE

  ! ***************************************************************************
 
!<subroutine>
  
  SUBROUTINE geom_calcSignedDistanceArray (rgeomObject, Dcoords, Ddistance)

!<description>
  
!</description>

!<input>
  ! The geometry object to calculate the distance from.
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
  
  ! An array holding the coordinates of the points that are to be tested.
  REAL(DP), DIMENSION(:,:), INTENT(IN)  :: Dcoords
  
!</input>

!<output>
  ! An array holding the calculated signed distances.
  ! The lower and upper bounds of the array are assumed to be the same as the ones
  ! for the coordinate array.
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Ddistance
  
!</output>

!</subroutine>

  INTEGER :: i, lb,ub

    ! Until now, this routine is a simple DO-loop
    lb = LBOUND(Dcoords, 2)
    ub = UBOUND(Dcoords, 2)

    DO i = lb, ub

      ! Call the geom_isInGeometry routine
      CALL geom_calcSignedDistance(rgeomObject, Dcoords(:,i), Ddistance(i))

    END DO
    
    ! That's it
   
 END SUBROUTINE

  ! ***************************************************************************
 
!<subroutine>
  
  SUBROUTINE geom_polygonise (rgeomObject, hpolyHandle, bconvertToWorld, &
                              ndesiredVerticeCount)
  
!<description>
  ! This routine converts a 2D geometry object to a polygon, so that it can
  ! be printed to an output file via ucd_addPolygon (see ucd.f90 for more
  ! details).
  !
  ! This routine does not work for composed geometry objects and 3D objects,
  ! and therefore it will set hpolyHandle to ST_NOHANDLE.
  !
  ! For all valid geometry objects this routine returns a handle to a 2D
  ! double array holding the vertices of the polygon in world coordinates.
  ! The caller of the routine is responsible for releasing the allocated
  ! storage memory after using it.
  !
  ! Calling this routine for a polygon type geometry object with
  ! bconvertToWorld = .FALSE. is quite senseless as the routine just copies
  ! the polygons vertices into a new array.
  !
  ! This routine is (more or less) just a wrapper for the geom_XXXX_polygonise
  ! routines.
!</description>

!<input>
  ! The geometry object to calculate the distance from.
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
  
  ! OPTIONAL: Decides whether the coordinates of the polygon should be world
  ! coordinates or coordinates relative to the geometry object's local
  ! coordinate system. If not given, the coordinates are given in world
  ! coordinates.
  LOGICAL, OPTIONAL, INTENT(IN) :: bconvertToWorld
  
  ! OPTIONAL: The desired number of vertices for the generated polygon.
  ! Is only used for circles and ellipses, and is ignored for all other
  ! geometry objects.
  ! If not given, 64 vertices are generated.
  INTEGER, OPTIONAL, INTENT(IN) :: ndesiredVerticeCount
  
!</input>

!<output>
  ! Handle to a 2D array holding the vertices of the polygon.
  INTEGER, INTENT(OUT) :: hpolyHandle
  
!</output>

!</subroutine>

  ! Desired vertice count
  INTEGER :: ndVC, i, ub, lb
  
  ! The polygon's vertices
  REAL(DP), DIMENSION(:,:), POINTER :: p_Dvertices
  
  ! A temporary vertice
  REAL(DP), DIMENSION(1:2) :: Dtemp
  
    ! Check the optional parameter
    IF (PRESENT(ndesiredVerticeCount)) THEN
      IF (ndesiredVerticeCount .GE. 3) THEN
        ndVC = ndesiredVerticeCount
      ELSE
        ndVC = 64
      END IF
    ELSE
      ndVC = 64
    END IF
    
    ! Set handle to ST_NOHANDLE
    hpolyHandle = ST_NOHANDLE
    
    ! Call the corresponding subroutine of the geometry object
    SELECT CASE(rgeomObject%ctype)
    CASE (GEOM_CIRCLE)
      CALL geom_circle_polygonise(rgeomObject, hpolyHandle, ndVC)
    CASE (GEOM_ELLIPSE)
      CALL geom_ellipse_polygonise(rgeomObject, hpolyHandle, ndVC)
    CASE (GEOM_SQUARE)
      CALL geom_square_polygonise(rgeomObject, hpolyHandle)
    CASE (GEOM_RECT)
      CALL geom_rect_polygonise(rgeomObject, hpolyHandle)
    CASE (GEOM_POLYGON)
      CALL geom_polygon_polygonise(rgeomObject, hpolyHandle)
    END SELECT
    
    ! Maybe the subroutine failed?
    IF (hpolyHandle .EQ. ST_NOHANDLE) THEN
      RETURN
    END IF
    
    ! Don't we need to convert them to world coordinates?
    IF (PRESENT(bconvertToWorld)) THEN
      IF (.NOT. bconvertToWorld) RETURN
    END IF
    
    ! Get the vertices
    CALL storage_getbase_double2D(hpolyHandle, p_Dvertices)
    
    ! Get bounds
    ub = UBOUND(p_Dvertices, 2)
    lb = LBOUND(p_Dvertices, 2)
    
    ! Go through all of them
    DO i=lb, ub

      ! Call transform routine
      CALL bgeom_transformPoint2D(rgeomObject%rcoord2D, p_Dvertices(1:2, i), &
                                  Dtemp)
      ! Copy back to buffer
      p_Dvertices(1:2,i) = Dtemp

    END DO
    
    ! That's it!
  
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  RECURSIVE SUBROUTINE geom_done(rgeomObject)

!<description>
  ! This routine releases allocated buffer from a geometry object
  
!</description>

!<inputoutput>
  ! The geometry object that is to be destroyed
  TYPE(t_geometryObject), INTENT(INOUT) :: rgeomObject

!</inputoutput>

!</subroutine>

  INTEGER :: i

    ! If the object is composed, we have some work to do...
    IF(rgeomObject%ctype .EQ. GEOM_COMPOSED) THEN
    
      ! Go through all sub-objects and destroy them
      DO i=1, rgeomObject%rcomposed%nsubObjects
      
        ! Destroy sub-object
        CALL geom_done(rgeomObject%rcomposed%p_RsubObjects(i))
      
      END DO
      
      ! Deallocate object array
      DEALLOCATE(rgeomObject%rcomposed%p_RsubObjects)
    
    END IF

    ! Reset the object's type
    rgeomObject%ctype = GEOM_NONE
    
    ! That's it
  
  END SUBROUTINE
  
END MODULE
