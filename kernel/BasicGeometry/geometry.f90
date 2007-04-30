!##############################################################################
!# ****************************************************************************
!# <name> geometry </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module realises some basic 2D/3D geometry objects...
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
!#  5.) geom_moveto
!#      -> Moves the origin of a 2D/3D object to a specified point.
!#
!#  6.) geom_rotate2D
!#      -> Overwrites the rotation and scaling factor of a 2D object.
!#
!#  7.) geom_isInGeometry
!#      -> Check whether a given point is inside a geometry object.
!#
!#  8.) geom_projectToBoundary
!#      -> Projects a point onto the boundary of a geometry object.
!#
!#  9.) geom_calcSignedDistance
!#      -> Calculates the shortest signed distance between a given point
!#         and a geometry object.
!#
!# </purpose>
!##############################################################################
 
MODULE geometry

  USE basicgeometry

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
  INTEGER, PARAMETER :: GEOM_RECTANGLE = 4

!</constantblock>

! *****************************************************************************
! *****************************************************************************
! *****************************************************************************

!<types>
!<typeblock>

  ! This structure realises the subnode for the composed geometry object.
  
  TYPE t_geometryComposed
    
    ! Number of sub-objects in the composed object.
    INTEGER :: nsubobjects
    
    ! An array of sub-objects with are composed.
    TYPE(t_geometryObject), DIMENSION(:), POINTER :: rsubobjects
    
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
  
CONTAINS

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
      ! Maybe it's inverted?
      IF (rgeomObject%binverted) THEN
        iisInObject = -1
      ELSE
        iisInObject = 1
      END IF
      
    ELSE
      ! We are not inside the circle
      iisInObject = 0
    END IF
        
    ! That's it
    
  END SUBROUTINE

  ! ***************************************************************************
  
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
      
      ! Is the square inverted?
      IF (rgeomObject%binverted) THEN
        iisInObject = -1
      ELSE
        iisInObject = 1
      END IF
      
    ELSE
      ! We are outside the square
      iisInObject = 0
      
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
      ! Is the ellipse inverted?
      IF (rgeomObject%binverted) THEN
        iisInObject = -1
      ELSE
        iisInObject = 1
      END IF
      
    ELSE
      ! We are outside the ellipse
      iisInObject = 0
      
    END IF

    ! That's it
    
  END SUBROUTINE

  ! ***************************************************************************
      
!<subroutine>
  
  SUBROUTINE geom_ellipse_prjToBoundary (rgeomObject, Dcoords, Dproj)

!<description>
  ! This routine calculates the projection of a point onto an ellipse's
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
  
!</output>

!</subroutine>

  ! Two variables for the ellipses radiuses
  REAL(DP) :: dradX, dradY
  
  ! Some other temporary variables
  REAL(DP) :: dT, dF, dFDer, dXDivA, dYDivB, dradXSqr, dradYSqr, dprX, dprY
  REAL(DP) :: dratio, dXDivASqr, dYDivBSqr, dASqr, dBSqr
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
        
        dradXSqr = dradX**2
        IF (DcoordsRef(1) .LT. dradX - (dradYSqr / dradX)) THEN
        
          dradYSqr = dradY**2
          DcoordsRef(1) = dradXSqr * DcoordsRef(1) / (dradXSqr - dradYSqr)
          dXDivA = DcoordsRef(1) / dradX
          DcoordsRef(2) = SQRT(ABS(1.0 - dXDivA**2))
          
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
  REAL(DP) :: dreldist
  
  ! And an array for the transformed X- and Y-coordinates of our point
  REAL(DP), DIMENSION(2) :: DrelCoords
  
  ! And one array for the projection of our point onto the ellipses boundary
  REAL(DP), DIMENSION(2) :: Dprojection

    ! First transfrom the point's coordinates into the ellipse's local
    ! coordinate system
    CALL bgeom_transformBackPoint2D(rgeomObject%rcoord2D, Dcoords, DrelCoords)
    
    ! Calculate the distance to the ellipse's middle point
    dreldist = SQRT(DrelCoords(1)**2 + DrelCoords(2)**2)
    
    ! SPECIAL CASE: The point we are measuring distance is the middle point
    !               of our ellipse.
    IF (dreldist .EQ. 0.0_DP) THEN
      ! In this case we're inside the ellipse and the shortest distance
      ! to the boundary is the minumum of the X- and Y-radius
      ddistance = -MIN(rgeomObject%rellipse%Dradius(1), &
                       rgeomObject%rellipse%Dradius(1)) * &
                       rgeomObject%rcoord2D%dscalingFactor
      
      ! Maybe the ellipse is inverted?
      IF (rgeomObject%binverted) THEN
        ddistance = -ddistance
      END IF
      
      ! That's it for this special case
      RETURN
      
    END IF
    
    ! The point we are checking is not the middle point of our ellipse.
    ! That means that it has got an unique projection onto the ellipse's
    ! boundary. 
        
    ! Normalize the transformed point and scale it by the radiuses of our 
    ! ellipse - here's our projection.
    Dprojection(1) = (DrelCoords(1) * rgeomObject%rellipse%Dradius(1)) &
                     / dreldist
    Dprojection(2) = (DrelCoords(2) * rgeomObject%rellipse%Dradius(2)) &
                     / dreldist
    
    ! Now subtract the projection from our relative point
    DrelCoords(1) = DrelCoords(1) - Dprojection(1)
    DrelCoords(2) = DrelCoords(2) - Dprojection(2)
    
    ! Caluclate the new length
    ddistance = SQRT(DrelCoords(1)**2 + DrelCoords(2)**2)
    
    ! And scale by the ellipse's coordinate system's scaling factor
    ddistance = ddistance * rgeomObject%rcoord2D%dscalingFactor

    ! Now we still need to know whether we are inside the ellipse or not.
    ! This information is implicitly given by dreldist
    IF ((dreldist .LE. 1.0_DP) .XOR. (rgeomObject%binverted)) THEN
      ddistance = -ddistance
    END IF

    ! That's it
    
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
    rgeomObject%ctype = GEOM_RECTANGLE
    
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
    rgeomObject%ctype = GEOM_RECTANGLE
    
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

  SUBROUTINE geom_rectangle_isInGeometry (rgeomObject, Dcoords, iisInObject)

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

  ! We need one local variable for distance calculation
  REAL(DP) :: ddistance
  
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
      ! Is the rectangle inverted?
      IF (rgeomObject%binverted) THEN
        iisInObject = -1
      ELSE
        iisInObject = 1
      END IF
      
    ELSE
      ! We are outside the rectangle
      iisInObject = 0
      
    END IF

    ! That's it
    
  END SUBROUTINE

  ! ***************************************************************************
      
!<subroutine>
  
  SUBROUTINE geom_rectangle_prjToBoundary (rgeomObject, Dcoords, Dproj)

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

  SUBROUTINE geom_rectangle_calcSignedDistance(rgeomObject, Dcoords, ddistance)
  
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

!<output>
  ! The geometry object that is to be rotated and/or scaled.
  TYPE(t_geometryObject), INTENT(OUT) :: rgeomObject
  
!</output>

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

  SUBROUTINE geom_isInGeometry (rgeomObject, Dcoords, iisInObject)

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
    !CASE (GEOM_COMPOSED)
      !CALL geom_composed_isInGeometry(rgeomObject, Dcoords, iisInObject)
    CASE (GEOM_CIRCLE)
      CALL geom_circle_isInGeometry(rgeomObject, Dcoords, iisInObject)
    CASE (GEOM_SQUARE)
      CALL geom_square_isInGeometry(rgeomObject, Dcoords, iisInObject)
    CASE (GEOM_ELLIPSE)
      CALL geom_ellipse_isInGeometry(rgeomObject, Dcoords, iisInObject)
    CASE (GEOM_RECTANGLE)
      CALL geom_rectangle_isInGeometry(rgeomObject, Dcoords, iisInObject)
    CASE DEFAULT
      iisInObject = 0
    END SELECT
    
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
   CASE (GEOM_RECTANGLE)
     CALL geom_rectangle_prjToBoundary(rgeomObject, Dcoords, Dproj)
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
   CASE (GEOM_RECTANGLE)
     CALL geom_rectangle_calcSignedDistance(rgeomObject, Dcoords, ddistance)
   CASE DEFAULT
     ddistance = 0.0_DP
   END SELECT
   
   ! That's it!
   
 END SUBROUTINE
 
END MODULE
