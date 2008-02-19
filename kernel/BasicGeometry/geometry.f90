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
!#  6.) geom_init_composed
!#      -> Creates a composed object.
!#
!#  7.) geom_composed_addNode
!#      -> Adds an existing geometry object into a composed one.
!#
!#  8.) geom_composed_updateBndApprox
!#      -> Updates the boundary approximation for a composed object.
!#
!#  9.) geom_done
!#      -> Releases a geometry object (including its subobjects).
!#
!# 10.) geom_moveto
!#      -> Moves the origin of a 2D/3D object to a specified point.
!#
!# 11.) geom_rotate2D
!#      -> Overwrites the rotation and scaling factor of a 2D object.
!#
!# 12.) geom_isInGeometry
!#      -> Checks whether a given point is inside a geometry object.
!#
!# 13.) geom_isInGeometryArray
!#      -> Checks whether an array of given points is inside a geometry object.
!#
!# 14.) geom_projectToBoundary
!#      -> Projects a point onto the boundary of a geometry object.
!#
!# 15.) geom_calcSignedDistance
!#      -> Calculates the shortest signed distance between a given point
!#         and a geometry object.
!#
!# 16.) geom_calcSignedDistanceArray
!#      -> Calculates the shortest signed distance between an array of given
!#         points and a geometry object.
!#
!# 17.) geom_polygonise
!#      -> Converts a (non-composed) geometry object to a polygon and
!#         optionally converts the vertices into world coordinates.
!#
!# Remark:
!# This module contains a lot more routines, which however are not meant to
!# be called by the user, e.g. routines which are called by wrapper routines
!# listed above or routines which are needed for the boundary approximation
!# of a composed object.
!#
!# 
!# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!# Composed Geometry Objects
!# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!# I. General information
!# ----------------------
!# Composed geometry objects are build up as a tree structure here. Every inner
!# node of the tree is a quantor (union or intersection) and every leaf of the
!# tree is an analytic geometry object.
!#
!# Let's take a look at a simple example:
!#
!#        +------+    +------+
!#        |      |    |      |
!#        |   +--+    +--+   |
!#        |   |          |   |
!#        |   |          |   |
!#        |   |          |   |
!#        |   |          |   |
!#        |   +----------+   |
!#        |                  |
!#        +------------------+
!#
!# The figure above is made of 3 analytic objects:
!#                                                        R
!#                S_1                                   +----+
!#        +------------------+                          |    |
!#        |                  |           S_2            |    |
!#        |                  |       +----------+       |    |
!#        |                  |       |          |       +----+
!#        |                  |       |          |
!#        |                  |       |          |
!#        |                  |       |          |
!#        |                  |       +----------+
!#        |                  |
!#        +------------------+
!#
!# One possible tree for the figure above would look like this:
!#
!#      INTERSECTION
!#           |
!#     +-----+------------+
!#     |                  |
!#    S_1        COMPLEMENT of UNION
!#                        |
!#                  +-----+-----+
!#                  |           |
!#                 S_2          R
!#
!#
!#
!# II. Creating composed geometry objects
!# --------------------------------------
!# To create a geometry object, you need to build up the tree in bottom-up
!# order, i.e. you first create the analytic geometry objects and connect them
!# using composed objects using the geom_composed_addNode routine.
!#
!# Here's the code for the creation of the tree above:
!#
!# <code>
!#  1. TYPE(t_geometryObject) :: S, R, SR, S_SR
!#  2.
!#  3. ! create a small square with edge length of 0.4 around (0.5, 0.5)
!#  4. CALL geom_init_square(S, 0.4_DP, (/0.5_DP, 0.5_DP/))
!#  5. 
!#  6. ! create a rectangle with edge lengths of (0.2, 0.4) around (0.5, 0.775)
!#  7. CALL geom_init_rectangle(R, (/0.2_DP, 0.4_DP/), (/0.5_DP, 0.775_DP/))
!#  8. 
!#  9. ! create the COMPLEMENT of a UNION-node with 2 children
!# 10. CALL geom_init_composed(SR, 2, GEOM_COMP_TYPE_OR, .TRUE.)
!# 11. 
!# 12. ! add S and R to the UNION-node
!# 13. CALL geom_composed_addNode(SR, S, 1)
!# 14. CALL geom_composed_addNode(SR, R, 2)
!# 15. 
!# 16. ! create a big square with edge length of 0.7 around (0.5, 0.5)
!# 17. CALL geom_init_square(S, 0.7_DP, (/0.5_DP, 0.5_DP/))
!# 18.
!# 19. ! create an INTERSECTION-node with 2 children
!# 20. CALL geom_init_composed(S_SR, 2, GEOM_COMP_TYPE_AND)
!# 21.  
!# 22. ! add S and SR to the INTERSECTION-node
!# 23. CALL geom_composed_addNode(S_SR,  S, 1)
!# 24. CALL geom_composed_addNode(S_SR, SR, 2)
!# 25.
!# 26. ! do something with our composed object...
!# 27. CALL foobar(S_SR)
!# 28.
!# 29. ! IMPORTANT: destroy the composed object afterwards
!# 30. CALL geom_done(S_SR)
!#</code>
!#
!# There is something important to say about the geom_composed_addNode routine:
!# After calling geom_composed_addNode, the object that is passed as the
!# sub-object is set undefined.
!# In the code above, in line(s)
!#  1 -  3  S is undefined
!#       4  S is defined as the small square
!#  5 - 12  S represents the small square
!#      13  S is passed to geom_composed_addNode
!# 14 - 16  S is undefined
!#      17  S is defined as the big square
!# 18 - 22  S represents the big square
!#      23  S is passed to geom_composed_addNode
!# 24 - 30  S is undefined
!#
!#
!# III. Some hints about adding objects into the tree
!# --------------------------------------------------
!# If you plan to call the inside-outside-test or distance calculation routines
!# on a composed object, it is recommended to add the sub-objects for which the
!# inside-outside-test is cheap to calculate before the sub-objects for which
!# it is expensive to calculate.
!#
!# Here's a small list for all analytic geometry objects and a rating of how
!# expensive the inside-outside-test routines are (++ = cheap, -- = expensive): 
!# Circle..........: ++
!# Ellipse.........:  +
!# Square..........:  +
!# Rectangle.......:  +
!# Convex Polygon..:  -
!# General Polygon.: --
!# Composed Object.: worst of its sub-objects
!#
!#
!# IV. Boundary approximation
!# --------------------------
!# Unlike other analytic geometry objects (e.g. circles, squares, etc.)
!# the boundary of a composed object is described by a set of points, which are
!# an approximation of the composed object's boundary.
!# These points are used for the distance calculation and boundary projection
!# of a composed object.
!# The boundary approximation is calculated at the first call of a distance
!# calculation or boundary projection routine, or optionally if the user calls
!# the geom_composed_updateBndApprox routine manually.
!# The geom_composed_updateBndApprox routine MUST be called, if there have been
!# made changes in the composed object tree (except for rotation / scaling /
!# translation of the root object) AFTER a boundary approximation has been
!# calculated / updated.
!# When calling the geom_composed_updateBndApprox routine, you can pass a
!# tolerance parameter to the routine which controls the quality of the
!# boundary approximation. If no custom tolerance parameter is passed, 0.01 is
!# used as a default parameter. When using a custom tolerance, it is
!# recommended to set the tolerance parameter between 10e-5 and 10e-2.
!# Choosing a small tolerance will result in a better approximation of the
!# composed object's boundary, however, it will also increase the number of
!# points used for the approximation which results in higher memory cost and
!# longer calculation time for distance calculation and boundary projection.
!#
!# Remark:
!# The inside-outside-test does not use the boundary approximation.
!#
!# Remark:
!# The boundary approximation is automatically destroyed when the root object
!# is destroyed.
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

  ! No composed object
  INTEGER, PARAMETER :: GEOM_COMP_TYPE_NONE = 0

  ! Union of objects
  INTEGER, PARAMETER :: GEOM_COMP_TYPE_OR   = 1
  
  ! Intersection of objects
  INTEGER, PARAMETER :: GEOM_COMP_TYPE_AND  = 2
  
!</constantblock>

!</constants>

! *****************************************************************************
! *****************************************************************************
! *****************************************************************************

!<types>
!<typeblock>

  ! This stucture realises a boundary approximation, which is needed for
  ! distance calculation and boundary projection of composed geometry objects.
  TYPE t_boundaryApprox
  
    ! Approximation tolerance
    REAL(DP) :: dtolerance = 0.0_DP
    
    ! Number of allocated vertices
    INTEGER :: nvertsAllocated = 0
  
    ! Number of vertices used for the boundary approximation
    INTEGER :: nverts = 0
    
    ! A handle for the vertice vector
    INTEGER :: h_Dverts = ST_NOHANDLE
    
    ! An array for the vertice vector
    REAL(DP), DIMENSION(:,:), POINTER :: p_Dverts => NULL()
    
  END TYPE
  
!</typeblock>

! *****************************************************************************

!<typeblock>

  ! This structure realises the subnode for the composed geometry object.
  TYPE t_geometryComposed
  
    ! Type of composed object. One of the GEOM_COMP_TYPE_XXXX constants.
    INTEGER :: ccomposedType = GEOM_COMP_TYPE_NONE
    
    ! Number of sub-objects in the composed object.
    INTEGER :: nsubObjects = 0
    
    ! An array of sub-objects with are composed.
    TYPE(t_geometryObject), DIMENSION(:), POINTER :: p_RsubObjects => NULL()
    
    ! A pointer to a boundary-approximation structure
    TYPE(t_boundaryApprox), POINTER :: p_rboundaryApprox => NULL()
    
  END TYPE
  
!</typeblock>
  
! *****************************************************************************

!<typeblock>

  ! This structure realises the subnode for the circle object.
  TYPE t_geometryCircle
    
    ! Radius of the circle. Must be positive.
    REAL(DP) :: dradius = 0.0_DP
    
  END TYPE
  
!</typeblock>
  
! *****************************************************************************

!<typeblock>

  ! This structure realises the subnode for the square object.
  TYPE t_geometrySquare
  
    ! Length of each edge of the square. Must be positive.
    REAL(DP) :: dlength = 0.0_DP
    
  END TYPE

!</typeblock>

! *****************************************************************************

!<typeblock>

  ! This structure realises the subnode for the ellipse object.
  TYPE t_geometryEllipse
  
    ! X- and Y-radius of the ellipse. Must be positive.
    REAL(DP), DIMENSION(2) :: Dradii = (/ 0.0_DP, 0.0_DP /)
    
  END TYPE

!</typeblock>

! *****************************************************************************

!<typeblock>

  ! This structure realises the subnode for the rectangle object.
  TYPE t_geometryRectangle
  
    ! Lengths of the horizontal and vertical edges of the rectangle. Must be
    ! positive.
    REAL(DP), DIMENSION(2) :: Dlength = (/ 0.0_DP, 0.0_DP /)
    
  END TYPE

!</typeblock>

! *****************************************************************************

!<typeblock>

  ! This structure realises the subnode for the polygon object.
  TYPE t_geometryPolygon

    ! The polygon's type. One of the GEOM_POLYGON_XXXX constants defined above.
    INTEGER :: npolyType = GEOM_POLYGON_GENERAL

    ! A pointer to the polygon's vertices
    REAL(DP), DIMENSION(:,:), POINTER :: p_Dvertices => NULL()

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
    LOGICAL                    :: binverted = .FALSE.
    
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
  
  ! Allocate boundary approximation structure
  ALLOCATE(rgeomObject%rcomposed%p_rboundaryApprox)
  
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
  ! The object against that the point is to be tested.
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
  
  ! some other temporary variables
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

!<subroutine>

  SUBROUTINE geom_composed_prjToBoundary(rgeomObject, Dcoords, Dproj)
  
!<description>
  ! This routine projects a point onto the boundary of a composed object.
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
  
  ! Some temporary variables
  REAL(DP), DIMENSION(2) :: DcoordsRef, Dray
  REAL(DP), DIMENSION(:,:), POINTER :: p_Dverts
  REAL(DP) :: dminDist, ddist
  INTEGER :: iminDist, i
  TYPE(t_boundaryApprox), POINTER :: p_rbndApprox

    ! Get boundary approximation
    p_rbndApprox => rgeomObject%rcomposed%p_rboundaryApprox
    
    ! Check if we already have a boundary approximation vector
    IF(p_rbndApprox%h_Dverts .EQ. ST_NOHANDLE) THEN
    
      ! Update the boundary approximation for this object
      CALL geom_composed_updateBndApprox(rgeomObject)

    END IF
    
    ! Get vertice vector
    p_Dverts => p_rbndApprox%p_Dverts
    
    ! Transform point to relative coordinate system
    CALL bgeom_transformBackPoint2D(rgeomObject%rcoord2D, Dcoords, DcoordsRef)
    
    ! Get vector from point to boundary
    Dray = DcoordsRef - p_Dverts(1:2, 1)
    
    ! Calculate distance
    dminDist = Dray(1)**2 + Dray(2)**2
    iminDist = 1

    ! Go through all points of our boundary projection
    DO i = 2, p_rbndApprox%nverts
    
      ! Get vector from point to boundary
      Dray = DcoordsRef - p_Dverts(1:2, i)
      
      ! Calculate distance
      ddist = Dray(1)**2 + Dray(2)**2
      
      IF (ddist .LT. dminDist) THEN
        dminDist = ddist
        iminDist = i
      END IF
    
    END DO
    
    ! Transform projection to world coordinates
    CALL bgeom_transformPoint2D(rgeomObject%rcoord2D, p_Dverts(1:2, &
                                iminDist), Dproj)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE geom_composed_calcSignDist (rgeomObject, Dcoords, ddistance)

!<description>
  ! This routine calculates the signed distance of a given point and a 
  ! composed geometry object.
!</description>

!<input>
  ! The composed object
  TYPE(t_geometryObject), INTENT(IN)  :: rgeomObject
  
  ! The coordinates of the point that is to be tested.
  REAL(DP), DIMENSION(:), INTENT(IN)  :: Dcoords
  
!</input>

!<output>
  ! The shortest signed distance between the point and the object's boundary.
  REAL(DP),               INTENT(OUT) :: ddistance
!</output>

!</subroutine>

  ! Some temporary variables
  REAL(DP), DIMENSION(2) :: DcoordsRef, Dray
  REAL(DP), DIMENSION(:,:), POINTER :: p_Dverts
  REAL(DP) :: dminDist, ddist
  INTEGER :: i, iInside
  TYPE(t_boundaryApprox), POINTER :: p_rbndApprox

    ! Get boundary approximation
    p_rbndApprox => rgeomObject%rcomposed%p_rboundaryApprox

    ! Check if we already have a boundary approximation vector
    IF(p_rbndApprox%h_Dverts .EQ. ST_NOHANDLE) THEN
    
      ! Update the boundary approximation for this object
      CALL geom_composed_updateBndApprox(rgeomObject)
      
    END IF
    
    ! Get boundary approximation vector
    p_Dverts => p_rbndApprox%p_Dverts

    ! Transform point to relative coordinate system
    CALL bgeom_transformBackPoint2D(rgeomObject%rcoord2D, Dcoords, DcoordsRef)
    
    ! Get vector from point to boundary
    Dray = DcoordsRef - p_Dverts(1:2, 1)
    
    ! Calculate distance
    dminDist = Dray(1)**2 + Dray(2)**2

    ! Go through all points of our boundary projection
    DO i = 2, rgeomObject%rcomposed%p_rboundaryApprox%nverts
    
      ! Get vector from point to boundary
      Dray = DcoordsRef - p_Dverts(1:2, i)
      
      ! Calculate distance
      ddist = Dray(1)**2 + Dray(2)**2
      
      IF (ddist .LT. dminDist) THEN
        dminDist = ddist
      END IF
    
    END DO
    
    ! Get square root of distance
    ddistance = SQRT(dminDist)
    
    CALL geom_composed_isInGeometry(rgeomObject, Dcoords, iInside)
    
    ! Maybe the composed object is inverted?
    IF (iInside .EQ. 1) THEN
      ddistance = -ddistance
    END IF

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  RECURSIVE SUBROUTINE geom_composed_getNAV(rgeomObject, dtolerance, nverts)
  
!<description>
  ! Calculates the number of vertices needed to approximate the composed
  ! object's boundary with a desired tolerance.
!</description>

!<input>
  ! The composed object
  TYPE(t_geometryObject), INTENT(IN) :: rgeomObject
  
  ! The desired tolerance. Must be > EPS
  REAL(DP), INTENT(IN) :: dtolerance
!</input>

!<output>
  ! The number of vertices needed.
  INTEGER, INTENT(OUT) :: nverts
!</output>

!</subroutine>

  ! Some temporary variables
  INTEGER :: nsubVerts, i, h_Dverts
  REAL(DP) :: drelTol
  TYPE(t_geometryObject), POINTER :: p_rsubObject
  
    ! Maybe we already have an boundary approximation?
    h_Dverts = rgeomObject%rcomposed%p_rboundaryApprox%h_Dverts
    
    ! If yes, then we can return the number of allocated vertices
    IF (h_Dverts .NE. ST_NOHANDLE) THEN
    
      nverts = rgeomObject%rcomposed%p_rboundaryApprox%nverts

      RETURN

    END IF

    ! Initialise output
    nverts = 0
    
    ! Go through all sub-objects
    DO i = 1, rgeomObject%rcomposed%nsubObjects
      
      ! Get i-th sub-object
      p_rsubObject => rgeomObject%rcomposed%p_RsubObjects(i)
      
      ! Calculate tolerance
      !IF (p_rsubObject%ndimension .EQ. NDIM3D) THEN
      !  drelTol = dtolerance / p_rsubObject%rcoord3D%dscalingFactor
      !ELSE
        drelTol = dtolerance / p_rsubObject%rcoord2D%dscalingFactor
      !END IF
            
      ! Get number of vertices for sub-object
      CALL geom_getNumApproxVerts(p_rsubObject, drelTol, nsubVerts)
      
      ! Add vertices
      nverts = nverts + nsubVerts
    
    END DO
    
    ! That's it
    
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  RECURSIVE SUBROUTINE geom_composed_getBndApprox(rgeomObject, dtolerance, &
                                                  Dverts, nverts)

!<description>
  ! Calculates the boundary approximation vertices for the composed object.
!</description>

!<input>
  ! The composed
  TYPE(t_geometryObject), INTENT(IN) :: rgeomObject
  
  ! The desired tolerance. Must be > SYS_EPS.
  REAL(DP), INTENT(IN) :: dtolerance
!</input>

!<output>
  ! The vertice array.
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dverts
  
  ! Number of vertices created
  INTEGER, INTENT(OUT) :: nverts
!</output>

!</subroutine>

  ! Some temporary variables
  INTEGER :: i, j, k, io, nsubverts, nmaxverts, ctype
  INTEGER, DIMENSION(:), ALLOCATABLE :: Iinout
  REAL(DP) :: dsubtol
  REAL(DP), DIMENSION(2) :: DmyVert
  TYPE(t_geometryObject), POINTER :: p_rsubObject
  TYPE(t_boundaryApprox), POINTER :: p_rbndApprox
  
    ! Get boundary approximation
    p_rbndApprox => rgeomObject%rcomposed%p_rboundaryApprox
    IF (p_rbndApprox%h_Dverts .NE. ST_NOHANDLE) THEN
    
      ! Simply copy the vertices of the boundary approximation
      nverts = p_rbndApprox%nverts
      
      DO i = 1, nverts
        Dverts(:, i) = p_rbndApprox%p_Dverts(:, i)
      END DO
      
      ! And destroy our boundary approximation
      CALL storage_free(p_rbndApprox%h_Dverts)
      p_rbndApprox%nverts = 0
      p_rbndApprox%nvertsAllocated = 0
      p_rbndApprox%p_Dverts => NULL()
      
      ! Exit here
      RETURN
    
    END IF
  
    ! Get composed object's type
    ctype = rgeomObject%rcomposed%ccomposedType
  
    ! Get number of vertices
    CALL geom_composed_getNAV(rgeomObject, dtolerance, nmaxverts)
    
    ! Allocate integer array for inside-outside-test
    ALLOCATE(Iinout(nmaxverts))
    
    nverts = 0
    
    ! Now go through all our objects
    DO i = 1, rgeomObject%rcomposed%nsubObjects
    
      ! Get i-th sub-object
      p_rsubObject => rgeomObject%rcomposed%p_RsubObjects(i)
    
      ! Calculate tolerance for sub-object
      dsubtol = dtolerance / p_rsubObject%rcoord2D%dscalingFactor
    
      ! Get vertices of i-th sub-object
      CALL geom_getBoundaryApprox(p_rsubObject, dsubtol, &
                                  Dverts(:, nverts+1:), nsubverts)
      
      ! A-priori all vertices of this object belong to the composed object's
      ! boundary
      DO j = nverts + 1, nverts + nsubverts
        ! A-priori this vertice belongs to our boundary
        Iinout(j) = 1
        
        ! Transform vertice to our coordinate system
        CALL bgeom_transformPoint2D(p_rsubObject%rcoord2D, Dverts(1:2, j), &
                                    DmyVert)
        Dverts(1:2, j) = DmyVert
        
      END DO
      
      ! Now perform inside-outside test with all other sub-objects
      DO k = 1, rgeomObject%rcomposed%nsubObjects
      
        ! Skip i-th sub-object
        IF (k .NE. i) THEN
        
          ! Loop through all vertices of i-th sub-object
          DO j = nverts + 1, nverts + nsubverts
          
            ! If this point is already out, then skip it
            IF (Iinout(j) .NE. 0) THEN
          
              ! Perform inside-outside-test
              CALL geom_isInGeometry(rgeomObject%rcomposed%p_RsubObjects(k), &
                                   Dverts(1:2, j), io)

              ! If this is a UNION-quantor node and the point is inside the
              ! object, or if this node is a INTERSECT-quantor node and the
              ! point is outside, then we can kick it
              IF (((ctype .EQ. GEOM_COMP_TYPE_OR) .AND. (io .EQ. 1)) .OR. &
                  ((ctype .EQ. GEOM_COMP_TYPE_AND) .AND. (io .EQ. 0))) THEN
                Iinout(j) = 0
              END IF
            
            END IF
            
          END DO
        
        END IF
      
      END DO
      
      ! Increment vertice offset
      nverts = nverts + nsubverts
    
    END DO
    
    ! Now we need to compress the vertice vector, or in other words:
    ! We will kick out all vertices which do not belong to the boundary of
    ! this composed object.
    j = 0
    i = 1
    DO WHILE (i+j .LE. nverts)
      
      IF (Iinout(i+j) .NE. 0) THEN
      
        IF (j .GT. 0) THEN
          ! copy vertice
          Dverts(1:2,i) = Dverts(1:2,i+j)
        END IF
      
        ! increment destination index
        i = i+1
      ELSE
        ! vertice has been kicked
        ! increment source index
        j = j+1
      END IF
      
    END DO
    
    ! Print some debug info
    !PRINT *, 'Boundary Approx: ', nmaxverts, ' vertices allocated'
    !PRINT *, 'Boundary Approx: ', nverts, ' vertices created'
    !PRINT *, 'Boundary Approx: ', j, ' vertices kicked'
    !PRINT *, 'Boundary Approx: ', (i-1), ' vertices left'
    
    ! Remember how many vertices we have
    nverts = i-1
    
    ! Deallocate integer array
    DEALLOCATE(Iinout)
    
    ! That's it
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE geom_composed_updateBndApprox(rgeomObject, dtolerance)

!<description>
  ! This routine updates the boundary approximation of the composed object.
!</description>

!<input>
  ! OPTIONAL: The tolerance for the boundary approximation.
  ! If not present, 10e-2 is used.
  REAL(DP), OPTIONAL, INTENT(IN) :: dtolerance  
!</input>

!<input>
  ! The composed object
  TYPE(t_geometryObject), INTENT(IN) :: rgeomObject
!</input>

!</subroutine>

  ! Some temporary variables
  REAL(DP) :: dTol
  INTEGER :: nvertsNeeded, nvertsUsed, h_Dverts
  INTEGER, DIMENSION(2) :: DmyDim
  TYPE(t_boundaryApprox), POINTER :: p_rbndApprox

    ! Make sure this is a composed object
    IF (rgeomObject%ctype .NE. GEOM_COMPOSED) RETURN
    
    ! Get tolerance
    dTol = 1E-2_DP
    IF (PRESENT(dtolerance)) THEN
      IF ((dtolerance .GT. 1E-8_DP) .AND. (dtolerance .LE. 1E-1_DP)) THEN
        dTol = dtolerance
      END IF
    END IF

    ! Get a pointer to the boundary approximation
    p_rbndApprox => rgeomObject%rcomposed%p_rboundaryApprox

    ! Backup boundary approximation vertice handle
    h_Dverts = p_rbndApprox%h_Dverts
    
    ! And remove it from the structure
    p_rbndApprox%h_Dverts = ST_NOHANDLE

    ! Get number of vertices needed for tolerance
    CALL geom_composed_getNAV(rgeomObject, dTol, nvertsNeeded)
    
    ! Do we already have a boundary approximation?
    IF (h_Dverts .NE. ST_NOHANDLE) THEN
    
      ! Do we need to re-allocate the array?
      IF (p_rbndApprox%nvertsAllocated .LT. nvertsNeeded) THEN
      
        ! Reallocate array
        CALL storage_realloc('geom_composed_updateBndApprox', nvertsNeeded, &
                             h_Dverts, ST_NEWBLOCK_NOINIT, .FALSE.)
        
        ! Store allocated size
        p_rbndApprox%nvertsAllocated = nvertsNeeded
        
      END IF
      
    ELSE
        
      ! Allocate some space for us
      DmyDim(1) = 2
      DmyDim(2) = nvertsNeeded
      CALL storage_new2D('geom_composed_updateBndApprox', 'Dverts', &
                         DmyDim, ST_DOUBLE, h_Dverts, ST_NEWBLOCK_NOINIT)

      ! Store allocated size
      p_rbndApprox%nvertsAllocated = nvertsNeeded

    END IF
    
    ! Get the vertice vector
    CALL storage_getbase_double2D(h_Dverts, p_rbndApprox%p_Dverts)
    
    ! Get the object's boundary approximation
    CALL geom_composed_getBndApprox(rgeomObject, dTol, p_rbndApprox%p_Dverts, &
                                    nvertsUsed)
    
    ! Store tolerance and number of used vertices
    p_rbndApprox%dtolerance = dTol
    p_rbndApprox%nverts = nvertsUsed
    p_rbndApprox%h_Dverts = h_Dverts
    
    ! Check if we have allocated more than 200% of the vertices we use.
    IF ((2 * nvertsUsed) .LE. p_rbndApprox%nvertsAllocated) THEN
    
      ! Reallocate array
      CALL storage_realloc('geom_composed_updateBndApprox', nvertsUsed, &
                           p_rbndApprox%h_Dverts, ST_NEWBLOCK_NOINIT, .TRUE.)
      
      ! Remember we have reallocated the array
     p_rbndApprox%nvertsAllocated = nvertsUsed
      
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

!<subroutine>

  PURE SUBROUTINE geom_circle_getNAV(rgeomObject, dtolerance, nverts)
  
!<description>
  ! Calculates the number of vertices needed to approximate the circle's
  ! boundary with a desired tolerance.
!</description>

!<input>
  ! The circle
  TYPE(t_geometryObject), INTENT(IN) :: rgeomObject
  
  ! The desired tolerance. Must be > EPS
  REAL(DP), INTENT(IN) :: dtolerance
!</input>

!<output>
  ! The number of vertices needed.
  INTEGER, INTENT(OUT) :: nverts
!</output>

!</subroutine>

  ! The length of the boundary
  REAL(DP) :: dboundLen
  
    ! The length of the boundary is 2 * radius * PI
    dboundLen = 2.0_DP * rgeomObject%rcircle%dradius * SYS_PI
    
    ! The number of vertices is simply the boundary length divided by the
    ! tolerance
    nverts = INT(dboundLen / dtolerance)
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE geom_circle_getBndApprox(rgeomObject, dtolerance, Dverts, nverts)

!<description>
  ! Calculates the boundary approximation vertices for the circle.
!</description>

!<input>
  ! The circle
  TYPE(t_geometryObject), INTENT(IN) :: rgeomObject
  
  ! The desired tolerance. Must be > SYS_EPS.
  REAL(DP), INTENT(IN) :: dtolerance
!</input>

!<output>
  ! The vertice array.
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dverts
  
  ! Number of vertices created
  INTEGER, INTENT(OUT) :: nverts
!</output>

!</subroutine>

  ! Some temporary variables
  REAL(DP) :: dangle, drad, dstep
  INTEGER :: i
  
    ! Get circle radius
    drad = rgeomObject%rcircle%dradius
    
    ! Get number of vertices to allocate
    CALL geom_circle_getNAV(rgeomObject, dtolerance, nverts)

    ! Calculate angle delta    
    dstep = (2.0_DP * SYS_PI) / REAL(nverts, DP)
    
    ! Loop though all vertices
    DO i = 1, nverts
    
      ! Calculate angle
      dangle = REAL(i-1, DP) * dstep
    
      ! Calculate point
      Dverts(1, i) = drad * COS(dangle)
      Dverts(2, i) = drad * SIN(dangle)
      
    END DO
    
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

!<subroutine>

  PURE SUBROUTINE geom_square_getNAV(rgeomObject, dtolerance, nverts)
  
!<description>
  ! Calculates the number of vertices needed to approximate the square's
  ! boundary with a desired tolerance.
!</description>

!<input>
  ! The square
  TYPE(t_geometryObject), INTENT(IN) :: rgeomObject
  
  ! The desired tolerance. Must be > EPS
  REAL(DP), INTENT(IN) :: dtolerance
!</input>

!<output>
  ! The number of vertices needed.
  INTEGER, INTENT(OUT) :: nverts
!</output>

!</subroutine>

  ! The length of the boundary
  REAL(DP) :: dboundLen
  
    ! The length of the boundary is 4 * edge_length
    dboundLen = 4.0_DP * rgeomObject%rsquare%dlength
    
    ! The number of vertices is simply the boundary length divided by the
    ! tolerance
    nverts = INT(dboundLen / dtolerance) + 4
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE geom_square_getBndApprox(rgeomObject, dtolerance, Dverts, nverts)

!<description>
  ! Calculates the boundary approximation vertices for the square.
!</description>

!<input>
  ! The square
  TYPE(t_geometryObject), INTENT(IN) :: rgeomObject
  
  ! The desired tolerance. Must be > SYS_EPS.
  REAL(DP), INTENT(IN) :: dtolerance
!</input>

!<output>
  ! The vertice array.
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dverts
  
  ! Number of vertices created
  INTEGER, INTENT(OUT) :: nverts
!</output>

!</subroutine>

  ! Some temporary variables
  REAL(DP) :: dedge
  INTEGER :: i, nvpe
  
    ! Get square's edge length
    dedge = rgeomObject%rsquare%dlength
    
    ! Calculate number of vertices per edge
    nvpe = INT(dedge / dtolerance)
    
    ! Get half of square's edge length
    dedge = 0.5_DP * dedge

    ! Edge 1
    DO i = 1, nvpe
      Dverts(1, i) = dedge
      Dverts(2, i) = -dedge + (REAL(i-1, DP) * dtolerance)
    END DO
    
    ! Edge 2
    DO i = 1, nvpe
      Dverts(1, nvpe + i) = dedge - (REAL(i-1, DP) * dtolerance)
      Dverts(2, nvpe + i) = dedge
    END DO

    ! Edge 3
    DO i = 1, nvpe
      Dverts(1, 2*nvpe + i) = -dedge
      Dverts(2, 2*nvpe + i) = dedge - (REAL(i-1, DP) * dtolerance)
    END DO
    
    ! Edge 4
    DO i = 1, nvpe
      Dverts(1, 3*nvpe + i) = -dedge + (REAL(i-1, DP) * dtolerance)
      Dverts(2, 3*nvpe + i) = -dedge
    END DO
    
    ! Store number of vertices written
    nverts = 4 * nvpe
    
    ! That's it
  
  END SUBROUTINE
  
  ! ***************************************************************************
  ! *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*
  ! *= 2D Ellipse Routines                                                   =*
  ! *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE geom_init_ellipse_indirect(rgeomObject, rcoordSys, Dradii, &
                                        binverted)

!<description>
  ! Creates a t_geometryObject representing a 2D ellipse.
!</description>

!<input>
  ! A 2D coordinate system for the ellipse.
  TYPE(t_coordinateSystem2D),  INTENT(IN)  :: rcoordSys
  
  ! An array holding the X- and Y-radii for the ellipse.
  REAL(DP), DIMENSION(:),      INTENT(IN)  :: Dradii
  
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
    
    ! Store the X- and Y-radii of the ellipse
    rgeomObject%rellipse%Dradii = Dradii(1:2)
    
    ! That's it!
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE geom_init_ellipse_direct(rgeomObject, Dradii, Dorigin, &
                                      drotation, dscalingFactor, binverted)

!<description>
  ! Creates a t_geometryObject representing a 2D ellipse.
!</description>

!<input>
  ! A array holding the X- and Y-radii for the ellipse.
  REAL(DP), DIMENSION(:),            INTENT(IN)  :: Dradii
  
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
    rgeomObject%rellipse%Dradii = Dradii(1:2)
    
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
    DrelCoords(1) = DrelCoords(1) / rgeomObject%rellipse%Dradii(1)
    DrelCoords(2) = DrelCoords(2) / rgeomObject%rellipse%Dradii(2)
    
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

    ! Get the ellipse's X- and Y-radii
    dradX = rgeomObject%rellipse%Dradii(1)
    dradY = rgeomObject%rellipse%Dradii(2)
  
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
    dradiusX = rgeomObject%rellipse%Dradii(1)
    dradiusY = rgeomObject%rellipse%Dradii(2)
  
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

!<subroutine>

  PURE SUBROUTINE geom_ellipse_getNAV(rgeomObject, dtolerance, nverts, &
                                      dcircumference)
  
!<description>
  ! Calculates the number of vertices needed to approximate the ellipse's
  ! boundary with a desired tolerance.
!</description>

!<input>
  ! The ellipse
  TYPE(t_geometryObject), INTENT(IN) :: rgeomObject
  
  ! The desired tolerance. Must be > EPS
  REAL(DP), INTENT(IN) :: dtolerance
!</input>

!<output>
  ! The number of vertices needed.
  INTEGER, INTENT(OUT) :: nverts
  
  ! OPTIONAL: An approximation of the ellipse's circumference.
  REAL(DP), OPTIONAL, INTENT(OUT) :: dcircumference
!</output>

!</subroutine>

  ! Some temporary variables
  REAL(DP) :: dalpha, dbeta, dlambda, dtheta, domega, dr1, dr2

  ! The length of the boundary
  REAL(DP) :: dboundLen
  
    ! The boundary length of an ellipse cannot be calculated analytically,
    ! therefore, we will use the approximation of Ramanujan here.
  
    ! Get the radii
    dr1 = rgeomObject%rellipse%Dradii(1)
    dr2 = rgeomObject%rellipse%Dradii(2)

    ! Calculate alpha and beta    
    IF (dr1 .LT. dr2) THEN
      dalpha = dr2
      dbeta = dr1
    ELSE
      dalpha = dr1
      dbeta = dr2
    END IF
    
    ! Calculate lambda
    dlambda = (dalpha - dbeta) / (dalpha + dbeta)
    
    ! Calculate theta
    dtheta = 3.0_DP * dlambda**2
    
    ! Calculate omega
    domega = 1.0_DP + (dtheta / (10.0_DP + SQRT(4.0_DP - dtheta)))
  
    ! The approximative length of the boundary
    dboundLen = SYS_PI * (dalpha + dbeta) * domega
    
    ! The number of vertices is simply the boundary length divided by the
    ! tolerance
    nverts = INT(dboundLen / dtolerance)
    
    ! Return circumference if desired
    IF (PRESENT(dcircumference)) THEN
      dcircumference = dboundlen
    END IF
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE geom_ellipse_getBndApprox(rgeomObject, dtolerance, Dverts, nverts)

!<description>
  ! Calculates the boundary approximation vertices for the ellipse.
!</description>

!<input>
  ! The ellipse
  TYPE(t_geometryObject), INTENT(IN) :: rgeomObject
  
  ! The desired tolerance. Must be > SYS_EPS.
  REAL(DP), INTENT(IN) :: dtolerance
!</input>

!<output>
  ! The vertice array.
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dverts
  
  ! Number of vertices created
  INTEGER, INTENT(OUT) :: nverts
!</output>

!</subroutine>

  ! Some temporary variables
  REAL(DP) :: dangle, dradX, dradY, dcf, dSin, dCos
  INTEGER :: i
  
    ! Get elipse's radii
    dradX = rgeomObject%rellipse%Dradii(1)
    dradY = rgeomObject%rellipse%Dradii(2)
    
    ! Get number of vertices to allocate
    CALL geom_ellipse_getNAV(rgeomObject, dtolerance, nverts, dcf)
    
    ! Special thanks to Dr. Marcus Stiemer and Manuel Jaraczewski for their
    ! help with the following piece of code...
    
    ! Calculate angle delta    
    dangle = 0.0_DP
    
    ! Loop though all vertices
    DO i = 1, nverts
    
      ! Calculate SIN and COS
      dCos = COS(dangle)
      dSin = SIN(dangle)
    
      ! Calculate point
      Dverts(1, i) = dradX * dCos
      Dverts(2, i) = dradY * dSin
      
      ! Update angle
      dangle = dangle + dtolerance / SQRT((dradX*dSin)**2 + (dradY*dCos)**2)
    
    END DO

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

!<subroutine>

  PURE SUBROUTINE geom_rect_getNAV(rgeomObject, dtolerance, nverts)
  
!<description>
  ! Calculates the number of vertices needed to approximate the rectangle's
  ! boundary with a desired tolerance.
!</description>

!<input>
  ! The rectangle
  TYPE(t_geometryObject), INTENT(IN) :: rgeomObject
  
  ! The desired tolerance. Must be > EPS
  REAL(DP), INTENT(IN) :: dtolerance
!</input>

!<output>
  ! The number of vertices needed.
  INTEGER, INTENT(OUT) :: nverts
!</output>

!</subroutine>

  ! The length of the boundary
  REAL(DP) :: dboundLen
  
    ! The length of the boundary is 2 * (X_edge_length + Y_edge_length)
    dboundLen = 2.0_DP * (rgeomObject%rrectangle%Dlength(1) + &
                          rgeomObject%rrectangle%Dlength(2))
    
    ! The number of vertices is simply the boundary length divided by the
    ! tolerance
    nverts = INT(dboundLen / dtolerance) + 4
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE geom_rect_getBndApprox(rgeomObject, dtolerance, Dverts, nverts)

!<description>
  ! Calculates the boundary approximation vertices for the rectangle.
!</description>

!<input>
  ! The rectangle
  TYPE(t_geometryObject), INTENT(IN) :: rgeomObject
  
  ! The desired tolerance. Must be > SYS_EPS.
  REAL(DP), INTENT(IN) :: dtolerance
!</input>

!<output>
  ! The vertice array.
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dverts
  
  ! Number of vertices created
  INTEGER, INTENT(OUT) :: nverts
!</output>

!</subroutine>

  ! Some temporary variables
  REAL(DP) :: dedgeX, dedgeY
  INTEGER :: i, nvpeX, nvpeY, off
  
    ! Get rectangle's edge lengths
    dedgeX = rgeomObject%rrectangle%Dlength(1)
    dedgeY = rgeomObject%rrectangle%Dlength(2)
    
    ! Calculate number of vertices per edge
    nvpeX = INT(dedgeX / dtolerance)
    nvpeY = INT(dedgeY / dtolerance)
    
    ! Get half of rectangle's edge lengths
    dedgeX = 0.5_DP * dedgeX
    dedgeY = 0.5_DP * dedgeY

    ! Edge 1
    DO i = 1, nvpeY
      Dverts(1, i) = dedgeX
      Dverts(2, i) = -dedgeY + (REAL(i-1, DP) * dtolerance)
    END DO
    
    ! Edge 2
    DO i = 1, nvpeX
      Dverts(1, nvpeY + i) = dedgeX - (REAL(i-1, DP) * dtolerance)
      Dverts(2, nvpeY + i) = dedgeY
    END DO

    ! Edge 3
    off = nvpeX + nvpeY
    DO i = 1, nvpeY
      Dverts(1, off + i) = -dedgeX
      Dverts(2, off + i) = dedgeY - (REAL(i-1, DP) * dtolerance)
    END DO
    
    ! Edge 4
    off = nvpeX + 2*nvpeY
    DO i = 1, nvpeX
      Dverts(1, off + i) = -dedgeX + (REAL(i-1, DP) * dtolerance)
      Dverts(2, off + i) = -dedgeY
    END DO
    
    ! Store number of vertices written
    nverts = 2 * (nvpeX + nvpeY)
    
    ! That's it
  
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

  PURE SUBROUTINE geom_polygon_getNAV(rgeomObject, dtolerance, nverts)
  
!<description>
  ! Calculates the number of vertices needed to approximate the polygon's
  ! boundary with a desired tolerance.
!</description>

!<input>
  ! The rectangle
  TYPE(t_geometryObject), INTENT(IN) :: rgeomObject
  
  ! The desired tolerance. Must be > EPS
  REAL(DP), INTENT(IN) :: dtolerance
!</input>

!<output>
  ! The number of vertices needed.
  INTEGER, INTENT(OUT) :: nverts
!</output>

!</subroutine>

  ! The length of the boundary
  REAL(DP) :: dboundLen
  
  ! The bounds of the polygon
  INTEGER :: lb, ub, i
  
  ! The polygon's vertices
  REAL(DP), DIMENSION(2) :: Dedge
  
    ! Get the polygon's bounds  
    lb = LBOUND(rgeomObject%rpolygon%p_Dvertices, 2)
    ub = UBOUND(rgeomObject%rpolygon%p_Dvertices, 2);
    
    ! Get last edge
    Dedge = rgeomObject%rpolygon%p_Dvertices(1:2, lb) &
          - rgeomObject%rpolygon%p_Dvertices(1:2, ub)
    
    ! Add edge's length
    dboundLen = SQRT(Dedge(1)**2 + Dedge(2)**2)
    
    ! Go through all other edges
    DO i = lb, ub-1
    
      ! Get i-th edge
      Dedge = rgeomObject%rpolygon%p_Dvertices(1:2, i+1) &
            - rgeomObject%rpolygon%p_Dvertices(1:2, i)
            
      ! Add edge's length
      dboundLen = dboundLen + SQRT(Dedge(1)**2 + Dedge(2)**2)
      
    END DO
  
    
    ! The number of vertices is simply the boundary length divided by the
    ! tolerance
    nverts = INT(dboundLen / dtolerance) + ub - lb + 1
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE geom_polygon_getBndApprox(rgeomObject, dtolerance, Dverts, nverts)

!<description>
  ! Calculates the boundary approximation vertices for the polygon.
!</description>

!<input>
  ! The polygon
  TYPE(t_geometryObject), INTENT(IN) :: rgeomObject
  
  ! The desired tolerance. Must be > SYS_EPS.
  REAL(DP), INTENT(IN) :: dtolerance
!</input>

!<output>
  ! The vertice array.
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dverts
  
  ! Number of vertices created
  INTEGER, INTENT(OUT) :: nverts
!</output>

!</subroutine>

  ! Some temporary variables
  REAL(DP) :: dlen, dalpha
  REAL(DP), DIMENSION(2) :: Dedge, DlastVert
  INTEGER :: i, j, nvpe, lb, ub
  
    ! Get the polygon's bounds  
    lb = LBOUND(rgeomObject%rpolygon%p_Dvertices, 2)
    ub = UBOUND(rgeomObject%rpolygon%p_Dvertices, 2);
    
    ! Get last edge
    Dedge = rgeomObject%rpolygon%p_Dvertices(1:2, lb) &
          - rgeomObject%rpolygon%p_Dvertices(1:2, ub)
    
    ! Get last vertice
    DlastVert = rgeomObject%rpolygon%p_Dvertices(1:2, ub)

    ! Get edge's length
    dlen = SQRT(Dedge(1)**2 + Dedge(2)**2)
    
    ! Calculate number of vertices for this edge
    nvpe = INT(dlen / dtolerance)
    
    ! Loop through this edge's vertices
    DO j = 1, nvpe
    
      ! Calculate interpolation factor
      dalpha = REAL(j-1, DP) * dtolerance / dlen
      
      ! Store vertice
      Dverts(1:2, j) = DlastVert + dalpha * Dedge
    
    END DO
    
    ! Update number of vertices written
    nverts = nvpe
    
    ! Loop through all other edges
    DO i = lb, ub-1
    
      ! Get i-th edge
      Dedge = rgeomObject%rpolygon%p_Dvertices(1:2, i+1) &
            - rgeomObject%rpolygon%p_Dvertices(1:2, i)
      
      ! Get last vertice
      DlastVert = rgeomObject%rpolygon%p_Dvertices(1:2, i)

      ! Get edge's length
      dlen = SQRT(Dedge(1)**2 + Dedge(2)**2)
      
      ! Calculate number of vertices for this edge
      nvpe = INT(dlen / dtolerance)
      
      ! Loop through this edge's vertices
      DO j = 1, nvpe
      
        ! Calculate interpolation factor
        dalpha = REAL(j-1, DP) * dtolerance / dlen
        
        ! Store vertice
        Dverts(1:2, nverts + j) = DlastVert + dalpha * Dedge
      
      END DO
      
      ! Update number of vertices written
      nverts = nverts + nvpe
      
    END DO

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
   CASE (GEOM_COMPOSED)
     CALL geom_composed_prjToBoundary(rgeomObject, Dcoords, Dproj)
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
   CASE (GEOM_COMPOSED)
     CALL geom_composed_calcSignDist(rgeomObject, Dcoords, ddistance)
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
  ! The geometry object to polygonise.
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

  ! Desired vertice count and other temporary variables
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
  ! This routine releases allocated buffer from a geometry object.
  
!</description>

!<inputoutput>
  ! The geometry object that is to be destroyed
  TYPE(t_geometryObject), INTENT(INOUT) :: rgeomObject

!</inputoutput>

!</subroutine>

  ! Some temporary variables
  INTEGER :: i, h_Dverts

    ! If the object is composed, we have some work to do...
    IF(rgeomObject%ctype .EQ. GEOM_COMPOSED) THEN
    
      ! Go through all sub-objects and destroy them
      DO i=1, rgeomObject%rcomposed%nsubObjects
      
        ! Destroy sub-object
        CALL geom_done(rgeomObject%rcomposed%p_RsubObjects(i))
      
      END DO
      
      ! Deallocate object array
      DEALLOCATE(rgeomObject%rcomposed%p_RsubObjects)
      
      ! Destroy boundary approximation
      h_Dverts = rgeomObject%rcomposed%p_rboundaryApprox%h_Dverts
      IF (h_Dverts .NE. ST_NOHANDLE) THEN
        CALL storage_free(h_Dverts)
      END IF
      
      ! Deallocate boundary approximation structure
      DEALLOCATE(rgeomObject%rcomposed%p_rboundaryApprox)
    
    END IF

    ! Reset the object's type
    rgeomObject%ctype = GEOM_NONE
    
    ! That's it
  
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  RECURSIVE SUBROUTINE geom_getNumApproxVerts(rgeomObject, dtolerance, nverts)
  
!<description>
  ! This routine calculates the number of vertices needed to approximate the
  ! object's boundary at a desired tolerance.
  
  ! This routine is just a wrapper for the geom_XXXX_getNAV routines.
!</description>

!<input>
  ! The geometry object
  TYPE(t_geometryObject), INTENT(IN) :: rgeomObject
  
  ! The tolerance. Must be > EPS
  REAL(DP), INTENT(IN) :: dtolerance
  
!</input>

!<output>
  ! The number of vertices for the boundary approximation
  INTEGER, INTENT(OUT) :: nverts
  
!</output>

!</subroutine>

    SELECT CASE(rgeomObject%ctype)
    CASE (GEOM_COMPOSED)
      CALL geom_composed_getNAV(rgeomObject, dtolerance, nverts)
    CASE (GEOM_CIRCLE)
      CALL geom_circle_getNAV(rgeomObject, dtolerance, nverts)
    CASE (GEOM_SQUARE)
      CALL geom_square_getNAV(rgeomObject, dtolerance, nverts)
    CASE (GEOM_ELLIPSE)
      CALL geom_ellipse_getNAV(rgeomObject, dtolerance, nverts)
    CASE (GEOM_RECT)
      CALL geom_rect_getNAV(rgeomObject, dtolerance, nverts)
    CASE (GEOM_POLYGON)
      CALL geom_polygon_getNAV(rgeomObject, dtolerance, nverts)
    
    END SELECT
    
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  RECURSIVE SUBROUTINE geom_getBoundaryApprox(rgeomObject, dtolerance, &
                                              Dverts, nverts)
  
!<description>
  ! Calculates the boundary approximation vertices for a geometry object.
  ! This routine is just a wrapper for the geom_XXXX_getBndApprox routines.
!</description>

!<input>
  ! The composed object
  TYPE(t_geometryObject), INTENT(IN) :: rgeomObject
  
  ! The desired tolerance. Must be > SYS_EPS.
  REAL(DP), INTENT(IN) :: dtolerance
!</input>

!<output>
  ! The vertice array.
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dverts
  
  ! Number of vertices created
  INTEGER, INTENT(OUT) :: nverts
!</output>

!</subroutine>

    SELECT CASE(rgeomObject%ctype)
    CASE (GEOM_COMPOSED)
      CALL geom_composed_getBndApprox(rgeomObject, dtolerance, Dverts, nverts)
    CASE (GEOM_CIRCLE)
      CALL geom_circle_getBndApprox(rgeomObject, dtolerance, Dverts, nverts)
    CASE (GEOM_SQUARE)
      CALL geom_square_getBndApprox(rgeomObject, dtolerance, Dverts, nverts)
    CASE (GEOM_ELLIPSE)
      CALL geom_ellipse_getBndApprox(rgeomObject, dtolerance, Dverts, nverts)
    CASE (GEOM_RECT)
      CALL geom_rect_getBndApprox(rgeomObject, dtolerance, Dverts, nverts)
    CASE (GEOM_POLYGON)
      CALL geom_polygon_getBndApprox(rgeomObject, dtolerance, Dverts, nverts)
    END SELECT
    
  END SUBROUTINE
    
END MODULE
