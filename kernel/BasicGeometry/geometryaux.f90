!##############################################################################
!# ****************************************************************************
!# <name> geometryaux </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains low level auxiliary functions for geometric
!# objects like lines, triangles, quads or similar.
!#
!# The following routines can be found here:
!#
!#  1.) gaux_getAspectRatio_quad2D
!#      -> Calculate the aspect ratio of a 2D quadrilateral
!#
!#  2.) gaux_getArea_tria2D
!#      -> Calculate the signed area of a 2D triangle
!#
!#  3.) gaux_getArea_quad2D
!#      -> Calculate the area of a 2D quadrilateral
!#
!#  4.) gaux_getVolume_tetra3D
!#      -> Calculate the volume of a 3D tetrahedron
!#
!#  5.) gaux_getVolume_hexa3D
!#      -> Calculate the volume of a 3D hexahedron.
!#
!#  6.) gaux_isIntersection_line2D
!#      -> Checks if two (finite) lines in 2D intersect
!#
!#  7.) gaux_getIntersection_ray2D
!#      -> Calculates the intersection point of two rays (if they intersect)
!#
!#  8.) gaux_isInElement_quad2D
!#      -> Checks if a point is inside of a 2D quadrilateral
!#
!#  9.) gaux_isInElement_tri2D
!#      -> Checks if a point is inside of a 2D triangle
!#
!# 10.) gaux_getBarycentricCoords_tri2D
!#      -> Calculates the barycentric coordinates of a point relative
!#         to a specified triangle in 2D
!#
!# 11.) gaux_isFlipped_hexa3D
!#      -> Checks whether a hexahedron is flipped or not.
!#
!# 12.) gaux_isInElement_hexa
!#      -> Checks whether a point is in a hexa or not
!#
!# 13.) gaux_isIntersection_face
!#      -> test for face/line segment intersection
!#
!# 14.) gaux_isIntersection_triangle
!#      -> test for triangle/line segment intersection
!#
!# 15.) gaux_intersect_edgecircle
!#      -> calculates intersection of edge and circle
!#
!# 16.) gaux_intersect_quadcircle
!#      -> calculates intersection of quad and circle
!#
!# 17.) gaux_getDirectedExtentQuad2D
!#      -> Calculates the maximum length of a line in a given direction
!#         when embedded in a 2D quadrilateral
!# 
!# </purpose>
!##############################################################################

module geometryaux

!$ use omp_lib
  use fsystem
  use basicgeometry
  use genoutput

  implicit none

  private

  public :: gaux_getAspectRatio_quad2D
  public :: gaux_getArea_tria2D
  public :: gaux_getArea_quad2D
  public :: gaux_getVolume_tetra3D
  public :: gaux_getVolume_hexa3D
  public :: gaux_isIntersection_line2D
  public :: gaux_getIntersection_ray2D
  public :: gaux_isIntersection_face
  public :: gaux_isInElement_quad2D
  public :: gaux_isInElement_tri2D
  public :: gaux_getBarycentricCoords_tri2D
  public :: gaux_isFlipped_hexa3D
  public :: gaux_isInElement_tetra
  public :: gaux_isInElement_hexa
  public :: gaux_calcDistPEdg2D
  public :: gaux_projectPointPlane
  public :: gaux_isInElement_hexa_aligned
  public :: gaux_intersect_edgecircle
  public :: gaux_intersect_quadcircle
  public :: gaux_getDirectedExtentQuad2D
  
contains

  ! ***************************************************************************

!<function>

  pure real(DP) function gaux_getAspectRatio_quad2D (Dpoints)

!<description>
  ! This routine calculates the aspect ratio of a 2D quadrilateral
  ! polygon. The polygon is given by the coordinates of its four corners,
  ! counterclockwise.
!</description>

!<input>
  ! The coordinates of the four corners of the polygon, ordered
  ! counterclockwise.
  ! Dpoints(1,.) = x-coordinates,
  ! Dpoints(2,.) = y-coordinates
  real(DP), dimension(2,4), intent(in) :: Dpoints
!</input>

!<result>
  ! The aspect ratio of the polygon.
!</result>

!</function>

    real(DP) :: XM1,XM2,XM3,XM4,YM1,YM2,YM3,YM4,LE1,LE2

    ! Calculate the position of the edge mitpoints
    XM1 = (Dpoints(1,2)+Dpoints(1,1))*0.5_DP
    YM1 = (Dpoints(2,2)+Dpoints(2,1))*0.5_DP
    XM2 = (Dpoints(1,3)+Dpoints(1,2))*0.5_DP
    YM2 = (Dpoints(2,3)+Dpoints(2,2))*0.5_DP
    XM3 = (Dpoints(1,4)+Dpoints(1,3))*0.5_DP
    YM3 = (Dpoints(2,4)+Dpoints(2,3))*0.5_DP
    XM4 = (Dpoints(1,1)+Dpoints(1,4))*0.5_DP
    YM4 = (Dpoints(2,1)+Dpoints(2,4))*0.5_DP

    ! Calculate the length of the edge between opposite midpoints
    ! (squared)...
    LE1 = (XM3-XM1)**2+(YM3-YM1)**2
    LE2 = (XM4-XM2)**2+(YM4-YM2)**2

    ! ...and the aspect ratio.
    gaux_getAspectRatio_quad2D = sqrt(LE1/max(LE2,1E-10_DP))

  end function

  ! ***************************************************************************

!<function>

  pure real(DP) function gaux_getArea_tria2D (Dpoints)

!<description>
    ! This routine calculates the signed area of a 2D triangular
    ! polygon. The polygon is given by the coordinates of its three
    ! corners, counterclockwise.
!</description>

!<input>
  ! The coordinates of the three corners of the polygon, ordered
  ! counterclockwise.
  ! Dpoints(1,.) = x-coordinates,
  ! Dpoints(2,.) = y-coordinates
    real(DP), dimension(2,3), intent(in) :: Dpoints
!</input>

!<result>
  ! The signed area of the polygon.
!</result>
!</function>

    gaux_getArea_tria2D = 0.5_DP*( &
        (Dpoints(1,2)-Dpoints(1,1))*(Dpoints(2,3)-Dpoints(2,1))-&
        (Dpoints(1,3)-Dpoints(1,1))*(Dpoints(2,2)-Dpoints(2,1)) )
  end function gaux_getArea_tria2D

  ! ***************************************************************************

!<function>

  real(DP) function gaux_getArea_quad2D (Dpoints)

!<description>
    ! This routine calculates the area of a 2D quadrilateral
    ! polygon. The polygon is given by the coordinates of its three
    ! corners, counterclockwise.
!</description>

!<input>
  ! The coordinates of the four corners of the polygon, ordered
  ! counterclockwise.
  ! Dpoints(1,.) = x-coordinates,
  ! Dpoints(2,.) = y-coordinates
    real(DP), dimension(2,4), intent(in) :: Dpoints
!</input>

!<result>
  ! The area of the polygon.
!</result>
!</function>

    ! local variables
    !INTEGER, DIMENSION(3), PARAMETER :: IpointSubset = (/1,3,4/)
    !REAL(DP), DIMENSION(2,3) :: Dpoints2
    !REAL(DP) :: daux1,daux2

    ! First approach (slow):
    ! Add the area of two subtriangles to get the area of the quad.
    !
    !Dpoints2 = Dpoints(:,IpointSubset)
    !daux1 = ABS(gaux_getArea_tria2D(Dpoints(:,1:3))) + &
    !        ABS(gaux_getArea_tria2D(Dpoints2))

    ! For computing the area of the quad, we use a formula documented at
    !
    !   http://softsurfer.com/Archive/algorithm_0101/algorithm_0101.htm
    !
    ! The formula was discovered by Pierre Varignon and first published in 1731.
    !
    ! `In modern linear algebra, as already noted, the area of a planar parallelogram
    !  is the magnitude of the cross product of two adjacent edge vectors.  So, for
    !  any 3D planar parallelogram V0V1V2V3, we have:`
    !
    !     A(V0V1V2V3) = 2 A(M0M1M2M3)
    !                 = 2 | (V1-V0) X (V3-V0) |     (cross product)
    !                 = ...
    !                 = (x1-x0)(y3-y0) - (x3-x1)(y1-y0)
    !
    !                     V3                V2
    !                       +--------------+
    !                     /              /
    !                   /              /
    !                 /              /
    !               /              /
    !             /              /
    !           /              /
    !         /              /
    !       +--------------+
    !     V0                V1
    !
    !
    ! `Next, for an arbitrary quadrilateral, one can compute its area using a
    !  parallelogram discovered by Pierre Varignon (first published in 1731).
    !  It is amazing that the Greeks missed Varignon`s simple result which was
    !  discovered 2000 years after Euclid!  Given any quadrilateral, one can take
    !  the midpoints of its 4 edges to get 4 vertices which form a new quadrilateral.
    !  It is then easy to show that this midpoint quadrilateral is always a
    !  parallelogram, called the "Varignon parallelogram", and that its area is
    !  exactly one-half the area of the original quadrilateral.  So, for a
    !  quadrilateral Q=V0V1V2V3, let this parallelogram have midpoint vertices
    !  M0M1M2M3 as shown in the diagram:`
    !
    !     V3         M2          V2
    !       +---------X---------+
    !       |       /   \       |
    !       |    /         \    |
    !       | /               \ |
    !    M3 X                   X M1
    !       | \               / |
    !       |    \         /    |
    !       |       \   /       |
    !       +---------X---------+
    !     V0         M0          V1
    !
    ! `From elementary geometry, we know that in triangle V0V1V2 the midpoint
    !  line M0M1 is parallel to the base V0V2.  In triangle V0V1V2V3, the line M3M2
    !  is parallel to that same base V0V2.  Thus, M0M1 and M3M2 are parallel to each other.
    !  Similarly, M0M3 and M1M2 are parallel, which shows that M0M1M2M3 is a parallelogram.
    !  The area relation is also easy to demonstrate, and we can compute the quadrilateral`s
    !  area as:`
    !
    !     A(V0V1V2V3) = 2 A(M0M1M2M3)
    !                 = 2 | (M1-M0) X (M3-M0) |     (cross product)
    !                 = ...
    !                 = 1/2 | (V2-V0) X (V3-V1) |
    !                 = ...
    !                 = (x2-x0)(y3-y1) - (x3-x1)(y2-y0)
    !
    ! `This formula for an arbitrary quadrilateral is just as efficient as the one for an
    !  arbitrary triangle, using only 2 multiplications and 5 additions.  For simple
    !  quadrilaterals, the area is positive when the vertices are oriented counterclockwise,
    !  and negative when they are clockwise.  However, it also works for nonsimple
    !  quadrilaterals and is equal to the difference in area of the two regions the
    !  quadrilateral bounds.  For example, in the following diagram where I is the
    !  self-intersection point of a nonsimple quadrilateral V0V1V2V3, we have:`
    !
    !     A(V0V1V2V3) = A(V0V1I) + A(IV2V3)
    !                 = A(V0V1I) + A(IV2V2)
    !
    !     V2                V3
    !       +--------------+
    !         \           /
    !           \       /
    !             \   /
    !               + I
    !             /   \                 .
    !           /       \               .
    !         /           \             .
    !       +--------------+
    !     V0                 V1

    gaux_getArea_quad2D = &
            0.5_DP *  abs( (Dpoints(1,3) - Dpoints(1,1) ) * &
                           (Dpoints(2,4) - Dpoints(2,2) ) - &
                           (Dpoints(1,4) - Dpoints(1,2) ) * &
                           (Dpoints(2,3) - Dpoints(2,1) ) )

  end function gaux_getArea_quad2D

!************************************************************************

!<function>

  pure real(DP) function gaux_getVolume_tetra3D (Dpoints)

!<description>
  ! This function calculates the volume of a 3D tetrahedron. The
  ! tetrahedron is given by the coordinates of its four corners.
!</description>

!<input>
  ! The coordinates of the four corners of the tetrahedron.
  ! Dpoints(1,.) = x-coordinates,
  ! Dpoints(2,.) = y-coordinates,
  ! Dpoints(3,.) = z-coordinates
  real(DP), dimension(3,4), intent(in) :: Dpoints
!</input>

!<result>
  ! The volume of the tetrahedron.
!</result>
!</function>

    ! A temporary array for the edge lengths
    real(DP), dimension(3,3) :: Dv

    Dv(1:3,1) = Dpoints(1:3,1) - Dpoints(1:3,4)
    Dv(1:3,2) = Dpoints(1:3,2) - Dpoints(1:3,4)
    Dv(1:3,3) = Dpoints(1:3,3) - Dpoints(1:3,4)

    ! Return the absolute volume
    gaux_getVolume_tetra3D = abs(&
        Dv(1,1) * (Dv(2,2)*Dv(3,3) - Dv(3,2)*Dv(2,3)) + &
        Dv(2,1) * (Dv(3,2)*Dv(1,3) - Dv(1,2)*Dv(3,3)) + &
        Dv(3,1) * (Dv(1,2)*Dv(2,3) - Dv(2,2)*Dv(1,3))) / 6.0_DP

  end function gaux_getVolume_tetra3D

!************************************************************************

!<function>

  pure real(DP) function gaux_getVolume_hexa3D (Dv)

!<description>
  ! This function calculates the volume of a 3D hexahedron. The
  ! hexahedron is given by the coordinates of its eight corners.
!</description>

!<input>
  ! The coordinates of the eight corners of the hexahedron.
  ! Dv(1,.) = x-coordinates,
  ! Dv(2,.) = y-coordinates,
  ! Dv(3,.) = z-coordinates
  real(DP), dimension(3,8), intent(in) :: Dv
!</input>

!<result>
  ! The volume of the hexahedron.
!</result>
!</function>

    ! Return the absolute volume
    gaux_getVolume_hexa3D = (1.0_DP / 6.0_DP) * (&
        abs((Dv(1,4)-Dv(1,1))*(Dv(2,4)-Dv(2,3))*(Dv(3,4)-Dv(3,8))&
           +(Dv(2,4)-Dv(2,1))*(Dv(3,4)-Dv(3,3))*(Dv(1,4)-Dv(1,8))&
           +(Dv(3,4)-Dv(3,1))*(Dv(1,4)-Dv(1,3))*(Dv(2,4)-Dv(2,8))&
           -(Dv(1,4)-Dv(1,8))*(Dv(2,4)-Dv(2,3))*(Dv(3,4)-Dv(3,1))&
           -(Dv(2,4)-Dv(2,8))*(Dv(3,4)-Dv(3,3))*(Dv(1,4)-Dv(1,1))&
           -(Dv(3,4)-Dv(3,8))*(Dv(1,4)-Dv(1,3))*(Dv(2,4)-Dv(2,1)))+&
        abs((Dv(1,2)-Dv(1,3))*(Dv(2,2)-Dv(2,1))*(Dv(3,2)-Dv(3,6))&
           +(Dv(2,2)-Dv(2,3))*(Dv(3,2)-Dv(3,1))*(Dv(1,2)-Dv(1,6))&
           +(Dv(3,2)-Dv(3,3))*(Dv(1,2)-Dv(1,1))*(Dv(2,2)-Dv(2,6))&
           -(Dv(1,2)-Dv(1,6))*(Dv(2,2)-Dv(2,1))*(Dv(3,2)-Dv(3,3))&
           -(Dv(2,2)-Dv(2,6))*(Dv(3,2)-Dv(3,1))*(Dv(1,2)-Dv(1,3))&
           -(Dv(3,2)-Dv(3,6))*(Dv(1,2)-Dv(1,1))*(Dv(2,2)-Dv(2,3)))+&
        abs((Dv(1,5)-Dv(1,8))*(Dv(2,5)-Dv(2,6))*(Dv(3,5)-Dv(3,1))&
           +(Dv(2,5)-Dv(2,8))*(Dv(3,5)-Dv(3,6))*(Dv(1,5)-Dv(1,1))&
           +(Dv(3,5)-Dv(3,8))*(Dv(1,5)-Dv(1,6))*(Dv(2,5)-Dv(2,1))&
           -(Dv(1,5)-Dv(1,1))*(Dv(2,5)-Dv(2,6))*(Dv(3,5)-Dv(3,8))&
           -(Dv(2,5)-Dv(2,1))*(Dv(3,5)-Dv(3,6))*(Dv(1,5)-Dv(1,8))&
           -(Dv(3,5)-Dv(3,1))*(Dv(1,5)-Dv(1,6))*(Dv(2,5)-Dv(2,8)))+&
        abs((Dv(1,7)-Dv(1,6))*(Dv(2,7)-Dv(2,8))*(Dv(3,7)-Dv(3,3))&
           +(Dv(2,7)-Dv(2,6))*(Dv(3,7)-Dv(3,8))*(Dv(1,7)-Dv(1,3))&
           +(Dv(3,7)-Dv(3,6))*(Dv(1,7)-Dv(1,8))*(Dv(2,7)-Dv(2,3))&
           -(Dv(1,7)-Dv(1,3))*(Dv(2,7)-Dv(2,8))*(Dv(3,7)-Dv(3,6))&
           -(Dv(2,7)-Dv(2,3))*(Dv(3,7)-Dv(3,8))*(Dv(1,7)-Dv(1,6))&
           -(Dv(3,7)-Dv(3,3))*(Dv(1,7)-Dv(1,8))*(Dv(2,7)-Dv(2,6)))+&
        abs((Dv(1,1)-Dv(1,3))*(Dv(2,1)-Dv(2,8))*(Dv(3,1)-Dv(3,6))&
           +(Dv(2,1)-Dv(2,3))*(Dv(3,1)-Dv(3,8))*(Dv(1,1)-Dv(1,6))&
           +(Dv(3,1)-Dv(3,3))*(Dv(1,1)-Dv(1,8))*(Dv(2,1)-Dv(2,6))&
           -(Dv(1,1)-Dv(1,6))*(Dv(2,1)-Dv(2,8))*(Dv(3,1)-Dv(3,3))&
           -(Dv(2,1)-Dv(2,6))*(Dv(3,1)-Dv(3,8))*(Dv(1,1)-Dv(1,3))&
           -(Dv(3,1)-Dv(3,6))*(Dv(1,1)-Dv(1,8))*(Dv(2,1)-Dv(2,3))))

  end function gaux_getVolume_hexa3D

!************************************************************************

!<subroutine>

  elemental subroutine gaux_isIntersection_line2D(&
      dx1,dy1,dx2,dy2,dx3,dy3,dx4,dy4, bintersect)

!<description>
  ! Checks whether the two 2D lines given by the start/endpoints
  ! (x1,y1)->(x2,y2) and (x3,y3)->(x4,y4) intersect each other.
!</description>

!<input>
  ! First point on ray 1.
  real(DP), intent(in) :: dx1,dy1

  ! A second point on ray 1. Must be different to (dx1,dy1)
  real(DP), intent(in) :: dx2,dy2

  ! First point on ray 2.
  real(DP), intent(in) :: dx3,dy3

  ! A second point on ray 2. Must be different to (dx3,dy3)
  real(DP), intent(in) :: dx4,dy4
!</input>

!<result>
  ! TRUE if the two rays intersect. FALSE otherwise.
  logical, intent(out) :: bintersect
!</result>

!</subroutine>

    ! local variables: aux parameters
    real(DP) :: daux1, daux2, daux3, daux4

    ! position of point 3 with respect to line between 1 and 2

    daux3 = (dx2-dx1)*(dy3-dy1) - (dy2-dy1)*(dx3-dx1)

    ! position of point 4 with respect to line between 1 and 2

    daux4 = (dx2-dx1)*(dy4-dy1) - (dy2-dy1)*(dx4-dx1)

    ! position of point 1 with respect to line between 3 and 4

    daux1 = (dx4-dx3)*(dy1-dy3) - (dy4-dy3)*(dx1-dx3)

    ! position of point 2 with respect to line between 3 and 4

    daux2 = (dx4-dx3)*(dy2-dy3) - (dy4-dy3)*(dx2-dx3)

    ! Determine if the lines truely intersect by checking the sign

    bintersect = ((daux3*daux4 .le. 0.0_DP) .and. &
                  (daux1*daux2 .le. 0.0_DP))

  end subroutine

!************************************************************************

!<subroutine>

  elemental subroutine gaux_getIntersection_ray2D(&
      dx0,dy0,dx1,dy1,dx2,dy2,dx3,dy3, dx,dy, iintersect, da)

!<description>
  ! Calculates the intersection point of two 2D rays given by
  ! (x1,y1)->(x2,y2) and (x3,y3)->(x4,y4).
!</description>

!<input>
  ! First point on ray 1.
  real(DP), intent(in) :: dx0,dy0

  ! A second point on ray 1. Must be different to (dx1,dy1)
  real(DP), intent(in) :: dx1,dy1

  ! First point on ray 2.
  real(DP), intent(in) :: dx2,dy2

  ! A second point on ray 2. Must be different to (dx3,dy3)
  real(DP), intent(in) :: dx3,dy3
!</input>

!<result>
  ! Intersection point.
  ! If the two rays do not intersect or are identical, this is set to (0,0).
  real(DP), intent(out) :: dx,dy

  ! Returns the type of intersection between the rays.
  ! =-1: The rays are the same
  ! = 0: The rays do not intersect.
  ! = 1: The rays intersect in exactly one point.
  integer, intent(out) :: iintersect

  ! Parameter value of the intersection.
  ! The intersection point (dx,dy) can be found at position
  ! (dx,dy) = (dx0,dy0) + da*(dx1-dx0,dy1-dy0).
  ! If iintersect<>1, da is set to 0.
  real(DP), intent(out) :: da
!</result>

!</subroutine>

    ! local variables
    real(DP) :: ddet

    ! Initial setting of the destination point
    dx = 0.0_DP
    dy = 0.0_DP
    iintersect = 0

    ! We have (hopefully) the situation
    !
    !               (X1,Y1)
    !                  |
    !                  |
    !  (X2,Y2) --------+--------- (X3,Y3)
    !                  |
    !                  |
    !               (X0,Y0)
    !
    ! and want to calculate the intersection point. This means
    ! we have to solve the linear system
    !
    !  ( X1-X0  X2-X3 ) * (a) = ( X2-X0 )
    !  ( Y1-Y0  Y2-Y3 )   (b)   ( Y2-Y0 )
    !
    ! to get the "parameter" values a,b along the two lines where
    ! the intersection occurs.
    !
    ! The determinant of the system is:

    ddet = dx1*dy2-dx1*dy3-dx0*dy2+dx0*dy3-dy1*dx2+dy1*dx3+dy0*dx2-dy0*dx3

    ! If it is =0, the lines are the same or completely different...

    if (ddet .eq. 0.0_DP) then

      ! If the vector (X2,Y2)->(X0,Y0) is linear dependent to
      ! (X2,Y2)->(X3,Y3), the lines are the same.

      ddet = -dy0*dx2-dx3*dy2+dy0*dx3+dx2*dy3+dx0*dy2-dx0*dy3

      if (ddet .eq. 0.0_DP) then
        iintersect = -1
      end if

      da = 0

    else

      ! There is an intersection point. Calculate one of the
      ! "parameter" values along the two lines.

      da = (dy0*dx2+dx3*dy2-dy0*dx3-dx2*dy3-dx0*dy2+dx0*dy3) / ddet

      !  The intersection point is then

      dx = da*dx1 + (1.0_DP-da)*dx0
      dy = da*dy1 + (1.0_DP-da)*dy0

      iintersect = 1

    end if

  end subroutine

!************************************************************************

#if defined(DEBUG) || defined(_DEBUG)
  subroutine gaux_isInElement_quad2D(dx,dy,DcornerCoords,binside)
#else

!<subroutine>
  pure subroutine gaux_isInElement_quad2D(dx,dy,DcornerCoords,binside)
#endif

!<description>
  ! Checks if a point (dx,dy) is inside of a 2D quadrilateral element
  ! given by the corners DcornerCoords.
!</description>

!<input>
  ! Point to check
  real(DP), intent(in) :: dx,dy

  ! Array with coordinates of the four corner points of the element.
  ! The corners must be ordered in counterclockwise order.
  !
  ! Note: For performance reasons, this array is defined as
  !   explicit array of dimension (2,4). As this deactivates array
  !   checking in Fortran, the caller must take care to specify exactly
  !   this type of array here!
#if defined(DEBUG) || defined(_DEBUG)
  real(DP), dimension(2,4), intent(in) :: DcornerCoords
#endif
!</input>

#if !defined(DEBUG) && !defined(_DEBUG)
  real(DP), dimension(:,:), intent(in) :: DcornerCoords
#endif

!<result>
  ! TRUE if (dx,dy) is inside of the element. FALSE otherwise.
  logical, intent(out) :: binside
!</result>

!</subroutine>

    ! local variables

    integer, parameter :: NVE = 4
    integer, dimension(4), parameter :: Inext = (/2,3,4,1/)
    real(DP) :: dxmid,dymid,dxdist,dydist,dxnormal,dynormal
    integer :: ive,ive2
    real(DP) :: dsproduct

#if defined(DEBUG) || defined(_DEBUG)
    if (ubound(DcornerCoords,1) .ne. NDIM2D) then
      call output_line ('Dimension 1 is not =2!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gaux_isInElement_quad2D')
      call sys_halt()
    end if
#endif

    binside = .true.

    ! Compute edge-midpoints and normal vectors to the four
    ! edges on element IEL

    do ive=1,NVE

      ive2 = Inext(ive)   ! Use of Inext avoids a division by avoiding MOD!

      ! compute midpoints of element edges

      dxmid = 0.5_DP*(DcornerCoords(1,ive) + DcornerCoords(1,ive2))
      dymid = 0.5_DP*(DcornerCoords(2,ive) + DcornerCoords(2,ive2))

      ! compute normal vectors to element edges

      dxnormal =  DcornerCoords(2,ive2) - DcornerCoords(2,ive)
      dynormal = -DcornerCoords(1,ive2) + DcornerCoords(1,ive)

      ! compute vectors from edge midpoints to node 'ivt'

      dxdist = dx - dxmid
      dydist = dy - dymid

      ! Check whether 'ivt' belongs to element 'iel' or not by
      ! multiplying distance vectors with corresponding normal vectors.
      ! The sign of this scalar product determines whether we are
      ! 'left' or 'right' of the edge (because of the cosine formula).
      ! If the point is "righthand" of all four edges, it is inside
      ! of the element.

      dsproduct = dxdist*dxnormal + dydist*dynormal

      ! Actually we have to check against <=0, but it is more advisable
      ! to check against something that is 'near' 0 in terms
      ! of machine exactness...
      binside = binside .and. (dsproduct .le. SYS_EPSREAL_DP*100.0_DP)

    end do

  end subroutine

!************************************************************************

!<subroutine>

  pure subroutine gaux_getBarycentricCoords_tri2D(&
      DcornerCoords,dx,dy,dxi1,dxi2,dxi3)

!<description>
  ! Calculates the barycentric coordinates (dxi1,dxi2,dxi3) of a
  ! point (dx,dy) relative to a 2D triangular element specified by
  ! the coordinates of the three corners in DcornerCoords.
!</description>

!<input>
  ! Point in real coordinates
  real(DP), intent(in) :: dx,dy

  ! Array with coordinates of the three corner points of the element.
  ! The corners must be ordered in counterclockwise order.
  real(DP), dimension(:,:), intent(in) :: DcornerCoords
!</input>

!<result>
  ! The barycentric coordinates of (dx,dy) relative to the element
  ! specified by DcornerCoords.
  real(DP), intent(out) :: dxi1,dxi2,dxi3
!</result>

!</subroutine>

    ! local variables
    real(DP) :: DAX, DAY, DBX, DBY, DCX, DCY, DDET

    DAX = DcornerCoords(1, 1)
    DAY = DcornerCoords(2, 1)
    DBX = DcornerCoords(1, 2)
    DBY = DcornerCoords(2, 2)
    DCX = DcornerCoords(1, 3)
    DCY = DcornerCoords(2, 3)

    ! Example where to find this formula here:
    ! http://home.t-online.de/home/nagel.klaus/matdir/bary.htm

    DDET = 1.0_DP / ( DAX*(DBY-DCY) + DBX*(DCY-DAY) + DCX*(DAY-DBY) )
    dxi1 = (dx*(DBY-DCY)+DBX*(DCY-dy)+DCX*(dy-DBY)) * DDET
    dxi2 = (DAX*(dy-DCY)+dx*(DCY-DAY)+DCX*(DAY-dy)) * DDET
    dxi3 = (DAX*(DBY-dy)+DBX*(dy-DAY)+dx*(DAY-DBY)) * DDET

  end subroutine

!************************************************************************

!<subroutine>

  pure subroutine gaux_isInElement_tri2D(dx,dy,DcornerCoords,binside)

!<description>
  ! Checks if a point (dx,dy) is inside of a 2D triangular element
  ! given by the corners DcornerCoords.
!</description>

!<input>
  ! Point to check
  real(DP), intent(in) :: dx,dy

  ! Array with coordinates of the four corner points of the element.
  ! The corners must be ordered in counterclockwise order.
  !
  ! Note: For performance reasons, this array is defined as
  !   explicit array of dimension (2,4). As this deactivates array
  !   checking in Fortran, the caller must take care to specify exactly
  !   this type of array here!
  real(DP), dimension(2,4), intent(in) :: DcornerCoords
!</input>

!<result>
  ! TRUE if (dx,dy) is inside of the element. FALSE otherwise.
  logical, intent(out) :: binside
!</result>

!</subroutine>

    real(DP) :: dxi1,dxi2,dxi3

    ! We use barycentric coordinates for that task.
    ! Calculate the barycentric coordinates of the point relative
    ! to the element specified by DcornerCoords.
    call gaux_getBarycentricCoords_tri2D (DcornerCoords,dx,dy,&
        dxi1,dxi2,dxi3)

    ! If all barycentric coordinates are in the range [0..1],
    ! we are inside of the element
    binside = (dxi1 .ge. 0.0_DP) .and. (dxi1 .le. 1.0_DP) .and. &
              (dxi2 .ge. 0.0_DP) .and. (dxi2 .le. 1.0_DP) .and. &
              (dxi3 .ge. 0.0_DP) .and. (dxi3 .le. 1.0_DP)

  end subroutine


!************************************************************************

!<function>

  pure logical function gaux_isFlipped_hexa3D (Dpoints)

!<description>
  ! This function checks whether a 3D hexahedron is flipped.
!</description>

!<input>
  ! The coordinates of the eight corners of the hexahedron.
  ! Dpoints(1,.) = x-coordinates,
  ! Dpoints(2,.) = y-coordinates,
  ! Dpoints(3,.) = z-coordinates
  real(DP), dimension(3,8), intent(in) :: Dpoints
!</input>

!<result>
  ! .TRUE. if the hexahedron is flipped, otherwise .FALSE.
!</result>
!</function>

    ! Three vectors connecting two opposite faces of the hexahedron,
    ! and a normal vector
    real(DP), dimension(3) :: Du,Dv,Dw,Dn
    real(DP) :: dt

    Du(:) = 0.25_DP * (Dpoints(:,5)+Dpoints(:,6)+Dpoints(:,7)+Dpoints(:,8)&
                      -Dpoints(:,1)-Dpoints(:,2)-Dpoints(:,3)-Dpoints(:,4))
    Dv(:) = 0.25_DP * (Dpoints(:,3)+Dpoints(:,4)+Dpoints(:,7)+Dpoints(:,8)&
                      -Dpoints(:,1)-Dpoints(:,2)-Dpoints(:,5)-Dpoints(:,6))
    Dw(:) = 0.25_DP * (Dpoints(:,1)+Dpoints(:,4)+Dpoints(:,5)+Dpoints(:,8)&
                      -Dpoints(:,2)-Dpoints(:,3)-Dpoints(:,6)-Dpoints(:,7))

    ! Calculate normal n := u x v with 3d cross product
    Dn(1) = Du(2)*Dv(3) - Du(3)*Dv(2)
    Dn(2) = Du(3)*Dv(1) - Du(1)*Dv(3)
    Dn(3) = Du(1)*Dv(2) - Du(2)*Dv(1)

    ! Calculate scalar product t := < n, w >
    dt = Dn(1)*Dw(1) + Dn(2)*Dw(2) + Dn(3)*Dw(3)

    ! Now if dt < 0, then the hexahedron is flipped
    gaux_isFlipped_hexa3D = (dt .lt. 0.0_DP)

  end function gaux_isFlipped_hexa3D

!************************************************************************

!<subroutine>

  pure subroutine gaux_isInElement_tetra(dx,dy,dz,Dpoints,binside)

!<description>
  ! Checks if a point (dx,dy) is inside of a 2D triangular element
  ! given by the corners DcornerCoords.
!</description>

!<input>
  ! Point to check
  real(DP), intent(in) :: dx,dy,dz

  ! Array with coordinates of the four corner points of the element.
  ! The corners must be ordered in counterclockwise order.
  !
  ! Note: For performance reasons, this array is defined as
  !   explicit array of dimension (3,4). As this deactivates array
  !   checking in Fortran, the caller must take care to specify exactly
  !   this type of array here!
  real(DP), dimension(3,4), intent(in) :: Dpoints
!</input>

!<result>
  ! TRUE if (dx,dy) is inside of the element. FALSE otherwise.
  logical, intent(out) :: binside
!</result>

!</subroutine>
  real(DP), dimension(3,4) :: Dnormals
  real(DP), dimension(3,4) :: Dpv
  real(DP), dimension(3,1) :: DP1
  real(DP) :: ddot1,ddot2,ddot3,ddot4

  DP1(1,1)=dx
  DP1(2,1)=dy
  DP1(3,1)=dz

  ! compute the face normals
  Dnormals(1,1) = ((Dpoints(2,2) - Dpoints(2,4)) * (Dpoints(3,3) - Dpoints(3,4))) &
                - ((Dpoints(3,2) - Dpoints(3,4)) * (Dpoints(2,3) - Dpoints(2,4)))

  Dnormals(2,1) = ((Dpoints(3,2) - Dpoints(3,4)) * (Dpoints(1,3) - Dpoints(1,4))) &
                - ((Dpoints(1,2) - Dpoints(1,4)) * (Dpoints(3,3) - Dpoints(3,4)))

  Dnormals(3,1) = ((Dpoints(1,2) - Dpoints(1,4)) * (Dpoints(2,3) - Dpoints(2,4))) &
                - ((Dpoints(2,2) - Dpoints(2,4)) * (Dpoints(1,3) - Dpoints(1,4)))

  Dnormals(1,2) = ((Dpoints(2,1) - Dpoints(2,3)) * (Dpoints(3,4) - Dpoints(3,3))) &
                - ((Dpoints(3,1) - Dpoints(3,3)) * (Dpoints(2,4) - Dpoints(2,3)))

  Dnormals(2,2) = ((Dpoints(3,1) - Dpoints(3,3)) * (Dpoints(1,4) - Dpoints(1,3))) &
                - ((Dpoints(1,1) - Dpoints(1,3)) * (Dpoints(3,4) - Dpoints(3,3)))

  Dnormals(3,2) = ((Dpoints(1,1) - Dpoints(1,3)) * (Dpoints(2,4) - Dpoints(2,3))) &
                - ((Dpoints(2,1) - Dpoints(2,3)) * (Dpoints(1,4) - Dpoints(1,3)))

  Dnormals(1,3) = ((Dpoints(2,4) - Dpoints(2,2)) * (Dpoints(3,1) - Dpoints(3,2))) &
                - ((Dpoints(3,4) - Dpoints(3,2)) * (Dpoints(2,1) - Dpoints(2,2)))

  Dnormals(2,3) = ((Dpoints(3,4) - Dpoints(3,2)) * (Dpoints(1,1) - Dpoints(1,2))) &
                - ((Dpoints(1,4) - Dpoints(1,2)) * (Dpoints(3,1) - Dpoints(3,2)))

  Dnormals(3,3) = ((Dpoints(1,4) - Dpoints(1,2)) * (Dpoints(2,1) - Dpoints(2,2))) &
                - ((Dpoints(2,4) - Dpoints(2,2)) * (Dpoints(1,1) - Dpoints(1,2)))

  Dnormals(1,4) = ((Dpoints(2,3) - Dpoints(2,1)) * (Dpoints(3,2) - Dpoints(3,1))) &
                - ((Dpoints(3,3) - Dpoints(3,1)) * (Dpoints(2,2) - Dpoints(2,1)))

  Dnormals(2,4) = ((Dpoints(3,3) - Dpoints(3,1)) * (Dpoints(1,2) - Dpoints(1,1))) &
                - ((Dpoints(1,3) - Dpoints(1,1)) * (Dpoints(3,2) - Dpoints(3,1)))

  Dnormals(3,4) = ((Dpoints(1,3) - Dpoints(1,1)) * (Dpoints(2,2) - Dpoints(2,1))) &
                - ((Dpoints(2,3) - Dpoints(2,1)) * (Dpoints(1,2) - Dpoints(1,1)))

  ! calculate all p-v
  Dpv(1:3,1)= DP1(1:3,1) - Dpoints(1:3,4)

  Dpv(1:3,2)= DP1(1:3,1) - Dpoints(1:3,3)

  Dpv(1:3,3)= DP1(1:3,1) - Dpoints(1:3,2)

  Dpv(1:3,4)= DP1(1:3,1) - Dpoints(1:3,1)

  ddot1 = Dpv(1,1) * Dnormals(1,1) + Dpv(2,1) * Dnormals(2,1) + Dpv(3,1) * Dnormals(3,1)
  ddot2 = Dpv(1,2) * Dnormals(1,2) + Dpv(2,2) * Dnormals(2,2) + Dpv(3,2) * Dnormals(3,2)
  ddot3 = Dpv(1,3) * Dnormals(1,3) + Dpv(2,3) * Dnormals(2,3) + Dpv(3,3) * Dnormals(3,3)
  ddot4 = Dpv(1,4) * Dnormals(1,4) + Dpv(2,4) * Dnormals(2,4) + Dpv(3,4) * Dnormals(3,4)

  if((ddot1 .le. 0.0001_dp).and.(ddot2 .le. 0.0001_dp).and.(ddot3 .le. 0.0001_dp).and.(ddot4 .le. 0.0001_dp))then
    binside = .true.
  else
    binside = .false.
  end if

  end subroutine

!************************************************************************

!<subroutine>

  pure subroutine gaux_isInElement_hexa(dx,dy,dz,Dpoints,binside)

!<description>
  ! Checks whether a point (dx,dy,dz) is inside of a hexahedron
!</description>

!<input>
  ! Point to check
  real(DP), intent(in) :: dx,dy,dz

  ! Array with coordinates of the four corner points of the element.
  ! The corners must be ordered in counterclockwise order!
  !
  ! Note: For performance reasons, this array is defined as
  !   explicit array of dimension (3,8). As this deactivates array
  !   checking in Fortran, the caller must take care to specify exactly
  !   this type of array here!
  real(DP), dimension(3,8), intent(in) :: Dpoints
!</input>

!<result>
  ! TRUE if (dx,dy) is inside of the element. FALSE otherwise.
  logical, intent(out) :: binside
!</result>

!</subroutine>
  real(dp), dimension(3,6) :: dnormals
  real(DP), dimension(3,6) :: Dpv
  real(DP), dimension(3,1) :: DP1
  real(DP) :: ddot1,ddot2,ddot3,ddot4,ddot5,ddot6,length
  real(DP) :: eps
  integer :: i

  DP1(1,1)=dx
  DP1(2,1)=dy
  DP1(3,1)=dz

  eps=1e-10

  ! compute the face normals 1
  Dnormals(1,1) = ((Dpoints(2,4) - Dpoints(2,1)) * (Dpoints(3,2) - Dpoints(3,1))) &
                - ((Dpoints(3,4) - Dpoints(3,1)) * (Dpoints(2,2) - Dpoints(2,1)))

  Dnormals(2,1) = ((Dpoints(3,4) - Dpoints(3,1)) * (Dpoints(1,2) - Dpoints(1,1))) &
                - ((Dpoints(1,4) - Dpoints(1,1)) * (Dpoints(3,2) - Dpoints(3,1)))

  Dnormals(3,1) = ((Dpoints(1,4) - Dpoints(1,1)) * (Dpoints(2,2) - Dpoints(2,1))) &
                - ((Dpoints(2,4) - Dpoints(2,1)) * (Dpoints(1,2) - Dpoints(1,1)))

  ! compute the face normals 2
  Dnormals(1,2) = ((Dpoints(2,2) - Dpoints(2,1)) * (Dpoints(3,5) - Dpoints(3,1))) &
                - ((Dpoints(3,2) - Dpoints(3,1)) * (Dpoints(2,5) - Dpoints(2,1)))

  Dnormals(2,2) = ((Dpoints(3,2) - Dpoints(3,1)) * (Dpoints(1,5) - Dpoints(1,1))) &
                - ((Dpoints(1,2) - Dpoints(1,1)) * (Dpoints(3,5) - Dpoints(3,1)))

  Dnormals(3,2) = ((Dpoints(1,2) - Dpoints(1,1)) * (Dpoints(2,5) - Dpoints(2,1))) &
                - ((Dpoints(2,2) - Dpoints(2,1)) * (Dpoints(1,5) - Dpoints(1,1)))
  ! compute the face normals 3
  Dnormals(1,3) = ((Dpoints(2,7) - Dpoints(2,3)) * (Dpoints(3,2) - Dpoints(3,3))) &
                - ((Dpoints(3,7) - Dpoints(3,3)) * (Dpoints(2,2) - Dpoints(2,3)))

  Dnormals(2,3) = ((Dpoints(3,7) - Dpoints(3,3)) * (Dpoints(1,2) - Dpoints(1,3))) &
                - ((Dpoints(1,7) - Dpoints(1,3)) * (Dpoints(3,2) - Dpoints(3,3)))

  Dnormals(3,3) = ((Dpoints(1,7) - Dpoints(1,3)) * (Dpoints(2,2) - Dpoints(2,3))) &
                - ((Dpoints(2,7) - Dpoints(2,3)) * (Dpoints(1,2) - Dpoints(1,3)))
  ! compute the face normals 4
  Dnormals(1,4) = ((Dpoints(2,4) - Dpoints(2,3)) * (Dpoints(3,7) - Dpoints(3,3))) &
                - ((Dpoints(3,4) - Dpoints(3,3)) * (Dpoints(2,7) - Dpoints(2,3)))

  Dnormals(2,4) = ((Dpoints(3,4) - Dpoints(3,3)) * (Dpoints(1,7) - Dpoints(1,3))) &
                - ((Dpoints(1,4) - Dpoints(1,3)) * (Dpoints(3,7) - Dpoints(3,3)))

  Dnormals(3,4) = ((Dpoints(1,4) - Dpoints(1,3)) * (Dpoints(2,7) - Dpoints(2,3))) &
                - ((Dpoints(2,4) - Dpoints(2,3)) * (Dpoints(1,7) - Dpoints(1,3)))
  ! compute the face normals 5
  Dnormals(1,5) = ((Dpoints(2,4) - Dpoints(2,8)) * (Dpoints(3,5) - Dpoints(3,8))) &
                - ((Dpoints(3,4) - Dpoints(3,8)) * (Dpoints(2,5) - Dpoints(2,8)))

  Dnormals(2,5) = ((Dpoints(3,4) - Dpoints(3,8)) * (Dpoints(1,5) - Dpoints(1,8))) &
                - ((Dpoints(1,4) - Dpoints(1,8)) * (Dpoints(3,5) - Dpoints(3,8)))

  Dnormals(3,5) = ((Dpoints(1,4) - Dpoints(1,8)) * (Dpoints(2,5) - Dpoints(2,8))) &
                - ((Dpoints(2,4) - Dpoints(2,8)) * (Dpoints(1,5) - Dpoints(1,8)))

  ! compute the face normals 6
  Dnormals(1,6) = ((Dpoints(2,7) - Dpoints(2,6)) * (Dpoints(3,5) - Dpoints(3,6))) &
                - ((Dpoints(3,7) - Dpoints(3,6)) * (Dpoints(2,5) - Dpoints(2,6)))

  Dnormals(2,6) = ((Dpoints(3,7) - Dpoints(3,6)) * (Dpoints(1,5) - Dpoints(1,6))) &
                - ((Dpoints(1,7) - Dpoints(1,6)) * (Dpoints(3,5) - Dpoints(3,6)))

  Dnormals(3,6) = ((Dpoints(1,7) - Dpoints(1,6)) * (Dpoints(2,5) - Dpoints(2,6))) &
                - ((Dpoints(2,7) - Dpoints(2,6)) * (Dpoints(1,5) - Dpoints(1,6)))


  do i=1,6
    length=sqrt(Dnormals(1,i)**2+Dnormals(2,i)**2+Dnormals(3,i)**2)
    Dnormals(1,i)=Dnormals(1,i)/length
    Dnormals(2,i)=Dnormals(2,i)/length
    Dnormals(3,i)=Dnormals(3,i)/length
  end do


  ! calculate all p-v
  Dpv(1:3,1)= DP1(1:3,1) - Dpoints(1:3,1)

  Dpv(1:3,2)= DP1(1:3,1) - Dpoints(1:3,2)

  Dpv(1:3,3)= DP1(1:3,1) - Dpoints(1:3,3)

  Dpv(1:3,4)= DP1(1:3,1) - Dpoints(1:3,4)

  Dpv(1:3,5)= DP1(1:3,1) - Dpoints(1:3,5)

  Dpv(1:3,6)= DP1(1:3,1) - Dpoints(1:3,6)

  ddot1 = Dpv(1,1) * Dnormals(1,1) + Dpv(2,1) * Dnormals(2,1) + Dpv(3,1) * Dnormals(3,1)
  ddot2 = Dpv(1,2) * Dnormals(1,2) + Dpv(2,2) * Dnormals(2,2) + Dpv(3,2) * Dnormals(3,2)
  ddot3 = Dpv(1,3) * Dnormals(1,3) + Dpv(2,3) * Dnormals(2,3) + Dpv(3,3) * Dnormals(3,3)
  ddot4 = Dpv(1,4) * Dnormals(1,4) + Dpv(2,4) * Dnormals(2,4) + Dpv(3,4) * Dnormals(3,4)
  ddot5 = Dpv(1,5) * Dnormals(1,5) + Dpv(2,5) * Dnormals(2,5) + Dpv(3,5) * Dnormals(3,5)
  ddot6 = Dpv(1,6) * Dnormals(1,6) + Dpv(2,6) * Dnormals(2,6) + Dpv(3,6) * Dnormals(3,6)

  if((ddot1 .le. eps).and.(ddot2 .le. eps).and.(ddot3 .le. eps).and.&
     (ddot4 .le. eps).and.(ddot5 .le. eps).and.(ddot6 .le. eps))then
    binside = .true.
  else
    binside = .false.
  end if

  end subroutine

!************************************************************************

!<subroutine>

  pure subroutine gaux_isIntersection_face(&
      Dpoint1,Dpoint2,Dface, bintersect)

!<description>
! This routine calculates whether there is
! an intersection between the face Dface and the
! line defined by the vertices Dpoint1 and Dface
! the face is split into two triangles, the
! triangles are tested for intersection
!</description>

!<input>
  ! First point on ray
  real(DP), dimension(3), intent(in) :: Dpoint1

  ! second point on ray
  real(DP), dimension(3), intent(in) :: Dpoint2

  ! Vertices of the face
  real(DP), dimension(3,4), intent(in) :: Dface

  ! TRUE if the two rays intersect. FALSE otherwise.
  logical, intent(out) :: bintersect
!</result>

!</subroutine>

    ! we split up the face into two triangles
    real(DP), dimension(3,3) :: Dtri1
    real(DP), dimension(3,3) :: Dtri2

    ! build the first triangle
    Dtri1(1:3,1)=Dface(1:3,1)
    Dtri1(1:3,2)=Dface(1:3,2)
    Dtri1(1:3,3)=Dface(1:3,4)

    ! build the second triangle
    Dtri2(1:3,1)=Dface(1:3,2)
    Dtri2(1:3,2)=Dface(1:3,3)
    Dtri2(1:3,3)=Dface(1:3,4)

    ! test for intersection
    call gaux_isIntersection_triangle(Dpoint1,Dpoint2,Dtri1,bintersect)

    if (bintersect) then
      return
    end if

    ! test for intersection
    call gaux_isIntersection_triangle(Dpoint1,Dpoint2,Dtri2,bintersect)

    if (bintersect) then
      return
    end if

    ! no intersection
    bintersect = .false.

  end subroutine

!************************************************************************

!<subroutine>

  pure subroutine gaux_isIntersection_triangle(&
      Dpoint1,Dpoint2,Dtri, bintersect)

!<description>
! This routine calculates whether there is
! an intersection between the triangle Dtri and
! the line segment between Dpoint1 and Dpoint2
!</description>

!<input>
  ! start of the ray
  real(DP), dimension(3), intent(in) :: Dpoint1

  ! 2nd point of the ray
  real(DP), dimension(3), intent(in) :: Dpoint2

  ! Vertices of the face
  real(DP), dimension(3,3), intent(in) :: Dtri

  ! TRUE if the two rays intersect. FALSE otherwise.
  logical, intent(out) :: bintersect
!</result>

!</subroutine>
  real(dp) :: eps, dot, dot2, dot3, u, v, t
  real(dp), dimension(3) :: De1, De2, p, s, q, Ddirect

  ! our epsilon when we test for zero
  eps = 0.000001_dp
  ! the algorithm is as follows, it is called the Moeller & Trumbore algorithm
  !
  ! A * B with
  !
  ! A := 1/dot(Ddirect x (Dtri(:,2)-Dtri(:,0)),(Dtri(:,1)-Dtri(:,0)))
  !
  !      |dot((Dpoint2-Dtri(:,0)) x (Dtri(:,1)-Dtri(:,0)),(Dtri(:,2)-Dtri(:,0)))|
  ! B := |dot(Ddirect x (Dtri(:,2)-Dtri(:,0))            ,(Dpoint2-Dtri(:,0)))  |
  !      |dot((Dpoint2-Dtri(:,0)) x (Dtri(:,1)-Dtri(:,0)),Ddirect))             |

  Ddirect(1:3) = Dpoint2(1:3) - Dpoint1(1:3)

  De1(1:3) = Dtri(1:3,2) - Dtri(1:3,1)
  De2(1:3) = Dtri(1:3,3) - Dtri(1:3,1)

  p(1) = Ddirect(2)*De2(3) - Ddirect(3)*De2(2)
  p(2) = Ddirect(3)*De2(1) - Ddirect(1)*De2(3)
  p(3) = Ddirect(1)*De2(2) - Ddirect(2)*De2(1)

  dot = p(1) * De1(1) + p(2) * De1(2) + p(3) * De1(3)

  if(dot > -eps .and. dot < eps)then
    bintersect = .false.
    return
  end if

  dot = 1.0_dp/dot

  s(:) = Dpoint1(:) - Dtri(:,1)

  dot2 =  s(1) * p(1) + s(2) * p(2) + s(3) * p(3)

  u = dot * dot2

  if(u < 0.0_dp .or. u > 1.0_dp)then
    bintersect = .false.
    return
  end if

  q(1) = s(2)*De1(3) - s(3)*De1(2)
  q(2) = s(3)*De1(1) - s(1)*De1(3)
  q(3) = s(1)*De1(2) - s(2)*De1(1)

  dot3 =  q(1) * Ddirect(1) + q(2) * Ddirect(2) + q(3) * Ddirect(3)

  v = dot * dot3

  if(v < 0.0_dp .or. v > 1.0_dp)then
    bintersect = .false.
    return
  end if

  t = dot * (De2(1) * q(1) + De2(2) * q(2) + De2(3) * q(3))

  if(t > 0.0_dp .and. t < 1.0_dp)then
    bintersect = .true.
  end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine gaux_calcDistPEdg2D(DpointA,Dedge1,Dedge2,dist,t)

!<description>
    ! Calculated the distance between a point and an edge in 2d
!</description>

!<input>
  !
  ! The coordinates of the point
  ! DpointA(1) = x-coordinates,
  ! DpointA(2) = y-coordinates
  real(DP), dimension(2), intent(in) :: DpointA
  ! the first vertex of the edge in real coordinates
  real(DP), dimension(2), intent(in) :: Dedge1
  ! the 2nd vertex of the edge in real coordinates
  real(DP), dimension(2), intent(in) :: Dedge2
!</input>
  real(dp), intent(out)              :: dist
  real(dp), intent(out)              :: t
!<result>
  ! The distance
!</result>
!</subroutine>
  ! parameter in the edge equation: edge=Dedge1+t*(Dedge2-Dedge1), t in [0,1]
  real(dp) :: r2
  real(DP), dimension(2) :: r,YP

  r(:)=Dedge2(:)-Dedge1(:)

  r2=r(1)**2+r(2)**2
  YP(:)=DpointA(:)-Dedge1(:)
  t=r(1)*YP(1)+r(2)*YP(2)
  t=t/r2

  if(t .le. 0.0_dp)then
    dist = sqrt(YP(1)**2 + YP(2)**2)
    return
  else if(t .ge. 1.0_dp)then
    YP(:)=DpointA(:)-Dedge2(:)
    dist = sqrt(YP(1)**2 + YP(2)**2)
  return
  else
    r(:)=r(:)*t
    r(:)=r(:)+Dedge1(:)
    YP(:)=DpointA(:)-r(:)
    dist = sqrt(YP(1)**2 + YP(2)**2)
    return
  end if

  end subroutine gaux_calcDistPEdg2D

!************************************************************************

!<subroutine>

 pure subroutine gaux_projectPointPlane(DPoint,DR1,DR2,DQ,DQP)
!<description>
  ! This routine projects a point onto a plane.
  ! the plane must be entered in parameterform that means:E= p + alpha*r1 + beta*r2
  ! the code will transform this to point normal form and calculate the result
!</description>

!<input>
  ! the p in the equation for E above
  real(DP), dimension(3), intent(in) :: DPoint
  ! the r1 in the equation for E above
  real(DP), dimension(3), intent(in) :: DR1
  ! the r2 in the equation for E above
  real(DP), dimension(3), intent(in) :: DR2
  ! the point that should be projected
!</input>
!<inputoutput>
  real(DP), dimension(3), intent(inout) :: DQ
  ! the projection of DQ onto the plane
  real(DP), dimension(3), intent(inout) :: DQP
!</inputoutput>
!</subroutine>
  ! local variables
  real(DP) :: d,n,q
  real(DP), dimension(3) :: DN, DNscaled

  DN(1)=DR1(2)*DR2(3)-DR1(3)*DR2(2)
  DN(2)=DR1(3)*DR2(1)-DR1(1)*DR2(3)
  DN(3)=DR1(1)*DR2(2)-DR1(2)*DR2(1)

  n=sqrt(DN(1)**2 + DN(2)**2 + DN(3)**2)
  DN(:)=DN(:)/n

  d = -1.0_dp * (DN(1)*DPoint(1)+DN(2)*DPoint(2)+DN(3)*DPoint(3))

  q = (DQ(1)*DN(1)+DQ(2)*DN(2)+DQ(3)*DN(3))+d

  DNscaled(:)=q*DN(:)

  DQP(:)=DQ(:)-DNscaled(:)

  end subroutine

!************************************************************************

!<subroutine>

  pure subroutine gaux_isInElement_hexa_aligned(dx,dy,dz,Dpoints,binside)

!<description>
  ! Checks whether a point (dx,dy,dz) is inside of a hexahedron
  ! The subroutine is faster than the general subroutine
  ! but it does not give correct results if the sides of the hexa
  ! are not aligned with the coordinate axes
!</description>

!<input>
  ! Point to check
  real(DP), intent(in) :: dx,dy,dz

  ! Array with coordinates of the four corner points of the element.
  ! The corners must be ordered in counterclockwise order!
  !
  ! Note: For performance reasons, this array is defined as
  !   explicit array of dimension (3,8). As this deactivates array
  !   checking in Fortran, the caller must take care to specify exactly
  !   this type of array here!
  real(DP), dimension(3,8), intent(in) :: Dpoints
!</input>

!<result>
  ! TRUE if (dx,dy) is inside of the element. FALSE otherwise.
  logical, intent(out) :: binside
!</result>

!</subroutine>
  real(DP), dimension(3,1) :: DP1
  real(DP) :: xmin,xmax,ymin,ymax,zmin,zmax,eps
  integer :: i

  DP1(1,1)=dx
  DP1(2,1)=dy
  DP1(3,1)=dz

  eps=1e-10

  xmin=1e+32
  xmax=-1e+32
  ymin=1e+32
  ymax=-1e+32
  zmin=1e+32
  zmax=-1e+32

  do i=1,8

    if(Dpoints(1,i).lt.xmin)then
      xmin=Dpoints(1,i)
    end if
    if(Dpoints(1,i).gt.xmax)then
      xmax=Dpoints(1,i)
    end if

    if(Dpoints(2,i).lt.ymin)then
      ymin=Dpoints(2,i)
    end if
    if(Dpoints(2,i).gt.ymax)then
      ymax=Dpoints(2,i)
    end if

    if(Dpoints(3,i).lt.zmin)then
      zmin=Dpoints(3,i)
    end if
    if(Dpoints(3,i).gt.zmax)then
      zmax=Dpoints(3,i)
    end if

  end do

  if((dx .ge. xmin).and.(dx .le. xmax).and. &
     (dy .ge. ymin).and.(dy .le. ymax).and. &
     (dz .ge. zmin).and.(dz .le. zmax))then
     binside=.true.
  else
     binside=.false.
  end if

  end subroutine

!************************************************************************

!<subroutine>

  pure subroutine gaux_intersect_quadcircle(drad,Dc,Dquad,Iresult,Dintersection)

!<description>
  ! This routine calculates the number of intersections
  ! of a circle with a quad and returns the intersection points
!</description>

!<input>
  ! Radius of the circle
  real(DP), intent(in) :: drad
  ! the center of the circle
  real(DP), dimension(2),intent(in) :: Dc

  ! The four corner vertices of the quad
  real(DP), dimension(2,4), intent(in) :: Dquad
!</input>

!<result>
  ! There are at maximum 8 intersections between a circle and a quad
  ! if the quad is
  !
  !   1----2
  !   |    |
  !   4----3
  ! Iresult(1)=2 if there are 2 intersection on edge 12
  ! Iresult(1)=1 if there is 1 intersection on edge 12
  ! Iresult(1)=0 if there is no intersection on edge 12
  ! Iresult(2)=2 if there are 2 intersection on edge 23
  ! ... and so on
  integer, dimension(4), intent(out) :: Iresult
  ! the possible 8 points of intersection
  ! if Iresult(1)=2 then Dintersection(:,1-2) stores the
  ! 2 points of intersection on edge 12
  ! if Iresult(1)=1 then Dintersection(:,1) stores the
  ! point of intersection on edge 12 and Dintersection(:,2) is empty
  ! ... and so on
  real(DP), dimension(2,8), intent(out) :: Dintersection
!</result>

!</subroutine>
  real(DP), dimension(2,2) :: Dedge
  real(DP), dimension(2,2) :: Dintersec
  integer :: i,iintersections

  Dintersection=0
  Iresult=0
  ! loop over all edges of the quad
  do i=1,4
    iintersections=0
    Dintersec=0
    Dedge(:,1) = Dquad(:,i)
    Dedge(:,2) = Dquad(:,mod(i,4)+1)
    call gaux_intersect_edgecircle(drad,Dc,Dedge,iintersections,Dintersec)
    Iresult(i)=iintersections
    Dintersection(:,2*i-1)=Dintersec(:,1)
    Dintersection(:,2*i)=Dintersec(:,2)
  end do

  end subroutine

!************************************************************************

!<subroutine>

  pure subroutine gaux_intersect_edgecircle(drad,Dc,Dedge,iresult,Dintersec)

!<description>
  ! computes the number and the intersection points of a circle and an edge
!</description>

!<input>
  ! radius of the circle
  real(DP), intent(in) :: drad
  ! Center of the circle
  real(DP), dimension(2),intent(in) :: Dc
  ! The two points of the edge
  real(DP), dimension(2,2), intent(in) :: Dedge
!</input>

!<result>
  ! The coordinates of the intersection points
  real(DP), dimension(2,2),intent(inout) :: Dintersec
  ! the number of intersection points on the edge
  integer, intent(out) :: iresult
!</result>

!</subroutine>
  real(dp), dimension(2) :: Ddir,Dp1
  real(DP), dimension(2) :: Ddelta
  real(DP) :: t,root,diff,dir2,eps,dirdotdelta,dt1,dt2
  integer :: idresult

  eps=1e-10

  ! initialise by 0
  idresult = 0
  Dintersec=0

  ! get direction vector of the edge
  Ddir(:)=Dedge(:,2)-Dedge(:,1)

  Dp1(:)=Dedge(:,1)

  Ddelta(:)=Dp1(:)-Dc(:)

  dir2=Ddir(1)**2+Ddir(2)**2

  diff=(Ddelta(1)**2+Ddelta(2)**2)-drad**2

  dirdotdelta = Ddir(1)*Ddelta(1)+Ddir(2)*Ddelta(2)

  root = dirdotdelta**2 - dir2*diff

  ! the term under the root is negative, there is no intersection
  if(root < 0.0_dp)then
    return
  end if

  ! the term under the root is equal to 0, we have one intersection 't'
  if((0.0 .le. root).and.(root .le. eps))then
    root=0.0_dp
    t= -Ddir(1)*Ddelta(1)-Ddir(2)*Ddelta(2) /dir2
    ! we need point of intersection to be on the segment
    if((0.0_dp .le. t).and.(t .le. 1.0_dp))then
      iresult = 1
      Dintersec(:,1)=Dp1(:)+ t*Ddir(:)
    end if
  ! the term under the root is greater than zero
  ! we have to intersections dt1,dt2
  else
    dt1=(-dirdotdelta + root)/dir2
    dt2=(-dirdotdelta - root)/dir2
    if((0.0_dp .le. dt1).and.(dt1 .le. 1.0_dp))then
      iresult = iresult + 1
      Dintersec(:,1)=Dp1(:)+ dt1*Ddir(:)
    end if
    if((0.0_dp .le. dt2).and.(dt2 .le. 1.0_dp))then
      iresult = iresult + 1
      Dintersec(:,2)=Dp1(:)+ dt2*Ddir(:)
    end if
  end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine gaux_getDirectedExtentQuad2D (Dcorners,dbeta1,dbeta2,dlen)
      
  !<description>
    ! Calculates the maximum length of a line in direction (dbeta1,dbeta2)
    ! when embedded in a 2D quadrilateral with corners Dcoords.
    !
    ! Remarks: The cell is assumed to be convex. The routine "shoots" a
    ! ray from each corner in direction (dbeta1,dbeta2) and detects
    ! where it hits the opposite side of the cell.
  !</description>

  !<input>
    ! Coordinates of the corners defining the cell.
    real(DP), dimension(:,:), intent(in) :: Dcorners
    
    ! Direction of the line. It is assumed that the vector is normalised,
    ! i.e., ||(dbeta1,dbeta2)||_2 = 1.
    real(DP), intent(in) :: dbeta1, dbeta2
    
  !</input>

  !<output>
    ! Maximum length of the line.
    real(DP), intent(out) :: dlen
  !</output>

!</subroutine>

    ! local variables
    real(DP) :: x1,y1,x2,y2,x3,y3,x4,y4, dx,dy, dalpha
    integer :: iintersect

    ! Fetch the coordinates of these corners

    x1=Dcorners(1,1)
    y1=Dcorners(2,1)
    x2=Dcorners(1,2)
    y2=Dcorners(2,2)
    x3=Dcorners(1,3)
    y3=Dcorners(2,3)
    x4=Dcorners(1,4)
    y4=Dcorners(2,4)

    dlen=0.0_DP

    ! Calculate the `maximum possible cell width in direction (dbeta1,dbeta2)`;
    ! The picture in mind is the following:
    !
    !          G3
    !   +-------------X-------+
    !   |            /        |
    !   |           /         |
    !   |          /          |
    !   |         /           |
    !   |        /            |
    ! G4|       /             | G2
    !   |      ^ (beta1,beta2)|
    !   |     /               |
    !   |    /                |
    !   |   /                 |
    !   |  /                  |
    !   | /                   |
    !   |/                    |
    !   O---------------------+
    !            G1
    !
    ! The vector (beta1,beta2) indicates the direction of the flow.
    ! A particle starting in point O moves at most up to point X.
    ! The length of the line (O,X) is returned.
    !
    ! Loop through the four corners of cell and check for a line
    ! with slope BETA=(xbeta1,xbeta2) starting in this corner whether it
    ! really intersects with one of the edges of the element. Remark
    ! that we only have to check the two opposite edges to the current
    ! corner!

    ! -----------------------------------------------------------------
    ! Check the first corner.
    ! We obtain dalpha, which may be negative if direction (dbeta1,dbeta2)
    ! points outside of the element and we would have to use
    ! (-dbeta1,-dbeta2) instead.

    call gaux_getIntersection_ray2D(&
        x1,y1, x1+dbeta1,y1+dbeta2, x2,y2, x3,y3, dx,dy, iintersect, dalpha)
    dlen=max(abs(dalpha),dlen)

    call gaux_getIntersection_ray2D(&
        x1,y1, x1+dbeta1,y1+dbeta2, x3,y3, x4,y4, dx,dy, iintersect, dalpha)

    dlen=max(abs(dalpha),dlen)

    ! -----------------------------------------------------------------
    ! The second one...

    call gaux_getIntersection_ray2D(&
        x2,y2, x2+dbeta1,y2+dbeta2, x4,y4, x1,y1, dx,dy, iintersect, dalpha)
    dlen=max(abs(dalpha),dlen)

    call gaux_getIntersection_ray2D(&
        x2,y2, x2+dbeta1,y2+dbeta2, x4,y4, x3,y3, dx,dy, iintersect, dalpha)
    dlen=max(abs(dalpha),dlen)

    ! -----------------------------------------------------------------
    ! The third one...

    call gaux_getIntersection_ray2D(&
        x3,y3, x3+dbeta1,y3+dbeta2, x1,y1, x2,y2, dx,dy, iintersect, dalpha)
    dlen=max(abs(dalpha),dlen)

    call gaux_getIntersection_ray2D(&
        x3,y3, x3+dbeta1,y3+dbeta2, x1,y1, x4,y4, dx,dy, iintersect, dalpha)
    dlen=max(abs(dalpha),dlen)

    ! -----------------------------------------------------------------
    ! And the fourth=last one...

    call gaux_getIntersection_ray2D(&
        x4,y4, x4+dbeta1,y4+dbeta2, x2,y2, x1,y1, dx,dy, iintersect, dalpha)
    dlen=max(abs(dalpha),dlen)

    call gaux_getIntersection_ray2D(&
        x4,y4, x4+dbeta1,y4+dbeta2, x2,y2, x3,y3, dx,dy, iintersect, dalpha)
    dlen=max(abs(dalpha),dlen)

  end subroutine

end module

