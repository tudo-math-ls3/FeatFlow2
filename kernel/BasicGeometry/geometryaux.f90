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
!# </purpose>
!##############################################################################

module geometryaux

  use fsystem
  
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

!<subroutine>
  
  pure subroutine gaux_isInElement_quad2D(dx,dy,DcornerCoords,binside)
  
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
  real(DP), dimension(2,4), intent(in) :: DcornerCoords
!</input>

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
      binside = binside .and. (dsproduct .le. SYS_EPSREAL*100.0_DP)

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
  real(DP) :: ddot1,ddot2,ddot3,ddot4,ddot5,ddot6
  
  DP1(1,1)=dx
  DP1(2,1)=dy
  DP1(3,1)=dz
  
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
  
  if((ddot1 .le. 0.0001_dp).and.(ddot2 .le. 0.0001_dp).and.(ddot3 .le. 0.0001_dp).and.&
     (ddot4 .le. 0.0001_dp).and.(ddot5 .le. 0.0001_dp).and.(ddot6 .le. 0.0001_dp))then
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

    ! local variables: aux parameters
    real(DP) :: daux1, daux2, daux3, daux4
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
  !                              1                                |dot((Dpoint2-Dtri(:,0)) x (Dtri(:,1)-Dtri(:,0)),(Dtri(:,2)-Dtri(:,0)))|
  !----------------------------------------------------------- *  |dot(Ddirect x (Dtri(:,2)-Dtri(:,0)),(Dpoint2-Dtri(:,0)))              |
  ! dot(Ddirect x (Dtri(:,2)-Dtri(:,0)),(Dtri(:,1)-Dtri(:,0)))    |dot((Dpoint2-Dtri(:,0)) x (Dtri(:,1)-Dtri(:,0)),Ddirect))             |
  !                                                               
  
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
  
  s(:) = Dpoint(:) - Dtri(:,1)
  
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

end module
