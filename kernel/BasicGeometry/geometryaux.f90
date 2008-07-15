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
!# </purpose>
!##############################################################################

MODULE geometryaux

  USE fsystem
  
  IMPLICIT NONE

CONTAINS

  ! ***************************************************************************
  
!<function>

  PURE REAL(DP) FUNCTION gaux_getAspectRatio_quad2D (Dpoints)
  
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
  REAL(DP), DIMENSION(2,4), INTENT(IN) :: Dpoints
!</input>

!<result>
  ! The aspect ratio of the polygon.
!</result>

!</function>

    REAL(DP) :: XM1,XM2,XM3,XM4,YM1,YM2,YM3,YM4,LE1,LE2

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
    gaux_getAspectRatio_quad2D = SQRT(LE1/MAX(LE2,1E-10_DP))

  END FUNCTION

  ! ***************************************************************************

!<function>
  
  PURE REAL(DP) FUNCTION gaux_getArea_tria2D (Dpoints)

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
    REAL(DP), DIMENSION(2,3), INTENT(IN) :: Dpoints
!</input>

!<result>
  ! The signed area of the polygon.
!</result>
!</function>

    gaux_getArea_tria2D = 0.5_DP*( &
        (Dpoints(1,2)-Dpoints(1,1))*(Dpoints(2,3)-Dpoints(2,1))-&
        (Dpoints(1,3)-Dpoints(1,1))*(Dpoints(2,2)-Dpoints(2,1)) )
  END FUNCTION gaux_getArea_tria2D

  ! ***************************************************************************

!<function>

  REAL(DP) FUNCTION gaux_getArea_quad2D (Dpoints)

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
    REAL(DP), DIMENSION(2,4), INTENT(IN) :: Dpoints
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
    ! "In modern linear algebra, as already noted, the area of a planar parallelogram 
    !  is the magnitude of the cross product of two adjacent edge vectors.  So, for 
    !  any 3D planar parallelogram V0V1V2V3, we have:"
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
    ! "Next, for an arbitrary quadrilateral, one can compute its area using a 
    !  parallelogram discovered by Pierre Varignon (first published in 1731). 
    !  It is amazing that the Greeks missed Varignon's simple result which was 
    !  discovered 2000 years after Euclid!  Given any quadrilateral, one can take 
    !  the midpoints of its 4 edges to get 4 vertices which form a new quadrilateral.  
    !  It is then easy to show that this midpoint quadrilateral is always a 
    !  parallelogram, called the "Varignon parallelogram", and that its area is 
    !  exactly one-half the area of the original quadrilateral.  So, for a 
    !  quadrilateral Q=V0V1V2V3, let this parallelogram have midpoint vertices 
    !  M0M1M2M3 as shown in the diagram:"
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
    ! "From elementary geometry, we know that in triangle V0V1V2 the midpoint 
    !  line M0M1 is parallel to the base V0V2.  In triangle V0V1V2V3, the line M3M2 
    !  is parallel to that same base V0V2.  Thus, M0M1 and M3M2 are parallel to each other.  
    !  Similarly, M0M3 and M1M2 are parallel, which shows that M0M1M2M3 is a parallelogram.  
    !  The area relation is also easy to demonstrate, and we can compute the quadrilateral's 
    !  area as:"
    !
    !     A(V0V1V2V3) = 2 A(M0M1M2M3)
    !                 = 2 | (M1-M0) X (M3-M0) |     (cross product)
    !                 = ...
    !                 = 1/2 | (V2-V0) X (V3-V1) |
    !                 = ...
    !                 = (x2-x0)(y3-y1) - (x3-x1)(y2-y0) 
    !
    ! "This formula for an arbitrary quadrilateral is just as efficient as the one for an 
    !  arbitrary triangle, using only 2 multiplications and 5 additions.  For simple 
    !  quadrilaterals, the area is positive when the vertices are oriented counterclockwise, 
    !  and negative when they are clockwise.  However, it also works for nonsimple 
    !  quadrilaterals and is equal to the difference in area of the two regions the 
    !  quadrilateral bounds.  For example, in the following diagram where I is the 
    !  self-intersection point of a nonsimple quadrilateral V0V1V2V3, we have:"
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
    !             /   \          
    !           /       \        
    !         /           \      
    !       +--------------+
    !     V0                 V1
    
    gaux_getArea_quad2D = &
            0.5_DP *  ABS( (Dpoints(1,3) - Dpoints(1,1) ) * &
                           (Dpoints(2,4) - Dpoints(2,2) ) - &
                           (Dpoints(1,4) - Dpoints(1,2) ) * &
                           (Dpoints(2,3) - Dpoints(2,1) ) )
                           
  END FUNCTION gaux_getArea_quad2D

!************************************************************************

!<function>

  PURE REAL(DP) FUNCTION gaux_getVolume_tetra3D (Dpoints)

!<description>
  ! This function calculates the volume of a 3D tetrahedron. The
  ! tetrahedron is given by the coordinates of its four corners.
!</description>

!<input>
  ! The coordinates of the four corners of the tetrahedron.
  ! Dpoints(1,.) = x-coordinates,
  ! Dpoints(2,.) = y-coordinates,
  ! Dpoints(3,.) = z-coordinates
  REAL(DP), DIMENSION(3,4), INTENT(IN) :: Dpoints
!</input>

!<result>
  ! The volume of the tetrahedron.
!</result>
!</function>

    ! A temporary array for the edge lengths
    REAL(DP), DIMENSION(3,3) :: Dv
    
    Dv(1:3,1) = Dpoints(1:3,1) - Dpoints(1:3,4)
    Dv(1:3,2) = Dpoints(1:3,2) - Dpoints(1:3,4)
    Dv(1:3,3) = Dpoints(1:3,3) - Dpoints(1:3,4)
    
    ! Return the absolute volume
    gaux_getVolume_tetra3D = ABS(&
        Dv(1,1) * (Dv(2,2)*Dv(3,3) - Dv(3,2)*Dv(2,3)) + &
        Dv(2,1) * (Dv(3,2)*Dv(1,3) - Dv(1,2)*Dv(3,3)) + &
        Dv(3,1) * (Dv(1,2)*Dv(2,3) - Dv(2,2)*Dv(1,3))) / 6.0_DP
        
  END FUNCTION gaux_getVolume_tetra3D

!************************************************************************

!<function>

  PURE REAL(DP) FUNCTION gaux_getVolume_hexa3D (Dv)

!<description>
  ! This function calculates the volume of a 3D hexahedron. The
  ! hexahedron is given by the coordinates of its eight corners.
!</description>

!<input>
  ! The coordinates of the eight corners of the hexahedron.
  ! Dv(1,.) = x-coordinates,
  ! Dv(2,.) = y-coordinates,
  ! Dv(3,.) = z-coordinates
  REAL(DP), DIMENSION(3,8), INTENT(IN) :: Dv
!</input>

!<result>
  ! The volume of the hexahedron.
!</result>
!</function>

    ! Return the absolute volume
    gaux_getVolume_hexa3D = (1.0_DP / 6.0_DP) * (&
        ABS((Dv(1,4)-Dv(1,1))*(Dv(2,4)-Dv(2,3))*(Dv(3,4)-Dv(3,8))&
           +(Dv(2,4)-Dv(2,1))*(Dv(3,4)-Dv(3,3))*(Dv(1,4)-Dv(1,8))&
           +(Dv(3,4)-Dv(3,1))*(Dv(1,4)-Dv(1,3))*(Dv(2,4)-Dv(2,8))&
           -(Dv(1,4)-Dv(1,8))*(Dv(2,4)-Dv(2,3))*(Dv(3,4)-Dv(3,1))&
           -(Dv(2,4)-Dv(2,8))*(Dv(3,4)-Dv(3,3))*(Dv(1,4)-Dv(1,1))&
           -(Dv(3,4)-Dv(3,8))*(Dv(1,4)-Dv(1,3))*(Dv(2,4)-Dv(2,1)))+&
        ABS((Dv(1,2)-Dv(1,3))*(Dv(2,2)-Dv(2,1))*(Dv(3,2)-Dv(3,6))&
           +(Dv(2,2)-Dv(2,3))*(Dv(3,2)-Dv(3,1))*(Dv(1,2)-Dv(1,6))&
           +(Dv(3,2)-Dv(3,3))*(Dv(1,2)-Dv(1,1))*(Dv(2,2)-Dv(2,6))&
           -(Dv(1,2)-Dv(1,6))*(Dv(2,2)-Dv(2,1))*(Dv(3,2)-Dv(3,3))&
           -(Dv(2,2)-Dv(2,6))*(Dv(3,2)-Dv(3,1))*(Dv(1,2)-Dv(1,3))&
           -(Dv(3,2)-Dv(3,6))*(Dv(1,2)-Dv(1,1))*(Dv(2,2)-Dv(2,3)))+&
        ABS((Dv(1,5)-Dv(1,8))*(Dv(2,5)-Dv(2,6))*(Dv(3,5)-Dv(3,1))&
           +(Dv(2,5)-Dv(2,8))*(Dv(3,5)-Dv(3,6))*(Dv(1,5)-Dv(1,1))&
           +(Dv(3,5)-Dv(3,8))*(Dv(1,5)-Dv(1,6))*(Dv(2,5)-Dv(2,1))&
           -(Dv(1,5)-Dv(1,1))*(Dv(2,5)-Dv(2,6))*(Dv(3,5)-Dv(3,8))&
           -(Dv(2,5)-Dv(2,1))*(Dv(3,5)-Dv(3,6))*(Dv(1,5)-Dv(1,8))&
           -(Dv(3,5)-Dv(3,1))*(Dv(1,5)-Dv(1,6))*(Dv(2,5)-Dv(2,8)))+&
        ABS((Dv(1,7)-Dv(1,6))*(Dv(2,7)-Dv(2,8))*(Dv(3,7)-Dv(3,3))&
           +(Dv(2,7)-Dv(2,6))*(Dv(3,7)-Dv(3,8))*(Dv(1,7)-Dv(1,3))&
           +(Dv(3,7)-Dv(3,6))*(Dv(1,7)-Dv(1,8))*(Dv(2,7)-Dv(2,3))&
           -(Dv(1,7)-Dv(1,3))*(Dv(2,7)-Dv(2,8))*(Dv(3,7)-Dv(3,6))&
           -(Dv(2,7)-Dv(2,3))*(Dv(3,7)-Dv(3,8))*(Dv(1,7)-Dv(1,6))&
           -(Dv(3,7)-Dv(3,3))*(Dv(1,7)-Dv(1,8))*(Dv(2,7)-Dv(2,6)))+&
        ABS((Dv(1,1)-Dv(1,3))*(Dv(2,1)-Dv(2,8))*(Dv(3,1)-Dv(3,6))&
           +(Dv(2,1)-Dv(2,3))*(Dv(3,1)-Dv(3,8))*(Dv(1,1)-Dv(1,6))&
           +(Dv(3,1)-Dv(3,3))*(Dv(1,1)-Dv(1,8))*(Dv(2,1)-Dv(2,6))&
           -(Dv(1,1)-Dv(1,6))*(Dv(2,1)-Dv(2,8))*(Dv(3,1)-Dv(3,3))&
           -(Dv(2,1)-Dv(2,6))*(Dv(3,1)-Dv(3,8))*(Dv(1,1)-Dv(1,3))&
           -(Dv(3,1)-Dv(3,6))*(Dv(1,1)-Dv(1,8))*(Dv(2,1)-Dv(2,3))))
             
  END FUNCTION gaux_getVolume_hexa3D
    
!************************************************************************

!<subroutine>
  
  ELEMENTAL SUBROUTINE gaux_isIntersection_line2D(&
      dx1,dy1,dx2,dy2,dx3,dy3,dx4,dy4, bintersect)
  
!<description>
  ! Checks whether the two 2D lines given by the start/endpoints 
  ! (x1,y1)->(x2,y2) and (x3,y3)->(x4,y4) intersect each other.
!</description>

!<input>
  ! First point on ray 1.
  REAL(DP), INTENT(IN) :: dx1,dy1
  
  ! A second point on ray 1. Must be different to (dx1,dy1)
  REAL(DP), INTENT(IN) :: dx2,dy2
  
  ! First point on ray 2.
  REAL(DP), INTENT(IN) :: dx3,dy3
  
  ! A second point on ray 2. Must be different to (dx3,dy3)
  REAL(DP), INTENT(IN) :: dx4,dy4
!</input>

!<result>
  ! TRUE if the two rays intersect. FALSE otherwise.
  LOGICAL, INTENT(OUT) :: bintersect
!</result>

!</subroutine>

    ! local variables: aux parameters
    REAL(DP) :: daux1, daux2, daux3, daux4

    ! position of point 3 with respect to line between 1 and 2

    daux3 = (dx2-dx1)*(dy3-dy1) - (dy2-dy1)*(dx3-dx1)

    ! position of point 4 with respect to line between 1 and 2

    daux4 = (dx2-dx1)*(dy4-dy1) - (dy2-dy1)*(dx4-dx1)

    ! position of point 1 with respect to line between 3 and 4

    daux1 = (dx4-dx3)*(dy1-dy3) - (dy4-dy3)*(dx1-dx3)

    ! position of point 2 with respect to line between 3 and 4

    daux2 = (dx4-dx3)*(dy2-dy3) - (dy4-dy3)*(dx2-dx3)

    ! Determine if the lines truely intersect by checking the sign

    bintersect = ((daux3*daux4 .LE. 0.0_DP) .AND. &
                  (daux1*daux2 .LE. 0.0_DP)) 

  END SUBROUTINE

!************************************************************************

!<subroutine>
  
  ELEMENTAL SUBROUTINE gaux_getIntersection_ray2D(&
      dx0,dy0,dx1,dy1,dx2,dy2,dx3,dy3, dx,dy, iintersect)
  
!<description>
  ! Calculates the intersection point of two 2D rays given by 
  ! (x1,y1)->(x2,y2) and (x3,y3)->(x4,y4).
!</description>

!<input>
  ! First point on ray 1.
  REAL(DP), INTENT(IN) :: dx0,dy0
  
  ! A second point on ray 1. Must be different to (dx1,dy1)
  REAL(DP), INTENT(IN) :: dx1,dy1
  
  ! First point on ray 2.
  REAL(DP), INTENT(IN) :: dx2,dy2
  
  ! A second point on ray 2. Must be different to (dx3,dy3)
  REAL(DP), INTENT(IN) :: dx3,dy3
!</input>

!<result>
  ! Intersection point.
  ! If the two rays do not intersect or are identical, this is set to (0,0).
  REAL(DP), INTENT(OUT) :: dx,dy

  ! Returns the type of intersection between the rays.
  ! =-1: The rays are the same
  ! = 0: The rays don't intersect.
  ! = 1: The rays intersect in exactly one point.
  INTEGER, INTENT(OUT) :: iintersect
!</result>

!</subroutine>

    ! local variables
    REAL(DP) :: ddet,da

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
    ! the intersection occurres.
    !
    ! The determinant of the system is:

    ddet = dx1*dy2-dx1*dy3-dx0*dy2+dx0*dy3-dy1*dx2+dy1*dx3+dy0*dx2-dy0*dx3 
       
    ! If it's =0, the lines are the same or completely different...
        
    IF (ddet .EQ. 0.0_DP) THEN
       
      ! If the vector (X2,Y2)->(X0,Y0) is linear dependent to
      ! (X2,Y2)->(X3,Y3), the lines are the same.

      ddet = -dy0*dx2-dx3*dy2+dy0*dx3+dx2*dy3+dx0*dy2-dx0*dy3
       
      IF (ddet .EQ. 0.0_DP) THEN
        iintersect = -1
      END IF
     
    ELSE

      ! There is an intersection point. Calculate one of the 
      ! "parameter" values along the two lines.

      da = (dy0*dx2+dx3*dy2-dy0*dx3-dx2*dy3-dx0*dy2+dx0*dy3) / ddet
        
      !  The intersection point is then

      dx = da*dx1 + (1.0_DP-da)*dx
      dy = da*dy1 + (1.0_DP-da)*dy
      
      iintersect = 1
       
    END IF
      
  END SUBROUTINE

!************************************************************************

!<subroutine>
  
  PURE SUBROUTINE gaux_isInElement_quad2D(dx,dy,DcornerCoords,binside)
  
!<description>
  ! Checks if a point (dx,dy) is inside of a 2D quadrilateral element
  ! given by the corners DcornerCoords.
!</description>

!<input>
  ! Point to check
  REAL(DP), INTENT(IN) :: dx,dy
  
  ! Array with coordinates of the four corner points of the element.
  ! The corners must be ordered in counterclockwise order.
  !
  ! Note: For performance reasons, this array is defined as
  !   explicit array of dimension (2,4). As this deactivates array
  !   checking in Fortran, the caller must take care to specify exactly
  !   this type of array here!
  REAL(DP), DIMENSION(2,4), INTENT(IN) :: DcornerCoords
!</input>

!<result>
  ! TRUE if (dx,dy) is inside of the element. FALSE otherwise.
  LOGICAL, INTENT(OUT) :: binside
!</result>

!</subroutine>

    ! local variables

    INTEGER, PARAMETER :: NVE = 4
    INTEGER, DIMENSION(4), PARAMETER :: Inext = (/2,3,4,1/)
    REAL(DP) :: dxmid,dymid,dxdist,dydist,dxnormal,dynormal
    INTEGER :: ive,ive2
    REAL(DP) :: dsproduct
      
    binside = .TRUE.

    ! Compute edge-midpoints and normal vectors to the four
    ! edges on element IEL

    DO ive=1,NVE

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
      ! If the point is "righthand" of all four edges, it's inside 
      ! of the element.

      dsproduct = dxdist*dxnormal + dydist*dynormal

      ! Actually we have to check against <=0, but it's more advisable
      ! to check against something that is 'near' 0 in terms
      ! of machine exactness...
      binside = binside .AND. (dsproduct .LE. dsproduct*100.0_DP)

    END DO
    
  END SUBROUTINE

!************************************************************************

!<subroutine>
  
  PURE SUBROUTINE gaux_getBarycentricCoords_tri2D(&
      DcornerCoords,dx,dy,dxi1,dxi2,dxi3)
  
!<description>
  ! Calculates the barycentric coordinates (dxi1,dxi2,dxi3) of a
  ! point (dx,dy) relative to a 2D triangular element specified by
  ! the coordinates of the three corners in DcornerCoords.
!</description>

!<input>
  ! Point in real coordinates
  REAL(DP), INTENT(IN) :: dx,dy
  
  ! Array with coordinates of the three corner points of the element.
  ! The corners must be ordered in counterclockwise order.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: DcornerCoords
!</input>

!<result>
  ! The barycentric coordinates of (dx,dy) relative to the element
  ! specified by DcornerCoords.
  REAL(DP), INTENT(OUT) :: dxi1,dxi2,dxi3
!</result>

!</subroutine>

    ! local variables
    REAL(DP) :: DAX, DAY, DBX, DBY, DCX, DCY, DDET

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

  END SUBROUTINE

!************************************************************************

!<subroutine>
  
  PURE SUBROUTINE gaux_isInElement_tri2D(dx,dy,DcornerCoords,binside)
  
!<description>
  ! Checks if a point (dx,dy) is inside of a 2D triangular element
  ! given by the corners DcornerCoords.
!</description>

!<input>
  ! Point to check
  REAL(DP), INTENT(IN) :: dx,dy
  
  ! Array with coordinates of the four corner points of the element.
  ! The corners must be ordered in counterclockwise order.
  !
  ! Note: For performance reasons, this array is defined as
  !   explicit array of dimension (2,4). As this deactivates array
  !   checking in Fortran, the caller must take care to specify exactly
  !   this type of array here!
  REAL(DP), DIMENSION(2,4), INTENT(IN) :: DcornerCoords
!</input>

!<result>
  ! TRUE if (dx,dy) is inside of the element. FALSE otherwise.
  LOGICAL, INTENT(OUT) :: binside
!</result>

!</subroutine>

    REAL(DP) :: dxi1,dxi2,dxi3

    ! We use barycentric coordinates for that task.
    ! Calculate the barycentric coordinates of the point relative
    ! to the element specified by DcornerCoords.
    CALL gaux_getBarycentricCoords_tri2D (DcornerCoords,dx,dy,&
        dxi1,dxi2,dxi3)

    ! If all barycentric coordinates are in the range [0..1],
    ! we are inside of the element
    binside = (dxi1 .GE. 0.0_DP) .AND. (dxi1 .LE. 1.0_DP) .AND. &
              (dxi2 .GE. 0.0_DP) .AND. (dxi2 .LE. 1.0_DP) .AND. &
              (dxi3 .GE. 0.0_DP) .AND. (dxi3 .LE. 1.0_DP) 

  END SUBROUTINE


!************************************************************************

!<function>

  PURE LOGICAL FUNCTION gaux_isFlipped_hexa3D (Dpoints)

!<description>
  ! This function checks whether a 3D hexahedron is flipped.
!</description>

!<input>
  ! The coordinates of the eight corners of the hexahedron.
  ! Dpoints(1,.) = x-coordinates,
  ! Dpoints(2,.) = y-coordinates,
  ! Dpoints(3,.) = z-coordinates
  REAL(DP), DIMENSION(3,8), INTENT(IN) :: Dpoints
!</input>

!<result>
  ! .TRUE. if the hexahedron is flipped, otherwise .FALSE.
!</result>
!</function>

    ! Three vectors connecting two opposite faces of the hexahedron,
    ! and a normal vector
    REAL(DP), DIMENSION(3) :: Du,Dv,Dw,Dn
    REAL(DP) :: dt
    
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
    gaux_isFlipped_hexa3D = (dt .LT. 0.0_DP)

  END FUNCTION gaux_isFlipped_hexa3D

END MODULE
