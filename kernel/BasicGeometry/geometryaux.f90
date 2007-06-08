!##############################################################################
!# ****************************************************************************
!# <name> geometryaux </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains low level auxiliary functions for geometric
!# objects like triangles, quads or similar.
!#
!# The following routines can be found here:
!#
!# 1.) gaux_getAspectRatio_quad2D
!#     -> Calculate the aspect ratio of a 2D quadrilateral
!#
!# 2.) gaux_getArea_tria2D
!#     -> Calculate the signed area of a 2D triangle
!#
!# 3.) gaux_getArea_quad2D
!#     -> Calculate the area of a 2D quadrilateral
!#
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
END MODULE
