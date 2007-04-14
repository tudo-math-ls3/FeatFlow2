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

  PURE REAL(DP) FUNCTION gaux_getArea_quad2D (Dpoints)

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
    
    gaux_getArea_quad2D = &
        ABS(gaux_getArea_tria2D(Dpoints(1:2,1:3))) +&
        ABS(gaux_getArea_tria2D(Dpoints(1:2,(/1,3,4/))))

  END FUNCTION gaux_getArea_quad2D
END MODULE
