!#########################################################################
!# ***********************************************************************
!# <name> ttrafo </name>
!# ***********************************************************************
!#
!# <purpose>
!# This module contains routines for the transformation and
!# back-transformation of points between a reference triangle
!# and the "real" element.
!#
!# The following routine can be found in this module:
!#
!# 1.) ttrafo_calcBaryCoord
!#     -> find the barycentric coordinates for a given point (x,y)
!#        in real coordinates.
!#
!# </purpose>
!#########################################################################

MODULE ttrafo

  USE fsystem

  IMPLICIT NONE
  
  CONTAINS
  
!<subroutine>
  SUBROUTINE ttrafo_calcBaryCoord (Dcoord, dbary1, dbary2, dbary3, &
                                   ddet, dx, dy)
  
  !<description>
    ! This subroutine is to find the barycentric coordinates for a given 
    ! point (x,y) in real coordinates.
    
    ! Explaination of the transformation:
    !
    ! In contrast to quadrilaterals the linear transformation on triangles
    ! is rather easy when using barycentric coordinates. Each point
    ! (X,Y) in a triangle is identified by a 3-tuple of coordinates
    ! (X1,X2,X3) giving the "relative" position of the point inside the
    ! triangle:
    !
    !        P3                      (0,1)
    !       /  \         Phi           |  \
    !      /    \        <--           |    \
    !     /  p   \                     | R    \
    !    P1------P2                  (0,0)---(1,0)
    !
    ! here:    p = X1*P1 + X2*P2 + X3*P3
    !
    ! The barycentric coordinates are "independent" of the triangle:
    ! When a 3-tuple of barycentric coordinates for a point p is obtained,
    ! this holds for both, the "real" triangle as well as the reference one.
    ! Therefore the corresponting point R of p on the reference triangle
    ! can be obtained by setting:
    !
    !          R = Phi^{-1}(p) = X1*(0,0) + X2*(1,0) + X3*(0,1)
    !
    ! So when the barycentric coordinates are known, the transformation
    ! is trivial and can be included in the code directly. The only crucial
    ! point is to calculate the barycentric coordinates from a set of
    ! corners of a triangle.
  !</description>
    
  !<input>
    ! coordinates of the element vertices
    REAL(DP), DIMENSION(2,3), INTENT(IN) :: Dcoord
    
    ! coordinated of the evaluation point
    REAL(DP), INTENT(IN) :: dx, dy
  !</input>
    
  !<output>
    ! barycentric coordinate values of (x,y)
    REAL(DP), INTENT(OUT) :: dbary1, dbary2, dbary3
    
    ! determinant for triangle
    REAL(DP), INTENT(OUT) :: ddet
  !</output>
    
!</subroutine>
    
    !local variables
    REAL(DP) :: dax, day, dbx, dby, dcx, dcy, ddetInv
    
    dax = Dcoord(1, 1) 
    day = Dcoord(2, 1)
    dbx = Dcoord(1, 2)
    dby = Dcoord(2, 2)
    dcx = Dcoord(1, 3)
    dcy = Dcoord(2, 3)

    ! Example where to find this formula here:
    ! http://home.t-online.de/home/nagel.klaus/matdir/bary.htm 
	
    ddet = dax*(dby-dcy) + dbx*(dcy-day) + dcx*(day-dby)
    IF (ddet .EQ. 0.0_DP) THEN
      PRINT *, 'ttrafo_calcBaryCoord: Critical error! Determinant is zero.'
      STOP
    END IF
    
    ddetInv = 1.0_DP / ddet
    dbary1 = (dx *(dby-dcy)+dbx*(dcy-dx )+dcx*(dy -dby)) * ddetInv 
    dbary2 = (dax*(dy -dcy)+dx *(dcy-day)+dcx*(day-dy )) * ddetInv
    dbary3 = (dax*(dby-dy )+dbx*(dy -day)+dx *(day-dby)) * ddetInv

       
  END SUBROUTINE ttrafo_calcBaryCoord

END MODULE
