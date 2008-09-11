************************************************************************
* This file contains routines for the transformation and
* back-transformation of points between a reference triangle
* and the "real" element.
************************************************************************

************************************************************************
* Explaination of the transformation:
*
* In contrast to quadrilaterals the linear transformation on triangles
* is rather easy when using barycentric coordinates. Each point
* (X,Y) in a triangle is identified by a 3-tuple of coordinates
* (X1,X2,X3) giving the "relative" position of the point inside the
* triangle:
*
*        P3                      (0,1)
*       /  \         Phi           |  \
*      /    \        <--           |    \
*     /  p   \                     | R    \
*    P1------P2                  (0,0)---(1,0)
*
* here:    p = X1*P1 + X2*P2 + X3*P3
*
* The barycentric coordinates are "independent" of the triangle:
* When a 3-tuple of barycentric coordinates for a point p is obtained,
* this holds for both, the "real" triangle as well as the reference one.
* Therefore the corresponting point R of p on the reference triangle
* can be obtained by setting:
*
*          R = Phi^{-1}(p) = X1*(0,0) + X2*(1,0) + X3*(0,1)
*
* So when the barycentric coordinates are known, the transformation
* is trivial and can be included in the code directly. The only crucial
* point is to calculate the barycentric coordinates from a set of
* corners of a triangle. Routines for this purpose can be found in the
* following.
************************************************************************

************************************************************************
* Triangular back-transformation
*
* This subroutine is to find the barycentric coordinates for a given 
* point (x,y) in real coordinates.
* 
* In:
*  DCOORD - array [1..2,1..3] of double
*           Coordinates of the three corners of the real triangle.
*           DCOORD(1,.) saves the X-, DCOORD(2,.) the Y-coordinates.
*  (DXREAL,DYREAL) - coordinates of a point in the real quadrilateral
*
* Out:
*  (X1,X2,X3)   - Barycentric coordinates of the point
*  DET          - Determinant of the transformation
*
* There is no check whether the point is really inside the triangle or
* not. This can easily be checked as the point is inside
* the triangle if 0 <= X1,X2,X3 <= 1.
************************************************************************

      SUBROUTINE TBTRAF (DCOORD, X1, X2, X3, DET, DXREAL, DYREAL)
      
      IMPLICIT NONE
      
C parameters
      
C coordinates of the evaluation point
      DOUBLE PRECISION DXREAL,DYREAL

C coordinates of the element vertices
      DOUBLE PRECISION DCOORD (2,3)
      
C barycentric coordinates values of (x,y)
      DOUBLE PRECISION X1, X2, X3, DET

C local variables
      DOUBLE PRECISION DAX, DAY, DBX, DBY, DCX, DCY, DDET

	DAX = DCOORD(1, 1) 
	DAY = DCOORD(2, 1)
	DBX = DCOORD(1, 2)
	DBY = DCOORD(2, 2)
	DCX = DCOORD(1, 3)
	DCY = DCOORD(2, 3)
	
C Example where to find this formula here:
C http://home.t-online.de/home/nagel.klaus/matdir/bary.htm 
	
	DET = DAX*(DBY-DCY) + DBX*(DCY-DAY) + DCX*(DAY-DBY)
	DDET = 1D0 / DET
	X1 = (DXREAL*(DBY-DCY)+DBX*(DCY-DYREAL)+DCX*(DYREAL-DBY)) * DDET 
	X2 = (DAX*(DYREAL-DCY)+DXREAL*(DCY-DAY)+DCX*(DAY-DYREAL)) * DDET
	X3 = (DAX*(DBY-DYREAL)+DBX*(DYREAL-DAY)+DXREAL*(DAY-DBY)) * DDET

      END