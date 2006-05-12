************************************************************************
* This file contains routines for the transformation and
* back-transformation of points between a reference quadrilateral
* and the "real" element.
************************************************************************

C **********************************************************************
C Explaination of the transformation:
C
C We want to perform a transformation from the reference quadrilateral
C [-1,1]x[-1,1] onto a "real" quadrilaterl:
C
C   (-1,1) ___ (1,1)         (x4,y4)   ___  (x3,y3)
C         |  |                        /  \
C         |__|           =>          /____\
C  (-1,-1)    (1,-1)         (x1,y1)        (x2,y2)   
C
C By theory this can by done with a bilinear mapping, i.e. a mapping
C Psi:R^2->R^2  of the form:
C
C      Psi (xi1) = ( a1 + a2*xi1 + a3*xi2 + a4*xi1*xi2 )
C          (xi2)   ( b1 + b2*xi1 + b3*xi2 + b4*xi1*xi2 )
C
C Our special transformation has to map:
C
C  Psi(-1,-1) = (x1,y1)
C  Psi( 1,-1) = (x2,y2)
C  Psi( 1, 1) = (x3,y3)
C  Psi(-1, 1) = (x4,y4)
C
C This gives the linear system:
C
C  a1 - a2 - a3 - a4 = x1       b1 - b2 - b3 - b4 = y1
C  a1 + a2 - a3 - a4 = x2       b1 + b2 - b3 - b4 = y2
C  a1 + a2 + a3 + a4 = x3       b1 + b2 + b3 + b4 = y3
C  a1 - a2 + a3 - a4 = x4       b1 - b2 + b3 - b4 = y4
C
C Reorder this to calculate the ai:
C
C  a1 = 1/4 * ( x1 + x2 + x3 + x4)      b1 = 1/4 * ( y1 + y2 + y3 + y4)
C  a2 = 1/4 * (-x1 + x2 + x3 - x4)      b2 = 1/4 * (-y1 + y2 + y3 - y4)
C  a3 = 1/4 * (-x1 - x2 + x3 + x4)      b3 = 1/4 * (-y1 - y2 + y3 + y4)
C  a4 = 1/4 * ( x1 - x2 + x3 - x4)      b4 = 1/4 * ( y1 - y2 + y3 - y4)
C
C The factors in the brackets in these equations are only dependent on 
C the corners of the quadrilateral, not of the current point. So they
C are constant for all points we want to map from the reference
C element to the real one. We call them here "auxiliary Jacobian factors".
C They can be calculated with QINIJF in advance for all points that
C have to be mapped. To be more exact, QINIJF calculates:
C
C  J1 = 1/2 * (-x1 - x2 + x3 - x4)
C  J2 = 1/2 * ( x1 - x2 + x3 - x4)
C
C  J3 = 1/2 * (-y1 + y2 - y3 + y4)
C  J4 = 1/2 * (-y1 + y2 + y3 - y4)
C
C Using these factors, one can write:
C
C  a1  =  1/4 * ( x1 + x2 + x3 + x4)  =  1/2 * (x1 + x2 + J1)
C  a2  =  1/4 * (-x1 + x2 + x3 - x4)  =  1/2 * (x2 - x1 + J2)
C  a3  =  1/4 * (-x1 - x2 + x3 + x4)  =  1/2 * J1 
C  a4  =  1/4 * ( x1 - x2 + x3 - x4)  =  1/2 * J2
C
C  b1  =  1/4 * ( y1 + y2 + y3 + y4)  =  1/2 * (y1 + y3 + J3)
C  b2  =  1/4 * (-y1 + y2 + y3 - y4)  =  1/2 * J4
C  b3  =  1/4 * (-y1 - y2 + y3 + y4)  =  1/2 * (y3 - y1 - J4)
C  b4  =  1/4 * ( y1 - y2 + y3 - y4)  =  1/2 * J3
C
C The Jacobian matrix of the bilinear transformation is now 
C calculated as usual by partial differentiation. Thing above 
C coefficients ai and bi one can write:
C
C DPhi ( xi1 )  =  ( a2 + a4*xi2    a3 + a4*xi1 )
C      ( xi2 )     ( b2 + b4*xi2    b3 + b4*xi1 )
C
C               =  ( 1/2*(x2-x1+J2) + 1/2*J2*xi2              1/2*J1 + 1/2*J2*xi1 )
C                  (         1/2*J4 + 1/2*J3*xi2      1/2*(y3-y1-J4) + 1/2*J3*xi1 )
C
C which gives the Jacobian determinant of the 2x2-matrix:
C
C   det DPhi (xi1,xi2)  =  (DPhi[1,1]*DPhi[2,2] - DPhi[1,2]*DPhi[2,1]) (xi1,xi2)
C
C Using these information makes it possible to map a point (XI1,XI2)
C on the reference element to coordinates (XX,YY) on the real element by:
C
C   (XX,YY) := Psi(XI1,XI2) 
C **********************************************************************

************************************************************************
* Calculate auxiliary Jacobian factors
*
* This routine builds up constant factors that are later used during
* the transformation from the reference element to the real element.
* This is used for saving some computational time...
*
* In:
*  DCOORD - array [1..2,1..4] of double
*           Coordinates of the four corners of the real quadrilateral.
*           DCOORD(1,.) saves the X-, DCOORD(2,.) the Y-coordinates.
* 
* Out:
*  DJF    - array [1..2,1..2] of double
*           Auxiliary constant factors for calculation of
*           Jacobian matrix
************************************************************************

      SUBROUTINE QINIJF (DCOORD,DJF)
      
      IMPLICIT NONE

C parameters
      
      DOUBLE PRECISION DCOORD (2,4),DJF(2,2)
      
      DJF(1,1)=0.5D0*(-DCOORD(1,1)-DCOORD(1,2)+DCOORD(1,3)+DCOORD(1,4))
      DJF(1,2)=0.5D0*( DCOORD(1,1)-DCOORD(1,2)+DCOORD(1,3)-DCOORD(1,4))
      DJF(2,1)=0.5D0*(-DCOORD(2,1)+DCOORD(2,2)-DCOORD(2,3)+DCOORD(2,4))
      DJF(2,2)=0.5D0*(-DCOORD(2,1)+DCOORD(2,2)+DCOORD(2,3)-DCOORD(2,4))
      
      END

************************************************************************
* Quadrilateral transformation of a point from the reference element
* onto the real element
*
* This subroutine performs two tasks:
* -Initialisation of a a given 2x2 matrix with the
*  mapping information from the reference element to the "real"
*  quadrilateral. Calculation of the Jacobian determinant
* -Transformation of a given point on the reference element onto
*  the "real" element
* Both things are performed simultaneously because the jacobian
* determinant is dependent of the point.
*
* Before this routine can be called, the auxiliary factors DJF
* have to be calculated with QINIJF for the considered element.
*
* In:
*  DCOORD - array [1..2,1..4] of double
*           Coordinates of the four corners of the real quadrilateral.
*           DCOORD(1,.) saves the X-, DCOORD(2,.) the Y-coordinates.
*  DJF    - array [1..2,1..2] of double
*           Auxiliary constants for the considered element with
*           coordinates in DCOORD; have to be computed previously
*           by QINIJF.
*  (DXPAR,DYPAR) - Coordinates of a point on the reference element
*  
* Out:
*  DJAC   - array [1..2,1..2] of double
*           2x2-system that receives the Jacobian matrix of the
*           mapping. When calling this routine this variable should
*           point to the DJAC variable in the COMMON block for
*           later calculations.
*  DETJ   - Determinant of the jacobian matrix
*  (DXREAL,DYREAL) - Coordinates of the point on the real element
************************************************************************

      SUBROUTINE QTRAF (DCOORD,DJF,DJAC,DETJ,DPARX,DPARY,DXREAL,DYREAL)
      
      IMPLICIT NONE

C parameters
      
      DOUBLE PRECISION DCOORD (2,4),DJF(2,2),DJAC(2,2),DETJ
      DOUBLE PRECISION DPARX,DPARY,DXREAL,DYREAL
      
C Jacobian matrix
      
      DJAC(1,1)=0.5D0 * (DCOORD(1,2)-DCOORD(1,1)+DJF(1,2)) +
     *          0.5D0 * DJF(1,2)*DPARY
      DJAC(1,2)=0.5D0 * DJF(1,1) + 
     *          0.5D0 * DJF(1,2)*DPARX
      DJAC(2,1)=0.5D0 * DJF(2,2) - 
     *          0.5D0 * DJF(2,1)*DPARY
      DJAC(2,2)=0.5D0 * (DCOORD(2,3)-DCOORD(2,1)-DJF(2,2)) - 
     *          0.5D0 * DJF(2,1)*DPARX

C Determinant of the mapping

      DETJ = DJAC(1,1)*DJAC(2,2) - DJAC(1,2)*DJAC(2,1)
      
C Map the point to the real element

      DXREAL = 0.5D0*(DCOORD(1,1)+DCOORD(1,2)+DJF(1,1)) +
     *         0.5D0*(DCOORD(1,2)-DCOORD(1,1)+DJF(1,2))*DPARX +
     *         0.5D0*DJF(1,1)*DPARY +
     *         0.5D0*DJF(1,2)*DPARX*DPARY
      DYREAL = 0.5D0*(DCOORD(2,1)+DCOORD(2,3)+DJF(2,1)) +
     *         0.5D0*DJF(2,2)*DPARX +
     *         0.5D0*(DCOORD(2,3)-DCOORD(2,1)-DJF(2,2))*DPARY -
     *         0.5D0*DJF(2,1)*DPARX*DPARY

      END

************************************************************************
* Calculate Jacobian determinant of mapping from reference- to
* real element.
*
* This routine only calculates the Jacobian determinant. of the
* mapping. this is on contrast to QTRAF, which not only calculates
* this determinant but also maps the point. So this routine
* can be used to speed up the code if the coordinates of the
* mapped point already exist.
*
* Before this routine can be called, the auxiliary factors DJF
* have to be calculated with QINIJF for the considered element.
*
* In:
*  DCOORD - array [1..2,1..4] of double
*           Coordinates of the four corners of the real quadrilateral.
*           DCOORD(1,.) saves the X-, DCOORD(2,.) the Y-coordinates.
*  DJF    - array [1..2,1..2] of double
*           Auxiliary constants for the considered element with
*           coordinates in DCOORD; have to be computed previously
*           by QINIJF.
*  
* Out:
*  DJAC   - array [1..2,1..2] of double
*           2x2-system that receives the Jacobian matrix of the
*           mapping. When calling this routine this variable should
*           point to the DJAC variable in the COMMON block for
*           later calculations.
*  DETJ   - Determinant of the jacobian matrix
************************************************************************

      SUBROUTINE QTRDET (DCOORD,DJF,DJAC,DETJ,DPARX,DPARY)
      
      IMPLICIT NONE

C parameters
      
      DOUBLE PRECISION DCOORD (2,4),DJF(2,2),DJAC(2,2),DETJ
      DOUBLE PRECISION DPARX,DPARY
      
C Jacobian matrix
      
      DJAC(1,1)=0.5D0 * (DCOORD(1,2)-DCOORD(1,1)+DJF(1,2)) +
     *          0.5D0 * DJF(1,2)*DPARY
      DJAC(1,2)=0.5D0 * DJF(1,1) + 
     *          0.5D0 * DJF(1,2)*DPARX
      DJAC(2,1)=0.5D0 * DJF(2,2) - 
     *          0.5D0 * DJF(2,1)*DPARY
      DJAC(2,2)=0.5D0 * (DCOORD(2,3)-DCOORD(2,1)-DJF(2,2)) - 
     *          0.5D0 * DJF(2,1)*DPARX

C Determinant of the mapping

      DETJ = DJAC(1,1)*DJAC(2,2) - DJAC(1,2)*DJAC(2,1)
      
      END

************************************************************************
* Quadrilateral back-transformation
*
* This subroutine is to find the parameter values for a given point
* (x,y) in real coordinates.
* 
* Remark: This is a difficult task, as usually in FEM codes the
* parameter values are known and one wants to obtain the real
* coordinates.
* 
* inverting the bilinear trafo in a straightforward manner by using 
* pq-formula does not work very well, as it is numerically unstable. 
* For parallelogram-shaped elements, one would have to introduce a
* special treatment. For nearly parallelogram-shaped elements, this 
* can cause a crash as the argument of the square root can become 
* negative due to rounding errors. In the case of points near
* the element borders, we divide nearly 0/0.
* 
* Therefore, we have implemented the algorithm described in
* 
* Introduction to Finite Element Methods, Carlos Felippa,
* Department of Aerospace Engineering Sciences and Center for 
* Aerospace Structures, http://titan.colorado.edu/courses.d/IFEM.d/  
*
* In:
*  DCOORD - array [1..2,1..4] of double
*           Coordinates of the four corners of the real quadrilateral.
*           DCOORD(1,.) saves the X-, DCOORD(2,.) the Y-coordinates.
*  (DXREAL,DYREAL) - coordinates of a point in the real quadrilateral
*
* Out:
*  (DXPAR,DYPAR)   - Coordinates of the transformed point in the
*                    reference element.
************************************************************************

      SUBROUTINE QBTRAF (DCOORD, DXPAR, DYPAR, DXREAL, DYREAL)
      
      IMPLICIT NONE
      
C parameters
      
C coordinates of the evaluation point
      DOUBLE PRECISION DXREAL,DYREAL

C coordinates of the element vertices
      DOUBLE PRECISION DCOORD (2,4)
C parameter values of (x,y)
      DOUBLE PRECISION DXPAR,DYPAR

C local variables

      DOUBLE PRECISION X1,X2,X3,X4,Y1,Y2,Y3,Y4,XB,YB,XCX,YCX,XCE,YCE
      DOUBLE PRECISION A,J1,J2,X0,Y0,BXI,BETA,CXI,XP0,YP0,CETA
      DOUBLE PRECISION ROOT1,ROOT2,XIP1,XIP2,ETAP1,ETAP2,D1,D2

C Get nodal x-coordinates
      X1 = DCOORD(1,1)
      X2 = DCOORD(1,2)
      X3 = DCOORD(1,3)
      X4 = DCOORD(1,4)

C Get nodal y-coordinates
      Y1 = DCOORD(2,1)
      Y2 = DCOORD(2,2)
      Y3 = DCOORD(2,3)
      Y4 = DCOORD(2,4)

      XB = X1-X2+X3-X4
      YB = Y1-Y2+Y3-Y4

      XCX = X1+X2-X3-X4
      YCX = Y1+Y2-Y3-Y4

      XCE = X1-X2-X3+X4
      YCE = Y1-Y2-Y3+Y4

      A = 0.5*((X3-X1)*(Y4-Y2)-(X4-X2)*(Y3-Y1))

      J1 = (X3-X4)*(Y1-Y2)-(X1-X2)*(Y3-Y4)
      J2 = (X2-X3)*(Y1-Y4)-(X1-X4)*(Y2-Y3)

      X0 = 0.25*(X1+X2+X3+X4)
      Y0 = 0.25*(Y1+Y2+Y3+Y4)

      XP0 = DXREAL-X0
      YP0 = DYREAL-Y0

      BXI  =  A-XP0*YB+YP0*XB
      BETA = -A-XP0*YB+YP0*XB

      CXI  = XP0*YCX-YP0*XCX
      CETA = XP0*YCE-YP0*XCE

      ROOT1 = -DSQRT(BXI**2-2.0*J1*CXI)-BXI
      ROOT2 =  DSQRT(BXI**2-2.0*J1*CXI)-BXI
      IF (ROOT1.NE.0D0) THEN
        XIP1 = 2.0*CXI/ROOT1
      ELSE
        XIP1 = 1E15
      END IF
      IF (ROOT2.NE.0D0) THEN
        XIP2 = 2D0*CXI/ROOT2
      ELSE
        XIP2 = 1E15
      END IF

      ROOT1 =  DSQRT(BETA**2+2D0*J2*CETA)-BETA
      ROOT2 = -DSQRT(BETA**2+2D0*J2*CETA)-BETA
      IF (ROOT1.NE.0D0) THEN
        ETAP1 = 2D0*CETA/ROOT1
      ELSE
        ETAP1 = 1D15
      END IF
      IF (ROOT2.NE.0D0) THEN
        ETAP2 = 2D0*CETA/ROOT2
      ELSE
        ETAP2 = 1D15
      END IF

      D1 = DSQRT(XIP1**2+ETAP1**2)
      D2 = DSQRT(XIP2**2+ETAP2**2)

      IF (D1.LT.D2) THEN
        DXPAR = XIP1
        DYPAR = ETAP1
      ELSE
        DXPAR = XIP2
        DYPAR = ETAP2
      END IF

      END