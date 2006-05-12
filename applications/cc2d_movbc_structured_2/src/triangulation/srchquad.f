**********************************************************************
* This file contains routines for searching points in the grid.
* All the routines here are trimmed for quadrilateral elements.
*
* There are also some general auxiliary routines to be found here
* which are valid for more general purposes.
**********************************************************************

**********************************************************************
* Purpose: returns the element containing the point (x0,y0) by a
* raytracing method
* Constructs a ray from the old node position (whose element number
* is known) and traces it to the actual position
*
* In:
*  X0,
*  Y0     - Coordinates of the (moved) point to determine the element 
*           number where this point resides in
*  IVT0   - Number of the vertex with coordinates (X0,Y0);
*           Only for output purposes in case of an error; can be 
*           ignored if necessary!
*  IELOLD - Element in which the not-moved point (XALT,YALT) 
*           resides in
* 
* Out:
*  IEL    - Number of the element where (X0,Y0) stays in
*  XALT, 
*  YALT   - Coordinates of the last found point of the ray inside
*           of the domain,
*           a) before the point (X0,Y0), if (X0,Y0) is inside of
*              the domain
*           b) before leaving the domain in direction of (X0,Y0)
*              in case of an error (domain left)
**********************************************************************

      SUBROUTINE PSRCH1(X0,Y0,XALT,YALT,IVT0,DCORVG,KVERT,KADJ,IELOLD,
     *                  IEL)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      
      INTEGER MAXIT
      PARAMETER (MAXIT=100)
      
C parameters
      
      DOUBLE PRECISION X0,Y0,XALT,YALT
      DOUBLE PRECISION DCORVG(2,*)
      INTEGER IVT0,KVERT,KADJ
      INTEGER IELOLD,IEL
      
      DIMENSION KVERT(NNVE,*),KADJ(NNVE,*)

C local variables

      DOUBLE PRECISION X1EDGE,Y1EDGE,X2EDGE,Y2EDGE
      INTEGER IVE,IDX,NVE
      LOGICAL BAUX,ISINE2,BINTSE

      NVE=4
      IEL=IELOLD
   
C We restrict our raytracing search to 100 neighbour cells;
C let's hope a point is not moved more than 100 elements in one step...

      DO IDX=1,MAXIT
      
C *** finished - our node is in element iel
C In that case, jump out of this subroutine.

        IF (ISINE2(X0,Y0,IEL,DCORVG,KVERT)) GOTO 110
        
C *** calculate point of intersection with one of the four surrounding edges
C *** take the center of the starting element as starting point of the ray

        XALT=0.5D0*(DCORVG(1,KVERT(1,IEL))+DCORVG(1,KVERT(3,IEL)))
        YALT=0.5D0*(DCORVG(2,KVERT(1,IEL))+DCORVG(2,KVERT(3,IEL)))
        
        DO IVE=1,NVE
C *** don't jump back to the old element
          IF (KADJ(IVE,IEL).NE.IELOLD) THEN
            X1EDGE=DCORVG(1,KVERT(IVE,IEL))
            Y1EDGE=DCORVG(2,KVERT(IVE,IEL))
            X2EDGE=DCORVG(1,KVERT(MOD(IVE,NVE)+1,IEL))
            Y2EDGE=DCORVG(2,KVERT(MOD(IVE,NVE)+1,IEL))
            BAUX=BINTSE(XALT,YALT,X0,Y0,X1EDGE,Y1EDGE,X2EDGE,Y2EDGE)

C *** found a point of intersection

            IF (BAUX) THEN
C go on searching in the neighbour element
              IELOLD=IEL
C Stop here if we leave our domain. To test that use the fact that
C KADJ()=0 if there is no neighbour element!
              IF (KADJ(IVE,IELOLD).EQ.0) GOTO 310
              IEL=KADJ(IVE,IELOLD)
              GOTO 210
            END IF
          END IF
        END DO

210   END DO

C If the search is successful, we cancel the DO-loop and jump out
C of this routine; otherwise we halt the program because we have not
C been able to find the element where the point resides in...
      
310   WRITE (*,FMT='(A,I7)') 
     *        'Found no Element containing the node ',IVT0
      WRITE (*,FMT='(A,2F12.7)') 'Position (x/y)=',X0,Y0
      STOP

110   END 

**********************************************************************
* This subroutine tests, if the two line segments between (1,2)
* and (3,4) truely intersect.
*
* To achive this, we test the position of point 1 with respect to
* line (3,4) and do same for point 2. We repeat this with the points 
* 3 and 4 with respect to line (1,2). If point 1 is "under" the line 
* (3,4) and point 2 is "over" line 3,4 (or vice versa) and 
* if point 3 is "left" from line (1,2) and point 4 is "right" from 
* line (1,2) (or vice versa), then the two lines intersect.
*
* To determine the position e.g. of point 1 with respect to line 
* 3,4, we consider the orientation of the triangle (3,4,1) by 
* computing the oriented area of the triangle.
* If the oriented area of triangle (3,4,1) has a different sign
* that triangle (3,4,2), line (1,2) intersects line (3,4)!
*
* In:
*  DX1,DY1,
*  DX2,DY2  - The start- and endpoint of the first line segment
*  DX3,DY3,
*  DX4,DY4  - The start- and endpoint of the 2nd line segment
*  
* Result:
*  =true,  if the two line segments intersecvt
*  =false, otherwise
**********************************************************************

      LOGICAL FUNCTION BINTSE(dx1,dy1,dx2, dy2, dx3, dy3, dx4, dy4)

      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      
C Parameters: coordinates of the for points defining the two lines
      DOUBLE PRECISION dx1, dy1, dx2, dy2, dx3, dy3, dx4, dy4

C local variables: aux parameters
      DOUBLE PRECISION daux1, daux2, daux3, daux4

C position of point 3 with respect to line between 1 and 2

      DAUX3 = (DX2-DX1)*(DY3-DY1) - (DY2-DY1)*(DX3-DX1)

C position of point 4 with respect to line between 1 and 2

      DAUX4 = (DX2-DX1)*(DY4-DY1) - (DY2-DY1)*(DX4-DX1)

C position of point 1 with respect to line between 3 and 4

      DAUX1 = (DX4-DX3)*(DY1-DY3) - (DY4-DY3)*(DX1-DX3)

C position of point 2 with respect to line between 3 and 4

      DAUX2 = (DX4-DX3)*(DY2-DY3) - (DY4-DY3)*(DX2-DX3)

C Determine if the lines truely intersect by checking the sign

      BINTSE = ((DAUX3*DAUX4.LE.0D0).AND.(DAUX1*DAUX2.LE.0D0)) 

C      if ((daux3*daux4.le.-1.0D-28).and.(daux1*daux2.le.-1.0D-28)) then

      END 

**********************************************************************
* Test if a point is inside of an element.
*
* Returns .true. if the position (x0,y0) is in the element IEL.
*
* In:
*   X0,Y0  - Position of an arbitrary point
*   DCORVG - coordinates of all points
*   KVERT  - Vertex numbers on all elements
*   IEL    - Element number
*
* Result:
*   =true,  if (X0,Y0) is in element IEL
*   =false, otherwise
**********************************************************************

      LOGICAL FUNCTION ISINEL(X0,Y0,IEL,DCORVG,KVERT)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      
C parameters
      
      DOUBLE PRECISION X0,Y0
      INTEGER IEL
      DOUBLE PRECISION DCORVG(2,*)
      INTEGER KVERT(NNVE,*)
      
C local variables

      DOUBLE PRECISION DPROD
      DOUBLE PRECISION X(4),Y(4),XMID(4),YMID(4),XNRM(4),YNRM(4)
      DOUBLE PRECISION XDST(4),YDST(4)
      INTEGER IVE,IVE2,NVE
      LOGICAL BINSIDE
      
      NVE=4
      BINSIDE=.TRUE.

C     Fetch the corners of element IEL to local arrays X and Y:

      DO IVE=1,NVE
        X(IVE)=DCORVG(1,KVERT(IVE,IEL))
        Y(IVE)=DCORVG(2,KVERT(IVE,IEL))
        
C       If (X0,Y0) is a corner of the element, we can immediately
C       return TRUE.
        
        IF ((X(IVE).EQ.X0).AND.(Y(IVE).EQ.Y0)) THEN
          ISINEL = .TRUE.
          RETURN
        END IF
      END DO
      
C     Compute edge-midpoints and normal vectors to the four
C     edges on element IEL

      DO IVE=1,NVE

        IVE2=MOD(IVE,NVE)+1

C       compute midpoints of element edges

        XMID(IVE)=0.5D0*(X(IVE)+X(IVE2))
        YMID(IVE)=0.5D0*(Y(IVE)+Y(IVE2))

C       compute normal vectors to element edges

        XNRM(IVE)= Y(IVE2)-Y(IVE)
        YNRM(IVE)=-X(IVE2)+X(IVE)

C       compute vectors from edge midpoints to node 'ivt'

        XDST(IVE)=X0-XMID(IVE)
        YDST(IVE)=Y0-YMID(IVE)

      END DO

C     Check whether 'ivt' belongs to element 'iel' or not by
C     multiplying distance vectors with corresponding normal vectors.
C     The sign of this scalar product determines whether we are
C     'left' or 'right' of the edge (because of the cosine formula).
C     If the point is "righthand" of all four edges, it's inside 
C     of the element.

      DO IVE=1,NVE
        DPROD = XDST(IVE)*XNRM(IVE)+YDST(IVE)*YNRM(IVE)
        BINSIDE = BINSIDE.AND.(DPROD.LE.1D-14)
      END DO
      
      ISINEL = BINSIDE

      END 

**********************************************************************
* Returns .true. if the position (x0,y0) is in the element iel.
*
* This is similar to ISINEL, but uses a different approach.
* Calculates the area of a triangle [(x0,y0), corner i and corner i+1]
* of the quad. If the area is positive, the corners are arranged in a
* counterclockwise sense. If this is the case for all all four
* possible i's, the point is inside of the element.
*
* In:
*   X0,Y0  - Position of an arbitrary point
*   DCORVG - coordinates of all points
*   KVERT  - Vertex numbers on all elements
*   IEL    - Element number
*
* Result:
*   =true,  if (X0,Y0) is in element IEL
*   =false, otherwise
**********************************************************************

      LOGICAL FUNCTION ISINE2(X0,Y0,IEL,DCORVG,KVERT)

      IMPLICIT NONE

      INCLUDE 'cbasictria.inc'

C parameters

      DOUBLE PRECISION X0,Y0,DCORVG
      INTEGER IEL,KVERT
      DIMENSION DCORVG(2,*),KVERT(NNVE,*)
      
C local variables
      
      DOUBLE PRECISION X,Y,DAREA
      DIMENSION X(4),Y(4)
      INTEGER IVE,IVE2,NVE
      
      NVE=4
      ISINE2 = .FALSE.

C     Fetch the corners of element IEL to local arrays X and Y:

      DO IVE=1,NVE
        X(IVE)=DCORVG(1,KVERT(IVE,IEL))
        Y(IVE)=DCORVG(2,KVERT(IVE,IEL))

C       If (X0,Y0) is a corner of the element, we can immediately
C       return TRUE.

        IF ((X(IVE).EQ.X0).AND.(Y(IVE).EQ.Y0)) THEN
          ISINE2 = .TRUE.
          RETURN
        END IF
      END DO


      DO IVE=1,NVE
      
        IVE2=MOD(IVE,NVE)+1
        DAREA=(X0-X(IVE))*(Y0+Y(IVE))+(X(IVE)-X(IVE2))*(Y(IVE)+Y(IVE2))
     *       +(X(IVE2)-X0)*(Y(IVE2)+Y0)
        
C       Cancel if one of the triangles has negative volume.
        
        IF (DAREA.LT.0D0) RETURN
        
      END DO
      
      ISINE2 = .TRUE.
      
      END
      
************************************************************************
* Calculate the intersection between two lines
*
* This routine calculates the intersection point (X,Y) of two
* lines [X1,Y1]->[X2,Y2] and [X3,Y3]->[X4,Y4].
* 
* In:
*   X0,Y0
*   X1,Y1   - two different points on the first line
*   X2,Y2
*   Y3,Y3   - two different points on the second line
*
* Out:
*   X, Y    - Intersection point
*   IRES    - =0: There is an intersection point
*             =1: The two lines are identical; X,Y is undefined
*             =2: The two lines do not intersect; X,Y is undefined
************************************************************************

      SUBROUTINE LINSCT (X0,Y0,X1,Y1,X2,Y2,X3,Y3,X,Y,IRES)
     
      IMPLICIT NONE
      
      DOUBLE PRECISION X0,Y0,X1,Y1,X2,Y2,X3,Y3,X,Y
      INTEGER IRES

C     local variables

      DOUBLE PRECISION DET,A
      
      IRES = 0
      
C     We have (hopefully) the situation
C                 
C                    (X1,Y1)
C                       |
C                       |
C       (X2,Y2) --------+--------- (X3,Y3)
C                       |
C                       |
C                    (X0,Y0)
C
C      and want to calculate the intersection point. This means
C      we have to solve the linear system
C
C       ( X1-X0  X2-X3 ) * (a) = ( X2-X0 )
C       ( Y1-Y0  Y2-Y3 )   (b)   ( Y2-Y0 )
C
C      to get the "parameter" values a,b along the two lines where
C      the intersection occurres.
C
C      The determinant of the system is:

       DET = X1*Y2-X1*Y3-X0*Y2+X0*Y3-Y1*X2+Y1*X3+Y0*X2-Y0*X3 
       
C      If it's =0, the lines are the same or completely different...
        
       IF (DET.EQ.0D0) THEN
       
C        If the vector (X2,Y2)->(X0,Y0) is linear dependent to
C        (X2,Y2)->(X3,Y3), the lines are the same.

          DET = -Y0*X2-X3*Y2+Y0*X3+X2*Y3+X0*Y2-X0*Y3
          IF (DET.EQ.0D0) THEN
            IRES = 1
          ELSE
            IRES = 2
          END IF
       
       ELSE

C       There is an intersection point. Calculate one of the 
C       "parameter" values along the two lines.

        A = (Y0*X2+X3*Y2-Y0*X3-X2*Y3-X0*Y2+X0*Y3)/DET
        
C       The intersection point is then

        X = A*X1 + (1D0-A)*X0
        Y = A*Y1 + (1D0-A)*Y0
       
       END IF
      
      END      

**********************************************************************
* Search for an element containing the point (x0,y0)
*
* Does a linear search on all elements until the one containing
* the position (x0,y0) is found.
*
* In:
*   X0,Y0  - Position of an arbitrary point
*   DCORVG - coordinates of all points
*   KVERT  - Vertex numbers on all elements
*   NEL    - Number of elements in the triangulation
*
* Out:
*   IEL    - Number of the element containing the point.
*            =0, if no element contains the point.
**********************************************************************
  
      SUBROUTINE PSRCH2(X0,Y0,DCORVG,KVERT,NEL,IEL)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'cbasictria.inc'
                           
C parameters
      
      DOUBLE PRECISION X0,Y0
      INTEGER IEL,NEL

      DOUBLE PRECISION DCORVG(2,*)
      INTEGER KVERT(NNVE,*)

C local variables

      INTEGER JEL
      LOGICAL ISINEL

      DO JEL=1,NEL
        IF (ISINEL(X0,Y0,JEL,DCORVG,KVERT)) THEN
          IEL = JEL
          RETURN
        END IF
      END DO

      IEL = 0

      END

**********************************************************************
* Extended ray tracing element search
*
* Purpose: returns the element containing the point (x0,y0) by a
* raytracing method
* Constructs a ray from the old node position (whose element number
* is known) and traces it to the current position.
*
* In contrast to the standard version of the ray-tracing search,
* this routine does not completely stop the program if no element is
* found. Instead it will tell the caller what happened, i.e.:
* - if an element was found
* - if no element was found in the inner
* - if the point left the domain, through which boundary component
*   and which element
*
* In:
*  X0,
*  Y0     - Coordinates of the (moved) point to determine the element 
*           number where this point resides in
*  IELOLD - Element in which the not-moved point (XALT,YALT) 
*           resides in
* 
* Out:
*  Return value:
*   = 0, if the element was found successfully. In this case:
*        IEL - Number of the element that contains (X0,Y0)
*   = 1, if the raytracing search broke down inside of the domain. 
*        In this case:
*        IEL - last element that was analyzed
*   = 2, if the search broke down because the domain was left.
*        In this case:
*        IEL - last analyzed element where the ray left the domain
**********************************************************************

      INTEGER FUNCTION BSRCH4(X0,Y0,
     *                        DCORVG,KVERT,KADJ,
     *                        IELOLD,IEL)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INTEGER MAXIT
      PARAMETER (MAXIT=100)
      
C parameters
      
      DOUBLE PRECISION X0,Y0,XALT,YALT
      DOUBLE PRECISION DCORVG(2,*)
      INTEGER KVERT,KADJ
      INTEGER IELOLD,IEL
      
      DIMENSION KVERT(NNVE,*),KADJ(NNVE,*)

C local variables

      DOUBLE PRECISION X1EDGE,Y1EDGE,X2EDGE,Y2EDGE
      INTEGER IVE,IDX,IELPRV
      LOGICAL BAUX,ISINE2,BINTSE

      IEL=IELOLD
      IELPRV=IELOLD
   
C     First a quick-check: If the point is in element IELOLD,
C     we are finished.

      IF ((IELOLD.NE.0).AND.(ISINE2(X0,Y0,IEL,DCORVG,KVERT))) THEN
        BSRCH4 = 0
        RETURN
      END IF

C     We restrict our raytracing search to 100 neighbour cells;
C     let's hope a point is not moved more than 100 elements in 
C     one step...

      DO IDX=1,MAXIT
      
C       Calculate point of intersection with one of the four
C       surrounding edges, take the center of the starting element
C       as starting point  of the ray

        XALT=0.5D0*(DCORVG(1,KVERT(1,IEL))+DCORVG(1,KVERT(3,IEL)))
        YALT=0.5D0*(DCORVG(2,KVERT(1,IEL))+DCORVG(2,KVERT(3,IEL)))
        
C       Check all edges to find an intersection point:
        
        DO IVE=1,NNVE
        
C         don't jump back to the old element

          IF (KADJ(IVE,IEL).NE.IELPRV) THEN
            X1EDGE=DCORVG(1,KVERT(IVE,IEL))
            Y1EDGE=DCORVG(2,KVERT(IVE,IEL))
            X2EDGE=DCORVG(1,KVERT(MOD(IVE,NNVE)+1,IEL))
            Y2EDGE=DCORVG(2,KVERT(MOD(IVE,NNVE)+1,IEL))
            BAUX=BINTSE(XALT,YALT,X0,Y0,X1EDGE,Y1EDGE,X2EDGE,Y2EDGE)

C           BAUX is TRUE if we found a point of intersection.

            IF (BAUX) THEN
            
C             Go on searching in the neighbour element

              IELPRV=IEL
              
C             Stop here if we leave our domain. To test that use the
C             fact that KADJ()=0 if there is no neighbour element!
C             The caller can use IEL to find out the information
C             about the boundary where the domain was left...

              IF (KADJ(IVE,IELPRV).EQ.0) THEN
                BSRCH4 = 2
                RETURN
              END IF
              
              IEL=KADJ(IVE,IELPRV)
              GOTO 210
              
            END IF
            
          END IF
          
        END DO ! IVE
        
C       We can't find a point of intersection!
C       This means, the line segment from the element midpoint
C       to the target mitpoint does not leave out element -
C       or in other words: we found the element that contains
C       the point!
       
C        IF (ISINE2(X0,Y0,IEL,DCORVG,KVERT)) THEN
          BSRCH4 = 0
          RETURN
C        ELSE
C          WRITE (*,*) 'BSRCH4: Can''t find point in element!?!'
C        END IF
        
210     CONTINUE

      END DO

C If the search is successful, we cancel the DO-loop and jump out
C of this routine; otherwise we end up here because we have not
C been able to find the element that contains the point...
      
      BSRCH4 = 1

      END 

**********************************************************************
* Extended ray tracing element search
*
* Purpose: returns the element containing the point (x0,y0) by a
* raytracing method.
* Constructs a ray from the old node position (whose element number
* is known) and traces it to the current position.
*
* This routine works very similar to BSRCH4. It has the following
* differences:
*  - It's a subroutine, not a function
*  - When the domain is left, the routine reports
*    - through which element the domain was left
*    - through which edge the domain was left
*
* In:
*  X0,
*  Y0     - Coordinates of the (moved) point to determine the element 
*           number where this point resides in
*  IELOLD - Element in which the not-moved point (XALT,YALT) 
*           resides in
*  DCORVG,
*  KVERT,
*  KMID,
*  KADJ,
*  NVT    - Usual geometry quantities
* 
* Out:
*  IRES = 0, if the element was found successfully. In this case:
*            IEL  - Number of the element that contains (X0,Y0)
*            IEDG = 0
*            X1 = Y1 = X2 = Y2 = 0
*       = 1, if the raytracing search broke down inside of the domain. 
*            In this case:
*            IEL - last element that was analyzed
*            IEDG = 0
*            X1 = Y1 = X2 = Y2 = 0
*       = 2, if the search broke down because the domain was left.
*            In this case:
*            IEL  - last analyzed element where the ray left 
*                   the domain
*            IEDG - Global number of the edge through which the domain
*                   was left; the value is i the range 1..NMT
*            X1,Y1 - starting point of the edge through which the 
*                    domain was left
*            X2,Y2 - ending point of the edge through which the 
*                    domain was left
**********************************************************************

      SUBROUTINE PSRCH4(X0,Y0,
     *                  DCORVG,KVERT,KMID,KADJ,NVT,
     *                  IELOLD,IEL,
     *                  IRES,IEDG,X1,Y1,X2,Y2)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INTEGER MAXIT
      PARAMETER (MAXIT=100)
      
C parameters
      
      DOUBLE PRECISION X0,Y0,XALT,YALT,X1,Y1,X2,Y2
      DOUBLE PRECISION DCORVG(2,*)
      INTEGER IELOLD,IEL,IRES,IEDG
      
      INTEGER KVERT(NNVE,*),KADJ(NNVE,*),KMID(NNVE,*),NVT

C local variables

      DOUBLE PRECISION X1EDGE,Y1EDGE,X2EDGE,Y2EDGE
      INTEGER IVE,IDX,IELPRV
      LOGICAL BAUX,ISINE2,BINTSE

      IEL=IELOLD
      IELPRV=IELOLD
      
C     Initialise output variables
        
      IRES = 0
      IEDG = 0
      X1 = 0
      Y1 = 0
      X2 = 0
      Y2 = 0
   
C     First a quick-check: If the point is in element IELOLD,
C     we are finished.

      IF ((IELOLD.NE.0).AND.(ISINE2(X0,Y0,IEL,DCORVG,KVERT))) THEN
        RETURN
      END IF

C     We restrict our raytracing search to 100 neighbour cells;
C     let's hope a point is not moved more than 100 elements in 
C     one step...

      DO IDX=1,MAXIT
      
C       Calculate point of intersection with one of the four
C       surrounding edges, take the center of the starting element
C       as starting point  of the ray

        XALT=0.5D0*(DCORVG(1,KVERT(1,IEL))+DCORVG(1,KVERT(3,IEL)))
        YALT=0.5D0*(DCORVG(2,KVERT(1,IEL))+DCORVG(2,KVERT(3,IEL)))
        
C       Check all edges to find an intersection point:
        
        DO IVE=1,NNVE
        
C         don't jump back to the old element

          IF (KADJ(IVE,IEL).NE.IELPRV) THEN
            X1EDGE=DCORVG(1,KVERT(IVE,IEL))
            Y1EDGE=DCORVG(2,KVERT(IVE,IEL))
            X2EDGE=DCORVG(1,KVERT(MOD(IVE,NNVE)+1,IEL))
            Y2EDGE=DCORVG(2,KVERT(MOD(IVE,NNVE)+1,IEL))
            BAUX=BINTSE(XALT,YALT,X0,Y0,X1EDGE,Y1EDGE,X2EDGE,Y2EDGE)

C           BAUX is TRUE if we found a point of intersection.

            IF (BAUX) THEN
            
C             Go on searching in the neighbour element

              IELPRV=IEL
              
C             Stop here if we leave our domain. To test that use the
C             fact that KADJ()=0 if there is no neighbour element!
C             The caller can use IEL to find out the information
C             about the boundary where the domain was left...

              IF (KADJ(IVE,IELPRV).EQ.0) THEN
              
                IRES = 2
                
C               What's the edge where we left the domain?

                IEDG = KMID(IVE,IELPRV)-NVT
                
C               Which start/endpoint does this edge have?

                X1 = DCORVG(1,KVERT(IVE,IELPRV))
                Y1 = DCORVG(2,KVERT(IVE,IELPRV))
                X2 = DCORVG(1,KVERT(MOD(IVE,4)+1,IELPRV))
                Y2 = DCORVG(2,KVERT(MOD(IVE,4)+1,IELPRV))
                
                RETURN
                
              END IF
              
              IEL=KADJ(IVE,IELPRV)
              GOTO 210
              
            END IF
            
          END IF
          
        END DO ! IVE
        
C       We can't find a point of intersection!
C       This means, the line segment from the element midpoint
C       to the target mitpoint does not leave out element -
C       or in other words: we found the element that contains
C       the point!
       
C        IF (ISINE2(X0,Y0,IEL,DCORVG,KVERT)) THEN
          IRES = 0
          RETURN
C        ELSE
C          WRITE (*,*) 'BSRCH4: Can''t find point in element!?!'
C        END IF
        
210     CONTINUE

      END DO

C If the search is successful, we cancel the DO-loop and jump out
C of this routine; otherwise we end up here because we have not
C been able to find the element that contains the point...
      
      IRES = 1

      END 

**********************************************************************
* Hierarchical element search
*
* This routine uses a hierarchical search to find the element number
* that contains a point (X0,Y0). The caller can specify a start-
* and an end-level where to search. Inside of each level, a raytracing
* technique is used to find the element.
* The caller can specify an element on the start-level where to start
* the search. If no element is specified, a linear search will be done
* on all elements on the start-level to find the initial element.
*
* In:
*   TRIAS  - array [1..SZTRIA,1..*] of integer
*            Triangulation structures on all relevant levels.
*   ILVMIN - Minimum level where to start the search
*   ILVMAX - Maximum level where to find the element number
*   X0,Y0  - Coordinates of the point that is to be searched
*   IELOLD - Initial element number on level ILVMIN where to start
*            the search.
*            =0, if the initial element number is unknown.
*
* Out:
*   IEL    - Element number on level ILVMAX that contains the
*            point (X0,Y0);
*            =0, if the element could not be found.
*
* Remark, IELOLD=IEL on call of the routine is allowed.
**********************************************************************

      SUBROUTINE PSRCH5(TRIAS,ILVMIN,ILVMAX,X0,Y0,IELOLD,IEL)

      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'stria.inc'
      
C     parameters      
      
      INTEGER TRIAS(SZTRIA,*),ILVMIN,ILVMAX,IEL,IELOLD
      DOUBLE PRECISION X0,Y0
      
C     externals

      INTEGER BSRCH4
      EXTERNAL BSRCH4
      
C     local variables

      INTEGER ILEV,IELTMP
      
      IEL = IELOLD
      
C     If the initial element number is unknown...

      IF (IEL.LE.0) THEN

C       Try to make a raytracing search from element number 1.

        IELTMP = 1
        IF (BSRCH4(X0,Y0,
     *             DWORK(L(TRIAS(OLCORVG,ILVMIN))),
     *             KWORK(L(TRIAS(OLVERT,ILVMIN))),
     *             KWORK(L(TRIAS(OLADJ,ILVMIN))),
     *             IELTMP,IEL).NE.0) THEN
     
C         Oops, raytracing failed.
C         Set the element number back to 0. The loop
C         will then make a linear search through all
C         elements to find one that contains the point.
      
          IEL = 0
     
        END IF ! BSRCH4 <> 0
        
C       If an element is found, good - then we save a (costly)
C       linear search!
        
      END IF ! IEL = 0
      
      DO ILEV=ILVMIN,ILVMAX
      
C       Detetmine the initial element that contains the point.
C       If the element is given, ok. If not, make a linear search
C       on the coarsest level to find it.

        IF (IEL.EQ.0) THEN
          CALL PSRCH2(X0,Y0,DWORK(L(TRIAS(OLCORVG,ILEV))),
     *                KWORK(L(TRIAS(OLVERT,ILEV))),TRIAS(ONEL,ILEV),
     *                IEL)
        
          IF (IEL.EQ.0) THEN
          
C           Oops, linear search failed. Then we start at element 1.
C           Don't stop here, as in "circle" domains it may happen that
C           the point on the higher level is in the domain while
C           the point on the lower level is not!
C           Let's hope that the domain is somehow convex so that
C           raytracing works...

            IEL = 1
          
          END IF
        END IF
        
C       Start the raytracing technique to find a better element
      
        IELTMP = IEL
        IF (BSRCH4(X0,Y0,
     *             DWORK(L(TRIAS(OLCORVG,ILEV))),
     *             KWORK(L(TRIAS(OLVERT,ILEV))),
     *             KWORK(L(TRIAS(OLADJ,ILEV))),
     *             IELTMP,IEL).NE.0) THEN
     
C         Oops, raytracing failed.
C         Set the element number back to 0 and start on the
C         next higher level with a linear search
      
          IEL = 0
          
        END IF ! BSRCH4 <> 0
        
C       Otherwise: Element number of lower level found.
C       Take this as initial element on the higher level.
C       Because of the 2-level ordering, this element
C       number is usually near to the element on the finer
C       level containing that point!
        
      END DO ! ILEV
      
C     If IEL<>0, this is the last element number on the highest
C     level that contains (X0,Y0)!
      
      END
      