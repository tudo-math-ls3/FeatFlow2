************************************************************************
* This file contains routines for reconstruction of fictitious
* boundaries - i.e. construction of an approximation of the
* analytical boundary of a fictitious boundary component.
************************************************************************

************************************************************************
* Search fictitious boundary intersection point
*
* This routine takes the coordinates of two points, assuming one point
* to be in the fictitious boundary component and one point outside.
* It performs N bisection steps to approximate the point where
* line between the points intersects with the fictitious boundary.
*
* In:
*  X1,Y1 - 1st point 
*  X2,Y2 - 2nd point
*  N     - Number of bisection steps to perform
*  IFBC  - Number of fictitious boundary component that is tested.
*          =0: automatically determine correct component.
*  BIN   - true:  the calculated point will be inside of the fictitious
*                 boundary domain
*          false: the calculated point will be outside of the
*                 fictitious boundary domain
*  IGEOM  - array [1..*] of integer 
*  DGEOM  - array [1..*] of double 
*           Integer- and double-precision parameter blocks with
*           geometry information. Passed to fictitious boundary
*           routines. Not used in this routine.
*
* Out:
*  X , Y - Approximation of intersection point of fictitious boundary
*          with the line [(X1,Y1),(X2,Y2)] - inside or outside,
*          depending on BIN.
*  BFOUND - true, if the approximation was successful,
*           false, if no intersection pont was found; e.g. because
*                  both line endpoints have been on the same side of
*                  the boundary domain
*  IPTIN  - 1=(X1,Y1) is in the fictitious boundary component,
*             (X2,Y2) is outside
*           2=(X2,Y2) is in the fictitious boundary component,
*             (X1,Y1) is outside
************************************************************************

      SUBROUTINE SFBISC (X1,Y1,X2,Y2,N,IFBC,BIN,X,Y,BFOUND,IPTIN,
     *                   IGEOM,DGEOM,IINFO,DINFO)
      
      IMPLICIT NONE
      
C parameters
      
      DOUBLE PRECISION X1,Y1,X2,Y2,X,Y,DGEOM(*),DINFO(*)
      INTEGER N,IFBC,IPTIN,IGEOM(*),IINFO(*)
      LOGICAL BIN,BFOUND
      
C externals

      INTEGER ISFBDY
      EXTERNAL ISFBDY
      
C local variables

      DOUBLE PRECISION XN1,YN1,XN2,YN2,XM,YM
      INTEGER I
      
      IF ((ISFBDY(X1,Y1,IFBC,IGEOM,DGEOM,IINFO,DINFO).EQ.0).EQV.
     *    (ISFBDY(X2,Y2,IFBC,IGEOM,DGEOM,IINFO,DINFO).EQ.0)) THEN
C No, that way we surely don't find a point
        BFOUND = .FALSE.
        GOTO 99998
      END IF     
      
      BFOUND = .TRUE. 
      
C Rearrange such that (XN1,YN1) is outside

      IF (ISFBDY(X1,Y1,IFBC,IGEOM,DGEOM,IINFO,DINFO).EQ.0) THEN
        IPTIN = 2
        XN1 = X1
        YN1 = Y1
        XN2 = X2
        YN2 = Y2
      ELSE 
        IPTIN = 1
        XN1 = X2
        YN1 = Y2
        XN2 = X1
        YN2 = Y1
      END IF      
  
C Perform bisection steps

      DO I=1,N
      
C Find the midpoint

        XM = 0.5*(XN1+XN2)
        YM = 0.5*(YN1+YN2)
        
C Move the correct line endpoint, depending on the status of the 
C current midpoint

        IF (ISFBDY(XM,YM,IFBC,IGEOM,DGEOM,IINFO,DINFO).NE.0) THEN
          XN2 = XM
          YN2 = YM
        ELSE
          XN1 = XM
          YN1 = YM
        END IF
      
      END DO

C Return either the left or the right endpoint of the current line,
C depending on BIN

99998 IF (BIN) THEN
        X = XN2
        Y = YN2
      ELSE
        X = XN1
        Y = YN1
      END IF
      
99999 END

************************************************************************
* Search fictitious boundary intersection point - parametrised version
*
* This routine takes the coordinates of two points, assuming one point
* to be in the fictitious boundary component and one point outside.
* It performs N bisection steps to approximate the point where
* line between the points intersects with the fictitious boundary.
* the return value will be a parameter value LAMBDA in [0,1]
* identifying the intersection point on the line [(X1,Y1),(X2,Y2)]
*
* In:
*  X1,Y1 - 1st point 
*  X2,Y2 - 2nd point
*  N     - Number of bisection steps to perform
*  IFBC  - Number of fictitious boundary component that is tested.
*          =0: automatically determine correct component.
*  BIN   - true:  the calculated point will be inside of the fictitious
*                 boundary domain
*          false: the calculated point will be outside of the
*                 fictitious boundary domain
*  IGEOM  - array [1..*] of integer 
*  DGEOM  - array [1..*] of double 
*           Integer- and double-precision parameter blocks with
*           geometry information. Passed to boundary
*           routines. Not used in this routine.
*
* Out:
*  BIN   - true:  the calculated point will be inside of the fictitious
*                 boundary domain
*          false: the calculated point will be outside of the
*                 fictitious boundary domain
*
* Return value:
*    The parameter value LAMBDA of the intersection point or 
*    -1D0 if no intersection point was found; e.g. because
*                  both line endpoints have been on the same side of
*                  the boundary domain
*
* The intersection point itself can be calculated by:
*            (X,Y) = LAMBDA*(X1,Y1) + (1-LAMBDA)*(X2,Y2)
************************************************************************

      DOUBLE PRECISION FUNCTION SFBISP (X1,Y1,X2,Y2,N,IFBC,BIN,IPTIN,
     *                                  IGEOM,DGEOM,IINFO,DINFO)
      
      IMPLICIT NONE
      
C parameters
      
      DOUBLE PRECISION X1,Y1,X2,Y2,DGEOM(*),DINFO(*)
      INTEGER N,IFBC,IPTIN,IGEOM(*),IINFO(*)
      LOGICAL BIN
      
C local variables

      DOUBLE PRECISION X,Y,LEN1,LEN2
      LOGICAL BFOUND
      
C Calculate the intersection point

      CALL SFBISC (X1,Y1,X2,Y2,N,IFBC,BIN,X,Y,BFOUND,IPTIN,
     *             IGEOM,DGEOM,IINFO,DINFO)
      
C Calculated the parameter value of that

      IF (BFOUND) THEN
        LEN1 = DSQRT((X2-X1)**2+(Y2-Y1)**2)
        LEN2 = DSQRT((X-X1)**2+(Y-Y1)**2)
        
        SFBISP = 0D0
        IF (LEN1.NE.0D0) SFBISP = LEN2/LEN1
      ELSE
        SFBISP = -1D0
      END IF
      
      END
      
************************************************************************
* Reconstruct fictitious boundary linearly
*
* This routine creates a linear approximation of the fictitious
* boundary on the element IEL on the current parametrisation.
* It returns the endpoints of the line on the edges of the
* quadrilateral element IEL.
*
* For proper later handling in integration routines, the returned
* starting- and ending points are *oriented*. The vector
* (X2-X1,Y2-Y1) points into the tangential direction of the boundary
* interface. Its normal counterpart (-Y2+Y1,X2-X1) points to the
* outside of the fictitious boundary domain.
*
* In:
*  DCORVG  - Array with coordinates of the vertices
*  KVERT   - Array with numbers of vertices for each element
*  IEL     - The element where to reconstruct the fictitious boundary
*  IFBC    - Number of fictitious boundary component that is tested.
*            =0: automatically determine correct component.
*  BIN     - whether the linear approximation should be on the inside
*            of the fictitious boundary domain or not
*  N       - Number of approximation steps on each edge to find an
*            approximation for the line endpoints
*  IGEOM  - array [1..*] of integer 
*  DGEOM  - array [1..*] of double 
*           Integer- and double-precision parameter blocks with
*           geometry information. Passed to boundary
*           routines. Not used in this routine.
*
* Out:
*  X1,Y1   - Starting point of the line approximating the fictitious
*            boundary on IEL
*  IEDGE1  - Number of the edge (1..4) containing X1,Y1
*  X2,Y2   - Ending point of the line approximating the fictitious
*            boundary on IEL
*  IEDGE2  - Number of the edge (1..4) containing X2,Y2
*  BFOUND  - true, of the approximation was successful
*            false, if the element does not intersect the fictitious
*            boundary; (X1,Y1) and (X2,Y2) is undefined in this case!
*
* KVERT(IEDGE1,.) and KVERT(IEDGE2,.) are the number of the vertices
* preceding the appropriate point X1,Y1 and X2,Y2, respectively,
* on the appropriate edge.
************************************************************************

      SUBROUTINE RCFBLI (DCORVG, KVERT, IEL, N, IFBC, BIN,
     *                   X1,Y1, IEDGE1, X2,Y2, IEDGE2, BFOUND,
     *                   IGEOM,DGEOM,IINFO,DINFO)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      
C parameters
      
      INTEGER IEL,IFBC,N,IEDGE1,IEDGE2,IGEOM(*),IINFO(*)
      DOUBLE PRECISION DCORVG(2,*),X1,Y1, X2,Y2, DTMP,DGEOM(*),DINFO(*)
      INTEGER KVERT(NNVE,*)
      LOGICAL BFOUND,BIN
      
C local variables

      INTEGER I,J,K1,K2,IPTIN
      
C Loop through all edges to find the first intersection point

      DO I=1,NNVE
        K1 = KVERT(I,IEL)
        K2 = KVERT(MOD(I,NNVE)+1,IEL)
        CALL SFBISC (DCORVG(1,K1),DCORVG(2,K1),
     *       DCORVG(1,K2),DCORVG(2,K2),N,IFBC,BIN,
     *       X1,Y1,BFOUND,IPTIN,IGEOM,DGEOM,IINFO,DINFO)
        IF (BFOUND) GOTO 10
      END DO
      
C Oops, no intersection point was found!
C BFOUND is FALSE - stop searching here.

      GOTO 99999

C Store the (local) number of the edge.
C Go on searching for the second intersection point
      
10    IEDGE1 = I

      DO J=I+1,NNVE
        K1 = KVERT(J,IEL)
        K2 = KVERT(MOD(J,NNVE)+1,IEL)
        CALL SFBISC (DCORVG(1,K1),DCORVG(2,K1),
     *       DCORVG(1,K2),DCORVG(2,K2),N,IFBC,BIN,
     *       X2,Y2,BFOUND,IPTIN,IGEOM,DGEOM,IINFO,DINFO)
        IF (BFOUND) GOTO 99998
      END DO

      GOTO 99999

C Store the (local) number of the edge

99998 IEDGE2 = J

C Now we have found two intersection points. Great question: how are they
C oriented? Well, that's rather easy if we consider the last calculated
C IPTIN and the anti-clockwise orientation of the element itself...

C If (IPTIN.EQ.1), the edge IEDGE2 starts with a point in the fictitious 
C boundary and ends with a point outside of it:
C
C 4-----3
C |   --+-  IEDGE2
C |  /XX|
C 1-+XXX2
C   IEDGE1
C
C So (X2,Y2) must be the ending point. This is the case where no coordinate
C exchange has to be done. In the other case...
C
C 4-----3
C |XXX/-+-  IEDGE2
C |XX/  |
C 1-+---2
C   IEDGE1
C
C ... we have to exchange the coordinates

      IF (IPTIN.EQ.2) THEN
        DTMP = X1
        X1 = X2
        X2 = DTMP
        
        DTMP = Y1
        Y1 = Y2
        Y2 = DTMP
      END IF

C That's all!
C If the 2nd intersection-point is found, BFOUND=true and the 2nd coordinates
C are known - otherwise BFOUND=false and (X2,Y2) is undefined!

99999 END

************************************************************************
* RCBGMV - Write reconstructed fictitious boundary interface to GMV file
*
* This routine performs a loop about all cells to find the reconstructed
* fictitious boundary. It writes out the line approximating the
* interface to file unit MUNIT using GMV-material MMAT,...,
* MMAT+NFBDRY(). For every fictitious boundary component another
* material is used.
*
* The routine assumes that the POLYGONs section of the GMV-file is
* administrated by the caller, i.e. only the content of this section
* is written to the file.
*
* In:
*  MUNIT  - File unit where to write the GMV output to.
*  MMAT   - The first material number to use.
*  TRIA   - array [1..SZTRIA] of integer
*           Triangulation structure on where the interface should
*           be reconstructed
*  IGEOM  - array [1..*] of integer 
*  DGEOM  - array [1..*] of double 
*           Integer- and double-precision parameter blocks with
*           geometry information. Passed to boundary
*           routines. Not used in this routine.
************************************************************************

      SUBROUTINE RCBGMV (MUNIT,MMAT,TRIA,IGEOM,DGEOM,IINFO,DINFO)

      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      INCLUDE 'cmem.inc'

C parameters

      INTEGER MUNIT,MMAT,IGEOM(*),TRIA(SZTRIA),IINFO(*)
      DOUBLE PRECISION DGEOM(*),DINFO(*)
      
C externals

      INTEGER NFBDYC,ISFBDY
      EXTERNAL NFBDYC,ISFBDY
      
C local variables

      INTEGER I,NL
      
      INTEGER IEDGE1,IEDGE2
      DOUBLE PRECISION X1,Y1,X2,Y2
      LOGICAL BFOUND
      
      INTEGER LCORVG,LVERT

C     Abort if there's no fictitious boundary component

      IF (NFBDYC(IGEOM,DGEOM).EQ.0) RETURN
      
C     Fetch some triangulation information

      LCORVG = TRIA(OLCORVG)
      LVERT  = TRIA(OLVERT)

C Write reconstructed fictitious boundary interface to file.
C The POLYGON's section is already opened by the caller.

C Ok, this is a little bit time intensive...
C Loop through all elements, test if the line intersects with the
C interface, write the line out:

      NL = 0
      DO I=1,TRIA(ONEL)
C Approximate the fictitious boundary by 8 approximation steps.
C Force the approximate points to be in the fictitious boundary
C component.
        CALL RCFBLI (DWORK(L(LCORVG)), KWORK(L(LVERT)), I, 8, 0, .TRUE.,
     *               X1,Y1, IEDGE1, X2,Y2, IEDGE2, BFOUND,
     *               IGEOM,DGEOM,IINFO,DINFO)
     
        IF (BFOUND) THEN

C Material number MMT+number of boundary component; 
C 2 nodes (a line with two endpoints).
C We forced both endpoints to be in the fictitious boundary component,
C so we can make a quick-and-dirty test which component it ment

          WRITE(MUNIT,*) MMAT+
     *                   ABS(ISFBDY (X1,Y1,0,IGEOM,DGEOM,IINFO,DINFO))-1
     *                   -TRIA(ONBCT) , 2
          WRITE(MUNIT,*) X1
          WRITE(MUNIT,*) X2
          WRITE(MUNIT,*) Y1
          WRITE(MUNIT,*) Y2
          WRITE(MUNIT,*) 0D0
          WRITE(MUNIT,*) 0D0

        END IF
      END DO

      END
