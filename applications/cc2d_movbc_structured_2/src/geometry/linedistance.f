************************************************************************
* This file contains basic routines to support calculating
* distances from a couple of points to a polygon given by straight
* lines.
* The routines are very basic ones without special handling by
* advanced mathematical routines (like Fast Marching,...)
************************************************************************

************************************************************************
* Redistancing via line segments (brute force method)
*
* This routine uses a brute force approach to calculate for all given
* points the distance to a polygon, which is given by line segments.
* NVT describes the number of points, DCORVG the (X/Y)-coordinates
* of the points. The distances of the points are saved into an
* array DU.
*
* In:
*  NVT    - Number of vertices
*  DCORVG - array [1..2,1..NVT] of double.
*           (X/Y)-coordinates of the vertices.
*           DCORVG(1) = X-Coordinate
*           DCORVG(2) = Y-Coordinate
*  NPLNSG - Number of points in the line segment array
*  DCORSG - array [1..2,1..NLINSG] of double
*           (X/Y)-coordinates of all the points defining the lines.
*           DCORSG(1) = X-Coordinate
*           DCORSG(2) = Y-Coordinate
*           This array represents a list of endpoints that define the
*           (NLINSG-1) line segments. The coordinates are stores as 
*           pairs of coordinates
*                   [(x1,y1), (x2,y2), ... ,(xn,yn)]
*           with (x1,y1) the first (and thus starting) point of 
*           the polygon and (xn,yn) (n=NLINSG) the last point.
*           The polygon itself is given by the lines between the
*           points: (x1,y1)->(x2,y2), (x2,y2)->(x3,y3), ...
*
* Out:
*  DU      - array [1..NVT] of double
*            DU(I) receives the minimum distance from vertex I to
*            all line segments.
*  IPROJ   - array [1..NVT] of integer
*            IPROJ(I) is the number of the line segment that defines
*            the minimum distance to vertex I
*  DPROJ   - array [1..2,1..NVT] of double
*            (DPROJ(1,I),DPROJ(2,I)) contains the projection of node I
*            onto the polygon.
*  DNML    - array [1..2,1..NVT] of double
*            Array with normal vectors pointing from the polygon to
*            the nodes. (DNML(1,I),DNML(2,I)) is a vector pointing from
*            the projected point on line segment IPROJ(I) onto node I.
*            Remark: The normal vectors are not necessarily outer
*             normals!!!
************************************************************************

      SUBROUTINE RDPLSB(NVT,DCORVG,NPLNSG,DCORSG,DU,IPROJ,DPROJ,DNML)
      
      IMPLICIT NONE
      
      INTEGER NVT,NPLNSG
      
      DOUBLE PRECISION DCORVG(2,NVT),DCORSG(2,NPLNSG),DPROJ(2,NVT)
      DOUBLE PRECISION DU(NVT),DNML(2,NVT)
      INTEGER IPROJ(NVT)
      
C local variables
      
      INTEGER IVT,ISG
      DOUBLE PRECISION D
C     DOUBLE PRECISION DAREA
      
C externals

      DOUBLE PRECISION DPO2SG
      EXTERNAL DPO2SG
 
C Cancel if the number of nodes is incorrect
 
      IF (NPLNSG.LE.0) RETURN
 
C Loop through all nodes
 
      DO IVT=1,NVT
      
C For each node loop through all line segments
      
        DO ISG=1,NPLNSG-1
        
C If the distance is lower than the previously calculated...

          D = DPO2SG(DCORVG(1,IVT),DCORVG(2,IVT),
     *               DCORSG(1,ISG),DCORSG(1,ISG+1))
        
          IF ((ISG.EQ.1).OR.(D.LT.DU(IVT))) THEN
          
C save the new best line segment

            DU(IVT) = D
            IPROJ(IVT) = ISG
            
          END IF
          
        END DO
        
C Finally calculate the projection and the normal vector
C from the projected point on the line segment to the original node:
C This is on purpose NOT THE OUTER normal vector!!!

        ISG = IPROJ(IVT)
        CALL PPO2SG(DCORVG(1,IVT),DCORVG(2,IVT),
     *              DCORSG(1,ISG),DCORSG(1,ISG+1),D,DPROJ(1,IVT))
        DNML(1,IVT) = DCORVG(1,IVT)-DPROJ(1,IVT)
        DNML(2,IVT) = DCORVG(2,IVT)-DPROJ(2,IVT)

C Check the sign of the volume of the triangle (point on line 1,
C point on line 2, IVT). If it's positive, we have the outer normal
C vector. If it's negative, we have the "inner" normal vector and must
C change the sign to get the outer one...
C --> This check is disabled because of a different implementation
C     now!!! (see remark above)
C        DAREA=(DCORVG(1,IVT)-DCORSG(1,ISG+1))*
C     *          (DCORVG(2,IVT)+DCORSG(2,ISG+1))
C     *       +(DCORSG(1,ISG+1)-DCORSG(1,ISG))*
C     *          (DCORSG(2,ISG+1)+DCORSG(2,ISG))
C     *       +(DCORSG(1,ISG)-DCORVG(1,IVT))*
C     *          (DCORSG(2,ISG)+DCORVG(2,IVT))
C     
C        IF (DAREA.LT.0) THEN
C          DNML(1,IVT) = -DNML(1,IVT)
C          DNML(2,IVT) = -DNML(2,IVT)
C        END IF

      END DO

      END

