************************************************************************
* This file contains basic transformation routines for lines and
* line segments.
************************************************************************

************************************************************************
* Determine minimum (Euclidean) distance from a point (DPX,DPY) to 
* a straight line segment. 
*
* In:
*  DPX    - x-coodinate of point
*  DPY    - y-coodinate of point
*  PSTART - starting point of the line segment.
*           PSTART(1) = X-coordinate
*           PSTART(2) = Y-coordinate
*  PEND   - endpoint of the line segment.
*           PEND(1) = X-coordinate
*           PEND(2) = Y-coordinate
*
* Out:
*  Return value = minimum distance
************************************************************************

      DOUBLE PRECISION FUNCTION DPO2SG(DPX,DPY,PSTART,PEND)
      
      IMPLICIT NONE
      DOUBLE PRECISION DPX,DPY,PSTART(2),PEND(2)
      
C local variables

      DOUBLE PRECISION DSX1,DSY1,DSX2,DSY2,DVX,DVY,DWX,DWY
      DOUBLE PRECISION DBX,DBY,DB
      DOUBLE PRECISION DVL,DPL

      DSX1 = PSTART(1)
      DSY1 = PSTART(2)
      DSX2 = PEND(1)
      DSY2 = PEND(2)

C Vector PSTART -> PEND and its length

      DVX = DSX2-DSX1
      DVY = DSY2-DSY1
      DVL = SQRT(DVX*DVX+DVY*DVY)

C Vector PSTART -> (DPX,DPY) 

      DWX = DPX-DSX1
      DWY = DPY-DSY1

C Scalar product to calculate the length of the
C projected vector

      DPL = DWX*DVX+DWY*DVY
      
C Relative length of the projection:

      DB  = DPL/(DVL*DVL)
      
C Relative Length <= 0 
C  => Connection between PSTART and (DPX,DPY) is shortest distance
C                
C Relative Length >= 1
C  => Connection between PEND and (DPX,DPY) is shortest distance
      
      IF (DB.LE.0D0) THEN
        DPO2SG = SQRT((DPX-DSX1)**2+(DPY-DSY1)**2)
      ELSE IF (DB.GE.1D0) THEN
        DPO2SG = SQRT((DPX-DSX2)**2+(DPY-DSY2)**2)
      ELSE
C Calculate the projection and the distance to that point.
C Remember to take the square root of the distances as they 
C are squared...
        DBX = DSX1+DVX*DB
        DBY = DSY1+DVY*DB
        DPO2SG = SQRT((DPX-DBX)**2+(DPY-DBY)**2)
      ENDIF

      END

************************************************************************
* Project a point (DPX,DPY) onto a line segment.
*
* In:
*  DPX    - x-coodinate of point
*  DPY    - y-coodinate of point
*  PSTART - starting point of the line segment.
*           PSTART(1) = X-coordinate
*           PSTART(2) = Y-coordinate
*  PEND   - endpoint of the line segment.
*           PEND(1) = X-coordinate
*           PEND(2) = Y-coordinate
*
* Out:
*  DPROJ  - Projection of (DPX,DPY) onto the line segment.
*           DPROJ(1) = X-coordinate
*           DPROJ(2) = Y-coordinate
*  DPARM  - Parameter value of the projected point along the line
*           segment. Values in [0,1].
*           DPROJ = PSTART + DPROJ*(PEND-PSTART)
************************************************************************

      SUBROUTINE PPO2SG(DPX,DPY,PSTART,PEND,DPARM,DPROJ)
      
      IMPLICIT NONE
      DOUBLE PRECISION DPX,DPY,PSTART(2),PEND(2),DPROJ(2),DPARM
      
C local variables

      DOUBLE PRECISION DSX1,DSY1,DSX2,DSY2,DVX,DVY,DWX,DWY
      DOUBLE PRECISION DVL, DPL

      DSX1 = PSTART(1)
      DSY1 = PSTART(2)
      DSX2 = PEND(1)
      DSY2 = PEND(2)

C Vector PSTART -> PEND and its length

      DVX = DSX2-DSX1
      DVY = DSY2-DSY1
      DVL = SQRT(DVX*DVX+DVY*DVY)

C Vector PSTART -> (DPX,DPY) 

      DWX = DPX-DSX1
      DWY = DPY-DSY1
      
C Scalar product to calculate the length of the
C projected vector 

      DPL = DWX*DVX+DWY*DVY

C Relative length of our projected vector:
      
      DPARM  = DPL/(DVL*DVL)

C Length <= 0 
C  => Point is projected onto the starting point PSTART
C                
C Length >= 1
C  => Point is projected onto the endpoint PEND
      
      IF (DPARM.LE.0D0) THEN
        DPARM = 0D0
        DPROJ(1)=DSX1
        DPROJ(2)=DSY1
      ELSE IF (DPARM.GE.1D0) THEN
        DPARM = 1D0
        DPROJ(1)=DSX2
        DPROJ(2)=DSY2
      ELSE
C Calculate the projection and the distance to that point.
C Remember to take the square root of the distances as they 
C are squared...
        DPROJ(1) = DSX1+DVX*DPARM
        DPROJ(2) = DSY1+DVY*DPARM
      ENDIF

      END

************************************************************************
* Determine minimum (Euclidean) distance from a point (DPX,DPY) to 
* a straight line, given by two points on the line. 
*
* In contrast to DPO2SG the line has no starting/ending point!
*
* In:
*  DPX    - x-coodinate of point
*  DPY    - y-coodinate of point
*  PPT1   - one point of the line segment.
*           PPT1(1) = X-coordinate
*           PPT1(2) = Y-coordinate
*  PPT2   - another point of the line segment. PEND != PSTART !!!
*           PPT2(1) = X-coordinate
*           PPT2(2) = Y-coordinate
*
* Out:
*  Return value = minimum distance, >= 0.
*  < 0 indicates an error: PPT1=PPT2.
************************************************************************

      DOUBLE PRECISION FUNCTION DPO2SL(DPX,DPY,PPT1,PPT2)
      
      IMPLICIT NONE
      DOUBLE PRECISION DPX,DPY,PPT1(2),PPT2(2)
      
C local variables

      DOUBLE PRECISION DSX1,DSY1,DSX2,DSY2,DVX,DVY,DWX,DWY
      DOUBLE PRECISION DB,DBX,DBY
      DOUBLE PRECISION DVL,DPL

      DSX1 = PPT1(1)
      DSY1 = PPT1(2)
      DSX2 = PPT2(1)
      DSY2 = PPT2(2)

C Vector PSTART -> PEND and its length DC2

      DVX = DSX2-DSX1
      DVY = DSY2-DSY1
      DVL = SQRT(DVX*DVX+DVY*DVY)
      
C Stop calculation if both points are the same.
      
      IF (DVL.EQ.0D0) THEN
        DPO2SL = -1
        RETURN
      END IF

C Vector PSTART -> (DPX,DPY) 

      DWX = DPX-DSX1
      DWY = DPY-DSY1
      
C Calculate length of projected vector

      DPL = DWX*DVX+DWY*DVY

C Calculate the projection

      DB  = DPL/DVL
      DBX = DSX1+DVX*DB
      DBY = DSY1+DVY*DB
      DPO2SL = SQRT((DPX-DBY)**2+(DPY-DBY)**2)

      END

************************************************************************
* Project a point (DPX,DPY) onto a straight line.
*
* In contrast to PPO2SG the line has no starting/ending point!
*
* In:
*  DPX    - x-coodinate of point
*  DPY    - y-coodinate of point
*  PPT1   - one point of the line segment.
*           PPT1(1) = X-coordinate
*           PPT1(2) = Y-coordinate
*  PPT2   - another point of the line segment. PEND != PSTART !!!
*           PPT2(1) = X-coordinate
*           PPT2(2) = Y-coordinate
*
* Out:
*  DPROJ  - Projection of (DPX,DPY) onto the line.
*           DPROJ(1) = X-coordinate
*           DPROJ(2) = Y-coordinate
*  DPARM  - Parameter value of the projected point along the line
*           segment. 
*           DPROJ = PPT1 + DPROJ*(PPT2-PPT1)
*           DPARM=1D99 indicates an error: PSTART=PEND
************************************************************************

      SUBROUTINE PPO2SL(DPX,DPY,PPT1,PPT2,DPARM,DPROJ)
      
      IMPLICIT NONE
      DOUBLE PRECISION DPX,DPY,PPT1(2),PPT2(2),DPROJ(2),DPARM
      
C local variables

      DOUBLE PRECISION DSX1,DSY1,DSX2,DSY2,DVX,DVY,DWX,DWY
      DOUBLE PRECISION DVL,DPL

      DSX1 = PPT1(1)
      DSY1 = PPT1(2)
      DSX2 = PPT2(1)
      DSY2 = PPT2(2)

C Vector PSTART -> PEND and its length DC2

      DVX = DSX2-DSX1
      DVY = DSY2-DSY1
      DVL = SQRT(DVX*DVX+DVY*DVY)
      
C Stop calculation if both points are the same.
      
      IF (DVL.EQ.0D0) THEN
        DPARM = 1D99
        RETURN
      END IF

C Vector PSTART -> (DPX,DPY) 

      DWX = DPX-DSX1
      DWY = DPY-DSY1
      
C Calculate length of projected vector

      DPL = DWX*DVX+DWY*DVY

C Calculate the projection and the distance to that point

      DPARM  = DPL/(DVL*DVL)
      DPROJ(1) = DSX1+DVX*DPARM
      DPROJ(2) = DSY1+DVY*DPARM
      
      END

************************************************************************
* Test if a point is right-side of a line
*
* This routine tests, if a given point (DPX,DPY) is right-side or
* left-side of a line, given by two points.
*
* In:
*  DPX    - x-coodinate of point
*  DPY    - y-coodinate of point
*  PPT1   - one point of the line segment.
*           PPT1(1) = X-coordinate
*           PPT1(2) = Y-coordinate
*  PPT2   - another point of the line segment. PEND != PSTART !!!
*           PPT2(1) = X-coordinate
*           PPT2(2) = Y-coordinate
*
* Out:
*  Return value >  0, if the point is right-side of the line PPT1->PPT2
*               <  0, if the point if left-side of the line
*               =  0, if the point is on the line or if an error
*                     occurred (PPT1=PPT2)
*
* Remark: More specifically the return value is the squared scalar
*  product of the normal vector of the line with the vector
*  PPT1->point...
************************************************************************

      INTEGER FUNCTION PTROLN (DPX,DPY,PPT1,PPT2)
      
      IMPLICIT NONE
      DOUBLE PRECISION DPX,DPY,PPT1(2),PPT2(2)

      DOUBLE PRECISION DNX,DNY

C To distinguish where the point is, relative to the line, we check
C the vector PPT1->point in relation to the normal vector of the line.
C
C At first build the normal vector by rotating the tangential vector
C by -90 degrees; this points to the right side of the line.

      DNX = (PPT2(2)-PPT1(2))
      DNY = -(PPT2(1)-PPT1(1))

C If this normal vector points into the "same" direction like the
C vector PPT1->point, the point is on the right.
C So we build the scalar product. This is:
C  > 0, if the point is on the right
C  < 0, if the point is on the left
C  = 0, if the point is on the line
C and thus a candidate to serve as a return value of this function!

      PTROLN = DNX*(DPX-PPT1(1)) + DNY*(DPY-PPT1(2))
      
      END
      