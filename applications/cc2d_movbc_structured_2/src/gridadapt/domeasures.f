************************************************************************
* This file contains auxiliary routines used in the grid adaption for
* measuring the domain measure, calculationg average size of
* elements,...
************************************************************************

************************************************************************
* Determine the area of the whole domain
*
* In:
*  NEL    - Number of vertices in domain
*  DONE   - array [1..NEL] of double
*
* Out:
*  DMEASRE - Measure of the domain
************************************************************************

      SUBROUTINE MEASURE(NEL,DONE,DMEASRE)
      IMPLICIT NONE
      INTEGER NEL
      DOUBLE PRECISION DONE
      DIMENSION DONE(NEL)
      INTEGER IEL
      DOUBLE PRECISION DMEASRE
      
      DMEASRE=0D0
      DO IEL=1,NEL
        DMEASRE=DMEASRE+DONE(IEL)
      END DO
      
      END 


************************************************************************
* Purpose: Calculate the average size of elements surrounding a node
* and normalize it such that 1d0 is the max. over the whole vector dsize
*
* In:
*  DCORVG, 
*  KVERT,
*  KADJ,
*  NEL,
*  NVT     - usual information about the current grid
*
* Out:
*  KCOUNT  - array [1..NVT] of integer
*            Number of elements that meet of a vertex 1..NVT
*  DSIZEM  - array [1..NEL] of double
*            Area of each element
*  DSIZE   - array [1..NVT] of double
*            Average size of elements surrounding a node, normalised
*  DRATIO  - array [1..NEL] of double
*            Maximum ratio of the area of each element and its
*            surrounding elements
************************************************************************

      SUBROUTINE MEASR2(DCORVG,KVERT,KADJ,NEL,NVT,KCOUNT,DSIZE,DSIZEM,
     *                  DRATIO)
     
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      
      DOUBLE PRECISION DCORVG,DSIZE,DSIZEM,DRATIO
      INTEGER KCOUNT,KVERT,KADJ
      DIMENSION DCORVG(2,*),DSIZE(*),KCOUNT(*),KVERT(4,*)
      DIMENSION DSIZEM(*),DRATIO(*),KADJ(4,*)
      INTEGER IEL,IVT,NEL,NVT,IVE,JEL
      DOUBLE PRECISION DAX,DAY,DNX,DNY,DNRMN,DIST,DIST2,DPX,DPY,DMAXS
      DOUBLE PRECISION DMINS,DQUOT

C Initialise the arrays with 0

      CALL LCL1 (DSIZE,NVT)
      CALL LCL3 (KCOUNT,NVT)

      DO IEL=1,NEL
!     DAREA=0.5D0*((DCORVG(1,KVERT(1,IEL))-DCORVG(1,KVERT(2,IEL)))
!     *            *(DCORVG(2,KVERT(1,IEL))+DCORVG(2,KVERT(2,IEL)))
!     *            +(DCORVG(1,KVERT(2,IEL))-DCORVG(1,KVERT(3,IEL)))
!     *            *(DCORVG(2,KVERT(2,IEL))+DCORVG(2,KVERT(3,IEL)))
!     *            +(DCORVG(1,KVERT(3,IEL))-DCORVG(1,KVERT(4,IEL)))
!     *            *(DCORVG(2,KVERT(3,IEL))+DCORVG(2,KVERT(4,IEL)))
!     *            +(DCORVG(1,KVERT(4,IEL))-DCORVG(1,KVERT(1,IEL)))
!     *            *(DCORVG(2,KVERT(4,IEL))+DCORVG(2,KVERT(1,IEL)))

C Get first point in quad-element

        DAX=DCORVG(1,KVERT(1,IEL))
        DAY=DCORVG(2,KVERT(1,IEL))
        
C Get 3rd point, build normal vector (DNX,DNY) to the diagonal of the
C element, normalise it
        
        DNX= DCORVG(2,KVERT(3,IEL))-DAY
        DNY=-DCORVG(1,KVERT(3,IEL))+DAX
        DNRMN=DSQRT(DNX*DNX+DNY*DNY)
        DNX=DNX/DNRMN
        DNY=DNY/DNRMN
        
! triangle 1,2,3

        DPX=DCORVG(1,KVERT(2,IEL))
        DPY=DCORVG(2,KVERT(2,IEL))
        DIST=ABS(DNX*DPX+DNY*DPY-DNX*DAX-DNY*DAY)

! triangle 1,3,4

        DPX=DCORVG(1,KVERT(4,IEL))
        DPY=DCORVG(2,KVERT(4,IEL))
        DIST2=ABS(DNX*DPX+DNY*DPY-DNX*DAX-DNY*DAY)
        
C Calculate size of the element
        
        DSIZEM(IEL)=0.5D0*DNRMN*(DIST+DIST2)
        
C Add the area information to each vertex, building an array with
C the sum of areas meeting in a vertex
        
        DSIZE(KVERT(1,IEL))=DSIZE(KVERT(1,IEL))+0.5D0*DNRMN*(DIST+DIST2)
        DSIZE(KVERT(2,IEL))=DSIZE(KVERT(2,IEL))+0.5D0*DNRMN*(DIST+DIST2)
        DSIZE(KVERT(3,IEL))=DSIZE(KVERT(3,IEL))+0.5D0*DNRMN*(DIST+DIST2)
        DSIZE(KVERT(4,IEL))=DSIZE(KVERT(4,IEL))+0.5D0*DNRMN*(DIST+DIST2)
        KCOUNT(KVERT(1,IEL))=KCOUNT(KVERT(1,IEL))+1
        KCOUNT(KVERT(2,IEL))=KCOUNT(KVERT(2,IEL))+1
        KCOUNT(KVERT(3,IEL))=KCOUNT(KVERT(3,IEL))+1
        KCOUNT(KVERT(4,IEL))=KCOUNT(KVERT(4,IEL))+1
        
      END DO

C Build the average size of the elements around each vertex:
C
C Divide the summed area by the number of elements meeting in a vertex
C and normalise it with the maximum area.

      DMAXS=0D0
      DMINS=1D0
      DO IVT=1,NVT
        DSIZE(IVT)=DSIZE(IVT)/DBLE(KCOUNT(IVT))
        DMAXS=MAX(DMAXS,DSIZE(IVT))
        DMINS=MIN(DMINS,DSIZE(IVT))
      END DO

      DO IVT=1,NVT
        DSIZE(IVT)=DSIZE(IVT)/DMAXS
      END DO

C      IF (MT.GE.2) THEN
C        WRITE (MTERM,FMT='(A,2F12.8)') '  Smallest/largest element:',
C     *                                 DMINS,DMAXS
C        WRITE (MTERM,FMT='(A,F8.6)')
C     *        '  Ratio smallest/largest element: ',DMINS/DMAXS
C      END IF

C Build the vector with the ratio between the area of each element
C and its surrounding elements

      DO IEL=1,NEL
        DRATIO(IEL)=0D0
        DO IVE=1,4
          JEL=KADJ(IVE,IEL)
          IF (JEL.NE.0) THEN
C Ratio between area of current and its neighbour element
            DQUOT=DSIZEM(IEL)/DSIZEM(JEL)
            IF (DQUOT.LT.1D0) DQUOT=1D0/DQUOT
C Find the maximum ratio
            IF (DQUOT.GT.DRATIO(IEL)) DRATIO(IEL)=DQUOT
          END IF
        END DO
      END DO

      END 

