**********************************************************************
* This file defines all the Xxxx-routines that are involved with
* the monitor function. The monitor function itself can be found
* (and modified by the user) in the file MON.F.
**********************************************************************

**********************************************************************
* Purpose: Set the general monitor function by dividing
* fmon_i by area_i.
* Handle-based variant
*
* In: 
*  TRIA  - array [1..SZTRIA] of integer
*          Triangulation structure that describes the mesh
*  LSIZE - Handle to rray with mean area of elements around each node
*  NVT   - Number of vertices
*  LMON0 - Handle to array for monitor function, will be filled 
*          with data
*  IREL  - Use relative monoitor function
*  EPS0  - minimum monitor function values; values below that
*          will be truncated.
*  IGEOM  - array [1..*] of integer 
*  DGEOM  - array [1..*] of double 
*           Integer- and double-precision parameter blocks with
*           geometry information. Passed to boundary
*           routines. Not used in this routine.
*  IDATA  - array [1..*] of integer 
*  DDATA  - array [1..*] of double 
*           User defined integer- and double-precision parameter 
*           blocks. Passed to callback routines like monitor
*           function. Not used in this routine.
*
* Out:
*  Information in LMON0-array: new monitor function
**********************************************************************

      SUBROUTINE XMON0(TRIA,LSIZE,LMON0,IREL,EPS0GS,IGEOM,DGEOM,
     *                 IDATA,DDATA)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'stria.inc'
      INTEGER LSIZE,LMON0,IREL,TRIA(SZTRIA)
      DOUBLE PRECISION EPS0GS
      INTEGER IGEOM(*)
      DOUBLE PRECISION DGEOM(*)
      INTEGER IDATA(*)
      DOUBLE PRECISION DDATA(*)
      
      INTEGER LCORVG,LIARR,LDARR
      
      IF (MT.GE.3) 
     *  WRITE (*,FMT='(A)') 'Setting general monitor function'
     
C Get information about the geometry 
      
      LCORVG = TRIA(OLCORVG)
      
C Precalculate information based on the current grid. This will give
C two handles LIARR/LDARR we can pass to the monitor function call.
C All geometry information is capsuled in the fictitious boundary
C components, so we can directly obtain the handles:

      LIARR = 0
      LDARR = 0
      CALL FBDPRC (TRIA,LIARR,LDARR,IGEOM,DGEOM)
      
C Calculate the monitor functions, if possible using the precalculated
C information.
      
      CALL MON0(DWORK(L(LCORVG)),DWORK(L(LSIZE)),TRIA(ONVT),
     *          DWORK(L(LMON0)),EPS0GS,IREL,TRIA,LIARR,LDARR,
     *          IGEOM,DGEOM,IDATA,DDATA)
C      CALL MON0(DWORK(L(LCORVG)),DWORK(L(LSIZE)),TRIA(ONVT),
C     *          DWORK(L(LMON0)),IREL,TRIA,0,0,IGEOM,DGEOM)
     
C Monitor function is completely calculated - release the handles,
C if FBDPRC allocated them (they might be 0 if nothing was
C precalculated)
      
      IF (LDARR.NE.0) CALL ZDISP(0,LDARR,'DARR  ')
      IF (LIARR.NE.0) CALL ZDISP(0,LIARR,'KARR  ')
     
      END 

**********************************************************************
* Calculate the scaled monitor function:
*
* In every vertex i compute FMON_i / AREA_i with AREA_i the
* mean area surrounding the node i.
*
* In:
*   DCORVG - Coordinates of the grid points
*   DSIZE  - array [1..NVT] of double
*            Array with mean area of elements around each node
*   EPS0   - minimum monitor function values; values below that
*            will be truncated.
*   IREL   - FALSE=determine monitor function with AREA_i by the
*                  formula FMON_i / AREA_i, thus prescribing
*                  the ratio between smallest and largest cell size
*                  absolutely
*            TRUE =determine monitor function with AREA_i by the
*                  formula FMON_i without concerning AREA_i, thus
*                  only prescribing the relative motion of the grid
*   TRIA   - array [1..SZTRIA] of integer
*            Triangulation structure that corresponds to DCORVG
*   LIARR  - integer
*            Handle to integer array with precalculated information.
*            Must be set to 0 if not used.
*   LDARR  - integer
*            Handle to double array with precalculated information.
*            Must be set to 0 if not used.
*  IGEOM   - array [1..*] of integer 
*  DGEOM   - array [1..*] of double 
*            Integer- and double-precision parameter blocks with
*            geometry information. Passed to boundary
*            routines. Not used in this routine.
*  IDATA  - array [1..*] of integer 
*  DDATA  - array [1..*] of double 
*           User defined integer- and double-precision parameter 
*           blocks. Passed to callback routines like monitor
*           function. Not used in this routine.
* 
* Out:
*   DMON0  - array [1..NVT] of double
*            Scaled monitor function
**********************************************************************

      SUBROUTINE MON0(DCORVG,DSIZE,NVT,DMON0,EPS0GS,IREL,
     *                TRIA,LIARR,LDARR,IGEOM,DGEOM,IDATA,DDATA)
      IMPLICIT NONE
      INCLUDE 'stria.inc'
      
      DOUBLE PRECISION DCORVG,DSIZE,DMON0,EPS0GS
      INTEGER NVT,IVT,IREL,TRIA(SZTRIA)
      DIMENSION DCORVG(2,*),DMON0(NVT),DSIZE(*)
      INTEGER LIARR,LDARR
      INTEGER IGEOM(*)
      DOUBLE PRECISION DGEOM(*)
      INTEGER IDATA(*)
      DOUBLE PRECISION DDATA(*)

      DOUBLE PRECISION FMON
      EXTERNAL FMON
      
      IF (IREL.EQ.0) THEN
        DO IVT=1,NVT
          DMON0(IVT) =
     *        FMON(DCORVG(1,IVT),DCORVG(2,IVT),IVT,EPS0GS,
     *             TRIA,LIARR,LDARR,IGEOM,DGEOM,IDATA,DDATA) / 
     *             DSIZE(IVT)
        END DO
      ELSE
        DO IVT=1,NVT
          DMON0(IVT) =
     *        FMON(DCORVG(1,IVT),DCORVG(2,IVT),IVT,EPS0GS,
     *             TRIA,LIARR,LDARR,IGEOM,DGEOM,IDATA,DDATA)
        END DO
      END IF      
      END 


**********************************************************************
* Calculate time dependent monitor function f(t)=1-t+t*f0,
* Handle-based variant
**********************************************************************

      SUBROUTINE XMON(LMON0,DT,NVT,LMON)
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INTEGER LMON0,NVT,LMON
      DOUBLE PRECISION DT
      IF (MT.GE.3) 
     *  WRITE (*,FMT='(A,F9.6)') 'Setting monitor function for t = ',DT
      CALL MON(DWORK(L(LMON0)),DT,NVT,DWORK(L(LMON)))
      END SUBROUTINE


**********************************************************************
* Calculates the time-dependent scaled monitor function
* f(t)=1-t+t*f0
**********************************************************************

      SUBROUTINE MON(DMON0,T,NVT,DMON)
      IMPLICIT NONE
      DOUBLE PRECISION DMON0,DMON,T
      INTEGER NVT,IVT
      DIMENSION DMON0(NVT),DMON(NVT)
      DO 100 IVT=1,NVT
        DMON(IVT)=1D0-T+T*DMON0(IVT)
100   CONTINUE
      END SUBROUTINE

