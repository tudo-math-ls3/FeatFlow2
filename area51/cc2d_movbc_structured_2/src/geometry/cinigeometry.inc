************************************************************************
* This file contains the variable representing the user provided
* information in the DAT file of the application.
* More specifically the variables here define the behavior of
* the (moving-boundary-) geometry.
************************************************************************

      INTEGER ICIRTP, IGEOCN, IGEOMS,IFBMNT,ISCMNT,IOVMNT
      COMMON /GEOINT/ ICIRTP,IFBMNT,IGEOCN,ISCMNT,IOVMNT,IGEOMS(2,100)
      
      DOUBLE PRECISION DCXPOS,DCYPOS,DCRPOS,DCRADX,DCRADY,DCRAD,
     *                 DCROT, DCRSIN, DCRCOS,DMNISF,DMNOSF
      COMMON /GEODBL/ DCXPOS,DCYPOS,DCRPOS,DCRADX,DCRADY,DCRAD,
     *                DCROT, DCRSIN, DCRCOS,DMNISF,DMNOSF
      
      SAVE /GEOINT/, /GEODBL/
      
C ICIRTP - Type of obstacle (circle, square, ellipse, ...)
C DCRAD  - Radius of a circle / half length of an edge on a square
C DCRADX - X-Radius of an ellipse
C DCRADY - Y-Radius of an ellipse
C DCROT  - Roration of the obstacle, 0..360°.
C
C For saving computational time the sin/cos of the rotation is
C precalculated:
C DCRSIN - SIN(DCROT); must always be updated if DCROT is changed
C DCRCOS - COS(DCROT); must always be updated if DCROT is changed
C
C IGEOMS - Saves handles for up to 100 (composed) geometry objects.
C          IGEOMS(1,*) is a handle to an integer array of type
C          TGEOCOMIDX, while IGEOMS(2,*) saves the handle to the
C          appropriate TGEOCOM structure (double prec. array) 
C          on the heap.
C IGEOCN - Saves the number of present geometries, i.e.
C          the number of handles allocated in IGEOMS.
C
C DMNISF,
C DMNOSF - Steepness factors for monitor function inside/outside of
C          objects; is multiplied with the monitor function value
C          to correct improper large values through distance
C          calculation when large geometries are used.
