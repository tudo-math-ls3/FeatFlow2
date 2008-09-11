************************************************************************
* This file implements basic composed geometries as a set of
* simple geometries. It provides the basic functionality like
* in the case of simple geometries but tries to augment it in this
* more general case.
*
* Like in the case of simple geometries the routines provided here
* separate into the following types of routines:
*
* 1.) Test if a point is inside of the geometry
* 2.) Write an approximation of the boundary into a GMV-file.
* 3.) Calculate the volume of the geometry
* 4.) Calculate the length of the boundary interface ("arc length")
* 5.) Calculate the outer normal vector of the geometry in a point on
*     the boundary interface
* 6.) Project a point onto the boundary
* 7.) Calculate the minimum distance of an arbitrary point to the
*     object.
*
* Each routine accepts (like described in SCOMPGEOMETRIES.INC)
* an index-integer-array and a double-array describing the composed
* geometry, as well as a wrapper function that calls the appropriate
* simple-geometry-routine according to the object type.
* It's acceptable to specify the appropriate GxxWRP as wrapper routines 
* for a basic functionality of simple-geometry wrappers.
************************************************************************

************************************************************************
* GxxWRP Wrapper-routines
*
* The following wrapper routines are designed to map the TYPE
* specifier in the object structure to a call to the appropriate
* routine defined in GEOMETRIES.F. They share the same syntax
* as the "simple"-routines and are fully interchangeable with them.
************************************************************************

*-----------------------------------------------------------------------
* Test if the point (X,Y) is inside of the geometry xxx.
*-----------------------------------------------------------------------

      INTEGER FUNCTION GINWRP (GEOSIM,CORSYS, X, Y)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(*),X,Y
      DOUBLE PRECISION CORSYS(SZCORS)
      INTEGER TP
      
      INTEGER  GINCIR,GINELL,GINSQR,GINREC,GINSWV,GINPSL
      EXTERNAL GINCIR,GINELL,GINSQR,GINREC,GINSWV,GINPSL
      
      TP = GEOSIM(OTYPE)
      
      GINWRP = 0
      IF (TP.EQ.SGTCIR) GINWRP = GINCIR(GEOSIM,CORSYS,X,Y)
      IF (TP.EQ.SGTELL) GINWRP = GINELL(GEOSIM,CORSYS,X,Y)
      IF (TP.EQ.SGTSQR) GINWRP = GINSQR(GEOSIM,CORSYS,X,Y)
      IF (TP.EQ.SGTREC) GINWRP = GINREC(GEOSIM,CORSYS,X,Y)
      IF (TP.EQ.SGTSWV) GINWRP = GINSWV(GEOSIM,CORSYS,X,Y)
      IF (TP.EQ.SGTPSL) GINWRP = GINPSL(GEOSIM,CORSYS,X,Y)
      
      END 

*-----------------------------------------------------------------------
* Calculate the volume of a geometry object.
*-----------------------------------------------------------------------
      
      DOUBLE PRECISION FUNCTION GVLWRP (GEOSIM,CORSYS)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(*)
      DOUBLE PRECISION CORSYS(SZCORS)
      INTEGER TP
      
      DOUBLE PRECISION GVLCIR,GVLELL,GVLSQR,GVLREC,GVLSWV,GVLPSL
      EXTERNAL GVLCIR,GVLELL,GVLSQR,GVLREC,GVLSWV,GVLPSL
      
      TP = GEOSIM(OTYPE)
      
      GVLWRP = -1
      IF (TP.EQ.SGTCIR) GVLWRP = GVLCIR(GEOSIM,CORSYS)
      IF (TP.EQ.SGTELL) GVLWRP = GVLELL(GEOSIM,CORSYS)
      IF (TP.EQ.SGTSQR) GVLWRP = GVLSQR(GEOSIM,CORSYS)
      IF (TP.EQ.SGTREC) GVLWRP = GVLREC(GEOSIM,CORSYS)
      IF (TP.EQ.SGTSWV) GVLWRP = GVLSWV(GEOSIM,CORSYS)
      IF (TP.EQ.SGTPSL) GVLWRP = GVLPSL(GEOSIM,CORSYS)
      
      END

*-----------------------------------------------------------------------
* Calculate the length of the boundary interface
*-----------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION GALWRP (GEOSIM,CORSYS)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(*)
      DOUBLE PRECISION CORSYS(SZCORS)
      INTEGER TP
      
      DOUBLE PRECISION GALCIR,GALELL,GALSQR,GALREC,GALSWV,GALPSL
      EXTERNAL GALCIR,GALELL,GALSQR,GALREC,GALSWV,GALPSL
      
      TP = GEOSIM(OTYPE)
      
      GALWRP = -1
      IF (TP.EQ.SGTCIR) GALWRP = GALCIR(GEOSIM,CORSYS)
      IF (TP.EQ.SGTELL) GALWRP = GALELL(GEOSIM,CORSYS)
      IF (TP.EQ.SGTSQR) GALWRP = GALSQR(GEOSIM,CORSYS)
      IF (TP.EQ.SGTREC) GALWRP = GALREC(GEOSIM,CORSYS)
      IF (TP.EQ.SGTSWV) GALWRP = GALSWV(GEOSIM,CORSYS)
      IF (TP.EQ.SGTPSL) GALWRP = GALPSL(GEOSIM,CORSYS)
      
      END

*-----------------------------------------------------------------------
* Write an approximation of the boundary interface into a GMV-file
* identified by the file handle MUNIT.
* Use material number MMAT to identify the coloring in the GMV-file.
*-----------------------------------------------------------------------

      SUBROUTINE GGMWRP(GEOSIM,CORSYS, MUNIT, MMAT)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(*)
      DOUBLE PRECISION CORSYS(SZCORS)
      INTEGER MUNIT,MMAT
      INTEGER TP
      
      TP = GEOSIM(OTYPE)
      
      IF (TP.EQ.SGTCIR) CALL GGMCIR(GEOSIM,CORSYS, MUNIT, MMAT)
      IF (TP.EQ.SGTELL) CALL GGMELL(GEOSIM,CORSYS, MUNIT, MMAT)
      IF (TP.EQ.SGTSQR) CALL GGMSQR(GEOSIM,CORSYS, MUNIT, MMAT)
      IF (TP.EQ.SGTREC) CALL GGMREC(GEOSIM,CORSYS, MUNIT, MMAT)
      IF (TP.EQ.SGTSWV) CALL GGMSWV(GEOSIM,CORSYS, MUNIT, MMAT)
      IF (TP.EQ.SGTPSL) CALL GGMPSL(GEOSIM,CORSYS, MUNIT, MMAT)

      END

*-----------------------------------------------------------------------
* Calculate the normalized outer normal vector in the boundary-point
* (X,Y) on the interface of the geometry object GEOSIM.
*-----------------------------------------------------------------------
      
      SUBROUTINE GNMWRP (GEOSIM,CORSYS, X,Y, XN,YN)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(*),X,Y,XN,YN
      DOUBLE PRECISION CORSYS(SZCORS)
      INTEGER TP
      
      TP = GEOSIM(OTYPE)
      
      IF (TP.EQ.SGTCIR) CALL GNMCIR(GEOSIM,CORSYS, X,Y, XN,YN)
      IF (TP.EQ.SGTELL) CALL GNMELL(GEOSIM,CORSYS, X,Y, XN,YN)
      IF (TP.EQ.SGTSQR) CALL GNMSQR(GEOSIM,CORSYS, X,Y, XN,YN)
      IF (TP.EQ.SGTREC) CALL GNMREC(GEOSIM,CORSYS, X,Y, XN,YN)
      IF (TP.EQ.SGTSWV) CALL GNMSWV(GEOSIM,CORSYS, X,Y, XN,YN)
      IF (TP.EQ.SGTPSL) CALL GNMPSL(GEOSIM,CORSYS, X,Y, XN,YN)
      
      END
      
      
*-----------------------------------------------------------------------
* Calculate the distance of the point (X,Y) to the boundary interface
* of the object.
*-----------------------------------------------------------------------

      INTEGER FUNCTION GDSWRP (GEOSIM,CORSYS, X,Y, DIST)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(*),X,Y,DIST
      DOUBLE PRECISION CORSYS(SZCORS)
      INTEGER TP
      
      INTEGER  GDSCIR,GDSELL,GDSSQR,GDSREC,GDSSWV,GDSPSL
      EXTERNAL GDSCIR,GDSELL,GDSSQR,GDSREC,GDSSWV,GDSPSL
      
      TP = GEOSIM(OTYPE)
      
      GDSWRP = 0
      IF (TP.EQ.SGTCIR) GDSWRP = GDSCIR(GEOSIM,CORSYS, X,Y, DIST)
      IF (TP.EQ.SGTELL) GDSWRP = GDSELL(GEOSIM,CORSYS, X,Y, DIST)
      IF (TP.EQ.SGTSQR) GDSWRP = GDSSQR(GEOSIM,CORSYS, X,Y, DIST)
      IF (TP.EQ.SGTREC) GDSWRP = GDSREC(GEOSIM,CORSYS, X,Y, DIST)
      IF (TP.EQ.SGTSWV) GDSWRP = GDSSWV(GEOSIM,CORSYS, X,Y, DIST)
      IF (TP.EQ.SGTPSL) GDSWRP = GDSPSL(GEOSIM,CORSYS, X,Y, DIST)
      
      END 

*-----------------------------------------------------------------------
* Project a point onto the boundary interface of the geometry object.
*-----------------------------------------------------------------------

      INTEGER FUNCTION GPRWRP (GEOSIM,CORSYS, X,Y, XP,YP)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(*),X,Y,XP,YP
      DOUBLE PRECISION CORSYS(SZCORS)
      INTEGER TP
      
      INTEGER  GPRCIR,GPRELL,GPRSQR,GPRREC,GPRSWV,GPRPSL
      EXTERNAL GPRCIR,GPRELL,GPRSQR,GPRREC,GPRSWV,GPRPSL
      
      TP = GEOSIM(OTYPE)
      
      GPRWRP = 0
      IF (TP.EQ.SGTCIR) GPRWRP = GPRCIR(GEOSIM,CORSYS, X,Y, XP,YP)
      IF (TP.EQ.SGTELL) GPRWRP = GPRELL(GEOSIM,CORSYS, X,Y, XP,YP)
      IF (TP.EQ.SGTSQR) GPRWRP = GPRSQR(GEOSIM,CORSYS, X,Y, XP,YP)
      IF (TP.EQ.SGTREC) GPRWRP = GPRREC(GEOSIM,CORSYS, X,Y, XP,YP)
      IF (TP.EQ.SGTSWV) GPRWRP = GPRSWV(GEOSIM,CORSYS, X,Y, XP,YP)
      IF (TP.EQ.SGTPSL) GPRWRP = GPRPSL(GEOSIM,CORSYS, X,Y, XP,YP)
      
      END 

************************************************************************
* The following set of routines now define the functionality
* of the composed-geometry objects. The same set as in the case of
* simple geometries is provided, but the calling conventions are
* different.
************************************************************************

************************************************************************
* INTEGER FUNCTION GINCMP (GEOIDX, GEOCOM, GINFCT, X, Y)
*
* Test if the point (X,Y) is inside of the geometry xxx.
*
* In:
*  X,Y    - coordinate of the point
*  GEOIDX - the GEOCOMIDX integer array describing the index positions
*           of the objects in GEOCOM
*  GEOCOM - the GEOCOM double array describing the sub-objects
*  GINFCT - a wrapper function that provides calling the correct
*           GINxxx-routine. GINWRP is acceptable.
*
* Out:
*  Result > 0, if the point is inside of the geometry object.
*              If the group-object is not completely inverted (INVER=1),
*              the number describes the number of sub-objects
*              that contain this point.
*         = 0, if the point is not inside of the geometry object
*         < 0, if the point is inside of the geometry object and
*              the object is inverted
************************************************************************
      
      INTEGER FUNCTION GINCMP (GEOIDX, GEOCOM, GINFCT, X, Y)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      INCLUDE 'scompgeometries.inc'
      DOUBLE PRECISION GEOCOM(*),X,Y
      INTEGER GEOIDX(*),GINFCT
      EXTERNAL GINFCT
      
      INTEGER I,CNT

      CNT = 0
      
C Loop through all sub-objects

      DO I=1,GEOIDX(OGCCNT)
      
C Add the results of the GINxxx-routines. Negative values
C decrease the number until it's <= 0 - that indicates
C that there are so many objects substracted from the first
C object, that the point is considered to be "outside".

        CNT = CNT + GINFCT(GEOCOM(GEOIDX(OGCIDX+I-1)), 
     *      GEOCOM(OGCMAIN+OCSYS-1), X, Y)
     
      END DO
      
      GINCMP = CNT

C Switch the result if we are in "inverted" mode
      
      IF (GEOCOM(OINVER).NE.0) THEN
        GINCMP = -GINCMP
      END IF
      
      END
      
************************************************************************
* DOUBLE PRECISION FUNCTION GVLCMP (GEOIDX, GEOCOM, GVLFCT)
* 
* Calculate sum of all volume of a geometry object.
*
* In:
*  GEOIDX - the GEOCOMIDX integer array describing the index positions
*           of the objects in GEOCOM
*  GEOCOM - the GEOCOM double array describing the sub-objects
*  GVLFCT - a wrapper function that provides calling the correct
*           GVLxxx-routine. GVLWRP is acceptable.
*
* Out:
*  Result = sum of all volumes of the sub-objects
*         = -1, if the volume cannot be calculated
************************************************************************

      DOUBLE PRECISION FUNCTION GVLCMP (GEOIDX, GEOCOM, GVLFCT)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      INCLUDE 'scompgeometries.inc'
      DOUBLE PRECISION GEOCOM(*)
      INTEGER GEOIDX(*)
      DOUBLE PRECISION GVLFCT
      EXTERNAL GVLFCT
      
      DOUBLE PRECISION VOL
      INTEGER I

      VOL = 0D0
      
C Loop through all sub-objects

      DO I=1,GEOIDX(OGCCNT)
      
C sum up the volumes

        VOL = VOL + GVLFCT(GEOCOM(GEOIDX(OGCIDX+I-1)),
     *                     GEOCOM(OGCMAIN+OCSYS-1) )
      
      END DO
      
      GVLCMP = VOL
      
      END

************************************************************************
* DOUBLE PRECISION FUNCTION GALCMP (GEOIDX, GEOCOM, GALFCT)
* 
* Calculate the length of the boundary interface
*
* In:
*  GEOIDX - the GEOCOMIDX integer array describing the index positions
*           of the objects in GEOCOM
*  GEOCOM - the GEOCOM double array describing the sub-objects
*  GALFCT - a wrapper function that provides calling the correct
*           GALxxx-routine. GALWRP is acceptable.
*
* Out:
*  Result = sum of all lenths of the interfaces of the sub-objects
*         = -1, if the length cannot be calculated
************************************************************************

      DOUBLE PRECISION FUNCTION GALCMP (GEOIDX, GEOCOM, GALFCT)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      INCLUDE 'scompgeometries.inc'
      DOUBLE PRECISION GEOCOM(*)
      INTEGER GEOIDX(*)
      DOUBLE PRECISION GALFCT
      EXTERNAL GALFCT
      
      DOUBLE PRECISION ALE
      INTEGER I

      ALE = 0D0
      
C Loop through all sub-objects

      DO I=1,GEOIDX(OGCCNT)
      
C sum up the arc lengths

        ALE = ALE + GALFCT(GEOCOM(GEOIDX(OGCIDX+I-1)),
     *                     GEOCOM(OGCMAIN+OCSYS-1) )
      
      END DO
      
      GALCMP = ALE

      END

************************************************************************
* SUBROUTINE GGMCMP (GEOIDX, GEOCOM, GGMFCT, GDSFCT, MUNIT, MMAT)
*
* Write an approximation of the boundary interface into a GMV-file
* identified by the file handle MUNIT.
* Use material number MMAT to identify the coloring in the GMV-file.
*
* The caller of this routine is responsible for the administration
* of the sections in the GMV file. This routine only has to write out
* the content for the POLYGON's section to approximate the boundary
* components. The header/footer of this section is written to the
* file by the caller.
*
* In:
*  GEOIDX - the GEOCOMIDX integer array describing the index positions
*           of the objects in GEOCOM
*  GEOCOM - the GEOCOM double array describing the sub-objects
*  GGMFCT - a wrapper function that provides calling the correct
*           GGMxxx-routine. GGMWRP is acceptable.
*  MUNIT  - File unit where to write the GMV output to.
*  MMAT   - The first material number to use.
************************************************************************

      SUBROUTINE GGMCMP (GEOIDX, GEOCOM, GGMFCT, MUNIT, MMAT)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      INCLUDE 'scompgeometries.inc'
      DOUBLE PRECISION GEOCOM(*)
      INTEGER GEOIDX(*),MUNIT, MMAT
      EXTERNAL GGMFCT
      
      INTEGER I
      
C Loop through all sub-objects

      DO I=1,GEOIDX(OGCCNT)
      
C tell the sub-object to write itself to disc

        CALL GGMFCT(GEOCOM(GEOIDX(OGCIDX+I-1)),
     *              GEOCOM(OGCMAIN+OCSYS-1),
     *              MUNIT, MMAT)
      
      END DO
      
      END

************************************************************************
* SUBROUTINE GNMCMP (GEOIDX, GEOCOM, GNMFCT, X,Y, XN,YN)
* 
* Calculate the normalized outer normal vector in the boundary-point
* (X,Y) on the interface of the geometry object GEOSIM.
* 
* In:
*  GEOIDX - the GEOCOMIDX integer array describing the index positions
*           of the objects in GEOCOM
*  GEOCOM - the GEOCOM double array describing the sub-objects
*  GNMFCT - a wrapper function that provides calling the correct
*           GGMxxx-routine. GNMWRP is acceptable.
*  GDSFCT - a wrapper function that provides calling the correct
*           GDSxxx-routine for distance calculation. GDSWRP is
*           acceptable.
*  X,Y    - coordinates of a point on the interface where to calculate
*           the outer normal vector
* 
* Out:
*  XN,YN  - Outer normal vector in the point (X,Y), 
*           normalized to length = 1
*           =(0,0) if the calculation is not possible.
*
* This routine only works completely correct for points outside of
* the geometry. Points inside of the geometry may be projected
* wrong if they are lying in more than one sub-component!
************************************************************************

      SUBROUTINE GNMCMP (GEOIDX, GEOCOM, GNMFCT, GDSFCT, X,Y, XN,YN)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      INCLUDE 'scompgeometries.inc'
      DOUBLE PRECISION GEOCOM(*),X,Y,XN,YN
      INTEGER GEOIDX(*),GDSFCT
      EXTERNAL GNMFCT,GDSCMP,GDSFCT
      INTEGER GDSCMP
      
      DOUBLE PRECISION DIST
      INTEGER I

      XN = 0D0
      YN = 0D0

C Calculate the sub-object with the minimum distance

      I = GDSCMP (GEOIDX, GEOCOM, GDSFCT, X,Y, DIST)
      
C Calculate the normal vector by that calculation routine
      
      IF (I.NE.0) THEN
        CALL GNMFCT(GEOCOM(GEOIDX(OGCIDX+I-1)),
     *              GEOCOM(OGCMAIN+OCSYS-1),
     *              X,Y, XN,YN)
      END IF

C Switch the result if we are in "inverted" mode
      
      IF (GEOCOM(OINVER).NE.0) THEN
        XN = -XN
        YN = -YN
      END IF
      
      END

************************************************************************
* INTEGER FUNCTION GDSCMP (GEOIDX, GEOCOM, GDSFCT, X,Y, DIST)
*
* Calculate the distance of the point (X,Y) to the boundary interface
* of the object.
*
* In:
*  GEOIDX - the GEOCOMIDX integer array describing the index positions
*           of the objects in GEOCOM
*  GEOCOM - the GEOCOM double array describing the sub-objects
*  GDSFCT - a wrapper function that provides calling the correct
*           GGMxxx-routine. GDSWRP is acceptable.
*  X,Y    - the point whose distance should be calculated
*
* Out:
*  DIST   - min. distance of (X,Y) to the boundary interfaces; the sign
*           identifies the position of the point:
*           > 0: the point is outside of the geometry object
*           <=0: the point is inside of the geometry object
*
*  Result: > 0, if the distance was calculated successfully.
*               The number describes the number of the sub-component
*               that served for calculating the distance.
*          <=0, if the computation was not possible
************************************************************************

      INTEGER FUNCTION GDSCMP (GEOIDX, GEOCOM, GDSFCT, X,Y, DIST)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      INCLUDE 'scompgeometries.inc'
      DOUBLE PRECISION GEOCOM(*),X,Y,DIST
      INTEGER GEOIDX(*)
      INTEGER GDSFCT
      EXTERNAL GDSFCT
      
      DOUBLE PRECISION DST
      INTEGER I,J

      GDSCMP = 0

C Loop through all sub-objects

      DO I=1,GEOIDX(OGCCNT)
      
C Calculate the distance to all objects that allow calculating the
C distance.

        J = GDSFCT(GEOCOM(GEOIDX(OGCIDX+I-1)),
     *             GEOCOM(OGCMAIN+OCSYS-1),
     *             X,Y,DST)
        
        IF (J.GT.0) THEN

C Build the minimum of the distances.
C Note: For points outside of the geometry this is the distance of the 
C  point to the geometry. For points inside the geometry, the minimum
C  is build from the *negative* distances: This results in the maximum
C  distance to all objects containing the point. This is of course not
C  the real distance to the boundary interface of the composed geometry,
C  but somehow an approximation to it.

          IF ((I.EQ.1).OR.(DST.LT.DIST)) THEN
            GDSCMP = I
            DIST = DST
          END IF

        END IF
      
      END DO
      
C Switch the result if we are in "inverted" mode
      
      IF (GEOCOM(OINVER).NE.0) THEN
        DIST = -DIST
      END IF
      
      END
      
************************************************************************
* INTEGER FUNCTION GPRCMP (GEOIDX, GEOCOM, GPRFCT, X,Y, XP,YP)
*
* Project a point onto the boundary interface of the geometry object.
*
* In:
*  GEOIDX - the GEOCOMIDX integer array describing the index positions
*           of the objects in GEOCOM
*  GEOCOM - the GEOCOM double array describing the sub-objects
*  GPRFCT - a wrapper function that provides calling the correct
*           GPRxxx-routine. GPRWRP is acceptable.
*  GDSFCT - a wrapper function that provides calling the correct
*           GDSxxx-routine for distance calculation. GDSWRP is
*           acceptable.
*  X,Y    - the point that should be projected
*
* Out:
*  XP,YP  - the projected point on the interface of GEOSIM
*
*  Result: > 0, if the projection was calculated successfully.
*               The number represents the number of the sub-component
*               onto which the projection was done.
*          <=0, if the computation was not possible
*
* This routine only works completely correct for points outside of
* the geometry. Points inside of the geometry may be projected
* wrong if they are lying in more than one sub-component!
************************************************************************

      INTEGER FUNCTION GPRCMP (GEOIDX, GEOCOM, GPRFCT, GDSFCT, 
     *                         X,Y, XP,YP)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      INCLUDE 'scompgeometries.inc'
      DOUBLE PRECISION GEOCOM(*),X,Y,XP,YP
      INTEGER GEOIDX(*)
      INTEGER GDSCMP,GPRFCT,GDSFCT
      EXTERNAL GPRFCT,GDSCMP,GDSFCT

      DOUBLE PRECISION DIST
      INTEGER I

C Calculate the sub-object with the minimum distance

      I = GDSCMP (GEOIDX, GEOCOM, GDSFCT, X,Y, DIST)
      
C Project onto that boundary. 
      
      IF (I.NE.0) THEN
        I = GPRFCT(GEOCOM(GEOIDX(OGCIDX+I-1)), 
     *             GEOCOM(OGCMAIN+OCSYS-1),
     *             X,Y, XP,YP)        
      END IF

      GPRCMP = I

      END
