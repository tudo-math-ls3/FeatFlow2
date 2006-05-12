************************************************************************
* This file contains a collection of simple, analytical geometries.
* These can be used e.g. in fictitious boundary routines to test
* whether a point is in a geometry or not. Furthermore this allowes
* to construct complex geometries out of simple ones.
*
* Each geometry object is specified by a geometry structure that defines
* its parameters. The structure itself is realised as a simple double-
* array that is passed to the corresponding geometry routines.
*
* The geometry structure is explained in the include file
* SGEOMETRIES.INC. The geometries themself consist of a set of different
* functions that provide the actual handling. In particular the
* following set of functions is provided by each geometry object:
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
* The function-type 1.) is absolutely mandatory for each object.
* The function-types 2.)-7.) should normally be possible for simple-
* type objects but maybe impossible for more complex objects.
* In this case the appropriate routines return values that indicate
* if they are abstract, i.e. if they have no function.
************************************************************************

************************************************************************
* In the following we describe the general interface to all types
* of geometry routines:
************************************************************************
* INTEGER FUNCTION GCRxxx (GEOSIM, X, Y, ROT, SCALE, INVER, ...)
*
* Initialise geometry structure
*
* This routine can be used for easy and direct initialisation
* of the GEOSIM-structure for the appropriate geometry xxx.
* The routine transfers the parameters into the structure GEOSIM.
* The return value is the size of the appropriate geometry structure.
*
* In:
*  GEOSIM - the geometry object that corresponds to geometry type xxx;
*           to be filled with data
*  X,Y    - position of the object 
*  ROT    - rotation of the object - in degrees (0..360)
*  SCALE  - scaling factor; 1=no scaling=standard
*  INVER  - whether the object is "inverted"
*  ...    - additional parameters that are transferred into the
*           structure; geometry dependent
*
* Out:
*  GEOSIM - the geometry object structure, filled with data
*
************************************************************************
* INTEGER FUNCTION GINxxx (GEOSIM,CORSYS, X, Y)
*
* Test if the point (X,Y) is inside of the geometry xxx.
*
* In:
*  X,Y    - coordinate of the point
*  GEOSIM - the geometry object that corresponds to geometry type xxx
*  CORSYS - the coordinate system of the geometry
*
* Out:
*  Result = 1, if the point is inside of the geometry object
*         = 0, if the point is not inside of the geometry object
*         = -1, if the point is inside of the geometry object
*               and the object is "inverted"
************************************************************************
* DOUBLE PRECISION FUNCTION GVLxxx (GEOSIM,CORSYS)
* 
* Calculate the volume of a geometry object.
*
* In:
*  GEOSIM - the geometry object that corresponds to geometry type xxx
*  CORSYS - the coordinate system of the geometry
*
* Out:
*  Result = volume of the geometry object
*         = -volume, if the object is inverted
*         = 0, if the object has no volume or
*              if volume cannot be calculated
*  The Result is negative in the case that the object is
*  inverted. In that case the volume has to be added/subtracted
*  to the volume of the whole geometry to get the actual volume.
************************************************************************
* DOUBLE PRECISION FUNCTION GALxxx (GEOSIM,CORSYS)
* 
* Calculate the length of the boundary interface
*
* In:
*  GEOSIM - the geometry object that corresponds to geometry type xxx
*  CORSYS - the coordinate system of the geometry
*
* Out:
*  Result = lenth of the interface of the boundary object
*         = -1, if the length cannot be calculated
************************************************************************
* SUBROUTINE GGMxxx (GEOSIM,CORSYS, MUNIT, MMAT)
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
*  GEOSIM - the geometry object that corresponds to geometry type xxx
*  CORSYS - the coordinate system of the geometry
*  MUNIT  - File unit where to write the GMV output to.
*  MMAT   - The first material number to use.
************************************************************************
* SUBROUTINE GNMxxx (GEOSIM,CORSYS, X,Y, XN,YN)
* 
* Calculate the normalized outer normal vector in the boundary-point
* (X,Y) on the interface of the geometry object GEOSIM.
* 
* In:
*  GEOSIM - the geometry object that corresponds to geometry type xxx
*  X,Y    - coordinates of a point on the interface where to calculate
*           the outer normal vector
* 
* Out:
*  XN,YN  - Outer normal vector in the point (X,Y), 
*           normalized to length = 1
*           =(0,0) if the calculation is not possible.
************************************************************************
* INTEGER FUNCTION GDSxxx (GEOSIM,CORSYS, X,Y, DIST)
*
* Calculate the distance of the point (X,Y) to the boundary interface
* of the object.
*
* In:
*  GEOSIM - the geometry object that corresponds to geometry type xxx
*  CORSYS - the coordinate system of the geometry
*  X,Y    - the point whose distance should be calculated
*
* Out:
*  DIST   - distance of (X,Y) to the boundary interface; the sign
*           identifies the position of the point:
*           > 0: the point is outside of the geometry object
*           <=0: the point is inside of the geometry object
*
*  Result: > 0, if the distance was calculated successfully
*          <=0, if the computation was not possible
************************************************************************
* INTEGER FUNCTION GPRxxx (GEOSIM,CORSYS, X,Y, XP,YP)
*
* Project a point onto the boundary interface of the geometry object.
*
* In:
*  GEOSIM - the geometry object that corresponds to geometry type xxx
*  CORSYS - the coordinate system of the geometry
*  X,Y    - the point that should be projected
*
* Out:
*  XP,YP  - the projected point on the interface of GEOSIM
*
*  Result: > 0, if the distance was calculated successfully
*          <=0, if the computation was not possible
************************************************************************

************************************************************************
* Auxiliary routine: Rotate and transform
*
* This routine performs a coordinate transformation from a
* "reference" geometry into the "real" grometry.
* It accepts a point (X,Y) which is assumed to be around the
* origin (0,0). This point is rotated, scaled and shifted
* according to the rotation angle / reference point of the object. 
*
* In:
*  GEOSIM - the geometry object that corresponds to geometry type xxx
*  CORSYS - the coordinate system of the geometry
*  X,Y    - the point that should be transformed
*
* Out:
*  XT,YT  - the transformed point
************************************************************************

      SUBROUTINE GRTTRF (GEOSIM,CORSYS, X,Y, XP,YP)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGSIM),X,Y,XP,YP,XP1,YP1
      DOUBLE PRECISION CORSYS(SZCORS)
      DOUBLE PRECISION C,S
      
C Rotate around the origin:
      
      IF (GEOSIM(OROT).NE.0D0) THEN
C        C = COS(GEOSIM(OROT)*PI/180D0)
C        S = SIN(GEOSIM(OROT)*PI/180D0)
        C = GEOSIM(OROTCOS)
        S = GEOSIM(OROTSIN)
        XP = C*X-S*Y
        YP = S*X+C*Y
      ELSE
        XP = X
        YP = Y
      END IF

C Scale and shift to reference point

      XP1 = GEOSIM(OSCALE)*XP+GEOSIM(OREFX)
      YP1 = GEOSIM(OSCALE)*YP+GEOSIM(OREFY)
      
C Rotate again around the origin for the rotation of the
C coordinate system
      
      IF (CORSYS(ORIGRT).NE.0D0) THEN
C        C = COS(CORSYS(ORIGRT)*PI/180D0)
C        S = SIN(CORSYS(ORIGRT)*PI/180D0)
        C = CORSYS(ORICOS)
        S = CORSYS(ORISIN)
        XP = C*XP1-S*YP1
        YP = S*XP1+C*YP1
      ELSE
        XP = XP1
        YP = YP1
      END IF

C Shift to the coordinate system origin

      XP = CORSYS(ORISCL)*XP+CORSYS(ORIGX)
      YP = CORSYS(ORISCL)*YP+CORSYS(ORIGY)
      
      END 

************************************************************************
* Auxiliary routine: Rotate and transform back
*
* This routine performs a coordinate transformation onto a
* "reference" geometry. It accepts a point (X,Y), shifts it to
* the origin (0,0) according to the reference point of the object
* and performs a back-rotation according to the
* (negative) angle in the geometry object. The caller can then work
* and test with the non-rotated coordinate system at the origin.
*
* In:
*  GEOSIM - the geometry object that corresponds to geometry type xxx
*  CORSYS - the coordinate system of the geometry
*  X,Y    - the point that should be transformed
*
* Out:
*  XT,YT  - the transformed point
************************************************************************

      SUBROUTINE GRTBCK (GEOSIM,CORSYS, X,Y, XP,YP)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGSIM),X,Y,XP,YP
      DOUBLE PRECISION CORSYS(SZCORS)
      DOUBLE PRECISION XP1,YP1,C,S,SC1,SC2

      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)

C get the scaling factors

      IF (CORSYS(ORISCL).NE.0D0) THEN
        SC1 = 1D0/CORSYS(ORISCL)
      ELSE 
        SC1 = 0D0
      END IF

      IF (GEOSIM(OSCALE).NE.0D0) THEN
        SC2 = 1D0/GEOSIM(OSCALE)
      ELSE 
        SC2 = 0D0
      END IF

C Shift to origin of the coordinate system and scale

      XP = SC1*(X-CORSYS(ORIGX))
      YP = SC1*(Y-CORSYS(ORIGY))

C Rotate according to the rotation of the coordinate sysrem

      IF (CORSYS(ORIGRT).NE.0D0) THEN
C        C = COS(CORSYS(ORIGRT)*PI/180D0)
C        S = SIN(CORSYS(ORIGRT)*PI/180D0)
        C = CORSYS(ORICOS)
        S = CORSYS(ORISIN)
        XP1 = C*XP+S*YP
        YP1 = -S*XP+C*YP
      ELSE 
        XP1 = XP
        YP1 = YP
      END IF

C Shift to reference point and scale

      XP1 = SC2*(XP1-GEOSIM(OREFX))
      YP1 = SC2*(YP1-GEOSIM(OREFY))

C Rotate according to the rotation of the coordinate sysrem

      IF (GEOSIM(OROT).NE.0D0) THEN
C        C = COS(GEOSIM(OROT)*PI/180D0)
C        S = SIN(GEOSIM(OROT)*PI/180D0)
        C = GEOSIM(OROTCOS)
        S = GEOSIM(OROTSIN)
        XP = C*XP1+S*YP1
        YP = -S*XP1+C*YP1
      ELSE 
        XP = XP1
        YP = YP1
      END IF

      END 

************************************************************************
* General routines for all geometries
************************************************************************

*-----------------------------------------------------------------------
* Create coordinate system
*
* This routine is designed to perform a standard initialization for
* a general coordinate system. It accepts a parameter CORSYS pointing
* to a coordinate system TCORSYS, and will initialize that according
* to the other parameters.
*-----------------------------------------------------------------------

      SUBROUTINE GCRSYS (CORSYS, X, Y, ROT, SCALE)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION CORSYS(SZCORS),X,Y,ROT,SCALE

      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)

      CALL LCL1(CORSYS,SZCORS)
      
      CORSYS(ORIGX)  = X
      CORSYS(ORIGY)  = Y
      CORSYS(ORIGRT) = ROT
      CORSYS(ORISIN) = SIN(ROT*PI/180D0)
      CORSYS(ORICOS) = COS(ROT*PI/180D0)
      CORSYS(ORISCL) = SCALE
      
      END

*-----------------------------------------------------------------------
* Create structure
*
* Type-identifier of the structure will be set to SGTNON.
* Information in X,Y,ROT,SCALE are transferred to geometry object.
*-----------------------------------------------------------------------
      INTEGER FUNCTION GCRSIM (GEOSIM, X, Y, ROT, SCALE, INVER)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGSIM),X,Y,ROT,SCALE
      INTEGER INVER

      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)
      
      GEOSIM(OTYPE)  = SGTNON
      GEOSIM(OINVER) = DBLE(INVER)
      GEOSIM(OREFX)  = X
      GEOSIM(OREFY)  = Y
      GEOSIM(OROT)   = ROT
      GEOSIM(OROTSIN) = SIN(ROT*PI/180D0)
      GEOSIM(OROTCOS) = COS(ROT*PI/180D0)
      GEOSIM(OSCALE)  = SCALE
      
      GCRSIM = SZGSIM
      
      END

*-----------------------------------------------------------------------
* Update coordinate system
*
* This routine allowes direct update of the coordinate system of
* an object. It will directly transfer all the information from the
* parameters into the coordinate system of the given object without
* modifying anything else.
*-----------------------------------------------------------------------
      SUBROUTINE GUPCOR (GEOSIM, X, Y, ROT, SCALE)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGSIM),X,Y,ROT,SCALE

      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)
      
      GEOSIM(OREFX) = X
      GEOSIM(OREFY) = Y
      GEOSIM(OROT)  = ROT
      GEOSIM(OROTSIN) = SIN(ROT*PI/180D0)
      GEOSIM(OROTCOS) = COS(ROT*PI/180D0)
      GEOSIM(OSCALE)  = SCALE
      
      END

************************************************************************
* Simple Geometry: Circle
*
* The "reference point" is the midpoint of the circle.
************************************************************************

*-----------------------------------------------------------------------
* Create structure
*
* Additional variables:
*  RAD    - radius of the circle
*-----------------------------------------------------------------------
      INTEGER FUNCTION GCRCIR (GEOSIM, X, Y, ROT, SCALE, INVER, RAD)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGCIR),X,Y,ROT, SCALE
      DOUBLE PRECISION RAD
      INTEGER INVER
      INTEGER GCRSIM
      EXTERNAL GCRSIM
C Dismiss return value of standard initialisation
      GCRCIR = GCRSIM (GEOSIM, X, Y, ROT, SCALE, INVER)
      GEOSIM(OTYPE)  = SGTCIR
      GEOSIM(OCRAD)  = RAD
      GCRCIR = SZGCIR
      END
*-----------------------------------------------------------------------
* Test if a point is inside
*-----------------------------------------------------------------------
      INTEGER FUNCTION GINCIR (GEOSIM,CORSYS, X, Y)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGCIR),X,Y
      DOUBLE PRECISION CORSYS(SZCORS)
      DOUBLE PRECISION XP,YP,RAD,DIST

C Transform the point back to the standard coordinate system

      CALL GRTBCK(GEOSIM,CORSYS, X,Y, XP,YP)
      
      RAD=GEOSIM(OCRAD)
C Using a circle we don't have to take care about rotation...
C Calculate the (squared) distance
      DIST=XP**2+YP**2
      IF (DIST.LE.RAD**2) THEN
        GINCIR = 1
      ELSE
        GINCIR = 0
      END IF
      
C Switch the result if we are in "inverted" mode
      
      IF (GEOSIM(OINVER).NE.0) THEN
        GINCIR = -GINCIR
      END IF
      
      END
*-----------------------------------------------------------------------
* Calculate volume
*-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION GVLCIR (GEOSIM,CORSYS)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGCIR)
      DOUBLE PRECISION CORSYS(SZCORS)
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)
      GVLCIR = PI*GEOSIM(OCRAD)**2
C Switch the sign if the object is inverted
      IF (GEOSIM(OINVER).NE.0) THEN
        GVLCIR = -GVLCIR
      END IF
      END
*-----------------------------------------------------------------------
* Calculate length of interface
*-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION GALCIR (GEOSIM,CORSYS)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGCIR)
      DOUBLE PRECISION CORSYS(SZCORS)
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)
      GALCIR = 2D0*PI*GEOSIM(OCRAD)
      END
*-----------------------------------------------------------------------
* Write GMV-Data
*-----------------------------------------------------------------------
      SUBROUTINE GGMCIR (GEOSIM,CORSYS, MUNIT, MMAT)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGCIR)
      DOUBLE PRECISION CORSYS(SZCORS)
      INTEGER MUNIT,MMAT
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)

      INTEGER I,NL
      DOUBLE PRECISION RAD,DEG,X,Y,ROT
      
C      return
      
C Calculate the position of the midpoint by transformation

      CALL GRTTRF (GEOSIM,CORSYS,0D0,0D0, X,Y)

C Calculate the radius with the help of the scaling factors of the two
C coordinate systems      
      RAD = GEOSIM(OCRAD)*GEOSIM(OSCALE)*CORSYS(ORISCL)
C Number of circle segments is dependent of radius to get fine resolution 
      NL = MIN(10000,MAX(1,NINT(RAD*2*PI*1000)))
C Material number MMAT, NL+1 nodes
      WRITE(MUNIT,*) MMAT,NL+1
C X-coordinates
      DO I=0,NL
        DEG = DBLE(I)*2D0*PI/DBLE(NL)
        WRITE(MUNIT,'(E20.10)') X+RAD*COS(DEG)
      END DO
C Y-coordinates
      DO I=0,NL
        DEG = DBLE(I)*2D0*PI/DBLE(NL)
        WRITE(MUNIT,'(E20.10)') Y+RAD*SIN(DEG)
      END DO
C Z-coordinates
      DO I=0,NL
        WRITE(MUNIT,'(E20.10)') 0E0
      END DO

C Line from midpoint to topmost point      
      
      WRITE(MUNIT,*) MMAT,2
      
      ROT = GEOSIM(OROT)*PI/180D0

      WRITE(MUNIT,'(E20.10)') X
      WRITE(MUNIT,*) X-RAD*SIN(ROT)

      WRITE(MUNIT,'(E20.10)') Y
      WRITE(MUNIT,'(E20.10)') Y+RAD*COS(ROT)

      WRITE(MUNIT,'(E20.10)') 0E0
      WRITE(MUNIT,'(E20.10)') 0E0
      
C Horizontal line

      WRITE(MUNIT,*) MMAT,2

      WRITE(MUNIT,'(E20.10)') X+RAD*COS(ROT)
      WRITE(MUNIT,'(E20.10)') X-RAD*COS(ROT)

      WRITE(MUNIT,'(E20.10)') Y-RAD*SIN(ROT)
      WRITE(MUNIT,'(E20.10)') Y+RAD*SIN(ROT)

      WRITE(MUNIT,'(E20.10)') 0E0
      WRITE(MUNIT,'(E20.10)') 0E0
      
C     To show the polygons filled with a color: Switch on
C     CELLS/SHADED and POLYGONS/SHADED in GMV.
C     To show the orientation: Switch on
C     POLYGONS/LINES in GMV. If POLYGONS/SHADED is also
C     activated, the orientation is shown by white lines!
      
      END
*-----------------------------------------------------------------------
* Calculate normal vector
*-----------------------------------------------------------------------
      SUBROUTINE GNMCIR (GEOSIM,CORSYS, X,Y, XN,YN)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGCIR),X,Y,XN,YN
      DOUBLE PRECISION CORSYS(SZCORS)
      
      DOUBLE PRECISION A
      INTEGER GINCIR
      EXTERNAL GINCIR

C Get the midpoint of the circle

      CALL GRTTRF(GEOSIM,CORSYS, 0D0,0D0, XN,YN)

C Calculate the difference vector to the actual point to get the normal.
C If the point is inside of the circle, take the inverse vector

      IF (GINCIR (GEOSIM,CORSYS, X, Y).LE.0) THEN
        XN = X-XN
        YN = Y-YN
      ELSE
        XN = XN-X
        YN = YN-Y
      END IF  

C Normalize
      IF ((XN.NE.0D0).OR.(YN.NE.0D0)) THEN
        A = 1D0/DSQRT(XN*XN+YN*YN)
        XN=A*XN
        YN=A*YN
      END IF
      
C Switch the result if we are in "inverted" mode
      
      IF (GEOSIM(OINVER).NE.0) THEN
        XN = -XN
        YN = -YN
      END IF
      
      END
*-----------------------------------------------------------------------
* Calculate distance
*-----------------------------------------------------------------------
      INTEGER FUNCTION GDSCIR (GEOSIM,CORSYS, X,Y, DIST)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGCIR),X,Y,DIST
      DOUBLE PRECISION CORSYS(SZCORS)
      DOUBLE PRECISION XP,YP,RAD

C Transform the point back to the standard coordinate system

      CALL GRTBCK(GEOSIM,CORSYS, X,Y, XP,YP)

      DIST = DSQRT(XP**2+YP**2) - GEOSIM(OCRAD)
      GDSCIR = 1
      
C Switch the result if we are in "inverted" mode
      
      IF (GEOSIM(OINVER).NE.0) THEN
        DIST = -DIST
      END IF
      
C Scale the result according to the scaling factor in the coordinate system

      DIST = DIST*CORSYS(ORISCL)
      
      END
*-----------------------------------------------------------------------
* Project onto boundary interface
*-----------------------------------------------------------------------
      INTEGER FUNCTION GPRCIR (GEOSIM,CORSYS, X,Y, XP,YP)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGCIR),X,Y,XP,YP
      DOUBLE PRECISION CORSYS(SZCORS)
      DOUBLE PRECISION XN,YN,RAD,DIST
      INTEGER I,GDSCIR
      EXTERNAL GDSCIR
      
C Transform the point back to the standard coordinate system

      CALL GRTBCK(GEOSIM,CORSYS, X,Y, XN,YN)
      
C Scale the vector (XP,YP) to norm 1 - this gives the projection

      DIST = 1D0/SQRT(XN**2+YN**2)
      XN=DIST*XN
      YN=DIST*YN
      
C Transform the point to the real coordinate system, finish

      CALL GRTTRF(GEOSIM,CORSYS, XN,YN, XP,YP)

      GPRCIR = 1
      
      END

************************************************************************
* Simple Geometry: Ellipse
*
* The "reference point" is the midpoint of the ellipse.
************************************************************************

*-----------------------------------------------------------------------
* Create structure
*
* Additional variables:
*  RADX   - X-radius of the circle
*  RADY   - Y-radius of the circle
*-----------------------------------------------------------------------
      INTEGER FUNCTION GCRELL (GEOSIM, X, Y, ROT, SCALE, INVER, 
     *                         RADX, RADY)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGELL),X,Y,ROT, SCALE
      DOUBLE PRECISION RADX,RADY
      INTEGER INVER
      INTEGER GCRSIM
      EXTERNAL GCRSIM
C Dismiss return value of standard initialisation
      GCRELL = GCRSIM (GEOSIM, X, Y, ROT, SCALE, INVER)
      GEOSIM(OTYPE)  = SGTELL
      GEOSIM(OERADX) = RADX
      GEOSIM(OERADY) = RADY
      GCRELL = SZGELL
      END
*-----------------------------------------------------------------------
* Test if a point is inside
*-----------------------------------------------------------------------
      INTEGER FUNCTION GINELL (GEOSIM,CORSYS, X, Y)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGELL),X,Y
      DOUBLE PRECISION CORSYS(SZCORS)
      DOUBLE PRECISION XM,YM,RADX,RADY,XP,YP
      
      DOUBLE PRECISION DIST
      
      XM=GEOSIM(OREFX)
      YM=GEOSIM(OREFY)
      RADX=GEOSIM(OERADX)
      RADY=GEOSIM(OERADY)

C Rotate the point to the "reference coordinate system" at the origin

      CALL GRTBCK (GEOSIM,CORSYS, X,Y, XP,YP)
      
C Get the distance to the origin

      DIST=( (XP**2)/(RADX**2) ) +
     *     ( (YP**2)/(RADY**2) )

C Building the square root is not necessary as we test against 1...
      
      IF (DIST.LE.1D0) THEN
        GINELL = 1
      ELSE
        GINELL = 0
      END IF

C Switch the result if we are in "inverted" mode
      
      IF (GEOSIM(OINVER).NE.0) THEN
        GINELL = -GINELL
      END IF
      
      END
*-----------------------------------------------------------------------
* Calculate volume
*-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION GVLELL (GEOSIM,CORSYS)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGELL)
      DOUBLE PRECISION CORSYS(SZCORS)
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)
      GVLELL = PI*GEOSIM(OERADX)*GEOSIM(OERADY)
C Switch the sign if the object is inverted
      IF (GEOSIM(OINVER).NE.0) THEN
        GVLELL = -GVLELL
      END IF
      END
*-----------------------------------------------------------------------
* Calculate length of interface
*-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION GALELL (GEOSIM,CORSYS)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGELL)
      DOUBLE PRECISION CORSYS(SZCORS)
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)
      DOUBLE PRECISION XR,YR
      XR = GEOSIM(OERADX)
      YR = GEOSIM(OERADY)
      GALELL = 0.5D0*PI*DSQRT(XR**2+YR**2)+2D0*PI*XR*YR
      END
*-----------------------------------------------------------------------
* Write GMV-Data
*-----------------------------------------------------------------------
      SUBROUTINE GGMELL (GEOSIM,CORSYS, MUNIT, MMAT)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGELL)
      DOUBLE PRECISION CORSYS(SZCORS)
      INTEGER MUNIT, MMAT

      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)

      INTEGER I,NL
      DOUBLE PRECISION RADX,RADY,DEG,XP,YP

C Calculate the position of the midpoint by transformation.
C Calculate the radii with the scaling factors in the coordinate systems

      RADX = GEOSIM(OERADX)
      RADY = GEOSIM(OERADY)
C Number of circle segments is dependent of radius to get fine resolution
      NL = MIN(10000,MAX(1,NINT(MAX(RADX*GEOSIM(OSCALE)*CORSYS(ORISCL),
     *              RADY*GEOSIM(OSCALE)*CORSYS(ORISCL))*2*PI*1000)))
C Material number MMAT, NL+1 nodes
      WRITE(MUNIT,*) MMAT,NL+1
C X-coordinates
      DO I=0,NL
        DEG = DBLE(I)*2D0*PI/DBLE(NL)
C Rotate and translate to obtain XP
        CALL GRTTRF (GEOSIM,CORSYS, RADX*COS(DEG),RADY*SIN(DEG),
     *               XP,YP)
        WRITE(MUNIT,'(E20.10)') XP
      END DO
C Y-coordinates
      DO I=0,NL
        DEG = DBLE(I)*2D0*PI/DBLE(NL)
C Rotate and translate to obtain YP
        CALL GRTTRF (GEOSIM,CORSYS, RADX*COS(DEG),RADY*SIN(DEG),
     *               XP,YP)
        WRITE(MUNIT,'(E20.10)') YP
      END DO
C Z-coordinates
      DO I=0,NL
        WRITE(MUNIT,'(E20.10)') 0E0
      END DO

C Line from midpoint to topmost point      
      
      WRITE(MUNIT,*) MMAT,2
      
      CALL GRTTRF (GEOSIM,CORSYS, 0D0,0D0,
     *             XP,YP)
      WRITE(MUNIT,'(E20.10)') XP
      CALL GRTTRF (GEOSIM,CORSYS, 0D0,RADY,
     *             XP,YP)
      WRITE(MUNIT,'(E20.10)') XP

      CALL GRTTRF (GEOSIM,CORSYS, 0D0,0D0,
     *             XP,YP)
      WRITE(MUNIT,'(E20.10)') YP
      CALL GRTTRF (GEOSIM,CORSYS, 0D0,RADY,
     *             XP,YP)
      WRITE(MUNIT,'(E20.10)') YP

      WRITE(MUNIT,'(E20.10)') 0E0
      WRITE(MUNIT,'(E20.10)') 0E0
      
C Horizontal line

      WRITE(MUNIT,*) MMAT,2

      CALL GRTTRF (GEOSIM,CORSYS, RADX,0D0,
     *             XP,YP)
      WRITE(MUNIT,'(E20.10)') XP
      CALL GRTTRF (GEOSIM,CORSYS, -RADX,0D0,
     *             XP,YP)
      WRITE(MUNIT,'(E20.10)') XP

      CALL GRTTRF (GEOSIM,CORSYS, RADX,0D0,
     *             XP,YP)
      WRITE(MUNIT,'(E20.10)') YP
      CALL GRTTRF (GEOSIM,CORSYS, -RADX,0D0,
     *             XP,YP)
      WRITE(MUNIT,'(E20.10)') YP

      WRITE(MUNIT,'(E20.10)') 0E0
      WRITE(MUNIT,'(E20.10)') 0E0
      

      END
*-----------------------------------------------------------------------
* Calculate normal vector
*-----------------------------------------------------------------------
      SUBROUTINE GNMELL (GEOSIM,CORSYS, X,Y, XN,YN)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGELL),X,Y,XN,YN
      DOUBLE PRECISION CORSYS(SZCORS)
C not implemented
      XN = 0D0
      YN = 0D0

C Switch the result if we are in "inverted" mode
      
      IF (GEOSIM(OINVER).NE.0) THEN
        XN = -XN
        YN = -YN
      END IF
      
      END
*-----------------------------------------------------------------------
* Calculate distance
*-----------------------------------------------------------------------
      INTEGER FUNCTION GDSELL (GEOSIM,CORSYS, X,Y, DIST)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGELL),X,Y,DIST
      DOUBLE PRECISION CORSYS(SZCORS)
      
      INTEGER GPRELL,GINELL
      EXTERNAL GPRELL,GINELL
      
      INTEGER I
      DOUBLE PRECISION XP,YP
      
      GDSELL = 1
      
C Project the point onto the boundary

      I = GPRELL (GEOSIM,CORSYS, X,Y, XP,YP)
      
C Measure the distance to that point

      DIST = DSQRT((XP-X)**2+(YP-Y)**2)

C Change the sign if we are in the square

      IF (GINELL (GEOSIM,CORSYS, X,Y).GT.0) DIST = -DIST
      
C Switch the result if we are in "inverted" mode
      
      IF (GEOSIM(OINVER).NE.0) THEN
        DIST = -DIST
      END IF

      END
*-----------------------------------------------------------------------
* Project onto boundary interface
*-----------------------------------------------------------------------
      INTEGER FUNCTION GPRELL (GEOSIM,CORSYS, X,Y, XP,YP)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGELL),X,Y,XP,YP
      DOUBLE PRECISION CORSYS(SZCORS)

      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)

      DOUBLE PRECISION PL,PR,T,D,C1,C2,C3,S,C,XPR,YPR
      DOUBLE PRECISION LENX,LENY,LEN,TMP1,TMP2
      INTEGER I
      
C first we define some auxiliary functions, see below

      DOUBLE PRECISION DS
      DOUBLE PRECISION DSS
      DS(C1,C2,C3,C,S) = (C1*S-C2*C-C3*S*C)
      DSS(C1,C2,C3,C,S) = (C1*C-C2*S+C3-C3*C**2)
      
      GPRELL = 1

      LENX = GEOSIM(OERADX)
      LENY = GEOSIM(OERADY)
      LEN = LENY

C Special treatment if the ellipse degenerates to a line:

      IF (LENX.EQ.0) THEN
      
C Rotate the point to the "reference coordinate system" at the origin.

        CALL GRTBCK (GEOSIM,CORSYS, X,Y, XP,YP)
      
C Only 3 cases: Above, below and in the middle of the line

        XPR = 0

        IF (YP.GT.LENY) THEN
          YPR = LENY
        ELSE IF (YP.LT.-LENY) THEN
          YPR = -LENY
        ELSE
          YPR = YP
        END IF

C Finally rotate the projected point to obtain the coordinates
C in the correct coordinate system.
      
        CALL GRTTRF (GEOSIM,CORSYS, XPR,YPR, XP,YP)
      
C Finish
      
        RETURN

      END IF

      IF (LENY.EQ.0) THEN
      
C Rotate the point to the "reference coordinate system" at the origin.

        CALL GRTBCK (GEOSIM,CORSYS, X,Y, XP,YP)
      
C Only 3 cases: Left, right and in the middle of the line

        YPR = 0

        IF (XP.GT.LENX) THEN
          XPR = LENX
        ELSE IF (XP.LT.-LENX) THEN
          XPR = -LENX
        ELSE
          XPR = XP
        END IF

C Finally rotate the projected point to obtain the coordinates
C in the correct coordinate system.
      
        CALL GRTTRF (GEOSIM,CORSYS, XPR,YPR, XP,YP)
      
C Finish
      
        RETURN

      END IF

C Rotate the point to the "reference coordinate system" at the origin.

      CALL GRTBCK (GEOSIM,CORSYS, X,Y, XP,YP)

C Calculating the real projection onto the interface of an ellipse
C is AWFUL! One has to solve an equation of 4th order, which is
C really hard and in no relationship to the efford one has to
C investigate!

C We make another approach here to only calculate an approximation
C to the distance: a small fix-point iteration. 
C We parametrise the circle in the parameter interval 0..2*PI
C and choose a left and a right "boundary" parameter according
C to the quadrant that contains (X,Y):
      
      IF (YP.GE.0) THEN
        IF (XP.GE.0) THEN
C top right quadrant
          PL = 0
          PR = 0.5D0*PI
        ELSE
C top left quadrant
          PL = 0.5D0*PI
          PR = PI
        END IF
      ELSE
        IF (XP.GE.0) THEN
C bottom right quadrant
          PL = 1.5D0*PI
          PR = 2D0*PI
        ELSE
C bottom left quadrant
          PL = PI
          PR = 1.5D0*PI
        END IF
      END IF
      
C Newtom-iteration to obtain a parameter value t in [0..2*PI] that
C gives us a point on the interface wih minimum distance.
C Let: (p,q)=(a*cos(t),b*sin(t)) the point on the interface we are searching.
C      d(p,q) the distance of (p,q) to (x,y)
C      d(t) = d(p(t),q(t)) the distance after parametrization (squared)
C           = ( x-a*cos(t) )^2 + ( y-b*sin(t) )^2
C           -> min
C  =>  d'(t) = 2ax sin(t) - 2by cos(t) + 2(b^2-a^2)sin(t)cos(t) 
C  with
C      d''(t) =   2ax cos(t) + 2by sin(t) 
C               + 2(a^2-b^2) + 4(b^2-a^2)cos^2(t)
C We want to solve:
C      d'(t) = 0
C Newton-iteration formula for that:
C      t_n+1 = t_n - d'(t)/d''(t)
C Our separation into the different quadrants guarantee that the
C function d'(t) is convex; should be enough for convergence.
      
C Calculate the coefficients of the derivatives in advance:

      C1 = 2D0*LENX*XP
      C2 = 2D0*LENY*YP
      C3 = 2D0*(LENX**2-LENY**2)
      
C First perform some bisection steps to come closer to our searched
C parameter:

      T = PL
      DO I=1,5
        T = 0.5D0*(PL+PR)
        IF (DS(C1,C2,C3,COS(PL),SIN(PL))*
     *      DS(C1,C2,C3,COS(T),SIN(T)) < 0D0) THEN
          PR = T
        ELSE
          PL = T
        END IF
      END DO
      
C Do some Newton-steps
      
      DO I=1,8
        S = SIN(T)
        C = COS(T)
        TMP1 = DS(C1,C2,C3,C,S)
        TMP2 = DSS(C1,C2,C3,C,S)
        
        IF (TMP2.LT.1D-40) THEN
C Cancel if the denominator is too low for Newton!
C Maybe the case for the midpoint of the ellipse!
          GOTO 10
        END IF
        
        T = T - TMP1 / TMP2
        
C Cancel if we can't see convergence. Maybe the case for the
C midpoint of the ellipse!
C In such a case we could do bisection...

        IF ((T.LT.PL).OR.(T.GT.PR)) THEN
          T = PL
          GOTO 10
        END IF
        
      END DO
      
C Calculate the point from the parameter value

10    XPR = LENX*COS(T)
      YPR = LENY*SIN(T)
      
C Finally rotate the projected point to obtain the coordinates
C in the correct coordinate system.
      
      CALL GRTTRF (GEOSIM,CORSYS, XPR,YPR, XP,YP)

      END

************************************************************************
* Simple Geometry: Square
*
* The "reference point" is the midpoint of the square.
************************************************************************

*-----------------------------------------------------------------------
* Create structure
*
* Additional variables:
*  LEN    - length of the edges
*-----------------------------------------------------------------------
      INTEGER FUNCTION GCRSQR (GEOSIM, X, Y, ROT, SCALE, INVER, 
     *                         LEN)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGCIR),X,Y,ROT, SCALE
      DOUBLE PRECISION LEN
      INTEGER INVER
      INTEGER GCRSIM
      EXTERNAL GCRSIM
C Dismiss return value of standard initialisation
      GCRSQR = GCRSIM (GEOSIM, X, Y, ROT, SCALE, INVER)
      GEOSIM(OTYPE)  = SGTSQR
      GEOSIM(OSLEN)  = LEN
      GCRSQR = SZGSQR
      END
*-----------------------------------------------------------------------
* Test if a point is inside
*-----------------------------------------------------------------------
      INTEGER FUNCTION GINSQR (GEOSIM,CORSYS, X, Y)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGSQR),X,Y
      DOUBLE PRECISION CORSYS(SZCORS)
      DOUBLE PRECISION XM,YM,LEN,XP,YP,DIST

      LEN = GEOSIM(OSLEN)

C Rotate the point to the "reference coordinate system" at the origin

      CALL GRTBCK (GEOSIM,CORSYS, X,Y, XP,YP)

      DIST=MAX(ABS(XP),ABS(YP))
      
C Test against the half of the edge-length
      
      IF (DIST.LE.0.5D0*LEN) THEN
        GINSQR = 1
      ELSE
        GINSQR = 0
      END IF

C Switch the result if we are in "inverted" mode
      
      IF (GEOSIM(OINVER).NE.0) THEN
        GINSQR = -GINSQR
      END IF

      END
*-----------------------------------------------------------------------
* Calculate volume
*-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION GVLSQR (GEOSIM,CORSYS)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGSQR)
      DOUBLE PRECISION CORSYS(SZCORS)
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)
      GVLSQR = GEOSIM(OSLEN)**2
C Switch the sign if the object is inverted
      IF (GEOSIM(OINVER).NE.0) THEN
        GVLSQR = -GVLSQR
      END IF
      END
*-----------------------------------------------------------------------
* Calculate length of interface
*-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION GALSQR (GEOSIM,CORSYS)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGSQR)
      DOUBLE PRECISION CORSYS(SZCORS)

      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)

      GALSQR = 4D0*GEOSIM(OSLEN)
      END
*-----------------------------------------------------------------------
* Write GMV-Data
*-----------------------------------------------------------------------
      SUBROUTINE GGMSQR (GEOSIM,CORSYS, MUNIT, MMAT)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGSQR)
      DOUBLE PRECISION CORSYS(SZCORS)
      INTEGER MUNIT,MMAT
      
      DOUBLE PRECISION XP(4),YP(4),X,Y
      INTEGER I
      
C Material number MMAT, 4+1 nodes
      
      WRITE(MUNIT,*) MMAT,5
      
C Calculate the coordinates of the corners

      CALL GRTTRF(GEOSIM,CORSYS, -0.5D0*GEOSIM(OSLEN),
     *            -0.5D0*GEOSIM(OSLEN),
     *            XP(1),YP(1))
      CALL GRTTRF(GEOSIM,CORSYS, 0.5D0*GEOSIM(OSLEN),
     *            -0.5D0*GEOSIM(OSLEN),
     *            XP(2),YP(2))
      CALL GRTTRF(GEOSIM,CORSYS, 0.5D0*GEOSIM(OSLEN),
     *            0.5D0*GEOSIM(OSLEN),
     *            XP(3),YP(3))
      CALL GRTTRF(GEOSIM,CORSYS, -0.5D0*GEOSIM(OSLEN),
     *            0.5D0*GEOSIM(OSLEN),
     *            XP(4),YP(4))
      
C X-coordinates

      DO I=1,5
        WRITE(MUNIT,'(E20.10)') XP(MOD(I,4)+1)
      END DO

C Y-coordinates

      DO I=1,5
        WRITE(MUNIT,'(E20.10)') YP(MOD(I,4)+1)
      END DO

C Z-coordinates

      DO I=1,5
        WRITE(MUNIT,'(E20.10)') 0E0
      END DO

      END
*-----------------------------------------------------------------------
* Calculate normal vector
*-----------------------------------------------------------------------
      SUBROUTINE GNMSQR (GEOSIM,CORSYS, X,Y, XN,YN)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGSQR),X,Y,XN,YN
      DOUBLE PRECISION CORSYS(SZCORS)
C not implemented
      XN = 0D0
      YN = 0D0

C Switch the result if we are in "inverted" mode
      
      IF (GEOSIM(OINVER).NE.0) THEN
        XN = -XN
        YN = -YN
      END IF
      
      END
*-----------------------------------------------------------------------
* Calculate distance
*-----------------------------------------------------------------------
      INTEGER FUNCTION GDSSQR (GEOSIM,CORSYS, X,Y, DIST)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGSQR),X,Y,DIST
      DOUBLE PRECISION CORSYS(SZCORS)
      
      INTEGER GPRSQR,GINSQR
      EXTERNAL GPRSQR,GINSQR
      
      INTEGER I
      DOUBLE PRECISION XP,YP
      
      GDSSQR = 1
      
C Project the point onto the boundary

      I = GPRSQR (GEOSIM,CORSYS, X,Y, XP,YP)
      
C Measure the distance to that point

      DIST = DSQRT((XP-X)**2+(YP-Y)**2)
      
C Change the sign if we are in the square

      IF (GINSQR (GEOSIM,CORSYS, X,Y).GT.0) DIST = -DIST
      
C Switch the result if we are in "inverted" mode
      
      IF (GEOSIM(OINVER).NE.0) THEN
        DIST = -DIST
      END IF
      
      END
*-----------------------------------------------------------------------
* Project onto boundary interface
*-----------------------------------------------------------------------
      INTEGER FUNCTION GPRSQR (GEOSIM,CORSYS, X,Y, XP,YP)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGSQR),X,Y,XP,YP
      DOUBLE PRECISION CORSYS(SZCORS)
      
      DOUBLE PRECISION XM,YM,LEN,XPR,YPR

      GPRSQR = 1

      LEN = 0.5D0*GEOSIM(OSLEN)

C Rotate the point to the "reference coordinate system" at the origin

      CALL GRTBCK (GEOSIM,CORSYS, X,Y, XP,YP)

C Decide where we are and choose a "projection" point.

      IF (XP.LE.-LEN) THEN
C Left side of the square. 
        XPR = -LEN
        IF (YP.GE.LEN) THEN
C Top left corner
          YPR = LEN
        ELSE IF (YP.LE.-LEN) THEN
C Bottom left corner
          YPR = -LEN
        ELSE
C a point inbetween
          YPR = YP
        END IF
      
      ELSE IF (XP.GE.LEN) THEN
C Right side of the square. 
        XPR = LEN
        IF (YP.GE.LEN) THEN
C Top left corner
          YPR = LEN
        ELSE IF (YP.LE.-LEN) THEN
C Bottom left corner
          YPR = -LEN
        ELSE
C a point inbetween
          YPR = YP
        END IF

      ELSE
      
C Horizontally in the square. The vertical position
C decides where our "projected point" has to be assumed.

        IF (YP.GE.LEN) THEN
C On the top edge
          XPR = XP
          YPR = LEN
        ELSE IF (YP.LE.-LEN) THEN
C On the bottom edge
          XPR = XP
          YPR = -LEN
        ELSE
        
C We are inside of the square. 

          IF (YP.GE.0D0) THEN
          
C Top half-plane

            IF (YP.GE.ABS(XP)) THEN
C Project to top edge
              XPR = XP
              YPR = LEN
            ELSE
C Project to left or right edge
              YPR = YP
              IF (XP.LT.0) THEN
                XPR = -LEN
              ELSE
                XPR = LEN
              END IF
            END IF
            
          ELSE
          
C Bottom half-plane
            IF (ABS(YP).GE.ABS(XP)) THEN
C Project to bottom edge
              XPR = XP
              YPR = -LEN
            ELSE
C Project to left or right edge
              YPR = YP
              IF (XP.LT.0) THEN
                XPR = -LEN
              ELSE
                XPR = LEN
              END IF
            END IF
            
          END IF
          
        END IF
        
      END IF

C Finally rotate the projected point to obtain the coordinates
C in the correct coordinate system
      
      CALL GRTTRF (GEOSIM,CORSYS, XPR,YPR, XP,YP)

      END

************************************************************************
* Simple Geometry: Rectangle
*
* The "reference point" is the midpoint of the rectangle.
************************************************************************

*-----------------------------------------------------------------------
* Create structure
*
* Additional variables:
*  LENX   - radius of the X-edges
*  LENY   - radius of the Y-edges
*-----------------------------------------------------------------------
      INTEGER FUNCTION GCRREC (GEOSIM, X, Y, ROT, SCALE, INVER, 
     *                         LENX, LENY)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGREC),X,Y,ROT, SCALE
      DOUBLE PRECISION LENX,LENY
      INTEGER INVER
      INTEGER GCRSIM
      EXTERNAL GCRSIM
C Dismiss return value of standard initialisation
      GCRREC = GCRSIM (GEOSIM, X, Y, ROT, SCALE, INVER)
      GEOSIM(OTYPE)  = SGTREC
      GEOSIM(ORLENX)  = LENX
      GEOSIM(ORLENY)  = LENY
      GCRREC = SZGREC
      END
*-----------------------------------------------------------------------
* Test if a point is inside
*-----------------------------------------------------------------------
      INTEGER FUNCTION GINREC (GEOSIM,CORSYS, X, Y)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGREC),X,Y
      DOUBLE PRECISION CORSYS(SZCORS)
      
      DOUBLE PRECISION DIST,XP,YP
      
C Rotate the point to the "reference coordinate system" at the origin

      CALL GRTBCK (GEOSIM,CORSYS, X,Y, XP,YP)
      
C Get the distance to the origin.
C Be careful with the transformation of the X/Y-coordinates as
C the origin is oriented in the midpoint of the rectangle
C => we have to divide by LEN/2 - equivalent to *2/LEN

      DIST=MAX(ABS( 2D0 * XP / GEOSIM(ORLENX) ),
     *         ABS( 2D0 * YP / GEOSIM(ORLENY) ) )

C Building the square root is not necessary as we test against 1...
      
      IF (DIST.LE.1D0) THEN
        GINREC = 1
      ELSE
        GINREC = 0
      END IF

C Switch the result if we are in "inverted" mode
      
      IF (GEOSIM(OINVER).NE.0) THEN
        GINREC = -GINREC
      END IF

      END
*-----------------------------------------------------------------------
* Calculate volume
*-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION GVLREC (GEOSIM,CORSYS)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGREC)
      DOUBLE PRECISION CORSYS(SZCORS)
      GVLREC = GEOSIM(ORLENX)*GEOSIM(ORLENY)
C Switch the sign if the object is inverted
      IF (GEOSIM(OINVER).NE.0) THEN
        GVLREC = -GVLREC
      END IF
      END
*-----------------------------------------------------------------------
* Calculate length of interface
*-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION GALREC (GEOSIM,CORSYS)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGREC)
      DOUBLE PRECISION CORSYS(SZCORS)
      GALREC = 2D0*(GEOSIM(ORLENX)+GEOSIM(ORLENY))
      END
*-----------------------------------------------------------------------
* Write GMV-Data
*-----------------------------------------------------------------------
      SUBROUTINE GGMREC (GEOSIM,CORSYS, MUNIT, MMAT)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGREC)
      DOUBLE PRECISION CORSYS(SZCORS)
      INTEGER MUNIT, MMAT
      
      DOUBLE PRECISION XP(4),YP(4)
      INTEGER I
      
C Material number MMAT, 4+1 nodes
      
      WRITE(MUNIT,*) MMAT,5
      
C Calculate the coordinates of the corners

      CALL GRTTRF(GEOSIM,CORSYS, 
     *            -0.5D0*GEOSIM(ORLENX),-0.5D0*GEOSIM(ORLENY),
     *            XP(1),YP(1))
      CALL GRTTRF(GEOSIM,CORSYS, 
     *            0.5D0*GEOSIM(ORLENX),-0.5D0*GEOSIM(ORLENY),
     *            XP(2),YP(2))
      CALL GRTTRF(GEOSIM,CORSYS, 
     *            0.5D0*GEOSIM(ORLENX),0.5D0*GEOSIM(ORLENY),
     *            XP(3),YP(3))
      CALL GRTTRF(GEOSIM,CORSYS, 
     *            -0.5D0*GEOSIM(ORLENX),0.5D0*GEOSIM(ORLENY),
     *            XP(4),YP(4))
      
C X-coordinates

      DO I=1,5
        WRITE(MUNIT,'(E20.10)') XP(MOD(I,4)+1)
      END DO

C Y-coordinates

      DO I=1,5
        WRITE(MUNIT,'(E20.10)') YP(MOD(I,4)+1)
      END DO

C Z-coordinates

      DO I=1,5
        WRITE(MUNIT,'(E20.10)') 0E0
      END DO
      END
*-----------------------------------------------------------------------
* Calculate normal vector
*-----------------------------------------------------------------------
      SUBROUTINE GNMREC (GEOSIM,CORSYS, X,Y, XN,YN)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGREC),X,Y,XN,YN
      DOUBLE PRECISION CORSYS(SZCORS)
C not implemented
      XN = 0D0
      YN = 0D0

C Switch the result if we are in "inverted" mode
      
      IF (GEOSIM(OINVER).NE.0) THEN
        XN = -XN
        YN = -YN
      END IF
      
      END
*-----------------------------------------------------------------------
* Calculate distance
*-----------------------------------------------------------------------
      INTEGER FUNCTION GDSREC (GEOSIM,CORSYS, X,Y, DIST)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGREC),X,Y,DIST
      DOUBLE PRECISION CORSYS(SZCORS)
      INTEGER I
      DOUBLE PRECISION XP,YP
      
      INTEGER GPRREC,GINREC
      EXTERNAL GPRREC,GINREC
      
      GDSREC = 1
      
C Project the point onto the boundary

      I = GPRREC (GEOSIM,CORSYS, X,Y, XP,YP)
      
C Measure the distance to that point

      DIST = DSQRT((XP-X)**2+(YP-Y)**2)
      
C Change the sign if we are in the square

      IF (GINREC (GEOSIM,CORSYS, X,Y).GT.0) DIST = -DIST

C Switch the result if we are in "inverted" mode
      
      IF (GEOSIM(OINVER).NE.0) THEN
        DIST = -DIST
      END IF
      
      END
*-----------------------------------------------------------------------
* Project onto boundary interface
*-----------------------------------------------------------------------
      INTEGER FUNCTION GPRREC (GEOSIM,CORSYS, X,Y, XP,YP)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGREC),X,Y,XP,YP
      DOUBLE PRECISION CORSYS(SZCORS)
      
      DOUBLE PRECISION LENX,LENY,LEN,XPR,YPR

      GPRREC = 1

      LENX = 0.5D0*GEOSIM(ORLENX)
      LENY = 0.5D0*GEOSIM(ORLENY)
      LEN = LENY

C Special treatment if the rectangle degenerates to a line:

      IF (LENX.EQ.0) THEN
      
C Rotate the point to the "reference coordinate system" at the origin.

        CALL GRTBCK (GEOSIM,CORSYS, X,Y, XP,YP)
      
C Only 3 cases: Above, below and in the middle of the line

        XPR = 0

        IF (YP.GT.LENY) THEN
          YPR = LENY
        ELSE IF (YP.LT.-LENY) THEN
          YPR = -LENY
        ELSE
          YPR = YP
        END IF

C Finally rotate the projected point to obtain the coordinates
C in the correct coordinate system.
      
        CALL GRTTRF (GEOSIM,CORSYS, XPR,YPR, XP,YP)
      
C Finish
      
        RETURN

      END IF

      IF (LENY.EQ.0) THEN
      
C Rotate the point to the "reference coordinate system" at the origin.

        CALL GRTBCK (GEOSIM,CORSYS, X,Y, XP,YP)
      
C Only 3 cases: Left, right and in the middle of the line

        YPR = 0

        IF (XP.GT.LENX) THEN
          XPR = LENX
        ELSE IF (XP.LT.-LENX) THEN
          XPR = -LENX
        ELSE
          XPR = XP
        END IF

C Finally rotate the projected point to obtain the coordinates
C in the correct coordinate system.
      
        CALL GRTTRF (GEOSIM,CORSYS, XPR,YPR, XP,YP)
      
C Finish
      
        RETURN

      END IF

C Scale the X- or Y-coordinate such that the rectangle deforms to a 
C square. Then we can use the routines that are known from the square 
C to calculate the projected point.
C Unfortunately this only works for rectangles, not for ellipses...

C Rotate the point to the "reference coordinate system" at the origin.

      CALL GRTBCK (GEOSIM,CORSYS, X,Y, XP,YP)

      IF (LENX.LT.LENY) THEN

        LEN = LENY
        XP = XP*LENY/LENX
        
      ELSE 
        
        LEN = LENX
        YP = YP*LENX/LENY
        
      END IF

C Decide where we are and choose a "projection" point.

      IF (XP.LE.-LEN) THEN
C Left side of the square. 
        XPR = -LEN
        IF (YP.GE.LEN) THEN
C Top left corner
          YPR = LEN
        ELSE IF (YP.LE.-LEN) THEN
C Bottom left corner
          YPR = -LEN
        ELSE
C a point inbetween
          YPR = YP
        END IF
      
      ELSE IF (XP.GE.LEN) THEN
C Right side of the square. 
        XPR = LEN
        IF (YP.GE.LEN) THEN
C Top left corner
          YPR = LEN
        ELSE IF (YP.LE.-LEN) THEN
C Bottom left corner
          YPR = -LEN
        ELSE
C a point inbetween
          YPR = YP
        END IF

      ELSE
      
C Horizontally in the square. The vertical position
C decides where our "projected point" has to be assumed.

        IF (YP.GE.LEN) THEN
C On the top edge
          XPR = XP
          YPR = LEN
        ELSE IF (YP.LE.-LEN) THEN
C On the bottom edge
          XPR = XP
          YPR = -LEN
        ELSE
        
C We are inside of the square. 

          IF (YP.GE.0D0) THEN
          
C Top half-plane

            IF (YP.GE.ABS(XP)) THEN
C Project to top edge
              XPR = XP
              YPR = LEN
            ELSE
C Project to left or right edge
              YPR = YP
              IF (XP.LT.0) THEN
                XPR = -LEN
              ELSE
                XPR = LEN
              END IF
            END IF
            
          ELSE
          
C Bottom half-plane
            IF (ABS(YP).GE.ABS(XP)) THEN
C Project to bottom edge
              XPR = XP
              YPR = -LEN
            ELSE
C Project to left or right edge
              YPR = YP
              IF (XP.LT.0) THEN
                XPR = -LEN
              ELSE
                XPR = LEN
              END IF
            END IF
            
          END IF
          
        END IF
        
      END IF

C Finally scale the X- or Y-coordinate back and 
C rotate the projected point to obtain the coordinates
C in the correct coordinate system.

      IF (LENX.LT.LENY) THEN
        CALL GRTTRF (GEOSIM,CORSYS, XPR * LENX/LENY,YPR, XP,YP)
      ELSE
        CALL GRTTRF (GEOSIM,CORSYS, XPR,YPR * LENY/LENX, XP,YP)
      END IF

      END

************************************************************************
* Simple Geometry: SIN-Wave
*
* The "reference point" is the midpoint of the Wave.
************************************************************************

*-----------------------------------------------------------------------
* Create structure
*
* Additional variables:
*  WID    - Width of the SIN-Wave
*  POS1   - Y-Position of the upper SIN-Wave; should be > 0
*  HEIGH1 - Height of the upper SIN-Wave
*  PHAS1  - Phase-Shift of the upper SIN-Wave
*  FREQ1  - Frequency-multiplicator of the upper SIN-Wave
*  POS2   - Y-Position of the lower SIN-Wave; should be < 0
*  HEIGH2 - Height of the lower SIN-Wave
*  PHAS2  - Phase-Shift of the lower SIN-Wave
*  FREQ2  - Frequency-multiplicator of the lower SIN-Wave
*-----------------------------------------------------------------------
      INTEGER FUNCTION GCRSWV (GEOSIM, X, Y, ROT, SCALE, INVER, 
     *        WID,POS1,HEIGH1,PHAS1,FREQ1,POS2,HEIGH2,PHAS2,FREQ2)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGSWV),X,Y,ROT, SCALE
      DOUBLE PRECISION WID,HEIGH1,PHAS1,FREQ1,HEIGH2,PHAS2,FREQ2
      DOUBLE PRECISION POS1,POS2
      INTEGER INVER
      INTEGER GCRSIM
      EXTERNAL GCRSIM
C Dismiss return value of standard initialisation
      GCRSWV = GCRSIM (GEOSIM, X, Y, ROT, SCALE, INVER)
      GEOSIM(OTYPE)  = SGTSWV
      GEOSIM(OWWID)  = WID
      GEOSIM(OWPOS1)   = POS1
      GEOSIM(OWHEI1)   = HEIGH1
      GEOSIM(OWPHAS1)  = PHAS1
      GEOSIM(OWFREQ1)  = FREQ1
      GEOSIM(OWPOS2)   = POS2
      GEOSIM(OWHEI2)   = HEIGH2
      GEOSIM(OWPHAS2)  = PHAS2
      GEOSIM(OWFREQ2)  = FREQ2
      GCRSWV = SZGSWV
      END
*-----------------------------------------------------------------------
* Auxiliary function
*
* At a defined X-position calculate the Y-positions of the upper-
* and lower SIN-Wave.
*
* In:
*  GEOSIM - The geometry
*  X      - The X-coordinate - in the reference coordinate system
*           around (0,0) with no rotation!
* Out
*  Y1, Y2 - The upper and lower Y-coordinate of the two SIN-Waves
*-----------------------------------------------------------------------
      SUBROUTINE GA1SWV (GEOSIM, X, Y1,Y2)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGSWV),X, Y1,Y2
      DOUBLE PRECISION HEIGH,PHAS,FREQ,POS,WID

      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)

      WID = GEOSIM(OWWID)

C Obtain geometry information from the array - for the upper
C SIN-Wave
      POS = GEOSIM(OWPOS1)
      HEIGH = GEOSIM(OWHEI1)
      PHAS  = GEOSIM(OWPHAS1)*PI/180D0
      FREQ  = GEOSIM(OWFREQ1)*2D0*PI
C Calculate the Y-coordinate
      Y1 = POS+HEIGH*SIN(X*FREQ/WID+PHAS)
C Obtain the information about the lower SIN-Wave
      POS = GEOSIM(OWPOS2)
      HEIGH = GEOSIM(OWHEI2)
      PHAS  = GEOSIM(OWPHAS2)*PI/180D0
      FREQ  = GEOSIM(OWFREQ2)*2D0*PI
C Calculate the Y-coordinate
      Y2 = POS+HEIGH*SIN(X*FREQ/WID+PHAS)
      
      END

*-----------------------------------------------------------------------
* Test if a point is inside
*-----------------------------------------------------------------------
      INTEGER FUNCTION GINSWV (GEOSIM,CORSYS, X, Y)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGSWV),X,Y
      DOUBLE PRECISION CORSYS(SZCORS)
      DOUBLE PRECISION XP,YP,Y1,Y2,YT,WID

      GINSWV = 0

      WID = 0.5D0*GEOSIM(OWWID)

C Transform the point back to the standard coordinate system

      CALL GRTBCK(GEOSIM,CORSYS, X,Y, XP,YP)
      
C Calculate the upper- and lower Y-position of the SIN-Wave

      CALL GA1SWV (GEOSIM, XP, Y1,Y2)
      
C Switch the Y-coordinates if they are in the wrong order

      IF (Y2.LT.Y1) THEN
        YT=Y2
        Y2=Y1
        Y1=YT
      END IF
      
C Is the X-coordinate in range?

      IF ((XP.LT.-WID).OR.(XP.GE.WID)) GOTO 10
      
C Is the Y-coordinate in range?
      
      IF ((YP.LT.Y1).OR.(YP.GT.Y2)) GOTO 10
      
C ok, we are inside

      GINSWV = 1
      
10    CONTINUE

C Switch the result if we are in "inverted" mode
      
      IF (GEOSIM(OINVER).NE.0) THEN
        GINSWV = -GINSWV
      END IF

      END
*-----------------------------------------------------------------------
* Calculate volume
*-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION GVLSWV (GEOSIM,CORSYS)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGCIR)
      DOUBLE PRECISION CORSYS(SZCORS)
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)
C not implemented
      GVLSWV = 0D0
C Switch the sign if the object is inverted
      IF (GEOSIM(OINVER).NE.0) THEN
        GVLSWV = -GVLSWV
      END IF
      END
*-----------------------------------------------------------------------
* Calculate length of interface
*-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION GALSWV (GEOSIM,CORSYS)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGCIR)
      DOUBLE PRECISION CORSYS(SZCORS)
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)
C not implemented
      GALSWV = -1D0
      END
*-----------------------------------------------------------------------
* Write GMV-Data
*-----------------------------------------------------------------------
      SUBROUTINE GGMSWV (GEOSIM,CORSYS, MUNIT, MMAT)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGSWV)
      DOUBLE PRECISION CORSYS(SZCORS)
      INTEGER MUNIT,MMAT
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)

      INTEGER I,NL1,NL2
      DOUBLE PRECISION WID,STP,Y1,Y2,X,XP,YP,FREQ1,STP1,FREQ2,STP2
      
      WID = 0.5D0*GEOSIM(OWWID)
      
C For every wave we reserve 100 polygon segments

C First consider the upper SIN-Wave

      FREQ1  = GEOSIM(OWFREQ1)
      NL1 = NINT(100D0*FREQ1)

C then the lower one

      FREQ2 = GEOSIM(OWFREQ2)
      NL2 = NINT(100D0*FREQ2)

C Material number MMAT, NL1+NL2+2 nodes
      WRITE(MUNIT,*) MMAT,NL1+NL2+2
      
C X-coordinates of the upper SIN-Wave
      DO I=0,NL1
        X = WID*DBLE(2*I-NL1)/DBLE(NL1)
C Calculate the Y-coordinate
        CALL GA1SWV (GEOSIM, X, Y1,Y2)
C Project the coordinates to obtain the correct X/Y-position
        CALL GRTTRF (GEOSIM,CORSYS, X,Y1, XP,YP)
        WRITE(MUNIT,'(E20.10)') XP
      END DO
C X-coordinates of the lower SIN-Wave
      DO I=NL2,0,-1
        X = WID*DBLE(2*I-NL2)/DBLE(NL2)
C Calculate the Y-coordinate
        CALL GA1SWV (GEOSIM, X, Y1,Y2)
C Project the coordinates to obtain the correct X/Y-position
        CALL GRTTRF (GEOSIM,CORSYS, X,Y2, XP,YP)
        WRITE(MUNIT,'(E20.10)') XP
      END DO

C Y-coordinates of the upper SIN-Wave
      DO I=0,NL1
        X = WID*DBLE(2*I-NL1)/DBLE(NL1)
C Calculate the Y-coordinate
        CALL GA1SWV (GEOSIM, X, Y1,Y2)
C Project the coordinates to obtain the correct X/Y-position
        CALL GRTTRF (GEOSIM,CORSYS, X,Y1, XP,YP)
        WRITE(MUNIT,'(E20.10)') YP
      END DO
C Y-coordinates of the lower SIN-Wave
      DO I=NL2,0,-1
        X = WID*DBLE(2*I-NL2)/DBLE(NL2)
C Calculate the Y-coordinate
        CALL GA1SWV (GEOSIM, X, Y1,Y2)
C Project the coordinates to obtain the correct X/Y-position
        CALL GRTTRF (GEOSIM,CORSYS, X,Y2, XP,YP)
        WRITE(MUNIT,'(E20.10)') YP
      END DO

C Z-coordinates = 0
      DO I=1,NL1+NL2+2
        WRITE(MUNIT,'(E20.10)') 0E0
      END DO

      END
*-----------------------------------------------------------------------
* Calculate normal vector
*-----------------------------------------------------------------------
      SUBROUTINE GNMSWV (GEOSIM,CORSYS, X,Y, XN,YN)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGSWV),X,Y,XN,YN
      DOUBLE PRECISION CORSYS(SZCORS)
      
      INTEGER GINCIR
      EXTERNAL GINCIR
C not implemented
      XN = 0D0
      YN = 0D0

C Switch the result if we are in "inverted" mode
      
      IF (GEOSIM(OINVER).NE.0) THEN
        XN = -XN
        YN = -YN
      END IF
      
      END
*-----------------------------------------------------------------------
* Calculate distance
*-----------------------------------------------------------------------
      INTEGER FUNCTION GDSSWV (GEOSIM,CORSYS, X,Y, DIST)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGCIR),X,Y,DIST
      DOUBLE PRECISION CORSYS(SZCORS)

      INTEGER GPRSWV,GINSWV
      EXTERNAL GPRSWV,GINSWV
      
      INTEGER I
      DOUBLE PRECISION XP,YP
      
      GDSSWV = 1
      
C Project the point onto the boundary

      I = GPRSWV (GEOSIM,CORSYS, X,Y, XP,YP)
      
C Measure the distance to that point

      DIST = DSQRT((XP-X)**2+(YP-Y)**2)
      
C Change the sign if we are in the object

      IF (GINSWV (GEOSIM,CORSYS, X,Y).GT.0) DIST = -DIST

C Switch the result if we are in "inverted" mode
      
      IF (GEOSIM(OINVER).NE.0) THEN
        DIST = -DIST
      END IF
      
C Scale the distance according to the coordinate system
      
      DIST = DIST*CORSYS(ORISCL)

      END
*-----------------------------------------------------------------------
* Project onto boundary interface
*-----------------------------------------------------------------------
      INTEGER FUNCTION GPRSWV (GEOSIM,CORSYS, X,Y, XP,YP)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGCIR),X,Y,XP,YP
      DOUBLE PRECISION CORSYS(SZCORS)
      DOUBLE PRECISION XN,YN
      DOUBLE PRECISION WID,Y1,Y2,YT

      WID = 0.5D0*GEOSIM(OWWID)

C Transform the point back to the standard coordinate system

      CALL GRTBCK(GEOSIM,CORSYS, X,Y, XP,YP)

C We project only in the Y-direction, if we are above/below the wave.

C Are we left/right from the waves?

      IF ((XP.LT.-WID).OR.(XP.GT.WID)) THEN
      
        IF (XP.LT.-WID) THEN
          XN = -WID
        ELSE
          XN = WID
        END IF
        
C Calculate the leftmost corners
        
        CALL GA1SWV (GEOSIM, XN, Y1,Y2)
        
        IF (Y2.LT.Y1) THEN
          YT=Y2
          Y2=Y1
          Y1=YT
        END IF

        IF (YP.LT.Y1) THEN
          YN = Y1
        ELSE IF (YP.GT.Y2) THEN
          YN = Y2
        ELSE
          YN = YP
        END IF
        
      ELSE
      
C We are above or below the wave. Project in Y-direction
C onto the curve.

        XN = XP

        CALL GA1SWV (GEOSIM, XN, Y1,Y2)

        IF (ABS(YP-Y2).LT.ABS(YP-Y1)) THEN
          YN = Y2
        ELSE
          YN = Y1
        END IF

      END IF

C Transform into the real coordinate system again

      CALL GRTTRF (GEOSIM,CORSYS, XN,YN, XP,YP)

      GPRSWV = 1

      END

************************************************************************
* Simple Geometry: User-defined pie-slice
*
* The "Pie slice" is a polygon given by a couple of line segments.
* Originally this is used for defining parts of a pointsymmetrical
* object (like a rotor). 
*
* The "reference point" is the midpoint of the object.
* See "sgeometries.inc" for details.
************************************************************************

*-----------------------------------------------------------------------
* Create structure
*
* Additional variables:
*  ICOUNT - number of points in boundary interface
*  ICLOS  - 0=pie slice is open
*           1=pie slice is closed
*  XBAS,
*  YBAS   - basis point of pie slice (midpoint of surrounding circle)
*  POINTS - array [1..2,1..COUNT] of double
*           X/Y-coordinates of all the points forming the bondary 
*           interface
*-----------------------------------------------------------------------
      INTEGER FUNCTION GCRPSL (GEOSIM, X, Y, ROT, SCALE, INVER,
     *                 ICOUNT, ICLOS, XBAS, YBAS, POINTS)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(SZGPSL),X,Y,ROT, SCALE
      INTEGER ICOUNT, ICLOS, INVER
      DOUBLE PRECISION XBAS, YBAS, POINTS(2,ICOUNT)
      INTEGER GCRSIM
      EXTERNAL GCRSIM
      
      INTEGER I
      
C Dismiss return value of standard initialisation
      GCRPSL = GCRSIM (GEOSIM, X, Y, ROT, SCALE, INVER)

      GCRPSL = SZGPSL+2*ICOUNT
      
      GEOSIM(OTYPE)  = SGTPSL
      GEOSIM(OPCNT)  = ICOUNT
      GEOSIM(OPCLOS) = ICLOS
      GEOSIM(OPXBAS) = XBAS
      GEOSIM(OPYBAS) = YBAS
      GEOSIM(OPILMP) = 0
      GEOSIM(OPIRMP) = 0

      DO I=1,ICOUNT

C Transfer all points. 

        GEOSIM(OPPOIN+2*(I-1)) = POINTS (1,I)
        GEOSIM(OPPOIN+2*(I-1)+1) = POINTS (2,I)

C Find the "leftmost" and "rightmost" point and store their indices.
C The basis point is the origin here.

C to be implemented...

      END DO
      
C Most easiest implementation: First point is rightmost, last point
C is leftmost. A general implementation can later be done...

      GEOSIM(OPIRMP) = 1
      GEOSIM(OPILMP) = ICOUNT
      
      END
*-----------------------------------------------------------------------
* Test if a point is inside
*-----------------------------------------------------------------------
      INTEGER FUNCTION GINPSL (GEOSIM,CORSYS, X, Y)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(*),X,Y
      DOUBLE PRECISION CORSYS(SZCORS)
      DOUBLE PRECISION XP,YP,RAD,DIST

C Transform the point back to the standard coordinate system

      CALL GRTBCK(GEOSIM,CORSYS, X,Y, XP,YP)
      
C Rest: not yet implemented.
      
      GINPSL = 0

C Switch the result if we are in "inverted" mode
      
      IF (GEOSIM(OINVER).NE.0) THEN
        GINPSL = -GINPSL
      END IF
      
      END
*-----------------------------------------------------------------------
* Calculate volume
*-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION GVLPSL (GEOSIM,CORSYS)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(*)
      DOUBLE PRECISION CORSYS(SZCORS)
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)
C not implemented
      GVLPSL = 0D0
C Switch the sign if the object is inverted
      IF (GEOSIM(OINVER).NE.0) THEN
        GVLPSL = -GVLPSL
      END IF
      END
*-----------------------------------------------------------------------
* Calculate length of interface
*-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION GALPSL (GEOSIM,CORSYS)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(*)
      DOUBLE PRECISION CORSYS(SZCORS)
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)
C not implemented
      GALPSL = -1D0
      END
*-----------------------------------------------------------------------
* Write GMV-Data
*-----------------------------------------------------------------------
      SUBROUTINE GGMPSL (GEOSIM,CORSYS, MUNIT, MMAT)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(*)
      DOUBLE PRECISION CORSYS(SZCORS)
      INTEGER MUNIT,MMAT
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)

      INTEGER I,CNT
      DOUBLE PRECISION XP,YP
      
C Material number MMAT, #nodes from the structure
      
      CNT = AINT(GEOSIM(OPCNT))
      
      IF (GEOSIM(OPCLOS).EQ.0) THEN
C an open pie slice has 2*CNT-2 nodes - 1x CNT nodes and 1x CNT-2 nodes
C for the way back.
        WRITE(MUNIT,*) MMAT,2*CNT-2
      ELSE
C A closed pie slice has exactly CNT nodes
        WRITE(MUNIT,*) MMAT,CNT
      END IF
      
C X-coordinates

      DO I=0,CNT-1
C Project the coordinates to obtain the correct X/Y-position
        CALL GRTTRF (GEOSIM,CORSYS, 
     *               GEOSIM(OPPOIN+2*I),GEOSIM(OPPOIN+2*I+1), XP,YP)
        
        WRITE(MUNIT,'(E20.10)') XP
      END DO

C When the pie-slice is open, just head the way back.
C So GMV paints the boundary in fact twice. This is to prevent some
C nasty lines connecting the first and last node!

      IF (GEOSIM(OPCLOS).EQ.0) THEN
        DO I=CNT-2,1,-1
C Project the coordinates to obtain the correct X/Y-position
          CALL GRTTRF (GEOSIM,CORSYS, 
     *                 GEOSIM(OPPOIN+2*I),GEOSIM(OPPOIN+2*I+1), XP,YP)
          
          WRITE(MUNIT,'(E20.10)') XP
        END DO
      END IF

C Y-coordinates

      DO I=0,CNT-1
C Project the coordinates to obtain the correct X/Y-position
        CALL GRTTRF (GEOSIM,CORSYS, 
     *               GEOSIM(OPPOIN+2*I),GEOSIM(OPPOIN+2*I+1), XP,YP)
        WRITE(MUNIT,'(E20.10)') YP
      END DO

C When the pie-slice is open, just head the way back.
C So GMV paints the boundary in fact twice. This is to prevent some
C nasty lines connecting the first and last node!

      IF (GEOSIM(OPCLOS).EQ.0) THEN
        DO I=CNT-2,1,-1
C Project the coordinates to obtain the correct X/Y-position
          CALL GRTTRF (GEOSIM,CORSYS, 
     *                 GEOSIM(OPPOIN+2*I),GEOSIM(OPPOIN+2*I+1), XP,YP)
          
          WRITE(MUNIT,'(E20.10)') YP
        END DO
      END IF

C Z-coordinates

      DO I=0,CNT-1
        WRITE(MUNIT,'(E20.10)') 0E0
      END DO

      IF (GEOSIM(OPCLOS).EQ.0) THEN
        DO I=CNT-2,1,-1
          WRITE(MUNIT,'(E20.10)') 0E0
        END DO
      END IF

      END
*-----------------------------------------------------------------------
* Calculate normal vector
*-----------------------------------------------------------------------
      SUBROUTINE GNMPSL (GEOSIM,CORSYS, X,Y, XN,YN)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(*),X,Y,XN,YN
      DOUBLE PRECISION CORSYS(SZCORS)
      
      DOUBLE PRECISION A

C not implemented

      XN = 0D0
      YN = 0D0

C Switch the result if we are in "inverted" mode
      
      IF (GEOSIM(OINVER).NE.0) THEN
        XN = -XN
        YN = -YN
      END IF
      
      END
*-----------------------------------------------------------------------
* Calculate distance
*-----------------------------------------------------------------------
      INTEGER FUNCTION GDSPSL (GEOSIM,CORSYS, X,Y, DIST)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(*),X,Y,DIST
      DOUBLE PRECISION CORSYS(SZCORS)
      DOUBLE PRECISION XP,YP,RAD

C Transform the point back to the standard coordinate system

      CALL GRTBCK(GEOSIM,CORSYS, X,Y, XP,YP)

C Rest: not implemented

      DIST = 0D0
      GDSPSL = 0
      
C Switch the result if we are in "inverted" mode
      
      IF (GEOSIM(OINVER).NE.0) THEN
        DIST = -DIST
      END IF
      
C Scale the distance according to the coordinate system
      
      DIST = DIST*CORSYS(ORISCL)

      END
*-----------------------------------------------------------------------
* Project onto boundary interface
*-----------------------------------------------------------------------
      INTEGER FUNCTION GPRPSL (GEOSIM,CORSYS, X,Y, XP,YP)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      DOUBLE PRECISION GEOSIM(*),X,Y,XP,YP
      DOUBLE PRECISION CORSYS(SZCORS)
      DOUBLE PRECISION XN,YN,RAD,DIST
      INTEGER I,GDSCIR
      EXTERNAL GDSCIR
      
C Transform the point back to the standard coordinate system

      CALL GRTBCK(GEOSIM,CORSYS, X,Y, XN,YN)

C rest: not implemented      

      GPRPSL = 0
      
      END
