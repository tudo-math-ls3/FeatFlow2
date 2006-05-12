************************************************************************
* FBDCRE - Create fictitious boundary structures
*
* This routine is called on start of the program after the parameters
* of the geometry are read into memory. It has to initialize
* fictitious boundary structures as necessary.
*
* For simple-type objects, normally nothing has to be done. For more
* complex objects (e.g. 100x100 balls, all with their own coordinate
* system or so), it might be necessary to allocate additional memory
* and build structures.
*
* In:
*  IGEOM - array [1..*] of integer
*  DGEOM - array [1..*] of double precision
*          Parameter blocks describing the geometry
*
* The routine might store necessary information in the IGEOUD/DGEOUD
* blocks of IGEOM/DGEOM. Allocated memory must be relased in FBDDIS!
************************************************************************

      SUBROUTINE FBDCRE (IGEOM,DGEOM)
      
      IMPLICIT NONE

      INCLUDE 'sinigeometry.inc'
      
      INTEGER IGEOM(*)
      DOUBLE PRECISION DGEOM(*)

      END
      
************************************************************************
* FBDDIS - Release fictitious boundary structures
*
* This routine is called to clean up fictitious boundary structures.
* It's called at the end of the program to clean up. Handles of
* allocated, which might be stored in IGEOUD/DGEOUD of the
* IGEOM/DGEOM data structures, must be released!
*
* In:
*  IGEOM - array [1..*] of integer
*  DGEOM - array [1..*] of double precision
*          Parameter blocks describing the geometry
************************************************************************

      SUBROUTINE FBDDIS (IGEOM,DGEOM)
      
      IMPLICIT NONE

      INCLUDE 'sinigeometry.inc'
      
      INTEGER IGEOM(*)
      DOUBLE PRECISION DGEOM(*)

      END
      
************************************************************************
* NFBDYC - Number of fictitious boundary components
*
* User provided function. Has to return the number of fictitious
* boundary components in the given geometry.
*
* Remark: If 0 is returned, fictitious boundary handling is disabled !!!
*
* In:
*  IGEOM - array [1..*] of integer
*  DGEOM - array [1..*] of double precision
*          Parameter blocks describing the geometry
************************************************************************

      INTEGER FUNCTION NFBDYC (IGEOM,DGEOM)
      
      IMPLICIT NONE
      
      INCLUDE 'sinigeometry.inc'
      
      INTEGER IGEOM(*)
      DOUBLE PRECISION DGEOM(*)
      
C     By default, we switch off fictitious boundary handling:

      NFBDYC = 0
      RETURN
      
C     Only one f.b.c.;
C     If we have a composed geometry, we have more perhaps. This can 
C     be seen in the IGEOCN variable.
      
      IF (IGEOM(OICIRTP).LE.50) THEN
        NFBDYC = 1
      ELSE
        NFBDYC = IGEOM(OIGEOCN)
      END IF
      
      END
      
************************************************************************
* FBDINF - Obtain information about a fictitious boundary object
*
* User provided function. Has to return information about a special
* fictitious boundary object.
*
* In:
*  IFBC  - Number of the fictitious boundary object, the routine
*          should return information about.
*          Fictitious boundary components are numbered
*          in the range NBCT+1..NBCT+NFBDYC with numbers behind
*          the normal boundary components.
*  IGEOM - array [1..*] of integer
*  DGEOM - array [1..*] of double precision
*          Parameter blocks describing the geometry
*
* Out:
*  IFIX  - =0: The object is a moving object which moves with
*              the flow depending on the object's characteristics.
*              Parameter DMASS must be defined.
*          =1: The object is a fixed object which does not move
*              with the flow. DMASS can be set to 0.0
*          =2: The object is at a fixed position but can rotate
*  XMASS - X-position of the mass center of the object
*  YMASS - Y-position of the mass center of the object
*  DMASS - Mass of the object; can be set to 0.0, if IFIX=1
*  DINRT - Moment of inertia for that object. This is
*                I = int_V |x|*rho dx
*          with V the domain of the object (ball,...)
*          For most objects there are explicit formulas,
*          see:
*             http://hyperphysics.phy-astr.gsu.edu/hbase/mi.html
*          Can be set to 0.0 if IFIX=1.
************************************************************************

      SUBROUTINE FBDINF (IFBC,
     *                   IFIX,XMASS,YMASS,DMASS,DINRT,
     *                   IGEOM,DGEOM)
      
      IMPLICIT NONE
      
      INCLUDE 'sinigeometry.inc'
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)
      
      INTEGER IGEOM(*)
      DOUBLE PRECISION DGEOM(*)
      INTEGER IFBC,IFIX
      DOUBLE PRECISION XMASS,YMASS,DMASS,DINRT
      
C     Return information from our geometry arrays.

      XMASS = DGEOM(ODCXPOS)
      YMASS = DGEOM(ODCYPOS)

C     The object is fixed and thus does not have a mass.

      IFIX = 1
      DMASS = 0.0
      DINRT = 0.0
      
C     Moment of inertia for a ball:
C       I = 2/5 * M * R^2
C     Mass of the ball:
C       M = int_V rho dx
C     For constant density:
C       M = rho int_V 1 dx = rho * V
C     with the volume
C       V = pi*R^2
C     of the ball. So the moment of inertia for a ball
C     with constant density is
C       I = 2/5 * rho * pi * R^4

      IF (IGEOM(OICIRTP).EQ.0) THEN

        IFIX = 0
        DMASS = 2D0*DGEOM(ODCRAD)*PI
        DINRT = (2D0*DMASS*DGEOM(ODCRAD)**2)/5D0
      
      END IF
      
      END
      
************************************************************************
* ISFBDY - Test for being in a fictitious boundary object
*
* This is a user provided routine. It has to check whether the given
* point (x,y) is in an object surrounded by a fictitious boundary.
* The return value must be = 0 if the point is a normal point and
* != 0 if the point is in an object with a fictitious boundary.
*
* in:
*  X     - X-coordinate of the point to test
*  Y     - Y-coordinate of the point to test
*  IFCB  - =0 or NBCT+1..NBCT+NFBDYC
*          =0: test, if (X,Y) is in any fictitious boundary component;
*              the return value is either 0 or the number of any
*              fictitious boundary component containing this point
*          >0: test, if (X,Y) is contained especially in boundary
*              component IFBC. The return value is either 0 or
*              +-IFBC, depending on whether IFBC contains (X,Y) or not.
*          Fictitious boundary components are numbered
*          in the range NBCT+1..NBCT+NFBDYC with numbers behind
*          the normal boundary components.
*  IGEOM - array [1..*] of integer
*  DGEOM - array [1..*] of double precision
*          Parameter blocks describing the geometry
*
* result:
*   0 - If (x,y) is not in a fictitious boundary object
*  >0 - If (x,y) is in a rigid fictitious boundary object.
*  <0 - If (x,y) if in a virtual fictitious boundary object
*
*  The absolute value defines the number of a/the fictitious boundary
*  component containing (X,Y). If IFCB=0, this is the number
*  of any of the fictitious boundary components containing (X,Y).
*  If IFCB>0, the return value is 0 or +-IFCB, depending on if
*  (X,Y) is contained especially in this component or not.
*
*  A positive number indicates a rigid fictitious boundary object (wall)
*  which must be treated as dirichlet-boundary in the framework.
*  A negative number indicates a "virtual" fictitious boundary
*  object that is treated as Neumann boundary. In such points the
*  framework does not perform any changes to the matrix/RHS/...
*  to implement any type of boundary conditions. This can be used
*  e.g. for defining a region where to refine the grid without
*  changing the type of geometry in that region.
*
* This routine has to consider the number of fixed boundary components,
* too. As all fixed boundary conponents take the numbers 1..NBCT,
* all moving boundary components take numbers NBCT+1..NBCT+NFBDYC.
* So the lowest number != 0 returned by this function in IFBC has to be
* +-(NBCT+1) !
************************************************************************

      INTEGER FUNCTION ISFBDY (X,Y,IFBC,IGEOM,DGEOM)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'sinigeometry.inc'
      
C     parameters

      DOUBLE PRECISION X,Y,DGEOM(*)
      INTEGER IFBC,IGEOM(*)

C     externals

      INTEGER NFBDYC,TNBC
      EXTERNAL NFBDYC,TNBC
      
      INTEGER GINCMP,GINWRP
      EXTERNAL GINCMP,GINWRP

C     local variables

      DOUBLE PRECISION XM1,YM1,DIST
      DOUBLE PRECISION XR,YR,XR1,YR1
      INTEGER I
      INTEGER L1,L2,L0
      
      DOUBLE PRECISION GEOTIM,DCRRLX,DCRRLY
      INTEGER ICIRTP,INSGEO
      DOUBLE PRECISION DCXPOS, DCYPOS, DCROT, DCRADX, DCRADY, DCRAD
      DOUBLE PRECISION DCRSIN,DCRCOS

      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)

C     quick parameter check

      IF ((IFBC.NE.0).AND.
     *  ((IFBC.LE.TNBC()).OR.(IFBC.GT.TNBC()+NFBDYC(IGEOM,DGEOM)))) THEN
        WRITE (*,'(A,I,A)') 'ERROR in call to ISFBDY()! IFBC=',IFBC,
     *                      ' Program halted!'
        STOP
      END IF

C     standard return value

      ISFBDY = 0
      
C     Get rotation, radius,... from the parameter blocks

      ICIRTP = IGEOM(OICIRTP)
      DCXPOS = DGEOM(ODCXPOS)
      DCYPOS = DGEOM(ODCYPOS)
      DCROT  = DGEOM(ODCROT )
      DCRAD  = DGEOM(ODCRAD )
      DCRADX = DGEOM(ODCRDX)
      DCRADY = DGEOM(ODCRDY)
      DCRSIN = DGEOM(ODCRSIN)
      DCRCOS = DGEOM(ODCRCOS)
      DCRRLX = DGEOM(ODCRLX)
      DCRRLY = DGEOM(ODCRLY)
      
      INSGEO = IGEOM(OINSGEO)
      
C     Get simulational time, i.e. time stamp how the geometry should
C     be acting
      
      GEOTIM = DGEOM(OGEOTIM)
      
C     In this implementation we only have one boundary component, so we
C     don't have to take care of IFBC for now...

      IF (ICIRTP.EQ.0) THEN

C ----------------------------------------------------------------------
C fictitious boundary component: circle at DCXPOS/DCYPOS, radius DCRAD
C ----------------------------------------------------------------------

        XM1=DCXPOS
        YM1=DCYPOS
 
C     Possibility:       
C     Pos = 1.1+sim(timens)
C     Vel = cos(timens)
C        XM1 = DCXPOS+0.5D0*SIN(0.8D0*GEOTIM)

C       Benchmark formula in nonsteady case.
C       Fixed point in steady case
      
        IF (INSGEO.EQ.1) THEN
          XM1 = DCXPOS+0.25*SIN(2D0*PI*0.25*GEOTIM)
        ELSE
          XM1 = DCXPOS
        END IF

        DIST=DSQRT((X-XM1)**2+(Y-YM1)**2)

        IF (DIST.LE.(DCRAD+DCRRLX)) THEN
          ISFBDY = (TNBC()+1)
        END IF

      ELSE IF (ICIRTP.EQ.1) THEN

C ----------------------------------------------------------------------
C fictitious boundary component: square with midpoint at DCXPOS/DCYPOS, 
C width/height=2xDCRAD
C ----------------------------------------------------------------------

        XM1=DCXPOS
        YM1=DCYPOS

        DIST=MAX(ABS(X-XM1),ABS(Y-YM1))

        IF (DIST.LE.(DCRAD+DCRRLX)) THEN
          ISFBDY = TNBC()+1
        END IF
        
      ELSE IF (ICIRTP.EQ.2) THEN

C       ---------------------------------------------------------------
C       fictitious boundary component: ellipse with midpoint 
C       at DCXPOS/DCYPOS, width=2*DCRADX, height=2*DCRADY.
C       The ellipse can be rotated by the angle DCROT!
C       ---------------------------------------------------------------

        XM1=DCXPOS
        YM1=DCYPOS

        IF (DCROT.EQ.0D0) THEN

C         Standard handling, no rotation

          DIST=( ((X-XM1)**2)/((DCRADX+DCRRLX)**2) ) +
     *         ( ((Y-YM1)**2)/((DCRADY+DCRRLY)**2) )

C         Building the square root is not necessary as we test against 1...

        ELSE
        
C         Extended handling: rotated ellipse

C         We expect DCRSIN/DCRCOS to be initialised properly in order to
C         calculate the rotation matrix.

C         Rotate the point "back" around our center:

C         Shift to 0

          XR1 = (X-XM1)
          YR1 = (Y-YM1)
          
C         rotate
          
          XR = DCRCOS*XR1+DCRSIN*YR1
          YR = -DCRSIN*XR1+DCRCOS*YR1

C         Get the distance
      
          DIST=( (XR**2)/((DCRADX+DCRRLX)**2) ) +
     *         ( (YR**2)/((DCRADY+DCRRLY)**2) )

        END IF

        IF (DIST.LE.1D0) THEN
          ISFBDY = TNBC()+1
        END IF

      ELSE IF (ICIRTP.EQ.3) THEN

C       ---------------------------------------------------------------
C       fictitious boundary component: rectangle with midpoint at 
C       DCXPOS/DCYPOS, width=2*DCRADX, height=2*DCRADY
C       ---------------------------------------------------------------

        XM1=DCXPOS
        YM1=DCYPOS

        DIST=MAX(ABS( (X-XM1) / (DCRADX+DCRRLX) ),
     *           ABS( (Y-YM1) / (DCRADY+DCRRLY) ) )

        IF (DIST.LE.1D0) THEN
          ISFBDY = TNBC()+1
        END IF

      ELSE IF (ICIRTP.GE.50) THEN

C       Multiple geometries.
C       Here we have to take into account whether to test for a special
C       geometry or generally for all geometries!

        IF (IFBC.NE.0) THEN

C         Get the handle of the integer- and double-prec. parameter
C         block from IGEOMS(1,IFBC-TNBC()) and
C         IGEOMS(2,IFBC-TNBC()), resp.; L0 gets
C         the type of the object.
C
C          L0 = IGEOMS(1,IFBC-TNBC())
C          L1 = IGEOMS(2,IFBC-TNBC())
C          L2 = IGEOMS(3,IFBC-TNBC())

          L0 = IGEOM(OIGEOMS+3*(IFBC-TNBC()-1))
          L1 = IGEOM(OIGEOMS+3*(IFBC-TNBC()-1)+1)
          L2 = IGEOM(OIGEOMS+3*(IFBC-TNBC()-1)+2)

C         Do we have a composed geometry?

          IF (L0.EQ.11) THEN

            IF (GINCMP (KWORK(L(L1)), 
     *                  DWORK(L(L2)), 
     *                  GINWRP, X, Y).GT.0) THEN
              ISFBDY = IFBC
            END IF
            
          ELSE
          
C           Not a composed geometry, but a group.
C           Not implemented. Could be a loop through all objects
C           and figure out the status for each object that belongs
C           to our main object. But that would be dear!

            ISFBDY = 0
          
          END IF ! L0=11
          
        ELSE
        
C         General test
          
          DO I=1,IGEOM(OIGEOCN)
          
C           Get the handle of the integer- and double-prec. parameter
C           block from IGEOMS(1,I) and IGEOMS(2,I), resp.; L0 gets
C           the type of the object.
C
C            L0 = IGEOMS(1,I)
C            L1 = IGEOMS(2,I)
C            L2 = IGEOMS(3,I)
C
            L0 = IGEOM(OIGEOMS+3*(I-1))
            L1 = IGEOM(OIGEOMS+3*(I-1)+1)
            L2 = IGEOM(OIGEOMS+3*(I-1)+2)

C           Do we have a composed geometry?

            IF (L0.EQ.11) THEN

              IF (GINCMP (KWORK(L(L1)), 
     *                    DWORK(L(L2)), 
     *                    GINWRP, X, Y).GT.0) THEN
     
C               Stop at the first found geometry containing the point

                ISFBDY = TNBC()+I
                GOTO 10
              END IF
              
            ELSE
            
C             Not a composed geometry, but a group.
C             Not implemented. Could be a loop through all objects
C             and figure out the status for each object that belongs
C             to our main object. But that would be dear!

            END IF
              
          END DO
        
10       CONTINUE        
        
        END IF

      ELSE IF (ICIRTP.GE.4) THEN

C       User-defined complex geometry.
C       nothing implemented here.

      END IF

C     Every other value for ICIRTP disables the fictitious boundary 
C     management, i.e. there is no fict. boundary component at all.

      END
      
************************************************************************
* FBDVOL - Get the (analytical) volume of the fictitious
*          boundary component
*
* This routine computes the reference volume of the fictitious boundary
* component IFBC. This is used as the reference volume in comparisons.
*
* In:
*  IFBC   - Number of fictitious boundary component to calculate
*           the volume for.
*           0 = calculate summed volume of all fictitious boundary
*               components in the domain
*           NBCT+1..NBCT+NFBDYC = calculate volume only for
*               fictitious boundary component IFBC
*  IGEOM - array [1..*] of integer
*  DGEOM - array [1..*] of double precision
*          Parameter blocks describing the geometry
*
* Return: calculated volume or -1D0, if the calculation is not supported
************************************************************************

      DOUBLE PRECISION FUNCTION FBDVOL (IFBC,IGEOM,DGEOM)

      IMPLICIT NONE
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)

      INCLUDE 'cbasictria.inc'
      INCLUDE 'sinigeometry.inc'
      
C     parameters

      INTEGER IFBC,IGEOM(*)
      DOUBLE PRECISION DGEOM(*)

C     local variables

      INTEGER ICIRTP
      DOUBLE PRECISION DCXPOS, DCYPOS, DCROT, DCRADX, DCRADY, DCRAD
      DOUBLE PRECISION DCRSIN,DCRCOS, DCRRLX, DCRRLY

C     Get rotation, radius,... from the parameter blocks

      ICIRTP = IGEOM(OICIRTP)
      DCXPOS = DGEOM(ODCXPOS)
      DCYPOS = DGEOM(ODCYPOS)
      DCROT  = DGEOM(ODCROT )
      DCRAD  = DGEOM(ODCRAD )
      DCRADX = DGEOM(ODCRDX)
      DCRADY = DGEOM(ODCRDY)
      DCRSIN = DGEOM(ODCRSIN)
      DCRCOS = DGEOM(ODCRCOS)
      DCRRLX = DGEOM(ODCRDX)
      DCRRLY = DGEOM(ODCRDY)

C     standard return value

      FBDVOL=-1D0

      IF ((IFBC.EQ.0).OR.
     *    ((FBDVOL.GT.0).AND.(IFBC.EQ.1))) THEN

        IF (ICIRTP.EQ.0) THEN
C circle
          FBDVOL = (DCRAD+DCRRLX)**2*PI
        ELSE IF (ICIRTP.EQ.1) THEN
C square
          FBDVOL = (2*(DCRAD+DCRRLX))**2
        ELSE IF (ICIRTP.EQ.2) THEN
C ellipse
          FBDVOL = (DCRADX+DCRRLX)*(DCRADY+DCRRLY)*PI
        ELSE IF (ICIRTP.EQ.3) THEN
C rectangle
          FBDVOL = 2*(DCRADX+DCRRLX)*2*(DCRADY+DCRRLY)
        END IF
        
      END IF
      
      END
      
      
************************************************************************
* FBDGMV - Write analytical fictitious boundaries to GMV file
*
* This routine is called by the framework in the postprocessing
* part. It should write out an analytical description of all
* fictitious boundary parts to the file unit FUNIT in syntax
* of GMV (as far as this is possible).
*
* The routine is allowed to use material number MMAT..MMAT+NFBDYC()
* for the fictitious bondary components.
*
* The caller of this routine is responsible for the administration
* of the sections in the GMV file. This routine only has to write out
* the content for the POLYGON's section to approximate the boundary
* components. The header/footer of this section is written to the
* file by the caller.
*
* In:
*  MUNIT  - File unit where to write the GMV output to.
*  MMAT   - The first material number to use.
*  IGEOM - array [1..*] of integer
*  DGEOM - array [1..*] of double precision
*          Parameter blocks describing the geometry
************************************************************************

      SUBROUTINE FBDGMV (MUNIT,MMAT,IGEOM,DGEOM)

      IMPLICIT NONE
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)

      INCLUDE 'cbasictria.inc'
      INCLUDE 'cmem.inc'

      INCLUDE 'sinigeometry.inc'
      
      INCLUDE 'sgeometries.inc'
            
C     parameters

      INTEGER MUNIT,MMAT,IGEOM(*)
      DOUBLE PRECISION DGEOM(*)
      
C     externals

      EXTERNAL GGMWRP
      
      INTEGER NFBDYC
      EXTERNAL NFBDYC
      
      INTEGER GCRELL,GCRCIR
      EXTERNAL GCRELL,GCRCIR
      
C     local variables

      INTEGER I,NL,L1,L2,L0
      DOUBLE PRECISION DEG
      DOUBLE PRECISION GEOSIM(SZGELL)
      DOUBLE PRECISION CORSYS(SZCORS)
      
      DOUBLE PRECISION GEOTIM
      INTEGER ICIRTP,INSGEO
      DOUBLE PRECISION DCXPOS, DCYPOS, DCROT, DCRADX, DCRADY, DCRAD
      DOUBLE PRECISION DCRSIN,DCRCOS

C     Get rotation, radius,... from the parameter blocks

      ICIRTP = IGEOM(OICIRTP)
      DCXPOS = DGEOM(ODCXPOS)
      DCYPOS = DGEOM(ODCYPOS)
      DCROT  = DGEOM(ODCROT )
      DCRAD  = DGEOM(ODCRAD )
      DCRADX = DGEOM(ODCRDX)
      DCRADY = DGEOM(ODCRDY)
      DCRSIN = DGEOM(ODCRSIN)
      DCRCOS = DGEOM(ODCRCOS)
      
      INSGEO = IGEOM(OINSGEO)

C     Get simulational time, i.e. time stamp how the geometry should
C     be acting
      
      GEOTIM = DGEOM(OGEOTIM)

      IF (NFBDYC(IGEOM,DGEOM).EQ.0) RETURN

C     GMV Header

      WRITE (MUNIT,'(A)') 'polygons'

C     One circle:

      IF (ICIRTP.EQ.0) THEN
        
C       As we only have one boundary component, this routine is only 
C       called once...
C
C       Create a geometry structure for a circle and write out 
      
        CALL GCRSYS (CORSYS, 0D0, 0D0, 0D0, 1D0)

C       Benchmark case in nonsteady simulation.
C       Fixed point in steady simulation

        IF (INSGEO.EQ.1) THEN
          I = GCRCIR (GEOSIM, DCXPOS+0.25*SIN(2D0*PI*0.25*GEOTIM), 
     *                DCYPOS, DCROT, 1.0D0, 0, DCRAD)
        ELSE
          I = GCRCIR (GEOSIM, DCXPOS, DCYPOS, DCROT, 1.0D0, 0, DCRAD)
        END IF
        CALL GGMCIR (GEOSIM,CORSYS, MUNIT, MMAT)
      
C       Emulate the circle by a number of lines...
C       Number of circle segments is dependent of radius
C
C        NL = MIN(10000,MAX(1,NINT(DCRAD*2*PI*1000)))
C
C       Material number MMAT=moving boundary, NL+1 nodes
C
C        WRITE(MUNIT,*) MMAT,NL+1
C
C       X-coordinates
C
C        DO I=0,NL
C          DEG = DBLE(I)*2D0*PI/DBLE(NL)
C          WRITE(MUNIT,*) DCXPOS+0.5D0*SIN(0.8D0*GEOTIM)+DCRAD*SIN(DEG)
C
C         Benchmark case in nonsteady simulation.
C         Fixed point in steady simulation
C          
C          IF (INSGEO.EQ.1) THEN
C            WRITE(MUNIT,*) DCXPOS+DCRAD*SIN(DEG)
C     *                     +0.25*SIN(2D0*PI*0.25*GEOTIM)
C          ELSE
C            WRITE(MUNIT,*) DCXPOS+DCRAD*SIN(DEG)
C          END IF
C        END DO
C
C       Y-coordinates
C
C        DO I=0,NL
C          DEG = DBLE(I)*2D0*PI/DBLE(NL)
C          WRITE(MUNIT,*) DCYPOS+DCRAD*COS(DEG)
C        END DO
C
C       Z-coordinates
C
C        DO I=0,NL
C          WRITE(MUNIT,*) 0D0
C        END DO
      
      ELSE IF (ICIRTP.EQ.2) THEN
      
C       Create a geometry structure for an ellipse and write out 
C       the ellipse
      
        CALL GCRSYS (CORSYS, 0D0, 0D0, 0D0, 1D0)
        I = GCRELL (GEOSIM, DCXPOS, DCYPOS, DCROT, 1.0D0, 0, 
     *              DCRADX, DCRADY)
        CALL GGMELL (GEOSIM,CORSYS, MUNIT, MMAT)
      
      ELSE IF (ICIRTP.GE.50) THEN

C       Multiple geometries.
C       Write out all sub-objects.
        
        DO I=1,IGEOM(OIGEOCN)
                                         
C         Get the handle of the integer- and double-prec. parameter
C         block from IGEOMS(1,I) and IGEOMS(2,I), resp.:
C
C          L0 = IGEOMS(1,I)
C          L1 = IGEOMS(2,I)
C          L2 = IGEOMS(3,I)

          L0 = IGEOM(OIGEOMS+3*(I-1))
          L1 = IGEOM(OIGEOMS+3*(I-1)+1)
          L2 = IGEOM(OIGEOMS+3*(I-1)+2)
          
C         Do we have a composed geometry?

          IF ((L0.EQ.11).OR.(L0.EQ.12)) THEN

            CALL GGMCMP (KWORK(L(L1)), 
     *                   DWORK(L(L2)), 
     *                   GGMWRP,MUNIT,MMAT+I) 
     
          ELSE IF (L0.EQ.0) THEN

C           Not a composed geometry but a standard simple-type 
C           geometry! 
C           Create an identity coordinate system for the current
C           object

            CALL GCRSYS (CORSYS, 0D0, 0D0, 0D0, 1D0)

C           The geometry data is stored in IGEODL at the starting
C           address given by L1!

            L2 = IGEOM(OIGEODL)
            IF (L2.EQ.0) THEN
              WRITE (*,*) 'Geometry data array IGEODL undefined!'
              STOP
            END IF
            L2 = L(L2)+L1-1

            CALL GGMWRP(DWORK(L2),CORSYS, MUNIT, MMAT+I)
          
          ELSE
          
C           Not a composed geometry, but a group.
C           Not implemented. Could be a loop through all objects
C           and figure out the status for each object that belongs
C           to our main object. But that would be dear!

          END IF
     
        END DO
        
      ELSE IF (ICIRTP.GE.4) THEN

C       User-defined geometry

      END IF
      
C     GMV Footer
      
      WRITE (MUNIT,'(A)') 'endpoly'
      
      END

************************************************************************
* Prescribe Dirichlet boundary conditions on fictitious boundary
*
* This routine is called by the framework when implementing Dirichlet 
* boundary conditions in a point of the fictitious boundary domain.
* It can be used to prescribe e.g. tangential velocity initiated
* by rotation.
* ITYP describes the information that has to be returned in the
* return value.
*
* In:
*  ITYP   - 1=X-velocity, 2=Y-velocity, 3=pressure
*  X,Y    - the point
*  TIMENS - Current time in instationary Navier-Stokes calculation
*           =0 for stationary simulation.
*  RE     - Reynolds number; from the DAT file
*  IPARAM - array [1..*] of integer 
*  DPARAM - array [1..*] of integer 
*           TIntAssembly/TDoubleAssembly assembly structures; gives
*           additional information about the discretization.
*           This is passed to user defined callback routines so that 
*           they can access more detailed information about the
*           problem. Not used in this routine.
*  IGEOM  - array [1..*] of integer 
*  DGEOM  - array [1..*] of double 
*           Integer- and double-precision parameter blocks with
*           geometry information. Passed to boundary
*           routines. Not used in this routine.
*
* "Optional" parameters:
*  INODE  - Number of the node corresponding to (X,Y). Can be 0.
*           If INODE=0, this function must create the necessary 
*              information only by the given (X,Y)-coordinates of 
*              the point
*           If INODE>0, the parameter TRIA must also be defined.
*              In this case INODE is a node number in the triangulation
*              TRIA with coordinates (X,Y). Then this routine can also
*              access additional information about the triangulation
*              (like precomputed values) to compute the result.
*  TRIA   - array [1..SZTRIA] of integer
*           If INODE>0, TRIA defines the triangulation INODE refers to.
*           If INODE=0, this parameter can be a dummy parameter.
* 
* Out:
*  Return value = the desired information
************************************************************************

      DOUBLE PRECISION FUNCTION FBINDT (ITYP,X,Y,TIMENS,RE,INODE,TRIA,
     *                                  IPARAM,DPARAM,IGEOM,DGEOM)
      
      IMPLICIT NONE
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)
      
      INCLUDE 'sinigeometry.inc'
      
      INCLUDE 'sassembly.inc'
      
C parameters

      INTEGER ITYP
      DOUBLE PRECISION X,Y
      DOUBLE PRECISION DPARAM(*),TIMENS,RE,DGEOM(*)
      LOGICAL BFIRST
      INTEGER INODE,TRIA(*),IPARAM(*),IGEOM(*)
      
C local variables

      DOUBLE PRECISION XM1,YM1,XR1,YR1,XN1,YN1,SPD
      
      DOUBLE PRECISION GEOTIM
      INTEGER ICIRTP,INSGEO
      DOUBLE PRECISION DCXPOS, DCYPOS, DCROT, DCRADX, DCRADY, DCRAD
      DOUBLE PRECISION DCRSIN,DCRCOS

C     Get rotation, radius,... from the parameter blocks

      ICIRTP = IGEOM(OICIRTP)
      DCXPOS = DGEOM(ODCXPOS)
      DCYPOS = DGEOM(ODCYPOS)
      DCROT  = DGEOM(ODCROT )
      DCRAD  = DGEOM(ODCRAD )
      DCRADX = DGEOM(ODCRDX)
      DCRADY = DGEOM(ODCRDY)
      DCRSIN = DGEOM(ODCRSIN)
      DCRCOS = DGEOM(ODCRCOS)
      
      INSGEO = IGEOM(OINSGEO)

C     Get simulational time, i.e. time stamp how the geometry should
C     be acting
      
      GEOTIM = DGEOM(OGEOTIM)

C Standard return value: no velocity.
      
      FBINDT = 0D0
      
      IF (ITYP.EQ.3) FBINDT=1D0
      return

C     Possibility:
C     Pos = 1.1+sim(timens)
C     Vel = cos(timens)
C      IF (ITYP.EQ.1) FBINDT = 0.5D0*COS(0.8D0*GEOTIM)
      
C     Benchmark formula in nonsteady simulation.
C     Fixed object in steady simulation
      
      IF (INSGEO.EQ.1) THEN
        IF (ITYP.EQ.1) FBINDT = 2D0*PI*0.25*0.25*COS(2D0*PI*0.25*GEOTIM)
      ELSE
        IF (ITYP.EQ.1) FBINDT = 0.0
      END IF
      
      GOTO 99999
      
C In case of a circle we might prescribe the rotation

      IF (ICIRTP.EQ.0) THEN

C fictitious boundary component: circle at DCXPOS/DCYPOS, radius DCRAD

C Midpoint:

10      XM1 = DCXPOS
        YM1 = DCYPOS
        
C difference vector

        XR1 = X-XM1
        YR1 = Y-YM1
        
C create the normal vector to that; this is at the same time the tangential
C vector on the circle

        XN1 = YR1
        YN1 = -XR1
        
C normalise it, such that points on the boundary have velocity=STD

        SPD = 0.05D0

        XN1 = SPD * XN1/DCRAD**2
        YN1 = SPD * YN1/DCRAD**2
        
        IF (ITYP.EQ.1) FBINDT=XN1
        IF (ITYP.EQ.2) FBINDT=YN1

      END IF
      
99999 END
      
************************************************************************
* FBDNML - Get normal vector of fictitious boundary
*
* This routine is called by the framework to compute analytical normal
* vectors of the fictitious boundary. To a given point (X,Y) -
* which may be inside the fictitious boundary or not, depending
* on the approximation - this routine has to calculate a normalised
* vector (XN,YN) perpendicular to the fictitious boundary domain IFBC.
*
* In:
*  X,Y    - Point where to calculate the normal vector
*  IFBC   - Number of fictitious boundary component to calculate
*           the volume for.
*           0 = automatically determine a fictitious boundary
*               component which boundary is taken for the calculation
*               of the normal vector (i.e. by analytically analysing
*               the distance)
*           NBCT+1..NBCT+NFBDYC = Number of fictitious boundary
*               component which is used for the calculation of the
*               normal vector.
*  IGEOM - array [1..*] of integer
*  DGEOM - array [1..*] of double precision
*          Parameter blocks describing the geometry
*
* Out:
*  XN,YN  - normal vector of the fictitious boundary component IFBC
*           in the point (X,Y); must be normalised to length=1D0 !
* 
* If calculation of (XN,YN) is possible, this routine should set
* (XN,YN)=(0D0,0D0) to indicate this to the framework.
************************************************************************

      SUBROUTINE FBDNML (X,Y,IFBC,XN,YN,IGEOM,DGEOM)

      IMPLICIT NONE
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)

      INCLUDE 'cbasictria.inc'

      INCLUDE 'sinigeometry.inc'
      
C     parameters

      DOUBLE PRECISION X,Y,XN,YN,A,DGEOM(*)
      INTEGER IFBC,IGEOM(*)

C     local variables

      INTEGER ICIRTP
      DOUBLE PRECISION DCXPOS, DCYPOS, DCROT, DCRADX, DCRADY, DCRAD
      DOUBLE PRECISION DCRSIN,DCRCOS

      XN = 0D0
      YN = 0D0

C     Get rotation, radius,... from the parameter blocks

      ICIRTP = IGEOM(OICIRTP)
      DCXPOS = DGEOM(ODCXPOS)
      DCYPOS = DGEOM(ODCYPOS)
      DCROT  = DGEOM(ODCROT )
      DCRAD  = DGEOM(ODCRAD )
      DCRADX = DGEOM(ODCRDX)
      DCRADY = DGEOM(ODCRDY)
      DCRSIN = DGEOM(ODCRSIN)
      DCRCOS = DGEOM(ODCRCOS)

C For the moment we only have one boundary component, so IFBC is unused.

C Circle case:

      IF (ICIRTP.EQ.0) THEN
C Calculate the normal by the taking the difference vector
C to the midpoint

        XN = X-DCXPOS
        YN = Y-DCYPOS
      END IF

C normalise the vector

      IF ((XN.NE.0D0).OR.(YN.NE.0D0)) THEN
        A = 1D0/DSQRT(XN*XN+YN*YN)
        XN=A*XN
        YN=A*YN
      END IF
      
      END
      
************************************************************************
* FBDDST - Get distance of point to fictitious boundary interface  
*
* This routine is called by the framework to compute the minimum 
* distance of a point (X,Y) to the interface of a fictitious
* boundary component. The parameter IFBC determines the number
* of the boundary component to determine the distance to.
* When called with IFBC=0, the routine should calculate the minimum
* distance to any of the fictitious boundary interfaces.
* When called with IFBC>0 (i.e. IFBC of NBCT+1..NBCT+NFBDYC), 
* the routine should calculate the distance to that specific
* fictitious boundary component.
*
* The distance is returned in the parameter DIST. For points inside 
* of a fictitious boundary component, this value is <0 whereas 
* points outside the fictitious boundary are described by a value >0.
*
* If the distance was successfully calculated, the routine must
* return a value > 0. If it's not possible to determine the 
* distance (either if it's not implemented or not possible 
* due to complexity) the routine must return a value <= 0 to indicate
* this. The framework will then take care of the fact that determining
* the distance to this fictitious boundary component is not
* possible.
*
* In:
*  X,Y    - Point to determine the distance to the fictitious
*           boundary interface
*  IFBC   - Number of fictitious boundary component to calculate
*           the distance to.
*           0 = Compute the minimum distance to all fictitious
*               boundaty interfaces
*           NBCT+1..NBCT+NFBDYC = Number of fictitious boundary
*               component where to calculate the distance to.
*  IGEOM  - array [1..*] of integer
*  DGEOM  - array [1..*] of double precision
*           Parameter blocks describing the geometry
*
* Out:
*  DIST   - Distance of (X,Y) to the interface of a/the fictitious
*           boundary
* 
* Return:
*  >  0   - If the distance was successfully computed
*  <= 0   - If the computation is not possible
* 
************************************************************************

      INTEGER FUNCTION FBDDST (X,Y,IFBC, DIST,IGEOM,DGEOM)

      IMPLICIT NONE

      INCLUDE 'sinigeometry.inc'
      
      DOUBLE PRECISION X,Y,DIST,DGEOM(*)
      INTEGER IFBC,IGEOM(*)
      
      INTEGER ICIRTP
      DOUBLE PRECISION DCXPOS, DCYPOS, DCROT, DCRADX, DCRADY, DCRAD
      DOUBLE PRECISION DCRSIN,DCRCOS

C     Get rotation, radius,... from the parameter blocks

      ICIRTP = IGEOM(OICIRTP)
      DCXPOS = DGEOM(ODCXPOS)
      DCYPOS = DGEOM(ODCYPOS)
      DCROT  = DGEOM(ODCROT )
      DCRAD  = DGEOM(ODCRAD )
      DCRADX = DGEOM(ODCRDX)
      DCRADY = DGEOM(ODCRDY)
      DCRSIN = DGEOM(ODCRSIN)
      DCRCOS = DGEOM(ODCRCOS)

C     standard return value

      FBDDST = 0

C We are only able to handle the circle case here:

      IF (ICIRTP.EQ.0) THEN
        DIST = DSQRT((X-DCXPOS)**2+(Y-DCYPOS)**2) - DCRAD
        FBDDST = 1
      END IF

      END

************************************************************************
* FBDMON - Get monitor function value of fictitious boundary domain
*
* This routine is called by the framework to compute a value for the
* monitor function in the point (X,Y) when performing grid adation. 
* It is called once with IFBC=0 to compute a general value for the
* monitor function and once for every fictitious boundary domain
* IFBC=NBCT+1..NBCT+NFBDYC. The monitor function itself is computed
* by taking the minimum value of all returned values in that point.
*
* This function has to return a positive real value. If the
* value is not in the interval [0,1], it will be restricted
* by the framework to this. (To be more exact, the value is
* restricted to the interval [EPS0>0,1] where EPS0 is given in the
* .DAT file of the grid adaption.)
* The variable EPS0 is given to inform this routine about the
* lower bound, the routine might ignore this value.
* The returned value should be ~ 0 if the point is near the interface
* of the fictitious boundary domain and ~ 1 if the point is far
* away from it. The grid adaption routine will then concentrate
* the distribution of points to the areas where this routine
* returned low values.
*
* In:
*  X,Y    - Point where to calculate the normal vector
*  IFBC   - Number of fictitious boundary component to calculate
*           the volume for.
*           0 = General computation of the monitor function value;
*               the routine can decide how to calculate
*           NBCT+1..NBCT+NFBDYC = Number of fictitious boundary
*               component whose weight should be calculated
*           If one case is not handled by this routine, 1D0 should
*           be returned.
*  EPS0   - minimum monitor function values; values below that
*           will be truncated.
*  IGEOM  - array [1..*] of integer
*  DGEOM  - array [1..*] of double precision
*           Parameter blocks describing the geometry
*
* Return:
*  Value of monitor function, either for all boundary components (if
*  IFBC=0) or for the fictitious boundary component IFBC.
* 
************************************************************************

      DOUBLE PRECISION FUNCTION FBDMON (X,Y,IFBC,EPS0,IGEOM,DGEOM)

      IMPLICIT NONE

      INCLUDE 'cmem.inc'

      INCLUDE 'sinigeometry.inc'
      INCLUDE 'sgeometries.inc'

C parameters

      DOUBLE PRECISION X,Y,EPS0,DGEOM(*)
      INTEGER IFBC,IGEOM(*)
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)

C local variables
      
      DOUBLE PRECISION DIST,WEIG,TMP,TMP2
      DOUBLE PRECISION XR,YR,XR1,YR1,XM1,YM1
      INTEGER I,L1,L2,L0
      DOUBLE PRECISION CORSYS(SZCORS)
      DOUBLE PRECISION REC (SZGREC)
      
      INTEGER GCRREC,GDSREC
      EXTERNAL GCRREC,GDSREC

C externals

      INTEGER GDSCMP,GDSWRP
      EXTERNAL GDSCMP,GDSWRP

      DOUBLE PRECISION GEOTIM
      INTEGER ICIRTP,IFBMNT
      DOUBLE PRECISION DCXPOS, DCYPOS, DCROT, DCRADX, DCRADY, DCRAD
      DOUBLE PRECISION DCRSIN,DCRCOS, DCRRLX, DCRRLY, DMNISF, DMNOSF

C     Get rotation, radius,... from the parameter blocks

      ICIRTP = IGEOM(OICIRTP)
      DCXPOS = DGEOM(ODCXPOS)
      DCYPOS = DGEOM(ODCYPOS)
      DCROT  = DGEOM(ODCROT )
      DCRAD  = DGEOM(ODCRAD )
      DCRADX = DGEOM(ODCRDX)
      DCRADY = DGEOM(ODCRDY)
      DCRSIN = DGEOM(ODCRSIN)
      DCRCOS = DGEOM(ODCRCOS)
      DCRRLX = DGEOM(ODCRLX)
      DCRRLY = DGEOM(ODCRLY)

      DMNISF = DGEOM(ODMNISF)
      DMNOSF = DGEOM(ODMNOSF)
      
      IFBMNT = IGEOM(OIFBMNT)

C     Get simulational time, i.e. time stamp how the geometry should
C     be acting
      
      GEOTIM = DGEOM(OGEOTIM)

C standard return value

      FBDMON = 1D0

      IF (IFBMNT.EQ.-1) THEN
      
C       Special monitor function: 
C       Equal distribution of cell size through the whole domain.
C       FBDMON is set to 1D0 (above), so return immediately.

        RETURN
      
      ELSE IF (IFBMNT.EQ.1) THEN

C       Special monitor function: polynomial of 4th order.
C       Concentrate everything leftmost and around the object.     
        
        FBDMON = EPS0 + x**2 * (x-DCXPOS)**2
        RETURN
        
      ELSE IF (IFBMNT.EQ.2) THEN
      
C       Special monitor function: concentrate everything in the inflow

        FBDMON = (X-DCXPOS-2*DCRAD)*1D0/(DCXPOS-2*DCRAD)
        RETURN
      
      ELSE IF (IFBMNT.EQ.3) THEN
        
C       Special monitor function: concentrate everything to the left
        
        FBDMON = 2*EPS0+X
        RETURN
        
      ELSE IF (IFBMNT.EQ.4) THEN

C       Modified monitor function: additional refinement on the inflow

        FBDMON = MAX(2*EPS0+(1D0-2*EPS0)*(0.2*X/(DCXPOS+2*DCRAD)),
     *               2*EPS0+0.1D0+0.5D0*(X-DCXPOS-2*DCRAD))
        
C       no RETURN here since this is a modification of the monitor
C       function defined below.

      ELSE IF (IFBMNT.EQ.5) THEN
       
C       Special monitor function: SIN-Wave

        FBDMON = EPS0 + 
     *          (1D0-EPS0)*(0.5D0*COS(PI+X*2D0*PI/DCXPOS)+0.5D0)**2
        IF (X.GT.3*0.5D0*DCXPOS) FBDMON = 1D0
        RETURN
      
      ELSE IF (IFBMNT.EQ.6) THEN

C       Modified monitor function: even more additional refinement 
C       on the inflow

        FBDMON = MAX(EPS0+(1D0-EPS0)*(0.2*X/(DCXPOS+2*DCRAD)),
     *               EPS0+0.1D0+0.5D0*(X-DCXPOS-2*DCRAD))
        
C       no RETURN here since this is a modification of the monitor
C       function defined below.

      ELSE IF (IFBMNT.EQ.7) THEN

C       Modified monitor function: 
C       Small cells from the beginning up to shortly past the cylinder,
C       then increasing cell size

        IF (X.LT.DCXPOS+2*DCRAD) THEN
          FBDMON = 5D0*EPS0
        ELSE 
          TMP = (X-DCXPOS-2*DCRAD) / (2.2D0-(DCXPOS+2*DCRAD))
          FBDMON = 5D0*EPS0 * (1D0-TMP) + 1D0*TMP
        END IF
        
C       no RETURN here since this is a modification of the monitor
C       function defined below.

      ELSE IF (IFBMNT.EQ.8) THEN

C       Modified monitor function: 
C       Small cells from the beginning up to shortly past the cylinder,
C       then increasing cell size.
C       Parabolic profile.

        IF (X.LT.DCXPOS+2*DCRAD) THEN
          FBDMON = 5D0*EPS0
        ELSE 
          TMP = (X-DCXPOS-2*DCRAD) / (2.2D0-(DCXPOS+2*DCRAD))
          FBDMON = 5D0*EPS0 * (1D0-TMP) + 1D0*TMP
        END IF
        
        TMP = Y/0.41D0
        FBDMON = ( 1D0+1.5D0-1.5D0*TMP*(1-TMP)/0.25D0 ) * FBDMON
        
C       no RETURN here since this is a modification of the monitor
C       function defined below.

      END IF
      
C We only have one global object, so return immediately if we are
C called for a special object
      
      IF (IFBC.NE.0) RETURN

C Calculate a monitor function value

C Circle case:

      IF (ICIRTP.EQ.0) THEN
      
        XM1 = DCXPOS
        YM1 = DCYPOS
C     Possibility:
C     Pos = 1.1+sim(timens)
C     Vel = cos(timens)

C        XM1 = DCXPOS+0.5D0*SIN(0.8D0*GEOTIM)

C       Benchmark formula
      
        XM1 = DCXPOS+0.25*SIN(2D0*PI*0.25*GEOTIM)

        DIST=SQRT((X-XM1)**2+(Y-YM1)**2)
C        DIST=ABS(X-XM1)
        
        IF (DIST.LT.DCRAD) THEN
C We use a slightly higher multiplicative factor in the inner to
C "press" some of the inner points to the boundary. But caution:
C A too large factor results in very anisotropic cells in the
C boundary!
          FBDMON = MIN(FBDMON,ABS(DIST-DCRAD)*DGEOM(ODMNISF))
        ELSE
          FBDMON = MIN(FBDMON,ABS(DIST-DCRAD)*DGEOM(ODMNOSF))
        END IF
        
C        FBDMON = EPS0 + 19.42793828 * x**2 * (x-1.1D0)**2

      END IF
      
C Ellipse case:

      IF (ICIRTP.EQ.2) THEN
      
C Shift to 0
        XR1 = (X-DCXPOS)
        YR1 = (Y-DCYPOS)
C rotate
        XR = DCRCOS*XR1+DCRSIN*YR1
        YR = -DCRSIN*XR1+DCRCOS*YR1
      
        DIST = ( (XR**2)/((DCRADX+DCRRLX)**2) ) +
     *         ( (YR**2)/((DCRADY+DCRRLY)**2) )
        
C Calculate the "distance" to the interface of the ellipse.
C For that purpose scale the Y-coordinate up/down such that
C the ellipse deforms to a circle (in another coordinate system)
C with radius DCXRAD.
C Then we can measure the distance as in the circle case.
        
        WEIG = DSQRT(XR**2 + (YR*(DCRADX+DCRRLX)/(DCRADY+DCRRLY))**2)
        
        IF (DIST.LT.1D0) THEN
C We use a slightly higher multiplicative factor in the inner to
C "press" some of the inner points to the boundary. But caution:
C A too large factor results in very anisotropic cells in the
C boundary!
          FBDMON = MIN(FBDMON,ABS(WEIG-(DCRADX+DCRRLX))*DGEOM(ODMNISF))
        ELSE
          FBDMON = MIN(FBDMON,ABS(WEIG-(DCRADX+DCRRLX))*DGEOM(ODMNOSF))
        END IF
      END IF

      IF (ICIRTP.GE.50) THEN

C Composed geometries. Loop through all geometries and
C determine the minimum monitor function value.

        TMP = 1D0

        DO I=1,IGEOM(OIGEOCN)

C         Get the handle of the integer- and double-prec. parameter
C         block from IGEOMS(1,IFBC-TNBC()) and
C         IGEOMS(2,IFBC-TNBC()), resp.; L0 gets
C         the type of the object.
C
C          L0 = IGEOMS(1,IFBC-TNBC())
C          L1 = IGEOMS(2,IFBC-TNBC())
C          L2 = IGEOMS(3,IFBC-TNBC())

          L0 = IGEOM(OIGEOMS+3*(I-1))
          L1 = IGEOM(OIGEOMS+3*(I-1)+1)
          L2 = IGEOM(OIGEOMS+3*(I-1)+2)

C         Do we have a composed geometry?

          IF (L0.EQ.11) THEN

            IF (GDSCMP (KWORK(L(L1)), 
     *                  DWORK(L(L2)),
     *                  GDSWRP, X, Y,DIST).GT.0) THEN
              TMP = MIN(TMP,ABS(DIST))
            END IF
            
          ELSE IF (L0.EQ.0) THEN

C           Not a composed geometry but a standard simple-type 
C           geometry! 
C           Create an identity coordinate system for the current
C           object

            CALL GCRSYS (CORSYS, 0D0, 0D0, 0D0, 1D0)

C           The geometry data is stored in IGEODL at the starting
C           address given by L1!

            L2 = IGEOM(OIGEODL)
            IF (L2.EQ.0) THEN
              WRITE (*,*) 'Geometry data array IGEODL undefined!'
              STOP
            END IF
            L2 = L(L2)+L1-1

            IF (GDSWRP(DWORK(L2),CORSYS, X,Y, DIST).GE.0) THEN
            
C             Correct distance depending on the sign with the help of the
C             multiplication factors. This enables us to absorb improper
C             geometry scaling manually.

              IF (DIST.LT.0D0) THEN
                DIST=DIST*DMNISF
              ELSE
                DIST=DIST*DMNOSF
              END IF
              
C             Save the distance in TMP

              DIST = MAX(EPS0,ABS(DIST))              
              IF (I.EQ.1) THEN
                TMP = ABS(DIST)
              ELSE
                IF (IGEOM(OIOVMNT).EQ.0) THEN
                  TMP = MIN(TMP,ABS(DIST))
                ELSE
                  TMP = TMP*ABS(DIST)
                END IF
              END IF
            END IF
          
          ELSE
          
C           Not a composed geometry, but a group.
C           Not implemented. Could be a loop through all objects
C           and figure out the status for each object that belongs
C           to our main object. But that would be dear!
          
          END IF
          
        END DO

        FBDMON = MAX(EPS0,TMP)

      ELSE IF (ICIRTP.GE.4) THEN

C User-defined complex geometry

      END IF

C     Special refinement in front of the circle

c      FBDMON = MAX(FBDMON,2D0*EPS0)
c
c      CALL GCRSYS (CORSYS, 0D0, 0D0, 0D0, 1D0)
c      
c      I = GCRREC (REC, 0.175D0, 0.155D0, 0D0, 1D0, 0, 
c     *            0.025D0, 0.025D0)
c      IF (GDSREC(REC,CORSYS,X,Y,TMP).GT.0) 
c     *  FBDMON = MAX(EPS0,MIN(4D0*TMP,FBDMON))
c
c      I = GCRREC (REC, 0.175D0, 0.245D0, 0D0, 1D0, 0, 
c     *            0.025D0, 0.025D0)
c      IF (GDSREC(REC,CORSYS,X,Y,TMP).GT.0) 
c     *  FBDMON = MAX(EPS0,MIN(4D0*TMP,FBDMON))
c
c      I = GCRREC (REC, 0.275D0, 0.251D0, 0D0, 1D0, 0, 
c     *            0.15D0, 0.03D0)
c      IF (GDSREC(REC,CORSYS,X,Y,TMP).GT.0) 
c     *  FBDMON = MIN(MAX(3D0*EPS0,4D0*TMP),FBDMON)
c
c      I = GCRREC (REC, 0.275D0, 0.149D0, 0D0, 1D0, 0, 
c     *            0.15D0, 0.03D0)
c      IF (GDSREC(REC,CORSYS,X,Y,TMP).GT.0) 
c     *  FBDMON = MIN(MAX(3D0*EPS0,4D0*TMP),FBDMON)

c       CALL LRSCLE (X,0.7D0,2.2D0,0.05D0,1D0,TMP)
c       TMP = MAX(0.05D0,TMP)
C       FBDMON = MIN(TMP,FBDMON)
c       FBDMON = MIN(1D0,TMP)
       
C Better resolution of drag/lift: by taking more points to the 
C inflow profile:
C      FBDMON = MIN(FBDMON,0.02*(1-X/1.1)+1*(X/1.1))
C      FBDMON = MIN(FBDMON, (X-0.7)**2+(Y-0.205)**2+EPS0)
C      FBDMON = MIN(FBDMON,(2.2-x)/2.2+0.1)
C      FBDMON = MIN(FBDMON,EPS0*6D0*(1-X/1.1)+1*(X/1.1))
C      FBDMON = MIN(FBDMON,x/2.2+0.3)
      
C Monitor function is between 0 and 1. Scale it to be
C between EPS0 and 1

      IF (IGEOM(OISCMNT).NE.0) THEN
        FBDMON = EPS0 + (1D0-EPS0)*FBDMON
      END IF

      END
      
**********************************************************************
* Cubature point projection method, analytical version
*
* This subroutine must project cubature points on the reference 
* interval [-1,1] to real cubature point coordinates on a given
* quadrilateral.
*
* The cubature points on the reference interval are given
* in the COMMON-block variable DXI(1..NCUBP,1), with /CUB/.NCUBP the
* number of cubature points. The points (X1,Y1) and (X2,Y2) are
* reconstructed intersection points of the fictitious boundary
* interface with the edges of the quadrilateral IEL. DCOORD
* denote the corner vertices of this element in counterclockwise
* order. If IEL=0 there is no element belonging to the coordinates
* DCOORD, but still DCOORD describes the four corner points of
* a quadrilateral.
*
* The routine must translates the parameter values in DXI(.,.) into
* real coordinates on the interface of the fictitious boundary
* interface, which intersects with the edges of the element in
* the coordinates (X1,Y1) and (X2,Y2). The coordinates of these 
* points have to be stored in 
*              (DCUBP(1,1..NCUBP),DCUBP(2,1..NCUBP)).
*
* DJAC will receive the Jacobian determinant of the mapping between
* the reference interval [-1,1] to the interface segment which
* intersects with the edges in the coordates (X1,Y1) / (X2,Y2).
* This coordinate pair is ordered in that way that the vector
* (X2-X1,Y2-Y1) denotes the tangential and (-(Y2-Y1),X2-X1)
* the normal vector of the fictitious boundary component pointing
* "to the outside".
*
* If the mapping between the interval [-1,1] to the
* interface segment does not have a constant Jacobian determinant,
* the value of DETJ has to be set to 1D0. Furthermore the
* values of the norm of the Jacobian determinant have to be
* included into the weights of the cubature formula, which are
* stored in the DOMEGA-variables of the COMMON block /CUB/ !
*
* Depending on the type of the fictitious boundary component this
* routine can perform an exact reconstruction of the quadrature
* points on the fictitious boundary interface. In the case that
* these points are too complicated to be reconstructed analytically,
* this routine can use the subroutine CBLPRQ, which creates a linear
* reconstruction of the interface to distribute the cubature points
* appropriately.
*
* In:
*  DCOORD - array [1..2,1..4] of double
*           Coordinates of the four corners of the real quadrilateral.
*           DCOORD(1,.) saves the X-, DCOORD(2,.) the Y-coordinates.
*  IEL    - The number of current element or 0, if there is no current
*           element belonging to the coordinates in DCOORD.
*  (X1,Y1),
*  (X2,Y2) - Start- and endpoint of the line segment intersecting
*            the edges of the element. These coordinates are given
*            in "real" coordinates on the "real" element, not on the 
*            reference element.
*  IFBC   - Number of fictitious boundary component currently being
*           treated
*           0 = no explicit handling of a special component.
*               Automatically determine the best approximation of
*               any fictitious boundary component that intersects
*               with the element
*           NBCT+1..NBCT+NFBDYC = Number of fictitious boundary
*               component whose interface is being meant. 
* 
* Out:
*  DJAC    - Norm of the Jacobi determinant of the mapping between 
*            [-1,1] and the interface segment in real coordinates, 
*            or 1D0 if the Jacobian determinant determinant is not
*            constant.
*  DCUBP   - array [1..2,1..NCUBP] of double
*            Coordinates of cubature points in real coordinates
*            on the quadrilateral given by DCOORD.
*  IGEOM   - array [1..*] of integer
*  DGEOM   - array [1..*] of double precision
*            Parameter blocks describing the geometry
*
* If the norm of the Jacobian determinant os not constant, the
* integration weights in /CUB/.DOMEGA have to be modified properly.
**********************************************************************
            
      SUBROUTINE FBCBPR (DCOORD,IEL,IFBC,X1,Y1,X2,Y2,DCUBP,DJAC,
     *                   IGEOM,DGEOM)
      
      IMPLICIT NONE
      
      INCLUDE 'ccub.inc'
      INCLUDE 'sinigeometry.inc'
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)

C parameters
      
      INTEGER IEL,IFBC,IGEOM(*)
      DOUBLE PRECISION X1,Y1,X2,Y2,DJAC,P,DGEOM(*)
      DOUBLE PRECISION DCOORD(2,4),DCUBP(2,NNCUBP)
      
C local variables
      
      DOUBLE PRECISION A1,A2,L1,L2,V1,V2,W1,W2,CURANG
      INTEGER I
      
      INTEGER ICIRTP
      DOUBLE PRECISION DCXPOS, DCYPOS, DCROT, DCRADX, DCRADY, DCRAD
      DOUBLE PRECISION DCRSIN,DCRCOS

C     Get rotation, radius,... from the parameter blocks

      ICIRTP = IGEOM(OICIRTP)
      DCXPOS = DGEOM(ODCXPOS)
      DCYPOS = DGEOM(ODCYPOS)
      DCROT  = DGEOM(ODCROT )
      DCRAD  = DGEOM(ODCRAD )
      DCRADX = DGEOM(ODCRDX)
      DCRADY = DGEOM(ODCRDY)
      DCRSIN = DGEOM(ODCRSIN)
      DCRCOS = DGEOM(ODCRCOS)

C For now we only have the case of exactly one fictitious bondary
C component, so we don't care about IFBC.
      
C In case of a circle we know the interface:

      IF (ICIRTP.EQ.0) THEN

C Ok, quadrilateral IEL intersects with our circle in two points.
C We must know the angle of these two points. Therefore we transform
C them back onto a "reference circle" and calculate the angles
C with trigonometric functions:

        L1 = 1D0/DSQRT((X1-DCXPOS)**2+(Y1-DCYPOS)**2)
        V1 = L1*(X1-DCXPOS)
        W1 = L1*(Y1-DCYPOS)

        L2 = 1D0/DSQRT((X2-DCXPOS)**2+(Y2-DCYPOS)**2)
        V2 = L2*(X2-DCXPOS)
        W2 = L2*(Y2-DCYPOS)
        
C Use the ARCCOS in the upper or lower half circle
        
        A1 = DACOS(V1)
        IF (W1.LT.0D0) A1 = 2*PI - A1 

        A2 = DACOS(V2)
        IF (W2.LT.0D0) A2 = 2*PI - A2 
        
C Make sure angle A2 is "before" A1 in the orientation
        
        IF (A2.GT.A1) A2 = A2 - 2*PI;
        
C Ok, we parametrize for the arc length. [-1,1] has to be mapped
C to the arc between the angles A1 and A2. Loop through the
C cubature points and map them to DCUBP. 

        DO I=1,NCUBP
        
C Calculate the angle of the cubature point
        
          P = 0.5D0+(0.5D0*DXI(I,1))
          CURANG = (1D0-P)*A1 + P*A2
          
C Use that to calculate the cubature point itself

          DCUBP (1,I) = DCXPOS + DCRAD*COS(CURANG)
          DCUBP (2,I) = DCYPOS + DCRAD*SIN(CURANG)

        END DO
        
C Fortran uses (as usual) the arc length parametrization 0..2*PI
C for trigonometrical functions. Therefore the difference between
C A1 and A2 is the arc length of our circle segment - up to
C a constant factor, depending on the radius. This is then the
C Jacobian "determinant" of the mapping. We have to divide it
C by 2 because our reference interval is [-1,1], not [0,1].

        DJAC = 0.5D0*DCRAD*ABS(A2-A1)
        
      ELSE
      
C General handling: linear reconstruction of the interface and
C distribution of the cubature points

        CALL CBLPRQ (DCOORD,IEL,X1,Y1,X2,Y2,DCUBP,DJAC)
        
      END IF
      
      END
      
