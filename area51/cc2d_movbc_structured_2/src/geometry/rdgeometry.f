***********************************************************************
* This file contains a procedure for reading the preferences
* for the geometry routines from a DAT file.
***********************************************************************

***********************************************************************
* Read geometry configuration DAT file
*
* Reads in the file CFNAME from hard disc into geometry
* COMMON blocks. It will open the file CFNAME, read the parameter and
* close the file.
* This will not initialise geometries that have to be loaded from
* disc, i.e. that use memory on the heap. To initialise objects
* that have to be loaded into memory, the routine INIGEO/DONGEO
* has to be called afterwards!
*
* In:
*   MDATA  - IO file handle to use.
*   MSHOW  - Level of output in the initialization phase
*   CFNAME - name of the file
***********************************************************************

      SUBROUTINE RDGEOM (MDATA,MSHOW,CFNAME)
      
      IMPLICIT NONE
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)
      
      INCLUDE 'cout.inc'
      INCLUDE 'cfileout.inc'
      INCLUDE 'cinigeometry.inc'
      
C parameters
      
      CHARACTER CFNAME*(*)
      INTEGER MSHOW,MDATA
      
C local variables

      INTEGER IFMTS,IGPOST
      CHARACTER CSTR*255
      
      WRITE (CSTR,*) 
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      WRITE (CSTR,'(A)') ' Geometry Parameters'
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      WRITE (CSTR,'(A)') '---------------------'
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        
C The file specifies the type of the geometry.
C
C Line  1: "100"            -> Version 1.0.0
C Line  2: -------------------------------------------------------------
C Line  3: "0"/"1"/"2"/"3"  -> 0=circle, 1=square, 
C                              2=ellipse, 3=rectangle,
C                              >50=complex fict. boundary geometry
C                              other values: no fict. boundary
C Line  4: X-Position
C Line  5: Y-Position
C Line  6: Size (Circle: radius; Square: half length of edge, 
C                similar to radius)
C          X-Size (when ellipse is used)
C Line  7: Y-Size (when ellipse is used; otherwise this line contains
C                  a dummy value)
C Line  8: Rotation angle
        
      IFMTS = 1
      CALL  OF0 (MDATA,CFNAME,IFMTS)

C Version number check

      CALL GETINT (MDATA,MSHOW,'Version number of DAT-file (geom.)  '//
     *            'IGPOST = ',IGPOST)
      
      IF (IGPOST.NE.102) THEN
        PRINT *,'Wrong version of geometry DAT-file!'
        CALL EXIT (1)
      END IF

C separation line

      WRITE (CSTR,9001) 
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*)
        
C Type of obstacle

      CALL GETINT (MDATA,MSHOW,'Type of obstacle:                   '//
     *              'ICIRTP = ',ICIRTP)

C position/dimension of the obstacle

      CALL GETDBL (MDATA,MSHOW,
     *              'X-Position of obstacle:             '//
     *              'DCXPOS = ',DCXPOS)
      CALL GETDBL (MDATA,MSHOW,
     *              'Y-Position of obstacle:             '//
     *              'DCYPOS = ',DCYPOS)
      CALL GETDBL (MDATA,MSHOW,
     *              'Radius of obstacle:                 '//
     *              'DCRAD  = ',DCRAD)
      CALL GETDBL (MDATA,MSHOW,
     *              'X-Radius of obstacle:               '//
     *              'DCRADX = ',DCRADX)
      CALL GETDBL (MDATA,MSHOW,
     *              'Y-Radius of obstacle:               '//
     *              'DCRADY = ',DCRADY)
      CALL GETDBL (MDATA,MSHOW,
     *              'Angle of rotation:                  '//
     *              'DCROT  = ',DCROT)
     
      WRITE (CSTR,9001) 
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*)

      CALL GETDBL (MDATA,MSHOW,
     *              'Steepness factor inside of objects: '//
     *              'DMNISF = ',DMNISF)

      CALL GETDBL (MDATA,MSHOW,
     *              'Steepness factor outside of objects:'//
     *              'DMNOSF = ',DMNOSF)

C separation line

        READ (MDATA,*)
        
C Type of problem dependent monitor function

        CALL GETINT (MDATA,MSHOW,
     *              'Type of moitor function:            '//
     *              'IFBMNT = ',IFBMNT)
        CALL GETINT (MDATA,MSHOW,
     *              'Rescale monitor function:           '//
     *              'ISCMNT = ',ISCMNT)
        CALL GETINT (MDATA,MSHOW,
     *              'Monitor function overlay method:    '//
     *              'IOVMNT = ',IOVMNT)
C and close the file, reading finished.
     
      CLOSE (MDATA)
      
C initialise sin/cos(rot)

      DCRSIN = SIN(DCROT*PI/180D0)
      DCRCOS = COS(DCROT*PI/180D0)

9000  FORMAT(79('-'))
9001  FORMAT(60('-'))

      END
      
***********************************************************************
* Initialise geometries
*
* This routine initialized the geometry parameter blocks
* IPARAM/DPARAM with standard values. If IC2PAR=1, the parameters
* from the COMMON block variables are transferred to the geometry
* data blocks and complex geometries are build up if configured.
* In this case, the geometry is loaded into memory and allocate space
* on the heap. The space can later be cleaned up by DONGEO.
*
* In: 
*   IGEODT : Integer parameter block of geometries
*            = TIGeometryData
*   DGEODT : Double precision parameter block of geometries
*            = TDGeometryData
*   IC2PAR : =0: initialize IGEODT/DGEODT only with standard values
*            =1: initialize IGEODT/DGEODT with standard values and
*                transfer values of COMMON-blocks into them.
*                Build up complex geometries if necessary, allocate
*                memory and load them from disc if necessary.
*   
* Out:
*  IGEODT,
*  DGEODT  : Initialized structures
*
* For nonsteady simulations, the caller should set the flag INSGEO to
* 1 to indicate that the simulation is noinsteady. INIGEO does not
* know anything anout the type of the simulation!
***********************************************************************
     
      SUBROUTINE INIGEO (IGEODT,DGEODT,IC2PAR)
      
      IMPLICIT NONE
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cout.inc'
      
      INCLUDE 'sgeometries.inc'
      INCLUDE 'scompgeometries.inc'
      
      INCLUDE 'sinigeometry.inc'
      INCLUDE 'cinigeometry.inc'
 
C     parameters

      INTEGER IGEODT(SZGEOI),IC2PAR
      DOUBLE PRECISION DGEODT(SZGEOD)
 
C     local variables
 
      INTEGER I,J,L1,L2
      
C     externals

      INTEGER  GCRCIR
      EXTERNAL GCRCIR

C     Hard-coded geometry data
      
      DOUBLE PRECISION RADII(2)
      DATA RADII /0.04D0, 0.02D0/
      
C     Transfer standard values into the parameter blocks:

      CALL LCL3(IGEODT,SZGEOI)
      CALL LCL1(DGEODT,SZGEOD)
      
      DGEODT(ODCRCOS) = 1D0
      DGEODT(ODMNISF) = 1D0
      DGEODT(ODMNOSF) = 1D0

C     Transfer COMMON block variables to parameter blocks if desired

      IF (IC2PAR.NE.0) THEN
      
        IGEODT(OICIRTP) = ICIRTP
        IGEODT(OIFBMNT) = IFBMNT
        IGEODT(OISCMNT) = ISCMNT
        IGEODT(OIOVMNT) = IOVMNT
        
C       Initialize the main geometry that "surrounds" everything

        DGEODT(ODCXPOS) = DCXPOS
        DGEODT(ODCYPOS) = DCYPOS

        DGEODT(ODCRAD ) = DCRAD 
        DGEODT(ODCRDX)  = DCRADX
        DGEODT(ODCRDY)  = DCRADY
        DGEODT(ODCROT ) = DCROT 
        DGEODT(ODCRSIN) = DCRSIN
        DGEODT(ODCRCOS) = DCRCOS
        DGEODT(ODMNISF) = DMNISF
        DGEODT(ODMNOSF) = DMNOSF
        
C       ---------------------------------------------------------------        
        
C       Initialise complex example geometry:

        IF (ICIRTP.GE.50) THEN
        
          IF (ICIRTP.EQ.60) THEN
          
C           Two independent circles

            IGEOCN = 2

C           Reserve a double-precision memory block
C           for the circles. Save the handle in the geometry
C           data block
            
            CALL ZNEW (100,1,L2,'IGEODL ')
            IGEODT(OIGEODL) = L2
            
C           The data block will receive two circle-structures.
C           In IGEOMS we have to store the starting addresses
C           of these blocks.
C           The first starting address is 1:

            J = 1
            
C           Loop through the circles and build their structures
C           in DWORK(L(L2)). Increment J on every new circle by
C           the memory consumption to get the new starting address. 

            DO I=1,IGEOCN

C             Save the starting address of that circle to 
C             IGEOMS(2,.). Set IGEOMS(1,.)=0 to indicate that we
C             are using a simple-type object with starting address

              IGEODT(OIGEOMS+3*(I-1))   = 0
              IGEODT(OIGEOMS+3*(I-1)+1) = J

C             Initialize the circle and calculate the starting
C             address of the next one.

              J = J + GCRCIR (DWORK(L(L2)+J-1), 0D0, 0D0, 0D0, 
     *                1D0, 0, RADII(I))
     
            END DO

C           Perform a second loop through the objects. This time
C           initialize the coordinate system to put all objects
C           to the position specified in the DAT file.

            DO I=1,IGEOCN

C             Get the starting address of that object and
C             update the coordinates.
      
              J = IGEODT(OIGEOMS+3*(I-1)+1)-1
              
              CALL GUPCOR (DWORK(L(L2)+J),
     *             DCXPOS, DCYPOS, 0D0, 1D0)

            END DO
              
            IGEODT(OIGEOCN)   = IGEOCN
        
          ELSE IF (ICIRTP.NE.56) THEN
        
C           Only one group; we reserve some memory in advance

            CALL ZNEW (500,3,L1,'IGEOM ')
            CALL ZNEW (500,1,L2,'DGEOM ')
          
C           Build the geometry in that 1st object

            CALL GINGEO (KWORK(L(L1)),
     *                   DWORK(L(L2)),
     *                   1,ICIRTP,DCRADX,DCRADY)
     
C           Initialise the coordinate system of the complete object
C           so that it is scaled/moved correctly.

            CALL GUPCOR (DWORK(L(L2)),
     *           DCXPOS, DCYPOS, DCROT, DCRADX)
     
C           Save the handles in the structure.
C           The handle of the integer block is saved to IGEOMS(2,1),
C           the handle of the double block in IGEOMS(3,1).
C           Set IGEOMS(1,1)=12 to indicate a group with
C           precomputed information.

            IGEODT(OIGEOCN)   = 1
            IGEODT(OIGEOMS  ) = 12
            IGEODT(OIGEOMS+1) = L1
            IGEODT(OIGEOMS+2) = L2
     
          ELSE
        
C           a group of two rotors in a case

            IGEOCN = 2
          
            DO I=1,IGEOCN
          
C             Reserve memory for that rotor

              IF (I.LE.2) THEN

                CALL ZNEW (30000,3,L1,'IGEOM ')
                CALL ZNEW (30000,1,L2,'DGEOM ')
                
              ELSE

                CALL ZNEW (MAXGEO,3,L1,'IGEOM ')
                CALL ZNEW (MAXGEO,1,L2,'DGEOM ')
              
              END IF

C             load both rotors into memory and create the cases
C             surrounding the rotors

              CALL GINGEO (KWORK(L(L1)),
     *                     DWORK(L(L2)),
     *                     I,ICIRTP,DCRADX,DCRADY)
        
C             Save the handles to the structure.
C             Set IGEOMS(1,1)=12 to indicate a group with
C             precomputed information.

              IGEODT(OIGEOMS+3*(I-1))   = 12
              IGEODT(OIGEOMS+3*(I-1)+1) = L1
              IGEODT(OIGEOMS+3*(I-1)+2) = L2
        
            END DO
            
            IGEODT(OIGEOCN)   = IGEOCN
          
C           Initialise the coordinate systems. DCRADY serves as the 
C           distance of the 2nd rotor to the 1st.

C           Update IGEOMS(2,2):
C            L2 = IGEODT(OIGEOMS+1) 
C            CALL GUPCOR (DWORK(L(L2)),
C      *          DCXPOS, DCYPOS, DCROT, DCRADX)
  
C            Update IGEOMS(2,2)
C            L1 = IGEODT(OIGEOMS+2+1) 
C            CALL GUPCOR (DWORK(L(L2)),
C      *          DCXPOS+DCRADY, DCYPOS, DCROT, DCRADX)

            CALL UPDGEO (IGEODT,DGEODT)

          END IF
          
        END IF
      
      END IF

      END

***********************************************************************
* Update composed geometries
*
* This routine must be called if the "geometry" consists of multiple
* "composed" geometries and the global coordinate system has changed.
* It will modify the "local" coordinate system of all objects
* according to the global coordinate system, which is saved
* in the IGEODT/DGEODT data blocks.
*
* If no composed geometry is used, this routine does nothing.
*
* In: 
*   IGEODT : Integer parameter block of geometries
*            = TIGeometryData
*   DGEODT : Double precision parameter block of geometries
*            = TDGeometryData
*
* Out:
*   Coordinate systems of composed objects in IGEODT/DGEODT
*   are updated.
***********************************************************************

      SUBROUTINE UPDGEO (IGEODT,DGEODT)
      
      IMPLICIT NONE
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'sinigeometry.inc'
      
      INCLUDE 'sgeometries.inc'
      INCLUDE 'scompgeometries.inc'
 
C     parameters

      INTEGER IGEODT(SZGEOI)
      DOUBLE PRECISION DGEODT(SZGEOD)

C     local variables
 
      INTEGER ICIRTP
      DOUBLE PRECISION T,DCXPOS, DCYPOS, DCROT, DCRADX, DCRADY
      INTEGER L1,L2
 
C     Do we have a composed geometry

      IF (IGEODT(OICIRTP).GT.50) THEN
      
C       Normally there is only one geometry.
C       Get the handle to the double-precision configurartion
C       of the first geometry (from IGEOMS(2,1)) and update.
      
        L2 = IGEODT(OIGEOMS+1)
        
C       Get rotation, radius,... from the parameter blocks

        ICIRTP = IGEODT(OICIRTP)
        DCXPOS = DGEODT(ODCXPOS)
        DCYPOS = DGEODT(ODCYPOS)
        DCROT  = DGEODT(ODCROT )
        DCRADX = DGEODT(ODCRDX )
        DCRADY = DGEODT(ODCRDY )
      
        CALL GUPCOR (DWORK(L(L2)),
     *       DCXPOS, DCYPOS, DCROT, DCRADX)
      
C       If we have a rotor geometry, there's more to update...
      
        IF (ICIRTP.EQ.56) THEN
        
C         Get the handles of both integer parameter blocks of 
C         the both rotors.

          L1 = IGEODT(OIGEOMS+1)
          L2 = IGEODT(OIGEOMS+2)
        
C         DCRADY serves as the distance of the 2nd rotor to the 1st.
C         The rotation of the 2nd rotor corresponds to the negative 
C         of the first one, modified by the relationship of the 
C         number of teeth between the rotors
          
          T = DBLE(KWORK(L(L1)+OGCCNT-1)) /
     *        DBLE(KWORK(L(L2)+OGCCNT-1))
          
C         Get the handles of the 2nd rotor

          L2 = IGEODT(OIGEOMS+3*1+1)
          
          CALL GUPCOR (DWORK(L(L2)),
     *         DCXPOS+DCRADY, DCYPOS, -T*DCROT+1D0/T*180D0+3.5D0,
     *         DCRADX)
C         -T*180D0
C
C         Correctly place the 3rd and 4th object on the positions 
C         of the 1st and 2nd rotor; they describe the case 
C         of the geometry.

          IF (IGEODT(OIGEOCN).GE.3) THEN
            L2 = IGEODT(OIGEOMS+3*2+1)
            CALL GUPCOR (DWORK(L(L2)),
     *           DCXPOS, DCYPOS, 0D0, DCRADX)
          END IF

          IF (IGEODT(OIGEOCN).GE.4) THEN
            L2 = IGEODT(OIGEOMS+3*3+1)
            CALL GUPCOR (DWORK(L(L2)),
     *           DCXPOS+DCRADY, DCYPOS, 0D0, DCRADX)
          END IF
        END IF
      END IF
      
      END
 
***********************************************************************
* Release complex geometries
*
* This deallocates all memory reserved by INIGEO for saving 
* geometries on the heap.
*
* In: 
*   IGEODT : Integer parameter block of geometries
*            = TIGeometryData
*   DGEODT : Double precision parameter block of geometries
*            = TDGeometryData
*
* Out:
*   All allocated memory is disposed.
***********************************************************************
     
      SUBROUTINE DONGEO (IGEODT)
      
      IMPLICIT NONE
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)
      
      INCLUDE 'cout.inc'
      INCLUDE 'sinigeometry.inc'
      
      INCLUDE 'sgeometries.inc'
      INCLUDE 'scompgeometries.inc'
 
C     parameters
 
      INTEGER IGEODT(SZGEOI)

C     local variables

      INTEGER I,L1,L2,L0
      
C     Free all memory; handles are found in IGEOMS.
      
      DO I=IGEODT(OIGEOCN),1,-1
        L0 = IGEODT(OIGEOMS+3*(I-1))
        L1 = IGEODT(OIGEOMS+3*(I-1)+1)
        L2 = IGEODT(OIGEOMS+3*(I-1)+2)
        
C       There's only memory allocated if L0=12
        
        IF (L0.EQ.12) THEN
          CALL ZDISP (0,L2,'DGEOM ')
          CALL ZDISP (0,L1,'IGEOM ')
        END IF
        
C       Otherwise: nothing to do.
        
      END DO
      
C     Release memory of IGEOIL and IGEODL, resp., if used.

      L1 = IGEODT(OIGEOIL)
      L2 = IGEODT(OIGEODL)
      IF (L2.NE.0) CALL ZDISP (0,L2,'IGEODL')
      IF (L1.NE.0) CALL ZDISP (0,L1,'IGEOIL')
 
      END
 
************************************************************************
* Initialise complex example-geometry
*
* The following routine provides some example geometries. It initialises
* the geometry-structures GEOIDX/GEOCOM according to the type
* ITYPE of the geometry. The parameter PAR1, PAR2,... allow to
* make some specifications to the geometry - their meaning is
* depending on the type of the example geometry.
*
* Some geometries consits of more than one object (e.g. two rotors).
* In that case this routine is called for every (sub-)object in that
* geometry (e.g. 1st rotor, 2nd rotor,...) with SBOBJ=number of
* current subobject to initialise.
*
* In:
*  GEOIDX - the GEOCOMIDX integer array that accepts the index-structure
*  GEOCOM - the GEOCOM double array that accepts the sub-objects
*           structure
*  SBOBJ  - Number of the current subobject to initialise.
*           normally = 1, but for geometries with more than one
*           subobject (2 rotors e.g.) this steps through 1..number of
*           subobjects
*  ITYPE  - type identifier; describes the type of example-geometry
*           GEOIDX/GEOCOM should be initialised to.
*           The TYPE-identifier in the MAIN structure of the geometry
*           receives this number, too.
*  PARx   - a set of parameters; example-specific
*
* Out:
*  GEOIDX,
*  GEOCOM - The initialised structure-arrays
*
* The arrays GEOIDX/GEOCOM must be large enough to hold the example
* geometry; there is no fixed array size assumed.
* The example-geometry is placed at position (0,0) with rotation=0.
* The user can move the geometry by modifying the parameters in the
* MAIN-part of the TGEOCOM-structure.
*
* EXAMPLES:
* ---------
*   ITYPE = 51 -- dump-bell geometry
*           0  -- a simple ball
*           1  -- Square
*           2  -- Ellipse
*           3  -- Rectangle
* Parameters:
*   PAR1  - length of the bar of the dump-bell
*   PAR2  - thickness of the bar, also modifying the size of the
*           circles/"weights"
* 
*   ITYPE = 52 -- enterprise geometry
* Parameters:
*   PAR1  - size of the geometry
*   PAR2  - not used
*
************************************************************************

      SUBROUTINE GINGEO (GEOIDX, GEOCOM, SBOBJ, ITYPE, PAR1, PAR2)
      
      IMPLICIT NONE

      INCLUDE 'sgeometries.inc'
      INCLUDE 'scompgeometries.inc'
      
      DOUBLE PRECISION GEOCOM(*), PAR1, PAR2
      INTEGER GEOIDX(*),ITYPE,SBOBJ
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)
      
      DOUBLE PRECISION PAR1SC
      INTEGER I,J
      CHARACTER*50 FNAME
      
      INTEGER GCRCIR,GCRSIM,GCRREC,GCRELL,GCRSWV,GCRSQR
      EXTERNAL GCRCIR,GCRSIM,GCRREC,GCRELL,GCRSWV,GCRSQR
      
      IF (ITYPE.EQ.51) THEN
        CALL LCL1(GEOCOM,SZGSIM+2*SZGCIR+SZGREC)
        GEOIDX(OGCCNT) = 3
C Guiding main-object        
        GEOCOM(OGCMAIN+OTYPE-1) = ITYPE
        GEOCOM(OGCMAIN+OREFX-1) = 0D0
        GEOCOM(OGCMAIN+OREFY-1) = 0D0
        GEOCOM(OGCMAIN+OROT-1) = 0D0
        GEOCOM(OGCMAIN+OSCALE-1) = 1D0
        GEOCOM(OGCMAIN+OROTSIN-1) = 0D0
        GEOCOM(OGCMAIN+OROTCOS-1) = 1D0
C circle
        GEOIDX(OGCIDX) = OGCOBJS
        GEOCOM(GEOIDX(OGCIDX) + OTYPE-1) = SGTCIR
        GEOCOM(GEOIDX(OGCIDX) + OREFX-1) = -0.5*PAR1
        GEOCOM(GEOIDX(OGCIDX) + OREFY-1) = 0D0
        GEOCOM(GEOIDX(OGCIDX) + OCRAD-1) = 1.5D0*PAR2
        GEOCOM(GEOIDX(OGCIDX) + OSCALE-1) = 1D0
C circle
        GEOIDX(OGCIDX+1) = OGCOBJS+SZGCIR
        GEOCOM(GEOIDX(OGCIDX+1) + OTYPE-1) = SGTCIR
        GEOCOM(GEOIDX(OGCIDX+1) + OREFX-1) = 0.5*PAR1
        GEOCOM(GEOIDX(OGCIDX+1) + OCRAD-1) = 1.5D0*PAR2
        GEOCOM(GEOIDX(OGCIDX+1) + OSCALE-1) = 1D0
C rectangle
        GEOIDX(OGCIDX+2) = GEOIDX(OGCIDX+1)+SZGCIR
        GEOCOM(GEOIDX(OGCIDX+2) + OTYPE-1) = SGTREC
        GEOCOM(GEOIDX(OGCIDX+2) + ORLENX-1) = PAR1
        GEOCOM(GEOIDX(OGCIDX+2) + ORLENY-1) = PAR2
        GEOCOM(GEOIDX(OGCIDX+2) + OSCALE-1) = 1D0
      END IF

C initialise enterprise geometry
      
      IF (ITYPE.EQ.52) THEN
        CALL LCL1(GEOCOM,SZGSIM+3*SZGCIR+5*SZGREC)
        PAR1SC = 1D0
        GEOIDX(OGCCNT) = 8
C Guiding main-object        
        GEOCOM(OGCMAIN+OTYPE-1) = ITYPE
        GEOCOM(OGCMAIN+OREFX-1) = 0
        GEOCOM(OGCMAIN+OREFY-1) = 0
        GEOCOM(OGCMAIN+OROT-1) = 0D0
        GEOCOM(OGCMAIN+OSCALE-1) = 1D0
        GEOCOM(OGCMAIN+OROTSIN-1) = 0D0
        GEOCOM(OGCMAIN+OROTCOS-1) = 1D0
C circle
        GEOIDX(OGCIDX) = OGCOBJS
        GEOCOM(GEOIDX(OGCIDX) + OTYPE-1) = SGTCIR
        GEOCOM(GEOIDX(OGCIDX) + OREFX-1) = -0.5*PAR1SC*15
        GEOCOM(GEOIDX(OGCIDX) + OREFY-1) = 0D0
        GEOCOM(GEOIDX(OGCIDX) + OCRAD-1) = 0.3D0*PAR1SC*15
        GEOCOM(GEOIDX(OGCIDX) + OSCALE-1) = 1D0
C circle
        GEOIDX(OGCIDX+1) = GEOIDX(OGCIDX)+SZGCIR
        GEOCOM(GEOIDX(OGCIDX+1) + OTYPE-1) = SGTCIR
        GEOCOM(GEOIDX(OGCIDX+1) + OREFX-1) = 0D0
        GEOCOM(GEOIDX(OGCIDX+1) + OREFY-1) = 4D0*PAR1SC
        GEOCOM(GEOIDX(OGCIDX+1) + OCRAD-1) = 0.5D0*PAR1SC
        GEOCOM(GEOIDX(OGCIDX+1) + OSCALE-1) = 1D0
C circle
        GEOIDX(OGCIDX+2) = GEOIDX(OGCIDX+1)+SZGCIR
        GEOCOM(GEOIDX(OGCIDX+2) + OTYPE-1) = SGTCIR
        GEOCOM(GEOIDX(OGCIDX+2) + OREFX-1) = 0D0
        GEOCOM(GEOIDX(OGCIDX+2) + OREFY-1) = -4D0*PAR1SC
        GEOCOM(GEOIDX(OGCIDX+2) + OCRAD-1) = 0.5D0*PAR1SC
        GEOCOM(GEOIDX(OGCIDX+2) + OSCALE-1) = 1D0
C rectangle
        GEOIDX(OGCIDX+3) = GEOIDX(OGCIDX+2)+SZGCIR
        GEOCOM(GEOIDX(OGCIDX+3) + OTYPE-1) = SGTREC
        GEOCOM(GEOIDX(OGCIDX+3) + OREFX-1) = 0D0
        GEOCOM(GEOIDX(OGCIDX+3) + OREFY-1) = 0D0
        GEOCOM(GEOIDX(OGCIDX+3) + ORLENX-1) = 0.9*PAR1SC*15
        GEOCOM(GEOIDX(OGCIDX+3) + ORLENY-1) = PAR1SC*2
        GEOCOM(GEOIDX(OGCIDX+3) + OSCALE-1) = 1D0
C rectangle
        GEOIDX(OGCIDX+4) = GEOIDX(OGCIDX+3)+SZGREC
        GEOCOM(GEOIDX(OGCIDX+4) + OTYPE-1) = SGTREC
        GEOCOM(GEOIDX(OGCIDX+4) + OREFX-1) = 0.25*PAR1SC*15
        GEOCOM(GEOIDX(OGCIDX+4) + OREFY-1) = 4D0*PAR1SC
        GEOCOM(GEOIDX(OGCIDX+4) + ORLENX-1) = 0.5*PAR1SC*15
        GEOCOM(GEOIDX(OGCIDX+4) + ORLENY-1) = PAR1SC
        GEOCOM(GEOIDX(OGCIDX+4) + OSCALE-1) = 1D0
C rectangle
        GEOIDX(OGCIDX+5) = GEOIDX(OGCIDX+4)+SZGREC
        GEOCOM(GEOIDX(OGCIDX+5) + OTYPE-1) = SGTREC
        GEOCOM(GEOIDX(OGCIDX+5) + OREFX-1) = 0.25*PAR1SC*15
        GEOCOM(GEOIDX(OGCIDX+5) + OREFY-1) = -4D0*PAR1SC
        GEOCOM(GEOIDX(OGCIDX+5) + ORLENX-1) = 0.5*PAR1SC*15
        GEOCOM(GEOIDX(OGCIDX+5) + ORLENY-1) = PAR1SC
        GEOCOM(GEOIDX(OGCIDX+5) + OSCALE-1) = 1D0
C rectangle
        GEOIDX(OGCIDX+6) = GEOIDX(OGCIDX+5)+SZGREC
        GEOCOM(GEOIDX(OGCIDX+6) + OTYPE-1) = SGTREC
        GEOCOM(GEOIDX(OGCIDX+6) + OREFX-1) = 0.5*PAR1SC*7
        GEOCOM(GEOIDX(OGCIDX+6) + OREFY-1) = 2D0*PAR1SC
        GEOCOM(GEOIDX(OGCIDX+6) + OROT-1) = -45D0
        GEOCOM(GEOIDX(OGCIDX+6) + OROTSIN-1) = 
     *    SIN(GEOCOM(GEOIDX(OGCIDX+6) + OROT-1)*PI/180D0)
        GEOCOM(GEOIDX(OGCIDX+6) + OROTCOS-1) = 
     *    COS(GEOCOM(GEOIDX(OGCIDX+6) + OROT-1)*PI/180D0)
        GEOCOM(GEOIDX(OGCIDX+6) + OSCALE-1) = 1D0
        GEOCOM(GEOIDX(OGCIDX+6) + ORLENX-1) = 0.5D0*PAR1SC*11
        GEOCOM(GEOIDX(OGCIDX+6) + ORLENY-1) = 0.85D0*PAR1SC
C rectangle
        GEOIDX(OGCIDX+7) = GEOIDX(OGCIDX+6)+SZGREC
        GEOCOM(GEOIDX(OGCIDX+7) + OTYPE-1) = SGTREC
        GEOCOM(GEOIDX(OGCIDX+7) + OREFX-1) = 0.5*PAR1SC*7
        GEOCOM(GEOIDX(OGCIDX+7) + OROT-1) = 45D0
        GEOCOM(GEOIDX(OGCIDX+7) + OROTSIN-1) = 
     *    SIN(GEOCOM(GEOIDX(OGCIDX+7) + OROT-1)*PI/180D0)
        GEOCOM(GEOIDX(OGCIDX+7) + OROTCOS-1) = 
     *    COS(GEOCOM(GEOIDX(OGCIDX+7) + OROT-1)*PI/180D0)
        GEOCOM(GEOIDX(OGCIDX+7) + OREFY-1) = -2D0*PAR1SC
        GEOCOM(GEOIDX(OGCIDX+7) + ORLENX-1) = 0.5D0*PAR1SC*11
        GEOCOM(GEOIDX(OGCIDX+7) + ORLENY-1) = 0.85D0*PAR1SC
        GEOCOM(GEOIDX(OGCIDX+7) + OSCALE-1) = 1D0
      END IF

C initialise simple test-geometry
      
      IF (ITYPE.EQ.53) THEN
        CALL LCL1(GEOCOM,SZGSIM+3*SZGCIR+5*SZGREC)
        GEOIDX(OGCCNT) = 2
C Guiding main-object        
        GEOCOM(OGCMAIN+OTYPE-1) = ITYPE
        GEOCOM(OGCMAIN+OREFX-1) = 0
        GEOCOM(OGCMAIN+OREFY-1) = 0
        GEOCOM(OGCMAIN+OSCALE-1) = 1D0
        GEOCOM(OGCMAIN+OROT-1) = 0D0
        GEOCOM(OGCMAIN+OROTSIN-1) = 0D0
        GEOCOM(OGCMAIN+OROTCOS-1) = 1D0
C rectangle
        GEOIDX(OGCIDX) = OGCOBJS
        GEOCOM(GEOIDX(OGCIDX) + OTYPE-1) = SGTELL
        GEOCOM(GEOIDX(OGCIDX) + OREFX-1) = 0D0
        GEOCOM(GEOIDX(OGCIDX) + OREFY-1) = 0D0
        GEOCOM(GEOIDX(OGCIDX) + ORLENX-1) = 1D0
        GEOCOM(GEOIDX(OGCIDX) + OSCALE-1) = 1D0
        GEOCOM(GEOIDX(OGCIDX) + ORLENY-1) = 1D0
C rectangle
        GEOIDX(OGCIDX+1) = GEOIDX(OGCIDX)+SZGELL
        GEOCOM(GEOIDX(OGCIDX+1) + OTYPE-1) = SGTREC
        GEOCOM(GEOIDX(OGCIDX+1) + OREFX-1) = 0
        GEOCOM(GEOIDX(OGCIDX+1) + OREFY-1) = 0
        GEOCOM(GEOIDX(OGCIDX+1) + ORLENX-1) = 1D0
        GEOCOM(GEOIDX(OGCIDX+1) + OSCALE-1) = 1D0
        GEOCOM(GEOIDX(OGCIDX+1) + ORLENY-1) = 1D0
      END IF

C     One simple object

      IF ((ITYPE.GE.0).AND.(ITYPE.LE.3)) THEN
        
C       ...with one circle as sub-object
        
        GEOIDX(OGCCNT) = 1
        
C       Guiding main-object        

        J = GCRSIM (GEOCOM(OGCMAIN), 0D0, 0D0, 0D0, 1D0, 0)
        GEOCOM(OGCMAIN+OTYPE-1) = ITYPE

C       initialise index-counter

        I = OGCOBJS

C       initialize the object

        GEOIDX(OGCIDX) = I
        
        IF (ITYPE.EQ.0) THEN
          I = I + GCRCIR (GEOCOM(I), 0D0, 0D0, 0D0, 1D0, 0, 
     *                    PAR1)
        ELSE IF (ITYPE.EQ.3) THEN
          I = I + GCRELL (GEOCOM(I), 0D0, 0D0, 0D0, 1D0, 0, 
     *                    PAR1,PAR2)
        ELSE IF (ITYPE.EQ.2) THEN
          I = I + GCRSQR (GEOCOM(I), 0D0, 0D0, 0D0, 1D0, 0, 
     *                    PAR1)
        ELSE IF (ITYPE.EQ.4) THEN
          I = I + GCRREC (GEOCOM(I), 0D0, 0D0, 0D0, 1D0, 0, 
     *                    PAR1,PAR2)
        END IF      
      END IF
      
C Dump-bell geometry with with initialisation by function calls:

      IF (ITYPE.EQ.54) THEN
        GEOIDX(OGCCNT) = 3
C Guiding main-object        
        J = GCRSIM (GEOCOM(OGCMAIN), 0D0, 0D0, 0D0, 1D0, 0)
        GEOCOM(OGCMAIN+OTYPE-1) = ITYPE
C initialise index-counter
        I = OGCOBJS
C circle
        GEOIDX(OGCIDX) = I
        I = I + GCRCIR (GEOCOM(I), -0.5*PAR1, 0D0, 0D0, 1D0, 0, 
     *                  1.5D0*PAR2)
C circle
        GEOIDX(OGCIDX+1) = I
        I = I + GCRCIR (GEOCOM(I), 0.5*PAR1, 0D0, 0D0, 1D0, 0,
     *                  1.5D0*PAR2)
C rectangle
        GEOIDX(OGCIDX+2) = I
        I = I + GCRREC (GEOCOM(I),0D0, 0D0, 0D0, 1D0, 0, PAR1, PAR2)
      END IF
      
C SIN-Wave

      IF (ITYPE.EQ.55) THEN
        GEOIDX(OGCCNT) = 2
C Guiding main-object        
        J = GCRSIM (GEOCOM(OGCMAIN), 0D0, 0D0, 0D0, 1D0, 0)
        GEOCOM(OGCMAIN+OTYPE-1) = ITYPE
C initialise index-counter
        I = OGCOBJS
C SIN-wave
        GEOIDX(OGCIDX) = I
        I = I + GCRSWV (GEOCOM(I), 0D0, PAR2, 0D0,  1D0, 0, PAR1,
     *        PAR2*1.1,1*PAR2/2D0,0D0,2D0, 0D0,PAR2/2D0,180D0,2D0)
C SIN-wave
        GEOIDX(OGCIDX+1) = I
        I = I + GCRSWV (GEOCOM(I), 0D0, -PAR2, 0D0,  1D0, 0, PAR1,
     *        PAR2*1.1,PAR2/2D0,180D0,2D0, 0D0,1*PAR2/2D0,0D0,2D0)
      END IF
      
C Rotor-geometry

      IF (ITYPE.EQ.56) THEN
      
C In this geometry we have 2 subobjects. Depending on the current
C index SBOBJ choose the correct filename and call the rotor
C initialisation routine to load that into memory.
C
C Apart from the rotors there is a case conisisting of two
C parts surrounding the rotors.
C The whole geometry has dimensions 250x170.
        
        IF (SBOBJ.EQ.1) THEN
        
C First rotor
        
          WRITE (FNAME,'(A)') 'pre/maleRotor_44_Punkte.txt'
C          WRITE (FNAME,'(A)') 'pre/maleRotor_400_Punkte.txt'
          CALL GINRTR (GEOIDX, GEOCOM, ITYPE, FNAME, PAR1)
          
        ELSE IF (SBOBJ.EQ.2) THEN
        
C 2nd rotor
        
          WRITE (FNAME,'(A)') 'pre/femaleRotor_47_Punkte.txt'
C          WRITE (FNAME,'(A)') 'pre/femaleRotor_399_Punkte.txt'
          CALL GINRTR (GEOIDX, GEOCOM, ITYPE, FNAME, PAR1)
          
        ELSE IF (SBOBJ.EQ.3) THEN
        
C left part of the case

          GEOIDX(OGCCNT) = 2
          
C Guiding main-object        
          J = GCRSIM (GEOCOM(OGCMAIN), 0D0, 0D0, 0D0, 1D0, 0)
          GEOCOM(OGCMAIN+OTYPE-1) = ITYPE
          
C initialise index-counter
          I = OGCOBJS
C rectangle
          GEOIDX(OGCIDX) = I
          I = I + GCRREC (GEOCOM(I), -32D0, 0D0, 0D0, 1D0, 0, 
     *                    100D0, 168D0)
C circle in the midpoint of the first rotor
          GEOIDX(OGCIDX+1) = I
          I = I + GCRCIR (GEOCOM(I), 0D0, 0D0, 0D0, 1D0, 1, 
     *                    51.1D0)


        ELSE IF (SBOBJ.EQ.4) THEN
        
C right part of the case        

          GEOIDX(OGCCNT) = 2
          
C Guiding main-object        
          J = GCRSIM (GEOCOM(OGCMAIN), 0D0, 0D0, 0D0, 1D0, 0)
          GEOCOM(OGCMAIN+OTYPE-1) = ITYPE
          
C initialise index-counter
          I = OGCOBJS
C rectangle
          GEOIDX(OGCIDX) = I
          I = I + GCRREC (GEOCOM(I), 32D0, 0D0, 0D0, 1D0, 0, 
     *                    100D0, 168D0)
C circle in the midpoint of the first rotor
          GEOIDX(OGCIDX+1) = I
          I = I + GCRCIR (GEOCOM(I), 0D0, 0D0, 0D0, 1D0, 1, 
     *                    51.1D0)
        
        END IF
        
        
      END IF
      
      END
      
************************************************************************
* Initialise rotor geometry.
*
* Reads and parses the file FNAME and build up the structure of a
* rotor according to it.
*
* In:
*  FNAME  - file name of the file to parse
*  ITYPE  - Type identifier of the object; saved in the structure
*  SCALE  - Scaling factor of the whole geometry
*
* Out:
*  GEOIDX,
*  GEOCOM - structure of a composed geometry describing the rotor
*           with pie slices
*
* The file must have the following structure:
*  Line 1: NSL  = Number of pie slices / teeth, the rotor is composed of
*  Line 2: NPT  = Number of points descibing the interface of one tooth
*  Line 3..3+NPT-1 = X/Y-coordinates of the points describing one tooth.
* Format of one line: 
*               "X ; Y" or "X Y" or "NUM: X Y"
*
* Each tooth is made of a list of lines between the points whose
* coordinates are specified starting in line 3.
*
* The routine will build NSL simple-geometry pie slices according to
* this data. The origin (0,0) is assumed to be the midpoint of the
* rotor.
************************************************************************
      
      SUBROUTINE GINRTR (GEOIDX, GEOCOM, ITYPE, FNAME, SCALE)
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      INCLUDE 'cerr.inc'

      INCLUDE 'sgeometries.inc'
      INCLUDE 'scompgeometries.inc'
      
      CHARACTER*(*) FNAME
      DOUBLE PRECISION GEOCOM(*),ROT,X,Y,RSIN,RCOS,SCALE
      INTEGER GEOIDX(*),ITYPE
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)

      INTEGER MFIL,IFM
      
      INTEGER NSL,NPT,I,J,AR,LSTR,IERR
      CHARACTER LIN*80
      
      INTEGER STCHR,GCRSIM,GCRPSL
      EXTERNAL STCHR,GCRSIM,GCRPSL
      
      IFM = 1
      MFIL = 71
      CALL OF0 (MFIL,FNAME,IFM)
      
C Number of teeth:

      READ (MFIL,*)  NSL
      
C DEBUG !!!
      
C      NSL = 1
      
      IF (NSL.LE.0) THEN
        WRITE (*,'(A)') 'Error reading geometry file!'
        WRITE (*,'(A)') 'NSL=0!'
        STOP
      END IF
      
      GEOIDX(OGCCNT) = NSL
      
C Number of points

      READ (MFIL,*) NPT

      IF (NPT.LE.0) THEN
        WRITE (*,'(A)') 'Error reading geometry file!'
        WRITE (*,'(A)') 'NPT=0!'
        STOP
      END IF
      
C Create an array for all the points

      CALL ZNEW (2*NPT+2,1,AR,'ARPTS ')
      
C Read in all the points

      DO I=0,NPT-1
        CALL STINL (LSTR, MFIL, IERR)
        
C Delete colons, replace semicolons and tabulators by spaces

        IF (LSTR.NE.0) THEN
          J = STCHR (LSTR,':')
          IF (J.NE.0) THEN
            CALL STDEL (LSTR,1,J)
          END IF
          
10        CONTINUE
          J = STCHR (LSTR,';')
          IF (J.NE.0) THEN
            CALL STSCH (LSTR,J,' ')
            GOTO 10  
          END IF

20        CONTINUE
          J = STCHR (LSTR,CHAR(9))
          IF (J.NE.0) THEN
            CALL STSCH (LSTR,J,' ')
            GOTO 20  
          END IF
          
C Trim and turn back into Fortran String
          
          CALL STTRM (LSTR)
          
          CALL STPUT (LSTR,LIN)
          CALL STDIS (LSTR)
          
C Interpret the line as 2 x double; X-coordinate and Y-coordinate

          READ (LIN,*) DWORK(L(AR)+2*I),DWORK(L(AR)+2*I+1)
          
C Scale the data

C          DWORK(L(AR)+2*I) = DWORK(L(AR)+2*I) 
C          DWORK(L(AR)+2*I+1) = DWORK(L(AR)+2*I+1) 
          
        ELSE
          WRITE (*,*) 'Error reading geometry file!'
          STOP
        END IF
      END DO

C That's all we read from that file
      
      CLOSE (MFIL)

C The last point in the list (number NPT+1) gets the same coordinates
C as the first point - rotated by one tooth. That way the pie slices
C overlap. 

C At first calculate the rotation angle depending on the 
C number of teeth:

      ROT = 2*PI/NSL
      RSIN = SIN(ROT)
      RCOS = COS(ROT)

      X = DWORK(L(AR))
      Y = DWORK(L(AR)+1)
      
      DWORK(L(AR)+2*NPT) = RCOS*X-RSIN*Y
      DWORK(L(AR)+2*NPT+1) = RSIN*X+RCOS*Y
      
C Initialise guiding main-object   
     
      J = GCRSIM (GEOCOM(OGCMAIN), 0D0, 0D0, 0D0, SCALE, 0)
      GEOCOM(OGCMAIN+OTYPE-1) = ITYPE
        
C initialise index-counter

      I = OGCOBJS
        
C Initialise all the teeth:

      DO J=1,NSL
        GEOIDX(OGCIDX+J-1) = I
        I = I + GCRPSL (GEOCOM(I), 0D0, 0D0, DBLE(J-1)*360.0/DBLE(NSL),
     *          1D0, 0, NPT+1,0,0D0,0D0,DWORK(L(AR)) )
C geht nicht mit ICLOS=1!
      END DO

C Free memory

      CALL ZDISP (0,AR,'ARPTS ')

      END
      