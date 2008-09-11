************************************************************************
* This file contains additional, non-basic routines for composed
* geometries. Like CPGEOMETRIES.F being an extension to GEOMETRIES.F,
* this file is an extension to COMPGEOMETRIES.F. It allowes to perform 
* the precalculation of values according to a complete given geometry,
* as defined in STRIA.INC.
*
* For this purpose there is a new set of wrapper routines introduced
* as well as a new set of computation routines that create a
* precomputed integer/double array from the sub-objects.
* 
* The node-dependent routines for fictitious boundary handling
* provided by this file normally receive a node-variable INODE
* with a node number to handle (determine distance to this node,...).
* The specification of the value INODE is generally mesh-dependent.
* Fortunately all routines here can pass the INODE value to the sub-
* objects, and so a different implementation of the node numbering
* will not alter this file. It must only be made sure that
* PCGEOMETRIES.F works correctly.
************************************************************************

************************************************************************
* GxPWRP Wrapper-routines
*
* The following wrapper routines are designed to map the TYPE
* specifier in the object structure to a call to the appropriate
* routine defined in GEOMETRIES.F. They share the same syntax
* as the "simple"-routines and are fully interchangeable with them.
*
* For simple geometries these routines call the calculation routines
* in GEOMETRIES.F directly. For more complex routines, precalculation
* can be handled here.
************************************************************************

*-----------------------------------------------------------------------
* Test if the node INODE is inside of the geometry xxx.
*-----------------------------------------------------------------------

      INTEGER FUNCTION GIPWRP (GEOSIM,CORSYS,TRIA,IINFO,DINFO, INODE)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      INCLUDE 'stria.inc'
      DOUBLE PRECISION GEOSIM(*),CORSYS(SZCORS),DINFO(*)
      INTEGER TRIA(SZTRIA)
      INTEGER IINFO(*),INODE
      DOUBLE PRECISION X,Y

      INTEGER TP
      
      INTEGER  GINCIR,GINELL,GINSQR,GINREC,GINSWV,GINPSL
      EXTERNAL GINCIR,GINELL,GINSQR,GINREC,GINSWV,GINPSL
      
      INTEGER GIPPSL
      EXTERNAL GIPPSL
      
C Standard implementation: Obtain coordinates and call the
C simple-geometry routine

      CALL NDE2XY (INODE,TRIA,X,Y)
      
      TP = GEOSIM(OTYPE)
      
      GIPWRP = 0
      IF (TP.EQ.SGTCIR) GIPWRP = GINCIR(GEOSIM,CORSYS,X,Y)
      IF (TP.EQ.SGTELL) GIPWRP = GINELL(GEOSIM,CORSYS,X,Y)
      IF (TP.EQ.SGTSQR) GIPWRP = GINSQR(GEOSIM,CORSYS,X,Y)
      IF (TP.EQ.SGTREC) GIPWRP = GINREC(GEOSIM,CORSYS,X,Y)
      IF (TP.EQ.SGTSWV) GIPWRP = GINSWV(GEOSIM,CORSYS,X,Y)
C      IF (TP.EQ.SGTPSL) GIPWRP = GINPSL(GEOSIM,CORSYS,X,Y)

C     For Pie-slices, the information can only be attained
C     by using precalculated information for now

      IF (TP.EQ.SGTPSL) GIPWRP = GIPPSL (GEOSIM,CORSYS,TRIA,
     *                                   IINFO,DINFO, INODE)
      
      END 

*-----------------------------------------------------------------------
* Calculate the normalized outer normal vector in the boundary-node
* INODE on the interface of the geometry object GEOSIM.
*-----------------------------------------------------------------------
      
      SUBROUTINE GNPWRP (GEOSIM,CORSYS,TRIA,IINFO,DINFO, INODE, XN,YN)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      INCLUDE 'stria.inc'
      DOUBLE PRECISION GEOSIM(*),CORSYS(SZCORS),DINFO(*)
      INTEGER TRIA(SZTRIA)
      INTEGER IINFO(*),INODE
      DOUBLE PRECISION X,Y,XN,YN
      INTEGER TP
      
C Standard implementation: Obtain coordinates and call the
C simple-geometry routine

      CALL NDE2XY (INODE,TRIA,X,Y)

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

      INTEGER FUNCTION GDPWRP (GEOSIM,CORSYS,TRIA,IINFO,DINFO, 
     *                         INODE, DIST)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      INCLUDE 'stria.inc'
      DOUBLE PRECISION GEOSIM(*),CORSYS(SZCORS),DINFO(*)
      INTEGER TRIA(SZTRIA)
      INTEGER IINFO(*),INODE
      DOUBLE PRECISION X,Y,DIST
      INTEGER TP
      
      INTEGER  GDSCIR,GDSELL,GDSSQR,GDSREC,GDSSWV,GDPPSL
      EXTERNAL GDSCIR,GDSELL,GDSSQR,GDSREC,GDSSWV,GDPPSL
      
C Standard implementation: Obtain coordinates and call the
C simple-geometry routine

      CALL NDE2XY (INODE,TRIA,X,Y)

      TP = GEOSIM(OTYPE)
      
      GDPWRP = 0
      IF (TP.EQ.SGTCIR) GDPWRP = GDSCIR(GEOSIM,CORSYS, X,Y, DIST)
      IF (TP.EQ.SGTELL) GDPWRP = GDSELL(GEOSIM,CORSYS, X,Y, DIST)
      IF (TP.EQ.SGTSQR) GDPWRP = GDSSQR(GEOSIM,CORSYS, X,Y, DIST)
      IF (TP.EQ.SGTREC) GDPWRP = GDSREC(GEOSIM,CORSYS, X,Y, DIST)
      IF (TP.EQ.SGTSWV) GDPWRP = GDSSWV(GEOSIM,CORSYS, X,Y, DIST)
      IF (TP.EQ.SGTPSL) 
     *    GDPWRP = GDPPSL (GEOSIM,CORSYS,TRIA,IINFO,DINFO,
     *                     INODE, DIST)
      
      END 

*-----------------------------------------------------------------------
* Project a point onto the boundary interface of the geometry object.
*-----------------------------------------------------------------------

      INTEGER FUNCTION GPPWRP (GEOSIM,CORSYS,TRIA,IINFO,DINFO, 
     *                         INODE, XP,YP)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      INCLUDE 'stria.inc'
      DOUBLE PRECISION GEOSIM(*),CORSYS(SZCORS),DINFO(*)
      INTEGER TRIA(SZTRIA)
      INTEGER IINFO(*),INODE
      DOUBLE PRECISION X,Y,XP,YP
      INTEGER TP
      
      INTEGER  GPRCIR,GPRELL,GPRSQR,GPRREC,GPRSWV,GPRPSL
      EXTERNAL GPRCIR,GPRELL,GPRSQR,GPRREC,GPRSWV,GPRPSL
      
C Standard implementation: Obtain coordinates and call the
C simple-geometry routine

      CALL NDE2XY (INODE,TRIA,X,Y)

      TP = GEOSIM(OTYPE)
      
      GPPWRP = 0
      IF (TP.EQ.SGTCIR) GPPWRP = GPRCIR(GEOSIM,CORSYS, X,Y, XP,YP)
      IF (TP.EQ.SGTELL) GPPWRP = GPRELL(GEOSIM,CORSYS, X,Y, XP,YP)
      IF (TP.EQ.SGTSQR) GPPWRP = GPRSQR(GEOSIM,CORSYS, X,Y, XP,YP)
      IF (TP.EQ.SGTREC) GPPWRP = GPRREC(GEOSIM,CORSYS, X,Y, XP,YP)
      IF (TP.EQ.SGTSWV) GPPWRP = GPRSWV(GEOSIM,CORSYS, X,Y, XP,YP)
      IF (TP.EQ.SGTPSL) GPPWRP = GPRPSL(GEOSIM,CORSYS, X,Y, XP,YP)
      
      END 

*-----------------------------------------------------------------------
* And finally - new in contrast to GEOMETRIES.F:
*
* Precalculate geometry information.
*-----------------------------------------------------------------------

      SUBROUTINE GPCWRP (GEOSIM,CORSYS, TRIA, ITYP, IINFO, DINFO)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      INCLUDE 'stria.inc'
      DOUBLE PRECISION GEOSIM(*),CORSYS(SZCORS),DINFO(*)
      INTEGER TRIA(SZTRIA)
      INTEGER IINFO(*),ITYP
      INTEGER TP
      
      EXTERNAL GPCPSL 
      
C Standard implementation: No precalculated values

      IF (ITYP.EQ.0) THEN
        IINFO(1) = 0
        IINFO(2) = 0
      END IF

      TP = GEOSIM(OTYPE)
      
      IF (TP.EQ.SGTPSL)
     *   CALL GPCPSL (GEOSIM,CORSYS, TRIA, ITYP, IINFO, DINFO)
      
      END 

************************************************************************
* The following set of routines now define the functionality
* of the precalculated composed-geometry objects. These routines
* are additional to the routines in COMPGEOMETRIES.F, the basic
* administration is perfomed there. What can be found here is:
*
* - calculation of an Integer/Double array containing gathered
*   information from the precalculated values in the sub-objects
* - Access to te precalculated values for distance calculation, 
*   normal vector calculation,...
************************************************************************

************************************************************************
* SUBROUTINE GPCCMP (GEOIDX,GEOCOM, 
*                    GPCFCT, GIPFCT, GNPFCT, GDPFCT, GPPFCT,
*                    TRIA, ITYP, IINFO, DINFO)
*
* Precalculate composed-geometry information.
*
* This routine pre-calculates information for composed geometries.
* The parameters and return values are similar to those described
* in GPCxxx in the file COMPGEOMETRIES.F - but some parameters are added
* to support composed geometries:
*
*  GEOIDX - the GEOCOMIDX integer array describing the index positions
*           of the objects in GEOCOM
*  GEOCOM - the GEOCOM double array describing the sub-objects
*  TRIA   - array [1..SZTRIA]
*           STRIA structure array; information about the geometry
*  ITYP   - =0: IINFO(0) returns the size necessary for IINFO when 
*               this routine is called with ITYP=1
*               IINFO(1) returns the size necessary for DINFO
*           =1: Calculate the information
*  GPCFCT - A wrapper function that provides calling the correct
*           GPCxxx-routine. GPCWRP is acceptable.
*  GIPFCT - A wrapper function that provides calling the correct
*           GIPxxx-routine. GIPWRP is acceptable.
*  GNPFCT - A wrapper function that provides calling the correct
*           GNPxxx-routine. GNPWRP is acceptable.
*  GDPFCT - A wrapper function that provides calling the correct
*           GDPxxx-routine. GDPWRP is acceptable.
*  GPPFCT - A wrapper function that provides calling the correct
*           GPPxxx-routine. GPPWRP is acceptable.
*
* Out:
*  IINFO  - array [1..*] of integer
*           If ITYP=0, IINFO(0) and IINFO(1) receive the necessary
*           size of IINFO/DINFO for storing precalculated values.
*           If ITYP=1, IINFO receives precalculated integer data
*  DINFO  - array [1..*] of double
*           If ITYP=1, IINFO receives precalculated double prec. data
*           Otherwise DINFO is not used.
************************************************************************
      
      SUBROUTINE GPCCMP (GEOIDX,GEOCOM, 
     *                   GPCFCT, GIPFCT, GNPFCT, GDPFCT, GPPFCT,
     *                   TRIA, ITYP, IINFO, DINFO)
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      INCLUDE 'stria.inc'
      INCLUDE 'sgeometries.inc'
      INCLUDE 'scompgeometries.inc'
      DOUBLE PRECISION GEOCOM(*),DINFO(*)
      INTEGER GEOIDX(*),IINFO(*),ITYP,TRIA(SZTRIA)
      
      INTEGER GIPFCT,GDPFCT,GPPFCT
      EXTERNAL GPCFCT, GIPFCT, GNPFCT, GDPFCT, GPPFCT 
      
      INTEGER GIPCMP
      EXTERNAL GIPCMP
      
C externals

      INTEGER TNDCNT
      EXTERNAL TNDCNT

C local variables      
      
      INTEGER IINFO2(2), DINFO2, H1, H2, I, INODE, K
      DOUBLE PRECISION D
      
C The easiest method for precalculation would be: Precalculate the
C information in every point for every sub-object. But this is crazy:
C Imagine an object composed of 100 sub-objects in a mesh with
C 1.000.000 unknowns - and your memory would be full ^^
C No, that makes no sense. Instead we save every information in one
C node, so if there are 100 sub-objects, every node receives a
C summary of the information about all the sub-objects.
C In detail, for the distance e.g. the minimum distance to all
C sub-objects is stored and so on...
C
C This of course makes it necessary to compute the information
C for every sub-object and then select the correct information to
C store. For this purpose we use the information whether a value
C was calculated by an object, or other geometrical information.
C If a value is not possible to calculate, it will not be used
C for the calculation of the "summary" value. If it is possible
C to calculate it, the type of the value decides about which
C geometrical approach should be used to find the "correct"
C final value.
C
C We assume the following structure of our integer/double
C arrays:
C
C  IINFO = record
C    INSOBJ : array [1..TNDCNT] of integer
C             number of sub-objects that contain a node
C  end
C
C  DINFO = record
C    DIST   : array [1..TNDCNT] of double
C             distances to all nodes
C  end
      
      IF (ITYP.EQ.0) THEN

C       Determine the size of the integer/double arrays

        IINFO(1) = 0
        IINFO(2) = 0
        
C       Ok, we have to store at least the distances 
C       -> for every node one double precision distance

        IINFO(2) = IINFO(2) + TNDCNT(TRIA)
        
C       and save how many sub-objects contain a node

        IINFO(1) = IINFO(1) + TNDCNT(TRIA)
        
      ELSE
      
C Precalculation routine. At first we want to calculate the minimum
C distance. We use precalculation routines for that. For this purpose
C we use the memory management, allocate a block of memory, store the
C precalculated information there, use this information to calculate
C our actual value, and release the block afterwards.
C 
C       Clear the distance memory block:

        CALL LCL1 (DINFO,TNDCNT(TRIA))

C       Clear the in-object memory block

        CALL LCL3 (IINFO,TNDCNT(TRIA))
        
C       Loop through all sub-objects 

        DO I=1,GEOIDX(OGCCNT)
        
C         How much memory does this object need for precalculated values
C         of the current geometry?
C         Use the wrapper-routine for general handling:

          CALL GPCFCT (GEOCOM(GEOIDX(OGCIDX+I-1)),
     *                 GEOCOM(OGCMAIN+OCSYS-1),
     *                 TRIA, 0, IINFO2, DINFO2);
     
C         Reserve some memory, large enough to store that information.
C         At least allocate 1 byte to prevent allocation routines from
C         allocating the whole memory.

          CALL ZNEW (MAX(1,IINFO2(1)),-3,H1,'H1    ')
          CALL ZNEW (MAX(1,IINFO2(2)),-1,H2,'H2    ')
          
C         Calculate the information of the sub-object

          CALL GPCFCT (GEOCOM(GEOIDX(OGCIDX+I-1)),
     *                 GEOCOM(OGCMAIN+OCSYS-1),
     *                 TRIA, 1, KWORK(L(H1)), DWORK(L(H2)))

C         Distance calculation: Loop through all sub-objects,
C         calculate the distance.
C         Our current distance is saved in DINFO(1..TNDCNT)

          DO INODE=1,TNDCNT(TRIA)
            
            K = GDPFCT (GEOCOM(GEOIDX(OGCIDX+I-1)),
     *                  GEOCOM(OGCMAIN+OCSYS-1),
     *                  TRIA,KWORK(L(H1)), DWORK(L(H2)),
     *                  INODE, D)
     
            IF (K.EQ.0) D=1D99
            
            IF (I.EQ.1) THEN
            
              DINFO (INODE) = D
              
            ELSE
            
C             Depending on whether the sub-object is inverted or not, 
C             take the minimum/maximum of the (perhaps negative) distance

              IF (GEOCOM(GEOIDX(OGCIDX+I-1)+OINVER-1).NE.0) THEN
                DINFO (INODE) = MAX (DINFO(INODE),D)
              ELSE
                DINFO (INODE) = MIN (DINFO(INODE),D)
              END IF
            END IF
      
          END DO

C         Make a second loop over the nodes. Count how many
C         objects contain each node.

          DO INODE=1,TNDCNT(TRIA)
            
C           Test whether a node is inside of a sub-object: count the 
C           number of objects containing a node.
            
            K = GIPFCT (GEOCOM(GEOIDX(OGCIDX+I-1)),
     *                  GEOCOM(OGCMAIN+OCSYS-1),
     *                  TRIA,KWORK(L(H1)), DWORK(L(H2)),INODE)
     
            IINFO (INODE) = IINFO (INODE) + K
            
          END DO
          
C         Release the temporary memory allocated for that object

          CALL ZDISP (0,H2,'H2    ')
          CALL ZDISP (0,H1,'H1    ')
     
        END DO
      
      END IF
      
      END

************************************************************************
* The next set of routines can be used to access the precomputed
* information.
************************************************************************
      
************************************************************************
* INTEGER FUNCTION GIPCMP (GEOIDX, GEOCOM, GINFCT, 
*                          TRIA,IINFO,DINFO, INODE)
*
* Test if the node INODE is inside of the composed geometry.
* Precalculated version.
*
* In:
*  GEOIDX - the GEOCOMIDX integer array describing the index positions
*           of the objects in GEOCOM
*  GEOCOM - the GEOCOM double array describing the sub-objects
*  GIPFCT - a wrapper function that provides calling the correct
*           GIPxxx-routine. GIPWRP is acceptable.
*  TRIA   - array [1..SZTRIA]
*           STRIA structure array; information about the geometry
*  IINFO  - array [1..*] of integer 
*           Precalculated integer information. Must correspond to TRIA.
*  DINFO  - array [1..*] of double
*           Precalculated double information. Must correspond to TRIA.
*  INODE  - number of the node 
*
* Out:
*  Result > 0, if the point is inside of the geometry object.
*              If the group-object is not completely inverted (INVER=1),
*              the number describes the number of sub-objects
*              that contain this point.
*         = 0, if the point is not inside of the geometry object
*         < 0, if the point is inside of the geometry object and the
*              object is inverted
************************************************************************
      
      INTEGER FUNCTION GIPCMP (GEOIDX, GEOCOM, GIPFCT, TRIA, 
     *                         IINFO, DINFO, INODE)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      INCLUDE 'scompgeometries.inc'
      INCLUDE 'stria.inc'
      DOUBLE PRECISION GEOCOM(*),DINFO(*)
      INTEGER TRIA(SZTRIA)
      INTEGER GEOIDX(*),IINFO(*),INODE
      INTEGER GIPFCT
      EXTERNAL GIPFCT
      
C     The number of nodes containing INODE is precalculated.
C     Obtain it from IINFO

      GIPCMP = IINFO(INODE)

C     Switch the result if we are in "inverted" mode
      
      IF (GEOCOM(OINVER).NE.0) THEN
        GIPCMP = -GIPCMP
      END IF
      
      END
      
************************************************************************
* SUBROUTINE GNPCMP (GEOIDX, GEOCOM, GNMFCT, TRIA, 
*                    IINFO, DINFO, INODE, XN,YN)
* 
* Calculate the normalized outer normal vector in the boundary-node
* INODE on the interface of the geometry object GEOSIM.
* Precalculated version.
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
*  TRIA   - array [1..SZTRIA]
*           STRIA structure array; information about the geometry
*  IINFO  - array [1..*] of integer 
*           Precalculated integer information. Must correspond to TRIA.
*  DINFO  - array [1..*] of double
*           Precalculated double information. Must correspond to TRIA.
*  INODE  - number of the node 
* 
* Out:
*  XN,YN  - Outer normal vector in the node INODE,
*           normalized to length = 1
*           =(0,0) if the calculation is not possible.
*
* This routine only works completely correct for points outside of
* the geometry. Points inside of the geometry may be projected
* wrong if they are lying in more than one sub-component!
************************************************************************

      SUBROUTINE GNPCMP (GEOIDX, GEOCOM, GNPFCT, GDPFCT, TRIA, 
     *                   IINFO, DINFO, INODE, XN,YN)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      INCLUDE 'scompgeometries.inc'
      INCLUDE 'stria.inc'
      DOUBLE PRECISION GEOCOM(*),DINFO(*),XN,YN
      INTEGER TRIA(SZTRIA)
      INTEGER GEOIDX(*),IINFO(*),INODE
      INTEGER GDPFCT
      EXTERNAL GNPFCT,GDPCMP,GDPFCT
      INTEGER GDPCMP
      
      DOUBLE PRECISION DIST
      INTEGER I

      DOUBLE PRECISION X,Y

      CALL NDE2XY (INODE,TRIA,X,Y)

C Standard implementation like in COMPGEOMETRIES.F - but for a different
C syntax of the functions!
      
      XN = 0D0
      YN = 0D0

C Calculate the sub-object with the minimum distance

      I = GDPCMP (GEOIDX, GEOCOM, GDPFCT, TRIA,IINFO,DINFO,INODE, DIST)
      
C Calculate the normal vector by that calculation routine.

C In the current implementation DINFO and IINFO are just passed as
C dummy parameters. They must not be used by GNPFCT since the definition
C of these arrays belong to the group object, not to the single object.
C This routine is therefore a matter of change: In later implementation
C this routine must act like other routines:
C - precalculate the normal vectors
C - obtain the normal vectors directly from the arrays without invoking
C   any other function
      
C      IF (I.NE.0) THEN
C        CALL GNPFCT(GEOCOM(GEOIDX(OGCIDX+I-1)),
C     *              GEOCOM(OGCMAIN+OCSYS-1),
C     *              TRIA,IINFO,DINFO,INODE, XN,YN)
C      END IF

C Switch the result if we are in "inverted" mode
      
      IF (GEOCOM(OINVER).NE.0) THEN
        XN = -XN
        YN = -YN
      END IF
      
      END

************************************************************************
* INTEGER FUNCTION GDPCMP (GEOIDX, GEOCOM, GDPFCT, TRIA, 
*                          IINFO, DINFO, INODE, DIST)
*
* Calculate the distance of the node INODE to the boundary interface
* of the object.
* Precalculated version.
*
* In:
*  GEOIDX - the GEOCOMIDX integer array describing the index positions
*           of the objects in GEOCOM
*  GEOCOM - the GEOCOM double array describing the sub-objects
*  GDPFCT - a wrapper function that provides calling the correct
*           GGMxxx-routine. GDPWRP is acceptable.
*  TRIA   - array [1..SZTRIA]
*           STRIA structure array; information about the geometry
*  IINFO  - array [1..*] of integer 
*           Precalculated integer information. Must correspond to TRIA.
*  DINFO  - array [1..*] of double
*           Precalculated double information. Must correspond to TRIA.
*  INODE  - number of the node 
*
* Out:
*  DIST   - min. distance of INODE to the boundary interfaces; the sign
*           identifies the position of the point:
*           > 0: the point is outside of the geometry object
*           <=0: the point is inside of the geometry object
*
*  Result: > 0, if the distance was calculated successfully.
*               The number describes the number of the sub-component
*               that served for calculating the distance.
*          <=0, if the computation was not possible
************************************************************************

      INTEGER FUNCTION GDPCMP (GEOIDX, GEOCOM, GDSFCT, TRIA, 
     *                         IINFO, DINFO, INODE, DIST)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      INCLUDE 'scompgeometries.inc'
      INCLUDE 'stria.inc'
      DOUBLE PRECISION GEOCOM(*),DINFO(*),DIST
      INTEGER TRIA(SZTRIA)
      INTEGER GEOIDX(*),IINFO(*),INODE
      INTEGER GDSFCT
      EXTERNAL GDSFCT

C The distance is precalculated!
C We can therefore directly take it from the array...

      DIST = DINFO(INODE)
      GDPCMP = 1

      END
      
************************************************************************
* INTEGER FUNCTION GPPCMP (GEOIDX, GEOCOM, GPRFCT, TRIA, 
*                          IINFO, DINFO, INODE, XP,YP)
*
* Project a point onto the boundary interface of the geometry object.
* Precalculated version.
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
*  TRIA   - array [1..SZTRIA]
*           STRIA structure array; information about the geometry
*  IINFO  - array [1..*] of integer 
*           Precalculated integer information. Must correspond to TRIA.
*  DINFO  - array [1..*] of double
*           Precalculated double information. Must correspond to TRIA.
*  INODE  - number of the node 
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

      INTEGER FUNCTION GPPCMP (GEOIDX, GEOCOM, GPPFCT, GDPFCT, 
     *                         TRIA, IINFO, DINFO, INODE, XP,YP)
      IMPLICIT NONE
      INCLUDE 'sgeometries.inc'
      INCLUDE 'scompgeometries.inc'
      INCLUDE 'stria.inc'
      DOUBLE PRECISION GEOCOM(*),DINFO(*),XP,YP
      INTEGER TRIA(SZTRIA)
      INTEGER GEOIDX(*),IINFO(*),INODE
      INTEGER GDPCMP,GPPFCT,GDPFCT
      EXTERNAL GDPCMP,GPPFCT,GDPFCT

      DOUBLE PRECISION DIST
      INTEGER I

C Calculate the sub-object with the minimum distance

      I = GDPCMP (GEOIDX, GEOCOM, GDPFCT,TRIA,IINFO, DINFO, INODE, DIST)
      
C Project onto that boundary. 
      
C In the current implementation DINFO and IINFO are just passed as
C dummy parameters. They must not be used by GPPFCT since the definition
C of these arrays belong to the group object, not to the single object.
C This routine is therefore a matter of change: In later implementation
C this routine must act like other routines:
C - precalculate the projection
C - obtain the points directly from the arrays without invoking
C   any other function

C      IF (I.NE.0) THEN
C        I = GPPFCT(GEOCOM(GEOIDX(OGCIDX+I-1)), 
C     *             GEOCOM(OGCMAIN+OCSYS-1),
C     *             TRIA, IINFO, DINFO, INODE, XP,YP)        
C      END IF

      GPPCMP = I

      END
