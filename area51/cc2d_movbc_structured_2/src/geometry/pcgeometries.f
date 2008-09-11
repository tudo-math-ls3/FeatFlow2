************************************************************************
* This file contains additional, non-basic routines for all simple
* geometries. The file GEOMETRIES.F and GEOMETRIES.INC only
* define the basic behaviour of simple geometries for single points. 
* In contrast the routines in this file allow to perform the
* precalculation of values according to a complete given geometry,
* as defined in STRIA.INC.
*
* The caller can use a geometry object in two ways:
*  1.) Directly use the geometry routines in GEOMETRIES.F or
*  2.) Pre-calculate the values using GPCxxx for the corresponding
*      geometry and then use GxPxxx to obtain the value of interest
*      for a specific node in the geometry.
*
* For simple geometries (like circle e.g.) it's not necessary to do
* precalculation. In this case the GxPxxx-routines should directly
* pass all queries from the framework to the appropriate routines in
* GEOMETRIES.F, i.e. all routines here should have a defined output.
* The framework can then generally use the precalculation routines
* without taking care about whether these routines really perform
* precalculation or not.
*
* On the other hand, all routines in GEOMETRIES.F should also be
* functional, as far as possible. But as the framework uses
* the precalculation routines, it's acceptable to let routines
* in GEOMETRIES.F abstract as long as the routines here work well --
* although this situalion should be avoided as long as possible!
*
* The node-dependent routines for fictitious boundary handling
* provided by this file normally receive a node-variable INODE
* with a node number to handle (determine distance to this node,...).
* The specification of the value INODE is mesh-dependent. The
* standard treatment of this value is:
*   INODE = 1..NVT                  - vertices
*         = NVT+1..NVT+NMT          - edge midpoints
*         = NVT+NMT+1..NVT+NMT+NEL  - element midpoints
* If the geometry is more complex (e.g. more than one node per
* edge or element), the user has to adapt the code according to that!
* The routines here have therefore to be synchronized to the node
* numbering defined/administrated in structria.f.
************************************************************************

************************************************************************
* SUBROUTINE GPCxxx (GEOSIM,CORSYS, TRIA, ITYP, IINFO, DINFO)
*
* Precalculate geometry information.
*
* This routine calculates depending on a given geometry TRIA
* (structure STRIA) information about the geometry. This infomation
* is then written into the IINFO and DINFO-array, which must be
* large enough to hold this information.
* By specifying ITYP=0 it's possible to obtain the size of the IINFO/
* DINFO array that's necessary by this routine.
*
* In:
*  GEOSIM - the geometry object that corresponds to geometry type xxx
*  CORSYS - the coordinate system of the geometry
*  TRIA   - array [1..SZTRIA]
*           STRIA structure array; information about the geometry
*  ITYP   - =0: IINFO(1) returns the size necessary for IINFO when 
*               this routine is called with ITYP=1,
*               IINFO(2) returns the size necessary for DINFO
*
* Out:
*  IINFO  - array [1..*] of integer
*           If ITYP=0:
*             IINFO(1) = necessary memory to store integer information
*             IINFO(2) = necessary memory to store double information
*           If ITYP=1: receives precalculated integer data
*  DINFO  - array [1..*] of double
*           If ITYP=1: receives precalculated double prec. data
*           Otherwise DINFO is not used.
*
* The structure of IINFO/DINFO is completely geometry-class dependent.
* If the caller wants to know a specific information from the structure,
* the GxPxxx-routines of that geometry have to be used to extract
* the information.
************************************************************************
* INTEGER FUNCTION GIPxxx (GEOSIM,CORSYS,TRIA,IINFO,DINFO, INODE)
*
* Test if the node INODE is inside of the geometry xxx.
* Precalculated version.
*
* In:
*  GEOSIM - the geometry object that corresponds to geometry type xxx
*  CORSYS - the coordinate system of the geometry
*  TRIA   - array [1..SZTRIA]
*           STRIA structure array; information about the geometry.
*  IINFO  - array [1..*} of integer 
*           Precalculated integer information. Must correspond to TRIA.
*  DINFO  - array [1..*} of double
*           Precalculated double information. Must correspond to TRIA.
*  INODE  - number of the node 
*
* Out:
*  Result = 1, if the node INODE is inside of the geometry object
*         = 0, if the node INODE is not inside of the geometry object
************************************************************************
* SUBROUTINE GNPxxx (GEOSIM,CORSYS,TRIA,IINFO,DINFO, INODE, XN,YN)
* 
* Calculate the normalized outer normal vector in the boundary-node
* INODE on the interface of the geometry object GEOSIM.
* Precalculated version.
* 
* In:
*  GEOSIM - the geometry object that corresponds to geometry type xxx
*  CORSYS - the coordinate system of the geometry
*  TRIA   - array [1..SZTRIA]
*           STRIA structure array; information about the geometry.
*  IINFO  - array [1..*} of integer 
*           Precalculated integer information. Must correspond to TRIA.
*  DINFO  - array [1..*} of double
*           Precalculated double information. Must correspond to TRIA.
*  INODE  - number of the node
* 
* Out:
*  XN,YN  - Outer normal vector in the point (X,Y), 
*           normalized to length = 1
*           =(0,0) if the calculation is not possible.
************************************************************************
* INTEGER FUNCTION GDPxxx (GEOSIM,CORSYS,TRIA,IINFO,DINFO, INODE, DIST)
*
* Calculate the distance of the point (X,Y) to the boundary interface
* of the object.
* Precalculated version.
*
* In:
*  GEOSIM - the geometry object that corresponds to geometry type xxx
*  CORSYS - the coordinate system of the geometry
*  TRIA   - array [1..SZTRIA]
*           STRIA structure array; information about the geometry.
*  IINFO  - array [1..*} of integer 
*           Precalculated integer information. Must correspond to TRIA.
*  DINFO  - array [1..*} of double
*           Precalculated double information. Must correspond to TRIA.
*  INODE  - number of the node 
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
* INTEGER FUNCTION GPPxxx (GEOSIM,CORSYS,TRIA,IINFO,DINFO, INODE, XP,YP)
*
* Project a point onto the boundary interface of the geometry object.
*
* In:
*  GEOSIM - the geometry object that corresponds to geometry type xxx
*  CORSYS - the coordinate system of the geometry
*  TRIA   - array [1..SZTRIA]
*           STRIA structure array; information about the geometry.
*  IINFO  - array [1..*} of integer 
*           Precalculated integer information. Must correspond to TRIA.
*  DINFO  - array [1..*} of double
*           Precalculated double information. Must correspond to TRIA.
*  INODE  - number of the node 
*
* Out:
*  XP,YP  - the projected point on the interface of GEOSIM
*
*  Result: > 0, if the distance was calculated successfully
*          <=0, if the computation was not possible
************************************************************************

************************************************************************
* Simple precomputed Geometry: User-defined pie-slice
************************************************************************

*-----------------------------------------------------------------------
* Precalculate values.
*
* The internal structure-ararys IINFO/DINFO are organized as follows:
*  DINFO(1..NNODES)        
*   -> Distances from the nodes to the polygon
*  DINFO(NNODES+1..NNODES+2*NNODES)
*   -> Projection of all nodes onto the polygon
*  DINFO(NNODES+2*NNODES+1..NNODES+4*NNODES)
*   -> Normal vectors of all nodes
*
*  IINFO(1..NNODES)
*   -> number of the line segment that is nearest to each node.
*-----------------------------------------------------------------------

      SUBROUTINE GPCPSL (GEOSIM,CORSYS, TRIA, ITYP, IINFO, DINFO)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'sgeometries.inc'
      INCLUDE 'stria.inc'

      DOUBLE PRECISION GEOSIM(SZGPSL+1),CORSYS(SZCORS),DINFO(*)
      INTEGER TRIA(SZTRIA),ITYP,IINFO(*)
      
      INTEGER NNODES,LCORVG,CNT,HND,I,LCORMG
      INTEGER TNDCNT,GIPPSL
      EXTERNAL TNDCNT,GIPPSL

C For the moment we only handle vertices; therefore we manually
C set NNODES to NVT instead of the total number of nodes!

      NNODES = TNDCNT(TRIA)
C      NNODES = TRIA(ONVT)

      IF (ITYP.EQ.0) THEN
        IINFO(1) = NNODES
        IINFO(2) = 5*NNODES
      ELSE

C Save the number of nodes in an integer variable. This type-cast
C is necessary because:
C - we can perform loops easier 
C - we must not call RDPLSB below with the address of the count-variable
C   pointing to a double - since RDPLSB expects an integer!

        CNT = AINT(GEOSIM(OPCNT))

C Before we can calculate anything, we have to transfer the coordinates
C of all our line-segment nodes into the correct ccordinate system.
C this is done dynamically:
C
C Allocate memory

        CALL ZNEW (2*CNT,-1,HND,'NDEPRJ')
        
C Transfer all nodes into the new coordinate system

        DO I=0,CNT-1
          CALL GRTTRF (GEOSIM,CORSYS, 
     *            GEOSIM(OPPOIN+2*I),GEOSIM(OPPOIN+2*I+1), 
     *            DWORK(L(HND)+I*2),DWORK(L(HND)+I*2+1))
        END DO

C Now we can calculate distances...

        LCORVG = TRIA(OLCORVG)
        
C At first we use the brute-force approach. These basic routines can 
C be found in the LINEDISTANCE.F-file.
C The routine RDPLSB does all the work for us. We plug in our
C transformed coordinates, so this routine can work with real
C coordinates.
C
C For the moment this routine does not support midpoints...
      
C Warning: RDPLSB accepts the number of vertices as integer, but we've 
C saved it as double - so don't forget to type-cast before !!!!!!!!!
C -> use CNT instead of GEOSIM(OPCNT)
      
        CALL RDPLSB(TRIA(ONVT),DWORK(L(LCORVG)),
     *              CNT,DWORK(L(HND)),
     *              DINFO(1),IINFO(1),DINFO(1+NNODES),DINFO(1+3*NNODES))
      
C Now handle all other nodes saved in DCORMG - if DCORMG exists

        LCORMG = TRIA(OLCORMG)
        
        IF (LCORMG.NE.0) THEN
      
C Warning: RDPLSB accepts the number of vertices as integer, but we've 
C saved it as double - so don't forget to type-cast before !!!!!!!!!
C -> use CNT instead of GEOSIM(OPCNT)
      
          CALL RDPLSB(NNODES-TRIA(ONVT),DWORK(L(LCORMG)),
     *                CNT,DWORK(L(HND)),
     *                DINFO(1+TRIA(ONVT)),IINFO(1+TRIA(ONVT)),
     *                DINFO(1+NNODES+2*TRIA(ONVT)),
     *                DINFO(1+3*NNODES+2*TRIA(ONVT)))
    
        END IF
      
C Release the memory of the node coordinates

        CALL ZDISP(0,HND,'NDEPRJ')
      
C So far so good, we have the distances and more... But the distances
C are specified to be negative insode of the geometries!
C So we have to perform another loop about the nodes to determine
C which node is inside.

        DO I=1,NNODES
          IF (GIPPSL (GEOSIM,CORSYS,TRIA,IINFO,DINFO, I).NE.0) THEN
            DINFO(I) = -DINFO(I)
          END IF
        END DO
      
      END IF
      
      END

*-----------------------------------------------------------------------
* Test for a point being inside of the object
*
* We assume the lines to be ordered counterclockwise, so we can use
* information about tangential/normal vector to calculate if the
* point is "left" or "right" of a line. Being "left" of the line
* is then equivalent to being "in" the object.
*
* Modified return value:
*  = -1 if an error occurred (Startpoint=Endpoint for a line segment)
*-----------------------------------------------------------------------
      
      INTEGER FUNCTION GIPPSL (GEOSIM,CORSYS,TRIA,IINFO,DINFO, INODE)
      
      IMPLICIT NONE
      
      INCLUDE 'sgeometries.inc'
      INCLUDE 'stria.inc'

      DOUBLE PRECISION GEOSIM(*),CORSYS(SZCORS),DINFO(*)
      INTEGER TRIA(SZTRIA),IINFO(*),INODE
      
      DOUBLE PRECISION TX,TY,T1X,T1Y,T2X,T2Y,LIN(2,2),XP,YP,X,Y,DNX,DNY
      INTEGER ISG,NNODES,LCORVG,I
      
      INTEGER TNDCNT
      EXTERNAL TNDCNT

      NNODES = TNDCNT(TRIA)
      
C Calculate the tangential vector using the precalculated line segment

      ISG = IINFO(INODE)-1
      
C Apply the coordinate system to the end- and to the starting point
C of the line to convert them to the current coordinate system.
C Subtract the coordinates, this gives the taggential vector in our
C current coordinate system.

      TX = GEOSIM(OPPOIN+2*ISG+2)
      TY = GEOSIM(OPPOIN+2*ISG+3)
      CALL GRTTRF (GEOSIM,CORSYS,TX,TY,T2X,T2Y)
      
      TX = GEOSIM(OPPOIN+2*ISG)
      TY = GEOSIM(OPPOIN+2*ISG+1)
      CALL GRTTRF (GEOSIM,CORSYS,TX,TY,T1X,T1Y)

      TX = T2X-T1X
      TY = T2Y-T1Y
      
      IF ((TX.EQ.0).AND.(TY.EQ.0)) THEN
C Startpoint=Endpoint is not allowed with line segments!
        GIPPSL = -1
        RETURN
      END IF
      
C Obtain the precalculated normal vector of the point,
C rotate it by 90 degrees.

      T2X = -DINFO(1+NNODES+2*NNODES+2*(INODE-1)+1)
      T2Y = DINFO(1+NNODES+2*NNODES+2*(INODE-1))
      
C If the "tangential" vector obtained by rotating of the normal vector
C points into the same direction as the "real" tangential vector,
C the normal vector points outside.

      IF ((TX*T2X+TY*T2Y).GE.0) THEN
      
C Point definitely outside
      
        GIPPSL = 0
        
      ELSE

C Ok, at a first glance the point is inside... but only at a 
C first glance!
C
C The problem is that the curve might not be closed!
      
        IF (GEOSIM (OPCLOS).NE.0) THEN
          
          GIPPSL = 1
          
        ELSE
        
C If the curve is not closed, we have a problem - what is inside
C and what outside?!?
C
C We define this with a picture of a pie-slice: the "leftmost"
C and "rightmost" point can be connected to the origin. If there
C are other "pie-slices" that stuck on our current one, this
C way we can determine exactly in which of these pie slices
C we stay at the moment - because they are all distinct
C (more specifically they meet only in one line).
C
C So we are "in" the pie slice, if the point is "left" of the
C polygon above, "left" of the line "origin->rightmost" and
C "right" of the line "origin->leftmost"

C First take the node and transform it back into our standard
C coordinate system:

          CALL NDE2XY (INODE,TRIA,X,Y)
          CALL GRTBCK (GEOSIM,CORSYS,X,Y,XP,YP)

C The vector "origin->(XP,YP)" is then to analyze...

C Calculate the tangential vector of the line
C      "origin->rightmost point"
C and rotate it be -90 degrees

          I = GEOSIM(OPIRMP)-1
          DNX = GEOSIM(OPPOIN+2*I+1)
          DNY = -GEOSIM(OPPOIN+2*I)
          
C This should not point into the same direction as the normal vector - 
C or we are "outside" of the pie-slice!

          IF ((DNX*XP+DNY*YP).GE.0D0) THEN
          
            GIPPSL = 0
            RETURN
          
          END IF

C Calculate the tangential vector of the line
C      "origin->leftmost point"
C and rotate it be 90 degrees:

          I = GEOSIM(OPILMP)-1
          DNX = -GEOSIM(OPPOIN+2*I+1)
          DNY = GEOSIM(OPPOIN+2*I)

C This should not point into the same direction as the normal vector - 
C or we are "outside" of the pie-slice!

          IF ((DNX*XP+DNY*YP).GE.0D0) THEN
          
            GIPPSL = 0
            RETURN
          
          END IF
          
C We are inside - finally :)
          
          GIPPSL = 1
        
        END IF
      END IF

      END
      
*-----------------------------------------------------------------------
* Calculate the normalized outer normal vector the geometry
* neat a node
*-----------------------------------------------------------------------

      SUBROUTINE GNPPSL (GEOSIM,CORSYS,TRIA,IINFO,DINFO, INODE, XN,YN)

      IMPLICIT NONE
      
      INCLUDE 'sgeometries.inc'
      INCLUDE 'stria.inc'

      DOUBLE PRECISION GEOSIM(*),CORSYS(SZCORS),DINFO(*)
      DOUBLE PRECISION XN,YN
      INTEGER TRIA(SZTRIA),IINFO(*),INODE

      INTEGER ISG
      DOUBLE PRECISION D,TX,TY
      
C Which line segment is the point projected to?

      ISG = IINFO(INODE)
      
C Use start/endpoint of that line segment and calculate the 
C normalized outer normal vector

      TX = GEOSIM(OPPOIN+2*ISG+2)-GEOSIM(OPPOIN+2*ISG)
      TY = GEOSIM(OPPOIN+2*ISG+3)-GEOSIM(OPPOIN+2*ISG+1)
      
      D = 1D0/SQRT(TX*TX+TY*TY)
      XN = -D*TY
      YN = D*TX
      
      END
      
*-----------------------------------------------------------------------
* Calculate the distance of a node
*-----------------------------------------------------------------------

      INTEGER FUNCTION GDPPSL (GEOSIM,CORSYS,TRIA,IINFO,DINFO,
     *                         INODE, DIST)
      
      IMPLICIT NONE
      
      INCLUDE 'sgeometries.inc'
      INCLUDE 'stria.inc'

      DOUBLE PRECISION GEOSIM(SZGPSL),CORSYS(SZCORS),DINFO(*),DIST
      INTEGER TRIA(SZTRIA),IINFO(*),INODE
      
C Can be directly obtained from the array with the precalculated values

      DIST = DINFO(INODE)
      GDPPSL = 1
      
      END
      
*-----------------------------------------------------------------------
* Calculate the projection of a node onto the geometry
*-----------------------------------------------------------------------

      INTEGER FUNCTION GPPPSL (GEOSIM,CORSYS,TRIA,IINFO,DINFO, INODE, 
     *                         XP,YP)

      IMPLICIT NONE
      
      INCLUDE 'sgeometries.inc'
      INCLUDE 'stria.inc'

      DOUBLE PRECISION GEOSIM(SZGPSL),CORSYS(SZCORS),DINFO(*)
      DOUBLE PRECISION XP,YP
      INTEGER TRIA(SZTRIA),IINFO(*),INODE

      INTEGER NNODES
      
      INTEGER TNDCNT
      EXTERNAL TNDCNT

      NNODES = TNDCNT(TRIA)
      
C Can be taken directly from the array with the precalculated values

      XP = DINFO(1+NNODES+2*(INODE-1))
      YP = DINFO(1+NNODES+2*(INODE-1)+1)
      
      GPPPSL = 1
      
      END
      