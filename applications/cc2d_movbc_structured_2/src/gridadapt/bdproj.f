************************************************************************
* This file contains additional routines for the correction of the
* geometry, thus improving those routines of rcgrid.f.
*
* A projection routine is introduced that projects boundary nodes back
* onto the boundary when the type of the boundary is a circle, not
* a line.
************************************************************************

************************************************************************
* Reconstruct parameter values of boundary nodes for lines and circles
* 
* This routine analyses the coordinates of vertices on the (real)
* boundary and reconstructs their parameter values.
*
* For now this works only if the grid fulfills the following
* conditions:
* - the line segments are [0,1]-parametrised with equally
*   distributed parametrisation (i.e. something like (1-t^2)*P1+t^2*P2
*   is not allowed as parametrisation of a line)
* The routine works for lines and circles as boundary components.
* Parameter values that are integer numbers are not corrected; it's
* assumed that these points are fixed anchors of the geometry.
* Apart of the parameters the coordinates of the boundary nodes
* are also corrected (if necessary, e.g. for circle boundary segments)
*
* Furthermore the routine assumes that the OMEGA-parametrization
* is used, thus the variables in the COMMON-block defined by the
* CPARQDATA.INC file must be initialised properly.
* 
* In:
*  DCORVG - Coordinates of the gridpoints
*  KVBD   - Numbers of vertices on the boundary
*  KBCT   - Pointer vector for vertices on the boundary
*  DVBDP  - Parameter values of boundary points
* From COMMON blocks:
*  NBCT   - Number of boundary components
*
* Out:
*  DVBDP  - Parameter values of boundary points, corrected
*  DCORVG - Coordinates of the boundary nodes are corrected
************************************************************************

      SUBROUTINE RXBNPC(DCORVG,KVBD,KBCT,DVBDP)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cbasictria.inc'
      INCLUDE 'cparqdata.inc'

      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)

C parameters
      
      INTEGER KVBD(*),KBCT(*)
      DOUBLE PRECISION DCORVG(2,*),DVBDP(*)
      
C externals

      DOUBLE PRECISION PARX,PARY,TMAX
      EXTERNAL PARX,PARY,TMAX
      
C local variables

      INTEGER IBCT,ISEGM,NSEGM,NVBS,ISI,IEI,I,ICS,ICE,ISV,IEV,IFIRST
      DOUBLE PRECISION DSSTRT,DSEND,SLEN,SPLEN,D
      
      INTEGER ICPTR,IICOMP,ITYP,IPPTR,KXPAR,KYPAR
      DOUBLE PRECISION DXM,DYM,DY1,DY2,DX1,DX2,DRAD,DPHI1,DPHI2,DPHI
      DOUBLE PRECISION DY1R, DY2R
      
C loop through the bondary components

      DO IBCT=1,NBCT
      
C For later use save the starting point of the current boundary
C inside the ITYP-array in ICPTR

        ICPTR = KWORK(L(LICPTR)+IBCT-1)
      
C In each boundary component loop through the boundary segments. One
C segment can contain arbitrary many intermediate segments of the
C refinement, so we first have to find the starting- and ending point
C of the segment. They are marked by the property that:
C - the first point contains parameter value 0D0
C - each corner point of a segment has an integer parameter value
C - the last point is equal to the first point
C The length of the segment we can obtain by using TMAX. Although this
C function returns a double, the returned number is indeed an integer
C because of the [0,1]-parametrization.

C We add a slight correction to TMAX to avoid probable rounding errors
C with DINT...

        NSEGM = DINT(TMAX(IBCT)+0.000001)

C Get the number of vertices on the current boundary segment

        NVBS = KBCT(IBCT+1)-KBCT(IBCT)

C Index of the first node on this boundary:

        IFIRST = KBCT(IBCT)

C ISI and IEI contain the index of the starting- and ending-point
C of the current segment. Initialise them with values that point
C "before" the current boundary component.

        ISI = KBCT(IBCT)-1
        IEI = KBCT(IBCT)
        
C DSSTRT and DSEND contain the parameter value of the starting-
C end ending point of the current segment

        DSSTRT = 0D0
        DSEND  = 0D0

C Loop through the segments

        DO ISEGM=1,NSEGM
        
C Find the starting- and ending point; start the search after the
C last found segment.
C Our last ending point is the new starting point
        
          ISI = IEI
          DSSTRT = DSEND

C Don't do anything if we are on the last segment and
C there's no node on that last boundary segment

          IF (ISI.LT.(NVBS+IFIRST-1)) THEN

C Search the endpoint of the segment

            DO IEI=ISI+1,NVBS+IFIRST-1
              IF (DVBDP(IEI).EQ.DINT(DVBDP(IEI))) GOTO 10
            END DO
10          CONTINUE

C Obtain the parameter value of the ending point.
C Obtain the starting- and ending-index of the points in the "inner"
C of the current boundary segment. Obtain the length of the
C current segment.

            IF (IEI.GT.(NVBS+IFIRST-1)) THEN

C We have left our allowed parameter interval since we have "behind"
C the "last point". The ending point is the starting point of the interval
C since the curve is closed. It has got exactly the maximum
C parameter value!

              DSEND = TMAX(IBCT)
            
C Vertex number of starting- and ending points

              ISV = KVBD (ISI)
              IEV = KVBD (IFIRST)

            ELSE
          
              DSEND = DVBDP(IEI)
            
C Vertex number of starting- and ending points

              ISV = KVBD (ISI)
              IEV = KVBD (IEI)

            END IF

C Starting- and ending index of "inner" nodes

            ICS = ISI+1
            ICE = IEI-1

C Now we are at the point where we have to decide whether the
C current boundary segment is a circle segment, a line segment
C or what else.
          
C Determine the type of the current boundary segment.
C Look in the ITYP-array that saves for every boundary segment
C the corresponding type. ICOMP provides the number of the
C boundary segment, ICPTR points to the beginning of the current
C boundary component inside of ITYP, then the type of the current 
C node is...

            IICOMP=ICPTR+ISEGM-1
            ITYP=KWORK(L(LITYP)+IICOMP)
            
C and the values describing the boundary segment are saved
C in the DWORK-Array at position...
            
            IPPTR=KWORK(L(LIPPTR)+IICOMP)
            KXPAR=L(LXPAR)+IPPTR
            KYPAR=L(LYPAR)+IPPTR

C Do we have a line?

            IF (ITYP.EQ.1) THEN
          
C Length of the line segment:

              SLEN = DSQRT((DCORVG(1,IEV)-DCORVG(1,ISV))**2+
     *                    (DCORVG(2,IEV)-DCORVG(2,ISV))**2)

              IF (SLEN.EQ.0D0) THEN
                WRITE (*,*) 'Error in RXMPTS: SLEN=0 !'
                STOP
              END IF
              
C Or do we have a circle?
              
            ELSE IF (ITYP.EQ.2) THEN
            
C Obtain midpoint, radius and angle of the starting point
C of the circle.
            
              DXM=DWORK(KXPAR)
              DYM=DWORK(KYPAR)
              DRAD=DWORK(KXPAR+1)
              DPHI1=DWORK(KXPAR+2)
              DPHI2=DWORK(KYPAR+2)
              
C Then obtain the vector pointing from the midpoint to the starting
C point of the circle:
              
C              DY1=DXM+DRAD*COS(DPHI1)
C              DY2=DYM+DRAD*SIN(DPHI1)
              DY1=DRAD*COS(DPHI1)
              DY2=DRAD*SIN(DPHI1)
              
C Rotate the vector by 90 degrees into the direction of the
C parametrization. This gives us the information where's the 
C upper and lower part of the circle
      
              DY1R = DRAD*COS(DPHI1+0.5*PI*SIGN(1D0,DPHI2-DPHI1))
              DY2R = DRAD*SIN(DPHI1+0.5*PI*SIGN(1D0,DPHI2-DPHI1))
            
            END IF

C Loop through all "inner" points of the segment (ICS..ICE) and
C correct their parameter values. This is done by analysing the
C distance from the point to the starting- and ending-point of
C the line.

            DO I=ICS,ICE
          
C Do we have a line segment?

              IF (ITYP.EQ.1) THEN
          
C Length of the part of the segment between current node and
C starting node of the segment

                SPLEN = DSQRT((DCORVG(1,KVBD(I))-DCORVG(1,ISV))**2+
     *                        (DCORVG(2,KVBD(I))-DCORVG(2,ISV))**2)
     
C Use this to rebuild the parameter value
     
                DVBDP (I) = DSSTRT + SPLEN / SLEN
                
C Or a circle?

              ELSE IF (ITYP.EQ.2) THEN 
              
C Our current boundary point is at X/Y-position...

                DX1=DCORVG(1,KVBD(I))-DXM
                DX2=DCORVG(2,KVBD(I))-DYM
              
C Determine the angle of the projected boundary node in relation
C to the starting point of the circle.
C Determine the angle
          
                DPHI=ACOS((DX1*DY1+DX2*DY2)/
     *                  (DSQRT(DX1*DX1+DX2*DX2)*DSQRT(DY1*DY1+DY2*DY2)))
     
C and if we are in "the other half" of the circle, correct it:
     
                IF (DY1R*DX1+DY2R*DX2.LT.0D0) DPHI=2*PI-DPHI

C and calculate its parameter value

                D = DPHI/ABS(DPHI2-DPHI1) + DBLE(ISEGM-1) 
              
                IF (ABS(DVBDP(I)-D).LE.0.5D0) THEN

C We must prevent the point from switching one side of the fixed 
C start/end point of the circle to the other side.
C If the point wants to switch, we must prevent this!

                  DVBDP(I) = D

                END IF
              
              END IF
              
C At last recalculate the coordinates of the boundary node.
C For line segments these are the same as before. For circle
C segments the coordinates have changed...

              DCORVG(1,KVBD(I)) = PARX(DVBDP(I),IBCT)
              DCORVG(2,KVBD(I)) = PARY(DVBDP(I),IBCT)
            
            END DO
              
          END IF

        END DO
       
      END DO
      
      END

**********************************************************************
* Boundary-projection, linear search
*
* This routine projects a point (X0,Y0) onto the nearest (real)
* boundary segment. For this purpose a linear serarch is performed
* along all boundary edges, either of all boundary components
* or on a specified boundary component.
*
* The routine returns 
* - the coordinates (XNEW,YNEW) of the projected point
* - the index IEDG (range: 1..NVBD) inside of the KMBD-array of the 
*   edge, that contains the point. 
* - if necessary the number of the boundary component IBC that
*   contains the point
* The routine projects on the current line-representation of the
* boundary, not onto the parametrized boundary!
*
* In:
*   TRIA   - array [1..SZTRIA] of integer
*            Triangulation structure containing inbformation about
*            the boundary
*   X0,
*   Y0     - coordinates of the point to project
*   IBC    - Boundary component where to search.
*            =0: search on all boundary components, return found
*                boundary component in IBC
*            >0: search only on boundary component IBC
*
* Out:
*   XNEW,
*   YNEW   - coordinates of the projected point
*   IEDG   - index inside of the KMBD-array of the boundary edge
*            where the point was projected to
*   IBC    - Number of the boundary component where the projected
*            point was found, if IBC=0 on call to this routine;
*            otherwise not changed
**********************************************************************

      SUBROUTINE BPRLSR(TRIA, X0, Y0, IBC, XNEW, YNEW, IEDG)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      INCLUDE 'cmem.inc'
      
C parameters
      
      DOUBLE PRECISION X0,Y0,XNEW,YNEW
      INTEGER TRIA(SZTRIA)
      INTEGER IBC,IEDG

C local variables

      INTEGER SBC,EBC,ICBC,IIED,IVT1,IVT2
      INTEGER KMBD,DCORVG,KBCT
      INTEGER IED,IBCNEW
      DOUBLE PRECISION DPPT(2),D,D1,D2
      
      KMBD = L(TRIA(OLMBD))
      KBCT = L(TRIA(OLBCT))
      DCORVG = L(TRIA(OLCORVG))
      
C either run through all BC's or only through one

      IF (IBC.EQ.0) THEN
        SBC = 1
        EBC = TRIA(ONBCT)
      ELSE
        SBC = IBC
        EBC = IBC
      END IF

      DO ICBC = SBC,EBC
      
C If we just started analyzing, remember the first edge - so that we have
C at least any edge to return

        IF (ICBC.EQ.SBC) THEN

          IEDG = KWORK(KBCT+ICBC-1)

C What are the adjacent nodes?

          CALL ED2NDS (TRIA,IEDG,IVT1,IVT2)
          
C Project the point to that edge
          
          CALL PPO2SG(X0,Y0,
     *         DWORK(DCORVG+2*(IVT1-1)),DWORK(DCORVG+2*(IVT2-1)),
     *         D,DPPT)
     
          XNEW = DPPT(1)
          YNEW = DPPT(2)
          D1 = (X0-DPPT(1))**2+(Y0-DPPT(2))**2
          IBCNEW = ICBC
     
        END IF
      
C Now loop through the boundary edges on this boundary component
      
        DO IIED = KWORK(KBCT+ICBC-1),KWORK(KBCT+ICBC)
        
C Current edge? (number, not index)

          IED = KWORK(KMBD+IIED-1)
        
C What are the adjacent nodes?

          CALL ED2NDS (TRIA,IED,IVT1,IVT2)
          
C Project the point to that edge
          
          CALL PPO2SG(X0,Y0,
     *         DWORK(DCORVG+2*(IVT1-1)),DWORK(DCORVG+2*(IVT2-1)),
     *         D,DPPT)
        
C Calculate the distance (squared) - are we nearer?

          D2 = (X0-DPPT(1))**2+(Y0-DPPT(2))**2
          
          IF (D2.LT.D1) THEN
C If we are nearer, save edge index, point, distance          
            IEDG = IIED
            XNEW = DPPT(1)
            YNEW = DPPT(2)
            IBCNEW = IBC
            D1 = D2
          END IF
        
        END DO
      
      END DO

C Return the boundary component if necessary

      IF (IBC.EQ.0) IBC=IBCNEW

      END

**********************************************************************
* Boundary-projection, linear search, enhanced
*
* This routine projects a point (X0,Y0) onto the nearest (real)
* boundary segment. The routine accepts a parameter IEL with the
* number of an element that contains the point (X0,Y0). If the element
* is a boundary element, the routine projects the point to the
* nearest boundary edge on that element. If IEL is not given or is
* not a boundary element, BPRLSR is called to perform a linear
* search to find the correct projection.
*
* In:
*   TRIA   - array [1..SZTRIA] of integer
*            Triangulation structure containing inbformation about
*            the boundary
*   X0,
*   Y0     - coordinates of the point to project
*   IEL    - element containing (X0,Y0); can be 0 if the element
*            is not known
*   IBC    - Boundary component where to search, if IEL=0 or
*            IEL is not a boundary element.
*            =0: search on all boundary components, return found
*                boundary component in IBC
*            >0: search only on boundary component IBC
*
* Out:
*   XNEW,
*   YNEW   - coordinates of the projected point
*   IEDG   - index inside of the KMBD-array of the boundary edge
*            where the point was projected to
*   IBC    - Number of the boundary component where the projected
*            point was found, if IBC=0 on call to this routine;
*            otherwise not changed
**********************************************************************

      SUBROUTINE BPRLSE(TRIA, X0, Y0, IEL, IBC, XNEW, YNEW, IEDG)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      INCLUDE 'cmem.inc'
      
C parameters
      
      DOUBLE PRECISION X0,Y0,XNEW,YNEW
      INTEGER TRIA(SZTRIA)
      INTEGER IBC,IEDG,IEL

C local variables

      INTEGER KMBD,DCORVG,KMEL,KMID
      INTEGER I,IED,IVT1,IVT2
      DOUBLE PRECISION DPPT(2),D,D1,D2

      KMID = L(TRIA(OLMID))
      KMBD = L(TRIA(OLMBD))
      KMEL = L(TRIA(OLMEL))
      DCORVG = L(TRIA(OLCORVG))

C Is there an element given to accelerate the search?

      IF (IEL.NE.0) THEN
      
        D1 = -1
      
C Search along the edges of that element if we find a boundary edge.
C Use the fact that KMEL(2,.)=0 for boundary edges.
      
        DO I=0,NNVE-1

          IED = KWORK(KMID+NNVE*(IEL-1)+I)

          IF (KWORK(KMEL+2*IED+1).EQ.0) THEN

C What are the adjacent nodes?
 
            CALL ED2NDS (TRIA,IED,IVT1,IVT2)
          
C Project the point to that edge
          
            CALL PPO2SG(X0,Y0,
     *         DWORK(DCORVG+2*(IVT1-1)),DWORK(DCORVG+2*(IVT2-1)),
     *         D,DPPT)
          
            D2 = (X0-DPPT(1))**2+(Y0-DPPT(2))**2

C Are we nearer than before?

            IF ((D1.LT.0).OR.(D2.LT.D1)) THEN
C If we are nearer, save edge number (we translate later), point, distance          
              IEDG = IED
              XNEW = DPPT(1)
              YNEW = DPPT(2)
              D1 = D2
            END IF
          
          END IF
        END DO
        
        IF (D1.GE.0D0) THEN

C We successfully found a boundary edge on our current element.
C In IEDG we only saved the number of the edge, not the index.
C So finally make a search inside of the KMBD-array to search for the index
C of our edge. Use the search routine to find it:

          CALL SRBDNI (TRIA, IEDG+TRIA(ONVT), IEDG, I)

        END IF
      
      END IF

C We haven't found the point - fall back to standard linear 
C search strategy

      CALL BPRLSR(TRIA, X0, Y0, IBC, XNEW, YNEW, IEDG)

      END 

************************************************************************
* Reconstruct parameter value of one boundary node
* 
* This routine analyses the coordinates of a given vertex on the (real)
* boundary and reconstructs its parameter values. It works
* independently of the current parametrization and might be slow
* in special cases.
*
* The routine accepts different parameters, characterizing the
* boundary node. It also accepts an optional parameter characterizing
* the cell that contains the node. If this parameter is given,
* the projection is done directly; otherwise there is a linear
* search on the boundary to find the segment closest to the point.
*
* For now this works only if the grid fulfills the following
* conditions:
* - the line segments are [0,1]-parametrised with equally
*   distributed parametrisation (i.e. something like (1-t^2)*P1+t^2*P2
*   is not allowed as parametrisation of a line)
*
* The routine does not use the internal structure of the
* parametrization, thus allowing lines, circles and arbitrary arcs
* to be used. 
*
* In:
*   TRIA   - array [1..SZTRIA] of integer
*            Structure specifying the triangulation and boundary
*            information.
*   IBND   - index in KVBD-array to the boundary node whose parameter
*            value should be computed from its coordinate
*            (this is an index, not the vertex number itself!)
*   IBCT   - number of the boundary component that contains KVBD(IBND).
*            Can be 0 if arbitrary, then the boundary component is
*            searched for automatically.
*   DTMAX  - Maximum parameter value on boundary component IBCT
*   IEL    - element number that contains node KVBD(IBND) or 0,
*            if the element is not known
*            
*   X,
*   Y      - Current coordinates of the boundary node IVT, which should
*            be corrected
*
* Out:
*   X,
*   Y      - Corrected coordinates of the boundary node
*   DPARM  - Corrected parameter value of the boundary node
*   IBCT   - If IBCT=0 on call of this routine, IBCT will be set to the
*            number of the boundary component that contains (X,Y) with
*            parameter value DPARM. Otherwise IBCT is not changed.
*   IELNEW - Number of the element that contains the corrected 
*            point (X,Y)
*
* The TRIA-structure itself is not changed. If there's an error in the
* reconstruction, X=Y=0, DPARM=-1 is returned.
*
* It's allowed to let IELNEW and IEL point to the same variable. in
* that case the element number will be updated if necessary.
************************************************************************

      SUBROUTINE RXBNPT(TRIA,IBND,IBCT,DTMAX,IEL, X,Y,DPARM,IELNEW)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'

      INTEGER TRIA(SZTRIA), IBND, IBCT, IEL, IELNEW
      DOUBLE PRECISION X,Y,DPARM,DTMAX

C local variables

      INTEGER DVBDP, DCORVG, KVBD, IVT, IED, IN1, IN2,I, IVT1, IVT2, NVT
      INTEGER IIED, J, IBC, IBC1, IBC2, IBCTMP
      INTEGER KMID, KMEL, KEBD,KBCT
      DOUBLE PRECISION MISG, MXSG, X0, Y0
      DOUBLE PRECISION DPPT(2),D,VD,DST,DP1,DP2
      
      DOUBLE PRECISION TMAX
      EXTERNAL TMAX
      
      DVBDP = L(TRIA(OLVBDP))
      DCORVG = L(TRIA(OLCORVG))
      
      KVBD = L(TRIA(OLVBD))
      KEBD = L(TRIA(OLEBD))

      KMID = L(TRIA(OLMID))
      KMEL = L(TRIA(OLMEL))
      KBCT = L(TRIA(OLBCT))
      DCORVG = L(TRIA(OLCORVG))
      DVBDP = L(TRIA(OLVBDP))
      NVT = TRIA(ONVT)
      
C Vertex and edge number of boundary node?

      IVT = KWORK(KVBD+IBND-1)
      IED = KWORK(KEBD+IBND-1)
      
C Remember the (non-projected) coordinates of our boundary point:

      X0 = X
      Y0 = Y
      
C Initialize the output variables
      
      X=0D0
      Y=0D0
      DPARM=-1D0
      
C What are the bounding parameter values of the parent segment that
C contains node IVT? We need this to ensure that the point does
C not leave the parent segment. The segments are guided by being
C integer values.

      IF (IBND.NE.0) THEN
        MISG = DINT(DWORK(DVBDP+IBND-1))
        MXSG = MISG+1D0
      ELSE
        MISG = 0D0
        MXSG = 0D0
      END IF
      
C Is there an element given to accelerate the search?

      IF (IEL.GT.0) THEN
      
C Search along the edges of that element if we find a boundary edge.
C Use the fact that KMEL(2,.)=0 for boundary edges.
C Furthermore we only search for boundary edges on the same
C boundary parent segment as the original point was. This is to prevent
C the point of "flipping around a guiding corner".
C We remark that there's at most one edge on every element belonging
C to such a boundary parent segment: If more that one edge of one element
C belongs to the boundary and the element is not degenerated,
C the adjacent note between the edges is integer-valued, so the
C edges belong to different boundary parent segments.
      
        DO I=0,NNVE-1

          IED = KWORK(KMID+NNVE*(IEL-1)+I)-NVT

          IF (KWORK(KMEL+2*(IED-1)+1).EQ.0) THEN

C What are the adjacent nodes?
 
            CALL ED2NDS (TRIA,IED,IVT1,IVT2)
          
C and the corresponding indices in the KVBD-array?

            CALL SRBDNI (TRIA, IVT1, IN1, J)
            CALL SRBDNI (TRIA, IVT2, IN2, J)

C Calculate the parameter values of these both nodes. The parameter
C value of IVT1 can be calculated directly:

            DP1 = DWORK(DVBDP+IN1-1)
            
C For the parameter value of IVT2 we check, whether IN2 >= IN1.
C In this case we can obtain DP2 directly. But if IN2 < IN1,
C IVT2 is the starting point of the boundary component, thus
C we have to calculate the parameter value by taking that of
C IVT1 and and round it up to the next higher integer

            IF (IN2.GE.IN1) THEN
              DP2 = DWORK(DVBDP+IN2-1)
            ELSE
              DP2 = DINT(DWORK(DVBDP+IN1-1)) + 1D0
            END IF
            
C Calculate a parameter value of any point between those
C two boundary points. This helps to get rid of routing
C errors when testing if this edge belongs to the current 
C parent segment.

            D = 0.5D0*(DP1+DP2)
            
C Is the edge on the allowed parent segment?

            IF ((D.GE.MISG).AND.(D.LE.MXSG)) THEN
            
C Very well. Calculate the projection onto that edge, or line, resp.:
            
              CALL PPO2SL(X0,Y0,
     *         DWORK(DCORVG+2*(IVT1-1)),DWORK(DCORVG+2*(IVT2-1)),
     *         D,DPPT)
     
C The value D (range 0..1) denotes the fractional part of the new
C parameter value, that interpolates the parameter values of IVT1
C and IVT2.
C
C D=0 or D=1 is allowed since this indicates that the point did
C not move.
C The only problem is if D is < 0 or > 1. In this case the
C projection leaves our current edge.
C This happens only on rather deformed elements, when the projection
C of an inner point to the boundary of the domain leaves that
C element:
C
C                ________________
C               /           ____/
C   |------------------>2_ /
C   |      _/      ___/ |
C   |     /       /     \/
C  -1----x-------x------3--------------
C
C In this case we have to fall back to linear search.
C Otherwise (i.e. if the projection of the inner point to the
C domain boundary hits the boundary of that element) we can
C directly compute the parameter-value and coordinates
C using D and MISG:

              IF ((D.GE.0D0).AND.(D.LE.1D0)) THEN
                DPARM = D*DP2 + (1D0-D)*DP1
                X = DPPT(1)
                Y = DPPT(2)
                IELNEW = IEL
                RETURN
              END IF

            END IF
          
          END IF
        END DO

      END IF
      
C We don't know anything about which edge the point could be
C projected to. We therefore have to search for the edge with
C the minimum distance. This is done by a linear search about
C all edges of the current boundary component. We only have to
C make sure that we don't leave our boundary parent segment.

      DST = -1D0

C If IBCT=0, we have to look into all boundary components, otherwise
C only into the IBCT'th:

      IF (IBCT.NE.0) THEN 
        IBC1 = IBCT
        IBC2 = IBCT
      ELSE
        IBC1 = 1
        IBC2 = TRIA(ONBCT)
        IBCTMP = 0
      END IF

C Loop over all boundary components we have to consider. Normally this
C is only one component, but maybe all in case of a brute force search.

      DO IBC = IBC1,IBC2

C If no boundary segment is given, we allow all boundary parameter
C values in the current boundary component to be used in the calculation
C of the minimum distance. MISG is already set to 0D0 above...
C --> Undocumented feature! Still has to be tested.

        IF (IBND.EQ.0) THEN

C         Determine the maximum parameter value of this boundary
C         component from the parametrization.

C          IIED = KWORK(KBCT+IBC)-1-1
C          MXSG = DINT(DWORK(DVBDP+IIED))+1D0

           MXSG = TMAX(IBC)
        END IF
 
C Loop about all boundary segments in our current boundary component.
C Whenever we find a better projection, take it. DST saves the current
C best distance.

        DO IIED = KWORK(KBCT+IBC-1),KWORK(KBCT+IBC)-1
        
C The current edge IIED starts at parameter value DVBCP(IIED)
C and ends at DVBCP(IIED+1) - except for if IIED is the last
C edge, that ends at DVBDP(first node of b.c.)!

          IF (IIED.NE.KWORK(KBCT+IBC)-1) THEN
              
            IVT1 = KWORK(KVBD+IIED-1)
            IVT2 = KWORK(KVBD+IIED)

C Save the parameter values of the endpoints

            DP1 = DWORK(DVBDP+IIED-1)
            DP2 = DWORK(DVBDP+IIED)
          
C Calculate a parameter value of any point between those
C two boundary points. this helps to get rid of routing
C errors when testing if this edge belongs to the current 
C parent segment.

          ELSE
          
C Last segment. The endpoint is the starting point.
          
            IVT1 = KWORK(KVBD+IIED-1)
            IVT2 = KWORK(KVBD+KWORK(KBCT+IBC-1)-1)

C Save the parameter values of the twi endpoints.
C To calculate the parameter value of the endpoint of the
C line segment, take TMAX and add the displacement from 0
C to that. Normally this gives TMAX (so the maximum parameter
C value of the boundary component), but in case that the 
C first point of the boundary was moved (should never happen),
C this will give the correct parameter value.

            DP1 = DWORK(DVBDP+IIED-1)
C            DP2 = DINT (DWORK(DVBDP+IIED-1)) + 1D0
            DP2 = DWORK(DVBDP+ KWORK(KBCT+IBC-1)-1) + TMAX(IBC)
            
C For line segments, DP2 should now be TMAX, so integer-based.
C For circles this value might be non-integer and > 1 !

          END IF

C Calculate the parameter value of the midpoint - or of any point
C of this edge that is not one of the endpoints

          D = DMOD(0.5D0*(DP1+DP2) , TMAX(IBC))
          
C Is the edge on the allowed parent segment?

          IF ((D.GE.MISG).AND.(D.LE.MXSG)) THEN
            
C Very well. Calculate the projection onto that edge:
            
            CALL PPO2SG(X0,Y0,
     *           DWORK(DCORVG+2*(IVT1-1)),DWORK(DCORVG+2*(IVT2-1)),
     *           D,DPPT)
     
C DPPT contains the projected point, D contains the parameter value 
C of the point relative to the line given by the endpoints IVT1
C and IVT2.
     
C What's the distance to our vertex?

            VD = (DPPT(1)-X0)**2+(DPPT(2)-Y0)**2
          
            IF ((DST.LT.0D0).OR.(VD.LT.DST)) THEN
          
C We found a better projection than the previous. Save the new
C data about the projected point:

              DPARM = D*DP2 + (1D0-D)*DP1
              
C The parameter value must be at most TMAX!
C If it's larger (might happen on circles if DP2>1), correct it
C to be < TMAX.

              DPARM = DMOD(DPARM,TMAX(IBC))
              
              X = DPPT(1)
              Y = DPPT(2)
              DST = VD
              IBCTMP = IBC
          
C Also determine the element that contains the new point

              IELNEW = KWORK(KEBD+IIED-1)
          
            END IF
     
          END IF
      
        END DO
        
      END DO
      
C Return the number of the boundary component if necessary

      IF (IBCT.EQ.0) IBCT=IBCTMP

C Make sure, DPARM is in allowed bounds
      
      IF (DPARM.GE.0D0) DPARM = DMOD(DPARM,DTMAX)
      
      END
