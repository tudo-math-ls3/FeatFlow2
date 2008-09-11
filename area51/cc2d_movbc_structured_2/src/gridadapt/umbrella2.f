************************************************************************
* Laplacian grid smoothing
*
* This routine performs Laplacian smoothing on the grid DCORVG,
* on inner nodes as well as on boundary nodes.
* The parameter values of boundary nodes are modified accordingly
* to the movement of boundary nodes. The parameters and coordinates
* of edge midpoints or any other nodes (DCORMG) are not touched,
* neighter of inner nodes nor on the boundary.
*
* The routine is designed to work with an extended grid, thus
* accepting the KXNPR-array instead of the KNPR-array for testing
* which nodes are boundary nodes and which not. 
*
* Boundary nodes with integer-valued parameters are not moved.
*
* In:
*   TRIA   - array [1..SZTRIA] of integer
*            Triangulation structure of the grid that should be
*            smoothed
*
* Out:
*   TRIA   - Smoothed grid. DCORVG and DBDPD are changed.
************************************************************************

      SUBROUTINE GSMTH5 (TRIA)

      IMPLICIT NONE

      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      
      INCLUDE 'stria.inc'
      
C parameters

      INTEGER TRIA(SZTRIA)
      
C local variables

      INTEGER LNCNT,LDCOPY,IED,IVT1,IVT2,KEAN,KCORVG,KDCOPY,KNCNT,IVT
      DOUBLE PRECISION DX1,DX2,DY1,DY2,D,W,DX,DY,TX,TY
      INTEGER LVBDP2,KVBDP,KVBDP2,KBCT,KVBD,IVBD,NVBD,IBCT,NPR,KXNPR
      INTEGER IEL,IELOLD,KVEL,NVEL,IVTPR
      DOUBLE PRECISION DPAR

C externals

      DOUBLE PRECISION PARX,PARY,TMAX
      EXTERNAL PARX,PARY,TMAX
      
      INTEGER BSRCH4
      EXTERNAL BSRCH4

C The Laplacian smoothing operator is defined as:
C
C       Pnew = (1-ALPHA)*P + ALPHA/#neighbours * SUM(Qj)
C
C with Qj being the neighbours of the point P. By setting
C ALPHA=0.5, we obtain the formula:
C
C       Pnew = 1/(2*#neighbours) * ( P*#neighbours +  SUM(Qj) )
C
C This is the formula implemented here.
      
C Small check:
      
      IF (TRIA(OTRIFL).EQ.0) RETURN
      
C Allocate a counter-array, counting how often a node was touched

      CALL ZNEW(TRIA(ONVT),3,LNCNT,'LNCNT ')
      IF (IER.NE.0) RETURN
      
C The new coordinates are written indo a backup-DCORVG; create this,
C filled with 0.
      
      CALL ZNEW(2*TRIA(ONVT),1,LDCOPY,'DCOPY ')
      IF (IER.NE.0) RETURN
      
C Reserve array for building new parameter values of boundary nodes.
      
      NVBD = TRIA(ONVBD)
      
      CALL ZNEW(NVBD,1,LVBDP2,'DVBD2')
      
C Some preparation:
      
      KXNPR  = L(TRIA(OLXNPR))
      KEAN  = L(TRIA(OLEAN))
      KCORVG = L(TRIA(OLCORVG))
      KDCOPY = L(LDCOPY)      
      KNCNT  = L(LNCNT)
      KVEL   = L(TRIA(OLVEL))
      NVEL   = TRIA(ONVEL)

      KVBDP2 = L(LVBDP2)
      KVBDP  = L(TRIA(OLVBDP))
      KVBD   = L(TRIA(OLVBD))
      KBCT   = L(TRIA(OLBCT))
      KCORVG = L(TRIA(OLCORVG))
      
C In the first step only smooth the inner-vertex coordinates.
C Make a loop about all edges:

      DO IED = 0,TRIA(ONMT)-1
      
C       On each edge find the two adjacent nodes
      
        IVT1 = KWORK(KEAN+2*IED)-1
        IVT2 = KWORK(KEAN+2*IED+1)-1
        
        DX1 = DWORK(KCORVG+2*IVT1)
        DY1 = DWORK(KCORVG+2*IVT1+1)
        DX2 = DWORK(KCORVG+2*IVT2)
        DY2 = DWORK(KCORVG+2*IVT2+1)
        
C       Add the coordinates of each node
C       - to itself and
C       - to the neighbour node
C       That way each node is weighted typically 4x to itself 
C       (or more, if there are more adjacent edges) and 1x to
C       the neighbour:

        DWORK(KDCOPY+2*IVT1)   = DWORK(KDCOPY+2*IVT1)   + DX1 + DX2 
        DWORK(KDCOPY+2*IVT1+1) = DWORK(KDCOPY+2*IVT1+1) + DY1 + DY2
        
        DWORK(KDCOPY+2*IVT2)   = DWORK(KDCOPY+2*IVT2)   + DX1 + DX2 
        DWORK(KDCOPY+2*IVT2+1) = DWORK(KDCOPY+2*IVT2+1) + DY1 + DY2
        
C       Increase the counter of both nodes:

        KWORK(KNCNT+IVT1) = KWORK(KNCNT+IVT1)+2
        KWORK(KNCNT+IVT2) = KWORK(KNCNT+IVT2)+2
        
C       Special case: If our vertex is a boundary vertex with only 3
C        neighbours, and if the current edge points to that single
C        neighbour inside the domain, we weight this edge 2x onto the
C        boundary node.
C        This way we artificially produce a 4th neighbour to that
C        boundary node. This prevents a "bellied" shape near the boundary
C        due to the lack of neighbours!
C      
C       At first make sure that the boundary point is IVT1
        IVT = -1
        IF ( ( IAND(KWORK(KXNPR+2*IVT1),2**8).NE.0 ).AND.
     *       ( IAND(KWORK(KXNPR+2*IVT2),2**8).EQ.0 ) ) THEN
         IVT = IVT1
        END IF
        IF ( ( IAND(KWORK(KXNPR+2*IVT1),2**8).EQ.0 ).AND.
     *       ( IAND(KWORK(KXNPR+2*IVT2),2**8).NE.0 ) ) THEN
         IVT = IVT2
         IVT2 = IVT1
         IVT1 = IVT
         DX = DX2
         DX2 = DX1
         DX1 = DX
         DY = DY2
         DY2 = DY1
         DY1 = DY
        END IF
        
C       Note that there are only three adjacent edges to the boundary
C       vertex IVT if there are exactly 2 adjacent elements!

        IF ((IVT.GE.0).AND.(KWORK(KVEL+IVT*NVEL+1).NE.0).AND.
     *      (KWORK(KVEL+IVT*NVEL+2).EQ.0)) THEN
        
C         Don't simply add the coordinates of the inner vertex again.
C         This works on simple domains, but on more complex domains
C         it will move the new point too much into the inner, so the
C         projection might be hard.
C         Instead mirror the inner point on the edge to the outside
C         of the domain and add the coordinates of that point:
C
C              IVTNEW
C                ^
C                |
C         O======|=IVT1=====IVTPR
C         |      | /          |
C         |      |/           |
C         |     IVT2          |
C        
C         At first we have to find an edge on the boundary where to
C         mirror IVT2. We take the edge between the predecessor of
C         IVT1 on the boundary and IVT1. We know that there's a pre-
C         decessor, because otherwise IVT1 would be a "corner" vertex
C         on the boundary - those are fixed and only have
C         one adjacent element!
C         Find IVT1 in the list of boundary points:

          CALL SRBDNI (TRIA, IVT1+1, IVBD, IBCT)
          
C         Get the number of the previous vertex
C         (we could also take the following vertex, or depending on
C          the position of the projection the following or the
C          previous vertex - but that would be much more complicated)

          IVTPR = KWORK(KVBD+IVBD-2)-1
          
C         For a description on how to create IVTNEW,
C         see http://mathworld.wolfram.com/Reflection.html.
C     
C         Build the tangential vector of that edge, normalized 
C         to length 1:

          TX = DX1-DWORK(KCORVG+2*IVTPR)
          TY = DY1-DWORK(KCORVG+2*IVTPR+1)
          W=1D0/SQRT(TX*TX+TY*TY)
          TX=W*TX
          TY=W*TY
        
C         build the vector pointing from IVT1 to the inner vertex
        
          DX = DX2-DX1
          DY = DY2-DY1
          
C         And build the reflection:

          W = DX*TX + DY*TY
          DX2 = -DX2 + 2D0*DX1 + 2D0*TX*W
          DY2 = -DY2 + 2D0*DY1 + 2D0*TY*W

C         Add that point to DCOPY as usual
        
          DWORK(KDCOPY+2*IVT)   = DWORK(KDCOPY+2*IVT)   + DX1 + DX2 
          DWORK(KDCOPY+2*IVT+1) = DWORK(KDCOPY+2*IVT+1) + DY1 + DY2
          KWORK(KNCNT+IVT)      = KWORK(KNCNT+IVT)+2
        
        END IF
      
      END DO
      
C Now correct the coordinates with the help of the counter-array KNCNT.
C Remember that each entry in DCOPY was calculated by summing up multiple
C coordinates of the surrounding nodes - and now we have to divide by
C the number of nodes we added together.
C Only handle those nodes here which are not (real) boundary nodes!

      DO IVT=0,TRIA(ONVT)-1

        NPR = KWORK(KXNPR+2*IVT)
        IF ( IAND(NPR,2**8).EQ.0 ) THEN

          D = 1D0/DBLE(KWORK(KNCNT+IVT))
          
          DWORK(KDCOPY+2*IVT)   = DWORK(KDCOPY+2*IVT)   * D
          DWORK(KDCOPY+2*IVT+1) = DWORK(KDCOPY+2*IVT+1) * D

        END IF
        
      END DO

C Now handle boundary nodes. Build the parameter values in the 
C temporary array and overwrite the original array later.
C      
C Loop over the boundary components
      
      DO IBCT=1,TRIA(ONBCT)
      
C Loop over the points in this boundary component

        DO IVBD=KWORK(KBCT+IBCT-1),KWORK(KBCT+IBCT)-1
        
C Don't move points with integer parameter value because these are "fixed".
        
          DPAR=DWORK(KVBDP+IVBD-1)
          
          IF (DPAR.NE.AINT(DPAR)) THEN
          
C           As in the case of the non-boundary nodes, calculate
C           the coordinates of the point. This might leave the
C           boundary!

            IVT = KWORK(KVBD+IVBD-1)-1

            D = 1D0/DBLE(KWORK(KNCNT+IVT))
            
            DX = DWORK(KDCOPY+2*IVT)   * D
            DY = DWORK(KDCOPY+2*IVT+1) * D

C           Use the search-routine to try to determine the element this
C           point has been moved to:

            IELOLD = KWORK(KVEL+IVT*NVEL)

            IF (BSRCH4(DX,DY,DWORK(KCORVG),
     *             KWORK(L(TRIA(OLVERT))),KWORK(L(TRIA(OLADJ))),
     *             IELOLD,IEL).NE.0) THEN
     
C             If we can't find the element, set IEL to 0 and perform later
C             linear search along the boundary for finding the element 
C             containing the vertex 
        
              IEL = 0
              
            END IF
            
C           Finally use the projection routine to project the point back
C           and calculate a new parameter value:
        
            CALL RXBNPT(TRIA,IVBD,IBCT,TMAX(IBCT),IEL,DX,DY,DPAR,IEL)
        
C           Save the parameter value in DVBDP2 and the coordinates
C           in the DCOPY array, overwriting the previously calculated
C           coordinates:

            DWORK(KVBDP2+IVBD-1) = DPAR
            
            DWORK(KDCOPY+2*IVT)   = DX
            DWORK(KDCOPY+2*IVT+1) = DY
            
          ELSE
          
            DWORK(KVBDP2+IVBD-1) = DPAR
            
          END IF
          
C DPAR now contains the current, perhaps corrected parameter value.
C Recalculate the X/Y-coordinates of that point, regardless of if it
C was moved or not. 
        
          IVT = KWORK(KVBD+IVBD-1) - 1
          
          DWORK(KDCOPY+2*IVT)   = PARX(DPAR,IBCT)
          DWORK(KDCOPY+2*IVT+1) = PARY(DPAR,IBCT)
        
        END DO
      
      END DO
      
C We are done; activate the new coordinate and parameter set

      CALL LCP1(DWORK(KDCOPY),DWORK(KCORVG),2*TRIA(ONVT))
      CALL LCP1(DWORK(KVBDP2),DWORK(KVBDP),TRIA(ONVBD))
      
C Release memory, finish
      
      CALL ZDISP(0,LVBDP2,'DVBD2')
      CALL ZDISP(0,LDCOPY,'DCOPY ')
      CALL ZDISP(0,LNCNT,'KCOUNT')
        
      END

************************************************************************
* Umbrella grid smoothing
*
* This routine performs Umbrella smoothing on the grid DCORVG,
* on inner nodes as well as on boundary nodes.
* The parameter values of boundary nodes are modified accordingly
* to the movement of boundary nodes. The parameters and coordinates
* of edge midpoints or any other nodes (DCORMG) are not touched,
* neighter of inner nodes nor on the boundary.
*
* The routine is designed to work with an extended grid, thus
* accepting the KXNPR-array instead of the KNPR-array for testing
* which nodes are boundary nodes and which not. 
*
* Boundary nodes with integer-valued parameters are not moved.
*
* In:
*   TRIA   - array [1..SZTRIA] of integer
*            Triangulation structure of the grid that should be
*            smoothed
*
* Out:
*   TRIA   - Smoothed grid. DCORVG and DBDPD are changed.
************************************************************************

      SUBROUTINE GSMTH6 (TRIA)

      IMPLICIT NONE

      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      
      INCLUDE 'stria.inc'
      
C parameters

      INTEGER TRIA(SZTRIA)
      
C local variables

      INTEGER LNLEN,LDCOPY,IED,IVT1,IVT2,KEAN,KCORVG,KDCOPY,KLEN,IVT
      DOUBLE PRECISION DX1,DX2,DY1,DY2,D,W,DX,DY,TX,TY
      INTEGER LVBDP2,KVBDP,KVBDP2,KBCT,KVBD,IVBD,NVBD,IBCT,NPR,KXNPR
      INTEGER KVEL,NVEL,IVTPR,IELOLD,IEL
      DOUBLE PRECISION DPAR

C externals

      DOUBLE PRECISION PARX,PARY,TMAX
      EXTERNAL PARX,PARY,TMAX
      
      INTEGER BSRCH4
      EXTERNAL BSRCH4
      
C The density weighted umbrella operator is defined as:
C
C  Pnew = (1-ALPHA) * P + ALPHA/SUM(|P-Qj|) * SUM(|P-Qj| * Qj)
C
C with Qj being the neighbours of the point P. Therefore the constant
C ALPHA defines the weighting of the midpoint P of the "patch".
C A typical choice is ALPHA=0.5, which leads to the formula:
C
C  Pnew = 1/(2*SUM(|P-Qj|)) * ( P*SUM(|P-Qj|) + SUM(|P-Qj| * Qj) )
C
C This is the formula implemented here.

C Small check:
      
      IF (TRIA(OTRIFL).EQ.0) RETURN
      
C Allocate a collector-array, collecting the sum of the length
C of all edges adjacent to each node:

      CALL ZNEW(TRIA(ONVT),1,LNLEN,'DLEN  ')
      IF (IER.NE.0) RETURN
      
C The new coordinates are written indo a backup-DCORVG; create this,
C filled with 0.
      
      CALL ZNEW(2*TRIA(ONVT),1,LDCOPY,'DCOPY ')
      IF (IER.NE.0) RETURN
      
C Some preparation:
      
      KXNPR  = L(TRIA(OLXNPR))
      KEAN  = L(TRIA(OLEAN))
      KCORVG = L(TRIA(OLCORVG))
      KDCOPY = L(LDCOPY)      
      KLEN  = L(LNLEN)
      KVEL   = L(TRIA(OLVEL))
      NVEL   = TRIA(ONVEL)

C Now handle boundary nodes. Build the new parameter values in a
C separate array and overwrite the old array later.
      
      NVBD = TRIA(ONVBD)
      
      CALL ZNEW(NVBD,1,LVBDP2,'DVBD2')
      
      KVBDP2 = L(LVBDP2)
      KVBDP  = L(TRIA(OLVBDP))
      KVBD   = L(TRIA(OLVBD))
      KBCT   = L(TRIA(OLBCT))
      KCORVG = L(TRIA(OLCORVG))
      
C In the first step only smooth the inner-vertex coordinates.
C Make a loop about all edges:

      DO IED = 0,TRIA(ONMT)-1
      
C       On each edge find the two adjacent nodes
      
        IVT1 = KWORK(KEAN+2*IED)-1
        IVT2 = KWORK(KEAN+2*IED+1)-1
        
        DX1 = DWORK(KCORVG+2*IVT1)
        DY1 = DWORK(KCORVG+2*IVT1+1)
        DX2 = DWORK(KCORVG+2*IVT2)
        DY2 = DWORK(KCORVG+2*IVT2+1)
        
C       Add the coordinates of each node to the neighbour, weighted by
C       the length of the edge.

        D = SQRT((DX2-DX1)**2+(DY2-DY1)**2)

        DWORK(KDCOPY+2*IVT1)   = DWORK(KDCOPY+2*IVT1)   + D*DX2 
        DWORK(KDCOPY+2*IVT1+1) = DWORK(KDCOPY+2*IVT1+1) + D*DY2
        
        DWORK(KDCOPY+2*IVT2)   = DWORK(KDCOPY+2*IVT2)   + D*DX1 
        DWORK(KDCOPY+2*IVT2+1) = DWORK(KDCOPY+2*IVT2+1) + D*DY1 
        
C       Add the length of the edge to the total length of all edges
C       adjacent to each node:

        DWORK(KLEN+IVT1) = DWORK(KLEN+IVT1) + D
        DWORK(KLEN+IVT2) = DWORK(KLEN+IVT2) + D
      
C       Special case: If our vertex is a boundary vertex with only 3
C        neighbours, and if the current edge points to that single
C        neighbour inside the domain, we weight this edge 2x onto the
C        boundary node.
C        This way we artificially produce a 4th neighbour to that
C        boundary node. This prevents a "bellied" shape near the boundary
C        due to the lack of neighbours!
C      
C       At first make sure that the boundary point is IVT1
        IVT = -1
        IF ( ( IAND(KWORK(KXNPR+2*IVT1),2**8).NE.0 ).AND.
     *       ( IAND(KWORK(KXNPR+2*IVT2),2**8).EQ.0 ) ) THEN
         IVT = IVT1
        END IF
        IF ( ( IAND(KWORK(KXNPR+2*IVT1),2**8).EQ.0 ).AND.
     *       ( IAND(KWORK(KXNPR+2*IVT2),2**8).NE.0 ) ) THEN
         IVT = IVT2
         IVT2 = IVT1
         IVT1 = IVT
         DX = DX2
         DX2 = DX1
         DX1 = DX
         DY = DY2
         DY2 = DY1
         DY1 = DY
        END IF
        
C       Note that there are only three adjacent edges to the boundary
C       vertex IVT if there are exactly 2 adjacent elements!

        IF ((IVT.GE.0).AND.(KWORK(KVEL+IVT*NVEL+1).NE.0).AND.
     *      (KWORK(KVEL+IVT*NVEL+2).EQ.0)) THEN
        
C         Don't simply add the coordinates of the inner vertex again.
C         This works on simple domains, but on more complex domains
C         it will move the new point too much into the inner, so the
C         projection might be hard.
C         Instead mirror the inner point on the edge to the outside
C         of the domain and add the coordinates of that point:
C
C              IVTNEW
C                ^
C                |
C         O======|=IVT1=====IVTPR
C         |      | /          |
C         |      |/           |
C         |     IVT2          |
C        
C         At first we have to find an edge on the boundary where to
C         mirror IVT2. We take the edge between the predecessor of
C         IVT1 on the boundary and IVT1. We know that there's a pre-
C         decessor, because otherwise IVT1 would be a "corner" vertex
C         on the boundary - those are fixed and only have
C         one adjacent element!
C         Find IVT1 in the list of boundary points:

          CALL SRBDNI (TRIA, IVT1+1, IVBD, IBCT)
          
C         Get the number of the previous vertex
C         (we could also take the following vertex, or depending on
C          the position of the projection the following or the
C          previous vertex - but that would be much more complicated)

          IVTPR = KWORK(KVBD+IVBD-2)-1
          
C         For a description on how to create IVTNEW,
C         see http://mathworld.wolfram.com/Reflection.html.
C     
C         Build the tangential vector of that edge, normalized
C         to length 1D0:

          TX = DX1-DWORK(KCORVG+2*IVTPR)
          TY = DY1-DWORK(KCORVG+2*IVTPR+1)
          W=1D0/SQRT(TX*TX+TY*TY)
          TX=W*TX
          TY=W*TY
        
C         build the vector pointing from IVT1 to the inner vertex
        
          DX = DX2-DX1
          DY = DY2-DY1
          
C         And build the reflection:

          W = DX*TX + DY*TY
          DX2 = -DX2 + 2D0*DX1 + 2D0*TX*W
          DY2 = -DY2 + 2D0*DY1 + 2D0*TY*W
        
C         Add that point to DCOPY as usual
        
          DWORK(KDCOPY+2*IVT1)   = DWORK(KDCOPY+2*IVT1)   + D*DX2 
          DWORK(KDCOPY+2*IVT1+1) = DWORK(KDCOPY+2*IVT1+1) + D*DY2
          DWORK(KLEN+IVT1)       = DWORK(KLEN+IVT1)+D
        
        END IF
      
      END DO
      
C Now correct the coordinates in DCOPY. Add the coordinates of the
C midpoint of each "patch" to the sum of the neighbours, weighted by
C the length of all adjacent edges. Afterwards divide everything by 
C 2xlength of all adjacent edges.
C This is done at first only for the inner nodes, bondary nodes are
C handled later!

      DO IVT=0,TRIA(ONVT)-1

        NPR = KWORK(KXNPR+2*IVT)
        IF ( IAND(NPR,2**8).EQ.0 ) THEN

          D = 0.5D0/DWORK(KLEN+IVT)
          
          DWORK(KDCOPY+2*IVT)   = D * 
     *      ( DWORK(KLEN+IVT)*DWORK(KCORVG+2*IVT) + 
     *        DWORK(KDCOPY+2*IVT) )
          DWORK(KDCOPY+2*IVT+1) = D * 
     *      ( DWORK(KLEN+IVT)*DWORK(KCORVG+2*IVT+1) + 
     *        DWORK(KDCOPY+2*IVT+1) )

        END IF
        
      END DO
   
C Now a special treatment for the boundary nodes.   
C Loop over the boundary components:
      
      DO IBCT=1,TRIA(ONBCT)
      
C Loop over the points in this boundary component

        DO IVBD=KWORK(KBCT+IBCT-1),KWORK(KBCT+IBCT)-1
        
C Don't move points with integer parameter value because these are "fixed".

          DPAR=DWORK(KVBDP+IVBD-1)
          
          IF (DPAR.NE.AINT(DPAR)) THEN
          
C           As in the case of the non-boundary nodes, calculate
C           the coordinates of the point. This might leave the
C           boundary!

            IVT = KWORK(KVBD+IVBD-1)-1

            D = 0.5D0/DWORK(KLEN+IVT)
            
            DX = D * 
     *        ( DWORK(KLEN+IVT)*DWORK(KCORVG+2*IVT) + 
     *          DWORK(KDCOPY+2*IVT) )
            DY = D * 
     *        ( DWORK(KLEN+IVT)*DWORK(KCORVG+2*IVT+1) + 
     *          DWORK(KDCOPY+2*IVT+1) )

C           Use the search-routine to try to determine the element this
C           point has been moved to:

            IELOLD = KWORK(KVEL+IVT*NVEL)

            IF (BSRCH4(DX,DY,DWORK(KCORVG),
     *             KWORK(L(TRIA(OLVERT))),KWORK(L(TRIA(OLADJ))),
     *             IELOLD,IEL).NE.0) THEN
     
C             If we can't find the element, set IEL to 0 and perform later
C             linear search along the boundary for finding the element 
C             containing the vertex 

              IEL = 0
              
            END IF
            
C           Finally use the projection routine to project the point back
C           and calculate a new parameter value:

            CALL RXBNPT(TRIA,IVBD,IBCT,TMAX(IBCT),IEL,DX,DY,DPAR,IEL)

C           Save the parameter value and the coordinates

            DWORK(KVBDP2+IVBD-1) = DPAR
            
            DWORK(KDCOPY+2*IVT)   = DX
            DWORK(KDCOPY+2*IVT+1) = DY
          
          ELSE
          
            DWORK(KVBDP2+IVBD-1) = DPAR
            
          END IF
          
C DPAR now contains the current, perhaps corrected parameter value.
C Recalculate the X/Y-coordinates of that point, regardless of if it
C was moved or not. 

          IVT = KWORK(KVBD+IVBD-1) - 1
          
          DWORK(KDCOPY+2*IVT)   = PARX(DPAR,IBCT)
          DWORK(KDCOPY+2*IVT+1) = PARY(DPAR,IBCT)

C !!!!!!!!!!!!! Routine does still not work for circles !!!!!!!!!!!
C Might be that a boundary node with parameter value 0.99 is moved
C to parameter 0.01 - this is ok. Unfortunately the KVBD/KEBD/...-
C arrays are yet not reordered according to the parameter value
C (the last node has to come to the front or vice versa), which
C quickly results in distorted grids because of wrong parameter
C values!

        END DO
      
      END DO
      
C We are done; activate the new coordinate and parameter set

      CALL LCP1(DWORK(KDCOPY),DWORK(KCORVG),2*TRIA(ONVT))
      CALL LCP1(DWORK(KVBDP2),DWORK(KVBDP),TRIA(ONVBD))
      
C Release memory, finish
      
      CALL ZDISP(0,LVBDP2,'DVBD2')
      CALL ZDISP(0,LDCOPY,'DCOPY ')
      CALL ZDISP(0,LNLEN,'DLEN  ')

      END

************************************************************************
* Laplacian smoothing of scalar function
*
* This routine performs Laplacian smoothing of a scalar function
* that is defined as function values in the vertices of a grid.
*
* In:
*   TRIA   - array [1..SZTRIA] of integer
*            Triangulation structure of the grid
*   NEQ    - Length of function vector; NEQ >= NVT!
*   DU     - array [1..NEQ] of double
*            Scalar function that should be smoothed.
*            The first NVT elements in DU are identified with the
*            NVT corner vertices of the grid. Only these elements
*            are smoothed.
*
* Out:
*   DU     - smoothed function
************************************************************************

      SUBROUTINE GSMTH9 (TRIA,NEQ,DU)

      IMPLICIT NONE

      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      
      INCLUDE 'stria.inc'
      
C parameters

      INTEGER TRIA(SZTRIA),NEQ
      DOUBLE PRECISION DU(NEQ)
      
C local variables

      INTEGER LNCNT,LDCOPY,IED,IVT1,IVT2,KDCOPY,KNCNT,IVT

C Allocate a counter-array, counting how often a node was touched

      CALL ZNEW(TRIA(ONVT),3,LNCNT,'LNCNT ')
      IF (IER.NE.0) RETURN
      
C The new coordinates are written indo a backup-DCORVG; create this,
C filled with 0.
      
      CALL ZNEW(TRIA(ONVT),1,LDCOPY,'DCOPY ')
      IF (IER.NE.0) RETURN
      
C Some preparation:
      
      KDCOPY = L(LDCOPY)      
      KNCNT  = L(LNCNT)
      
C Make a loop about all edges:

      DO IED = 1,TRIA(ONMT)
      
C       On each edge find the two adjacent nodes
      
        CALL ED2NDS (TRIA,IED+1,IVT1,IVT2)
        
C       Add the value of the function in each node
C       - to itself and
C       - to the neighbour node
C       That way each node is weighted typically 4x to itself 
C       (or more, if there are more adjacent edges) and 1x to
C       the neighbour:

        DWORK(KDCOPY+IVT1-1) = DWORK(KDCOPY+IVT1-1) + DU(IVT1)+DU(IVT2)
        
        DWORK(KDCOPY+IVT2-1) = DWORK(KDCOPY+IVT2-1) + DU(IVT1)+DU(IVT2)
        
C       Increase the counter of both nodes:

        KWORK(KNCNT+IVT1-1) = KWORK(KNCNT+IVT1-1)+2
        KWORK(KNCNT+IVT2-1) = KWORK(KNCNT+IVT2-1)+2
      
      END DO
      
C Now copy the calculated values back to DU.
C Divide the coordinates by the number of how often each node was added
C something to:

      DO IVT=0,TRIA(ONVT)-1

        DU(1+IVT) = DWORK(KDCOPY+IVT) / DBLE(KWORK(KNCNT+IVT))

      END DO
      
C Release the backup- and counter-array again. 

      CALL ZDISP(0,LDCOPY,'DCOPY ')
      CALL ZDISP(0,LNCNT,'KCOUNT')
      
      END
