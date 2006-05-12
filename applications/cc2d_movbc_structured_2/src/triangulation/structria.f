***********************************************************************
* This file contains general triangulation routines to support the
* STRIA structure. The routines here define how the node numbering
* in the grid is organized.
*
* The standard treatment of the node numbering when DCORMG is not
* given is:
*   INODE = 1..NVT                  - vertices
*         = NVT+1..NVT+NMT          - edges; number associated to
*                                     edge midpoints
*         = NVT+NMT+1..NVT+NMT+NEL  - element midpoints
*
* If DCORMG is given, the node numbering is extended in the following
* way:
*  - 1..NVT 
*    represent all corner vertices.
*  - NVT+1..NVT+NVEDT 
*    represent all regularly distributed vertices on edges.
*  - NVT+NVEDT+1..NVT+NVEDT+NIEVT 
*    represent all regularly distributed vertices in the inner of elements.
*  - NVT+NVEDT+NIEVT..NVT+NVEDT+NIEVT+NANT 
*    represent all additional nodes that don't belong to regular 
*    distributed vertices on edges or inside of the elements.
* Note that the original numbering can be obtained if NVEDT=NMT and
* NIEVT=NEL is given and all coordinates are realized in DCORMG.
*
* The coordinates of the vertices are saved in DCORVG. The coordinates
* of vertices on edges and inner-element vertices are saved/expected in
* DCORMG.
*
* Advanced node numbering
* -----------------------
* The above numbering is not the completely correct one. More
* specifically for compatibility reasons, the node numbering is depending 
* on the values in NVEDT, NIEVT and NANT and defines some kind of
* "mixing" between the old style and the new.
* The general picture of the underlying vertex numbering is:
*
*   1..NVT 
*   --> always accesses corners
*
*   NVT + 1 .. NVT + MAX(NMT,NVEDT)
*   --> always accesses regularly distributed vertices on edges; either 
*       only edge-midpoints or general points, realized in DCORMG
*
*   NVT + MAX(NMT,NVEDT) + 1 .. NVT + MAX(NMT,NVEDT) + MAX(NEL,NIELVT)
*   --> always accesses regularly distributed inner-element vertices; 
*       either element midpoints or general inner element points,
*       realized in DCORMG
*
*   NVT + MAX(NMT,NVEDT) + MAX(NEL,NIELVT) + 1 .. 
*         NVT + MAX(NMT,NVEDT) + MAX(NEL,NIELVT) + NANT
*   --> always accesses additional vertices, realized in DCORMG
*
* This vertex numbering is possible because NVEDT and NIELVT are
* only defined if there is a regular distribution of points
* realized in DCORMG. Therefore NVEDT=i*NMT and NIELVT=j*NEL !
* Any additional vertices that do not fit into this form must
* be realized by the caller "by hand" in the DCORMG-array as
* additional vertices.
*
* There are standard converter-routines to convert "old-style"
* FEAT triangulations into "new-style" and back; these try to
* make a closest match, thus not introducing any additional points.
* If additional points should be added to DCORMG, the caller must
* modify DCORMG and NANT directly.
*
* Whether the structure is converted can be determined with
* the TRIFL-flag in the array:
*  TRIFL = 0  - standard FEAT style
*  TRIFL = 1  - extended style
* The extended routines here only use Bit0 of TRIFL to determine
* whether the structure is converted or not. Bit1-Bit7 are reserved
* for future use. Bit8-Bit31 can be user-defined.
***********************************************************************

************************************************************************
* Node to X/Y-coordinate
*
* This routine accepts a mesh TRIA and a node number INODE.
* The returned values are the X/Y-coodinate of the node, depending
* on the current implementation of the TRIA-structure.
*
* In:
*  INODE  - Node number of the node to check
*  TRIA   - array [1..SZTRIA] of integer
*           Triangulation structure STRIA.
*
* Out:
*  X,Y    - X/Y-coordinate of the node
************************************************************************

      SUBROUTINE NDE2XY (INODE,TRIA,X,Y)
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      INCLUDE 'stria.inc'
      INCLUDE 'cbasictria.inc'
      INTEGER INODE, TRIA(SZTRIA)
      DOUBLE PRECISION X,Y
      
      INTEGER KCORVG, KCORMG, KVERT, IVT1, IVT2, IVT3, IVT4, IEL
      INTEGER NVT,NVEDT,NIEVT,NANT
      
C Functions to take the mean of the X- or Y-coordinates of vertices
C in DCORVG
      
      INTEGER IXY,V1,V2,V3,V4,KVG
      DOUBLE PRECISION MID2,MID4
      MID2(KVG,V1,V2,IXY) = 0.5D0 * (DWORK(KVG+2*(V1-1)+IXY)+
     *                               DWORK(KVG+2*(V2-1)+IXY) )
      MID4(KVG,V1,V2,V3,V4,IXY) = 0.25D0 * (DWORK(KVG+2*(V1-1)+IXY)+
     *                                      DWORK(KVG+2*(V2-1)+IXY)+
     *                                      DWORK(KVG+2*(V3-1)+IXY)+
     *                                      DWORK(KVG+2*(V4-1)+IXY) )
      
      
C Standard implementation
C
C   INODE = 1..NVT                  - corner vertices
C           NVT+1 .. NVT+MAX(NMT,NVEDT)
C                                   - vertices on edges
C           NVT+MAX(NMT,NVEDT)+1 .. NVT+MAX(NMT,NVEDT)+MAX(NEL,NIELVT)
C                                   - vertices in elements
C           all other numbers       - additional vertices
C
C Load some variables for easier access:

      NVT   = TRIA(ONVT)
      NANT  = TRIA(ONANT)
      NVEDT = MAX(TRIA(ONMT),TRIA(ONVEDT)) 
      NIEVT = MAX(TRIA(ONEL),TRIA(ONIEVT))
      
      IF (INODE < 1) THEN      
      
        PRINT *,'NDE2XY: INODE=',INODE,' < 1 !?!'
        STOP
        
      ELSE IF (INODE.LE.NVT) THEN
      
C       Corner vertex
      
        KCORVG = L(TRIA(OLCORVG))
        X = DWORK(KCORVG+2*(INODE-1))
        Y = DWORK(KCORVG+2*(INODE-1)+1)
        RETURN

      ELSE IF (INODE.LE.NVT+NVEDT) THEN
      
C       Vertex on an edge. Can we get it from DCORMG or do we have
C       to calculate a midpoint coordinate?

        IF ((TRIA(OLCORMG).EQ.0).OR.(TRIA(ONVEDT).EQ.0)) THEN
        
C         Calculate edge midpoint

          CALL ED2NDS (TRIA,INODE-TRIA(ONVT),IVT1,IVT2)
        
          X = MID2 (L(TRIA(OLCORVG)),IVT1,IVT2,0)
          Y = MID2 (L(TRIA(OLCORVG)),IVT1,IVT2,1)
        
        ELSE
        
C         Get the coordinate from DCORMG
        
          KCORMG = L(TRIA(OLCORMG))
          X = DWORK(KCORMG+2*(INODE-NVT-1))
          Y = DWORK(KCORMG+2*(INODE-NVT-1)+1)
        
        END IF
        
      ELSE IF (INODE.LE.NVT+NVEDT+NIEVT) THEN
      
C       Vertex inside of an element; check if the element midpoint
C       is ment

        IF ((TRIA(OLCORMG).EQ.0).OR.(TRIA(ONIEVT).EQ.0)) THEN

C         We have to manually calculate the element midpoint.
C         INODE represents the number of an element:

          KVERT = L(TRIA(OLVERT))
          IEL = INODE-NVT-NVEDT-1
          IVT1 = KWORK(KVERT+IEL*NNVE)
          IVT2 = KWORK(KVERT+IEL*NNVE+1)
          IVT3 = KWORK(KVERT+IEL*NNVE+2)
          IVT4 = KWORK(KVERT+IEL*NNVE+3)
        
          X = MID4 (L(TRIA(OLCORVG)),IVT1,IVT2,IVT3,IVT4,0)
          Y = MID4 (L(TRIA(OLCORVG)),IVT1,IVT2,IVT3,IVT4,1)

        ELSE
        
C         INODE is the number of an inner-element node, realized in 
C         DCORMG. Get it's coordinate from DCORMG directly:
        
          KCORMG = L(TRIA(OLCORMG))
        
C         Don't use NVEDT in the next commands, but use TRIA(ONVEDT)
C         directly - as NVEDT might be NMT in case of TRIA(ONVEDT)=0!
          
          X = DWORK(KCORMG+2*(INODE-NVT-TRIA(NVEDT)-1))
          Y = DWORK(KCORMG+2*(INODE-NVT-TRIA(NVEDT)-1)+1)
        
        END IF

      ELSE IF ((TRIA(OLCORMG).NE.0).AND.
     *         (INODE.LE.NVT+NVEDT+NIEVT+NANT)) THEN
      
C       Number of an additional node, must be realized in DCORMG.
C       Don't use NVEDT/NIEVT in the next commands, but use TRIA(ONVEDT)
C       and TRIA(ONIEVT) directly! This will directly compute the index
C       in the DCORMG-array, anso if TRIA(ONVEDT) or TRIA(ONIEVT)
C       is 0!

        KCORMG = L(TRIA(OLCORMG))
        X = DWORK(KCORMG+2*(INODE-NVT-TRIA(NVEDT)-TRIA(ONIEVT)-1))
        Y = DWORK(KCORMG+2*(INODE-NVT-TRIA(NVEDT)-TRIA(ONIEVT)-1)+1)
      
      ELSE

        PRINT *,'NDE2XY: INODE=',INODE,': unknown node number!?!'
        STOP
      
      END IF

      END
      
************************************************************************
* Calculate total number of nodes
*
* This routine calculates depending on the current implementation of the
* mesh information the total number of nodes in the mesh.
* This is calculated by:
* - number of vertices
* - number of vertices on edges
* - number of inner-element vertices
* - number of additional nodes
* Whether midpoints and additional nodes exist depends on whether the
* DCORMG-array exists.
************************************************************************

      INTEGER FUNCTION TNDCNT (TRIA)
      IMPLICIT NONE
      INCLUDE 'stria.inc'
      INTEGER TRIA(SZTRIA)
      
      IF (TRIA(OLCORMG).NE.0) THEN
C Vertices in corners, on edges, on elements and additional nodes
        TNDCNT = TRIA(ONVT)+TRIA(ONVEDT)+TRIA(ONIEVT)+TRIA(ONANT)
      ELSE
C only vertices, no DCORMG array existing
        TNDCNT = TRIA(ONVT)
      END IF
      
      END
      
***********************************************************************
* Generate extended triangulation information
*
* This routine adds extended triangulation information to an existing
* FEAT triangulation. This allowes to generate a STRIA-structure from
* the /TRIAx/-common blocks with C2TRIA and add missing information
* to it.
*
* In:
*   TRIA   - array [1..SZTRIA] of integer
*            Source triangulation structure, which should be extended.
*
* Out:
*   TRIA   - array [1..SZTRIA] of integer
*            The extended triangulation structure.
***********************************************************************

      SUBROUTINE GENETR (TRIA)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cout.inc'
      
      INTEGER TRIA(SZTRIA)
      
C Cancel if the structure is already converted

      IF (IAND(TRIA(OTRIFL),1).NE.0) RETURN

C Call the auxiliary routines to generate the extended arrays.

      CALL GENET1 (TRIA)
      CALL GENET2 (TRIA)
      CALL GENET3 (TRIA)
      CALL GENET4 (TRIA)
      CALL GENET5 (TRIA)
      CALL GENET6 (TRIA)
      CALL GENET7 (TRIA)

      IF (IER.NE.0) THEN
        WRITE (MTERM,'(A)') 'GENETR: Not enough memory'
        RETURN
      END IF

C Mark that the structure is converted.
C We don't have incorporated midpoints, so the TRIFL-flag is
C simply set to 1 = basic converted structure.

      TRIA(OTRIFL) = 1
      
C We build a 2D-structure - denote that.

      TRIA(ONDIM)  = 2
      
      END 

***********************************************************************
* Generate extended triangulation information
*
* Auxiliary routine 1: Generate KXNPR
*
* This routine accepts a TRIA-structure and creates the data in KXNPR
* to it. If the handle to KXNPR is not yet allocated, a new handle
* will be created. Otherwise the old KXNPR-array is simply
* overwritten.
* All boundary nodes are initialized to be Dirichlet nodes by default.
* All inner nodes (including fictitious boundary inner nodes) are
* initialized to be of Neumann type by default.
*
* In:
*   TRIA   - array [1..SZTRIA] of integer
*            Source triangulation structure, which should be extended.
*
* Out:
*   TRIA   - array [1..SZTRIA] of integer
*            The extended triangulation structure.
***********************************************************************

      SUBROUTINE GENET1 (TRIA)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cout.inc'
      
      INTEGER TRIA(SZTRIA),NVT,NMT,NEL
      
      INTEGER I,KNPR,KXNPR
      
      NVT = TRIA(ONVT)
      NMT = TRIA(ONMT)
      NEL = TRIA(ONEL)

C Create the array if necessary

      IF (TRIA(OLXNPR).EQ.0) THEN
        CALL ZNEW (2*(NVT+NMT),3,TRIA(OLXNPR),'KXNPR ')
        IF (IER.NE.0) THEN
          WRITE (MTERM,'(A)') 'GENET1: Not enough memory'
          RETURN
        END IF
      END IF
      
      KNPR = L(TRIA(OLNPR))
      KXNPR = L(TRIA(OLXNPR))
      
C Transform corner vertices, generate KXNPR
      
      DO I=0,NVT-1
        KWORK(KXNPR+2*I) = 0
        KWORK(KXNPR+2*I+1) = 0
        IF (KWORK(KNPR+I).GT.0) THEN
C Boundary node - to be handled as Dirichlet in standard setting
          KWORK(KXNPR+2*I+1) = KWORK(KNPR+I)
          KWORK(KXNPR+2*I) = IOR(KWORK(KXNPR+2*I),2**8+2**12)
        ELSE IF (KWORK(KNPR+I).LT.-TRIA(ONEL)) THEN
C Fictitious boundary node
          KWORK(KXNPR+2*I+1) = -KWORK(KNPR+I)-TRIA(ONEL)
          KWORK(KXNPR+2*I) = IOR(KWORK(KXNPR+2*I),2**13)
        ELSE IF (KWORK(KNPR+I).LT.0) THEN
C Irregular inner node
          KWORK(KXNPR+2*I+1) = -KWORK(KNPR+I)
          KWORK(KXNPR+2*I) = IOR(KWORK(KXNPR+2*I),2**9)
        END IF
      END DO
      
C Transform edges

      DO I=NVT,NVT+NMT-1
        IF (KWORK(KNPR+I).GT.0) THEN
C Boundary edge
          KWORK(KXNPR+2*I+1) = KWORK(KNPR+KWORK(KNPR+I)-1)
          KWORK(KXNPR+2*I) = IOR(KWORK(KXNPR+2*I),2**8)
        ELSE IF (KWORK(KNPR+I).LT.0) THEN
C Irregular inner edge
          KWORK(KXNPR+2*I+1) = -KWORK(KNPR+I)
          KWORK(KXNPR+2*I) = IOR(KWORK(KXNPR+2*I),2**9)
        END IF
      END DO
      
      END 
      
***********************************************************************
* Generate extended triangulation information
*
* Auxiliary routine 2: Generate KEAN
*
* This routine accepts a TRIA-structure and creates the data in KEAN
* to it. If the handle to KEAN is not yet allocated, a new handle
* will be created. Otherwise the old KEAN-array is simply
* overwritten.
*
* In:
*   TRIA   - array [1..SZTRIA] of integer
*            Source triangulation structure, which should be extended.
*
* Out:
*   TRIA   - array [1..SZTRIA] of integer
*            The extended triangulation structure.
***********************************************************************

      SUBROUTINE GENET2 (TRIA)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cout.inc'
      
      INTEGER TRIA(SZTRIA),NVT,NMT,NEL
      
      INTEGER I,J,KEAN,KVERT,KMID,INDE,IEDG
      
      NVT = TRIA(ONVT)
      NMT = TRIA(ONMT)
      NEL = TRIA(ONEL)

C Cancel if we don't have midpoints

      IF (NMT.LE.0) RETURN

C Reserve memory for KEAN if necessary

      IF (TRIA(OLEAN).EQ.0) THEN
        CALL ZNEW (2*NMT,3,TRIA(OLEAN),'KEAN  ')
        IF (IER.NE.0) THEN
          WRITE (MTERM,'(A)') 'GENET2: Not enough memory'
          RETURN
        END IF
      END IF

C Build KEAN-array

      KVERT = L(TRIA(OLVERT))
      KMID = L(TRIA(OLMID))
      KEAN = L(TRIA(OLEAN))
      
      DO I=0,NEL-1
        DO J=0,TRIA(ONVE)-1
C Current node?
          INDE = KWORK(KVERT+I*NNVE+J)
C Current edge
          IEDG = KWORK(KMID+I*NNVE+J)-TRIA(ONVT) - 1
C save the connection
          KWORK(KEAN+2*IEDG) = INDE
C Number of the following node
          INDE = KWORK(KVERT+I*NNVE+MOD(J+1,TRIA(ONVE)))
C Save that, too
          KWORK(KEAN+2*IEDG+1) = INDE
        END DO
      END DO
      
      END 

***********************************************************************
* Generate extended triangulation information
*
* Auxiliary routine 3: Generate KVBDI
*
* This routine accepts a TRIA-structure and creates the data in KVBDI
* to it. If the handle to KVBDI is not yet allocated, a new handle
* will be created. Otherwise the old KVBDI-array is simply
* overwritten.
*
* In:
*   TRIA   - array [1..SZTRIA] of integer
*            Source triangulation structure, which should be extended.
*
* Out:
*   TRIA   - array [1..SZTRIA] of integer
*            The extended triangulation structure.
***********************************************************************

      SUBROUTINE GENET3 (TRIA)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cout.inc'
      
      INTEGER TRIA(SZTRIA),NVT,NMT,NEL,NVBD
      
      INTEGER I,KVBDI,KVBD
      
      NVT = TRIA(ONVT)
      NMT = TRIA(ONMT)
      NEL = TRIA(ONEL)

C Reserve memory for KVBDI if necessary

      IF (TRIA(OLVBDI).EQ.0) THEN
        CALL ZNEW (2*TRIA(ONVBD),3,TRIA(OLVBDI),'KVBDI ')
        IF (IER.NE.0) THEN
          WRITE (MTERM,'(A)') 'GENET3: Not enough memory'
          RETURN
        END IF
      END IF

C Build KVBDI-array.
C
C At first copy the KVBD-structure and assign the index in the second
C tag at each entry:

      IF (TRIA(OLVBD).NE.0) THEN
      
        KVBD = L(TRIA(OLVBD))
        KVBDI = L(TRIA(OLVBDI))
        NVBD = TRIA(ONVBD)
      
        DO I=0,NVBD-1
          KWORK(KVBDI+2*I) = KWORK(KVBD+I)
          KWORK(KVBDI+2*I+1) = I+1
        END DO
        
C Then sort this array (using heap-sort) using the vertex number
C as key. Each sub-array has 2 entries: the vertex number and the
C index number.

        CALL HSRTIN (TRIA(ONVBD),2,KWORK(KVBDI),1)
        
      END IF
      
      END 

***********************************************************************
* Generate extended triangulation information
*
* Auxiliary routine 4: Generate KMBD
*
* This routine accepts a TRIA-structure and creates the data in KMBD
* to it. If the handle to KVBD is not yet allocated, a new handle
* will be created. Otherwise the old KMBD-array is simply
* overwritten.
*
* In:
*   TRIA   - array [1..SZTRIA] of integer
*            Source triangulation structure, which should be extended.
*
* Out:
*   TRIA   - array [1..SZTRIA] of integer
*            The extended triangulation structure.
***********************************************************************

      SUBROUTINE GENET4 (TRIA)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cout.inc'
      
      INTEGER TRIA(SZTRIA),NVT,NMT,NEL,NVBD
      
      INTEGER I,J,KVERT,KMID
      INTEGER KVBD,KMBD,KBCT,KEBD
      INTEGER IBCT,IEL,NVEL,IVT
      
C Cancel if the structure is already converted

      NVT = TRIA(ONVT)
      NMT = TRIA(ONMT)
      NEL = TRIA(ONEL)

C Reserve memory for KMBD if necessary

      IF (TRIA(OLMBD).EQ.0) THEN
        CALL ZNEW (TRIA(ONVBD),3,TRIA(OLMBD),'KMBD  ')
        IF (IER.NE.0) THEN
          WRITE (MTERM,'(A)') 'GENET4: Not enough memory'
          RETURN
        END IF
      END IF

C Build KMBD-array

      IF (TRIA(OLMBD).NE.0) THEN
        
C We have to make a loop over all boundary nodes. The edge
C that follows a node counterclockwise is a boundary edge.

        KVERT = L(TRIA(OLVERT))
        KVBD = L(TRIA(OLVBD))
        KBCT = L(TRIA(OLBCT))
        KEBD = L(TRIA(OLEBD))
        KMBD = L(TRIA(OLMBD))
        KMID = L(TRIA(OLMID))
        NVBD = TRIA(ONVBD)
        NVEL = TRIA(ONVEL)
        
C At first a loop about all boundary components:
        
        DO IBCT = 0,TRIA(ONBCT)-1
        
C         On every component a loop over all vertices - or more
C         precisely, we make a loop about all elements along
C         this boundary component

          DO I=KWORK(KBCT+IBCT),KWORK(KBCT+IBCT+1)-1
          
C           I is the index of the boundary node. The real node number is
      
            IVT = KWORK(KVBD+I-1)
          
C           The adjacent element (- 1) on the boundary is

            IEL = KWORK(KEBD+I-1)-1

C           Now loop through all vertices on that element until we find
C           the node, which is the starting point of our edge:

            DO J=0,TRIA(ONVE)-1
              IF (KWORK(KVERT+IEL*NNVE+J).EQ.IVT) GOTO 10
            END DO    
              
            IF (J.EQ.TRIA(ONVE)) THEN
              WRITE (MTERM,'(A)') 'GENET4: KVERT/KEBD destroyed'
              RETURN
            END IF
          
10          CONTINUE

C           Ok, we found the node IVT. The edge that follows this vertex
C           counterclockwise is then our boundary edge:

            KWORK(KMBD+I-1) = KWORK(KMID+IEL*NNVE+J)
          
          END DO
        
        END DO
        
      END IF

      END 

***********************************************************************
* Generate extended triangulation information
*
* Auxiliary routine 5: Generate KMBDI
*
* This routine accepts a TRIA-structure and creates the data in KMBDI
* to it. If the handle to KMBDI is not yet allocated, a new handle
* will be created. Otherwise the old KMBDI-array is simply
* overwritten.
*
* The KMBD-array must have been created prior to calling this routine.
*
* In:
*   TRIA   - array [1..SZTRIA] of integer
*            Source triangulation structure, which should be extended.
*
* Out:
*   TRIA   - array [1..SZTRIA] of integer
*            The extended triangulation structure.
***********************************************************************

      SUBROUTINE GENET5 (TRIA)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cout.inc'
      
      INTEGER TRIA(SZTRIA),NVT,NMT,NEL,NVBD
      
      INTEGER I,KMBDI,KMBD
      
C Cancel if KMBD does not exist

      IF (TRIA(OLMBD).EQ.0) RETURN

      NVT = TRIA(ONVT)
      NMT = TRIA(ONMT)
      NEL = TRIA(ONEL)

C Reserve memory for KMBDI if necessary

      IF (TRIA(OLMBDI).EQ.0) THEN
        CALL ZNEW (2*TRIA(ONVBD),3,TRIA(OLMBDI),'KMBDI ')
        IF (IER.NE.0) THEN
          WRITE (MTERM,'(A)') 'GENET5: Not enough memory'
          RETURN
        END IF
      END IF

C Build KMBDI-array like KVBDI
      
C At first copy the KVBD-structure and assign the index in the second
C tag at each entry:

      IF ((TRIA(OLMBD).NE.0).AND.(TRIA(OLMBDI).NE.0)) THEN
      
        KMBD = L(TRIA(OLMBD))
        KMBDI = L(TRIA(OLMBDI))
        NVBD = TRIA(ONVBD)
      
        DO I=0,NVBD-1
          KWORK(KMBDI+2*I) = KWORK(KMBD+I)
          KWORK(KMBDI+2*I+1) = I+1
        END DO
        
C Then sort this array (using heap-sort) using the vertex number
C as key. Each sub-array has 2 entries: the vertex number and the
C index number.

        CALL HSRTIN (TRIA(ONVBD),2,KWORK(KMBDI),1)
        
      END IF

C Mark that the structure is converted.
C We don't have incorporated midpoints, so the TRIFL-flag is
C simply set to 1 = basic converted structure.

      TRIA(OTRIFL) = 1
      
      END 

***********************************************************************
* Generate extended triangulation information
*
* Auxiliary routine 6: Generate DAREA
*
* This routine accepts a TRIA-structure and creates the data in DAREA
* to it. If the handle to DAREA is not yet allocated, a new handle
* will be created. Otherwise the old DAERA-array is simply
* overwritten.
*
* In:
*   TRIA   - array [1..SZTRIA] of integer
*            Source triangulation structure, which should be extended.
*
* Out:
*   TRIA   - array [1..SZTRIA] of integer
*            The extended triangulation structure.
***********************************************************************

      SUBROUTINE GENET6 (TRIA)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cout.inc'
      
C parameters

      INTEGER TRIA(SZTRIA)

C local variables

      DOUBLE PRECISION AAA,SUM
      INTEGER IEL,I1,I2,I3,I4
      DOUBLE PRECISION X1,X2,X3,X4,Y1,Y2,Y3,Y4
      INTEGER KAREA,KVERT,KCORVG
      
C Reserve memory for DAREA if necessary

      IF (TRIA(OLAREA).EQ.0) THEN
        CALL ZNEW (TRIA(ONEL)+1,1,TRIA(OLAREA),'DAREA ')
        IF (IER.NE.0) THEN
          WRITE (MTERM,'(A)') 'GENET6: Not enough memory'
          RETURN
        END IF
      END IF

C Transfer some variables for easier access

      KAREA = L(TRIA(OLAREA))
      KVERT = L(TRIA(OLVERT))
      KCORVG = L(TRIA(OLCORVG))

      SUM = 0D0

      DO  IEL=0,TRIA(ONEL)-1

        I1 = KWORK(KVERT+IEL*NNVE)-1
        I2 = KWORK(KVERT+IEL*NNVE+1)-1
        I3 = KWORK(KVERT+IEL*NNVE+2)-1
        I4 = KWORK(KVERT+IEL*NNVE+3)-1

        X1 = DWORK(KCORVG+2*I1)
        X2 = DWORK(KCORVG+2*I2)
        X3 = DWORK(KCORVG+2*I3)
        X4 = DWORK(KCORVG+2*I4)

        Y1 = DWORK(KCORVG+2*I1+1)
        Y2 = DWORK(KCORVG+2*I2+1)
        Y3 = DWORK(KCORVG+2*I3+1)
        Y4 = DWORK(KCORVG+2*I4+1)

C       Calculate the area with the formula for general polygons;
C       compare e.g. "http://mathworld.wolfram.com/PolygonArea.html"

        AAA = 0.5D0*(  DABS((X1-X2)*(Y3-Y2)-(Y1-Y2)*(X3-X2))
     *                +DABS((X1-X4)*(Y3-Y4)-(Y1-Y4)*(X3-X4)) )
        
        DWORK(KAREA+IEL)=AAA
        
        SUM = SUM+AAA
      END DO
 
      DWORK(KAREA+TRIA(ONEL)) = SUM
 
      END

***********************************************************************
* Generate extended triangulation information
*
* Auxiliary routine 7: Translate DCORMG
*
* This routine initializes the variables NVPED,NVEDT,NIELV,NIEVT
* and NANT with compatibility data for FEAT structures.
***********************************************************************

      SUBROUTINE GENET7 (TRIA)

      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cout.inc'
      
C parameters

      INTEGER TRIA(SZTRIA)

C Check whether we have a DCORMG-array

      IF (TRIA(OLCORMG).EQ.0) THEN
      
C       No, we don't have. Then neighter element nor edge mitpoints are
C       realized there:

        TRIA(ONVPED) = 0
        TRIA(ONVEDT) = 0
        TRIA(ONIELV) = 0
        TRIA(ONIEVT) = 0
        TRIA(ONANT)  = 0

      ELSE

C       Yes, we have. Then DCORMG realizes edge midpoints only. No
C       element midpoints and no additional vertices are realized there
C       at the moment!
      
        TRIA(ONVPED) = 1
        TRIA(ONVEDT) = TRIA(ONMT)
        TRIA(ONIELV) = 0
        TRIA(ONIEVT) = 0
        TRIA(ONANT)  = 0
      
      END IF

      END 

***********************************************************************
* ZDISPQ - Quick dispose; auxiliary routine
* Works like ZDISP, but does not throw an error message if
* the handle is 0.
***********************************************************************

      SUBROUTINE ZDISPQ(WORK,LHAND,CNAME)
      IMPLICIT NONE
      INTEGER WORK,LHAND
      CHARACTER*(*) CNAME
      
      IF (LHAND.NE.0) CALL ZDISP(WORK,LHAND,CNAME)
      
      END

***********************************************************************
* Release extended triangulation information
*
* This routine disposes the memory that was used for an extended
* triangulation. After applying this routine to a STRIA-structure,
* the structure is FEAT-conform again. This is done by simply
* deleting unnecessary information from it.
* This also deletes "additional nodes" from the DCORMG-array
* if they exist.
*
* In:
*   TRIA  - array [1..SZTRIA] of integer
*           Triangulation structure
*   ICBCK - whether to convert the structure back.
*           =0: don't convert, simply dispose extended structures
*               from heap.
*               DCORMG is not touched.
*           =1: convert back and release structure from heap;
*               recalculate midpoints for DCORMG and release unused
*               memory for that array.
*           =2: convert back but don't release from heap
*               recalculate midpoints for DCORMG and release unused
*               memory for that array.
***********************************************************************

      SUBROUTINE DISETR (TRIA,ICBCK)

      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cmem.inc'
      
      INTEGER TRIA(SZTRIA)
      INTEGER ICBCK
      
      INTEGER I,KNPR,KXNPR
      
C Cancel if this is not in converted style

      IF (IAND(TRIA(OTRIFL),1).EQ.0) RETURN
      
C Don't perform any reconstruction if not desired

      IF (ICBCK.NE.0) THEN
      
C Reconstruct KNPR-array from KXNPR

        print *,TRIA(OLNPR)
        KNPR = L(TRIA(OLNPR))
        KXNPR = L(TRIA(OLXNPR))

C Vertex information

        DO I=0,TRIA(ONVT)-1
C Inner node or real boundary node?
C          print *,'1: ',KWORK(KNPR+I)
          
          IF (IAND(KWORK(KXNPR+2*I),2**8).NE.0) THEN 
          
C           This is a real boundary node. Should it be treated
C           as Neumann?

            IF (IAND(KWORK(KXNPR+2*I),2**12).EQ.0) THEN

C             Yep, Neumann boundary node

              KWORK(KNPR+I) = 0         
              
            ELSE
            
C             Standard boundary node             
          
              KWORK(KNPR+I) = KWORK(KXNPR+2*I+1)
              
            END IF
            
          ELSE
          
C           This is an inner node - maybe a fictitious boundary
C           one, which should be treated as Dirichlet

            IF (IAND(KWORK(KXNPR+2*I),2**12+2**13).EQ.
     *          (2**12+2**13)) THEN
            
C             Ok, fictitious boundary Dirichlet node!

              KWORK(KNPR+I) = -(KWORK(KXNPR+2*I+1)+TRIA(ONVT))
              
C           Otherwise test, if we have an irregular inner node:              
              
            ELSE IF (IAND(KWORK(KXNPR+2*I),2**9).NE.0) THEN
            
              KWORK(KNPR+I) = -KWORK(KXNPR+2*I+1)
            
            ELSE
            
C             Standard inner node

              KWORK(KNPR+I) = 0
            
            END IF
            
          END IF
          
C          print *,'2: ',KWORK(KNPR+I)
        
        END DO

C Edge information

        DO I=TRIA(ONVT),TRIA(ONVT)+TRIA(ONMT)-1
        
C Inner node or real bondary node?

          IF (IAND(KWORK(KNPR+I),2**8).NE.0) THEN 
          
C inner vertex; fictitious boundary edges do not exist

            KWORK(KNPR+I) = 0
            
          ELSE
          
C boundary edge; save the preceding vertex number

            KWORK(KNPR+I) = KWORK(KXNPR+2*I-2*TRIA(ONVT)+1)
            
          END IF
          
        END DO

C Remove "Additional nodes" and "element midpoints", keep the
C first vertex of each edge as edge midpoint as converted earlier:

        IF (TRIA(OLCORMG).GT.0) THEN
        
          IF ( (TRIA(ONANT).GT.0) .OR.
     *         (TRIA(ONVEDT).GT.TRIA(ONMT)) .OR.
     *         (TRIA(ONIEVT).GT.TRIA(ONEL)) ) THEN
            CALL ZDISP(2*TRIA(ONMT),TRIA(OLCORMG),'DCORMG')
            TRIA(ONVEDT) = TRIA(ONMT)
            TRIA(ONVPED) = 1
            TRIA(ONIEVT) = 0
            TRIA(ONIELV) = 0
            TRIA(ONANT)  = 0
          END IF
          
C         Recalculate edge midpoints to DCORMG:

          CALL GXVTED (TRIA,1,0.5D0)

        END IF

      END IF
 
      IF (ICBCK.NE.2) THEN 
      
C       Remove extended boundary information

        CALL GENFBA (TRIA,0)
 
C       Remove unnecessary arrays - without any error message

        CALL ZDISPQ(0,TRIA(OLAREA),'DAREA ')
        CALL ZDISPQ(0,TRIA(OLMBDI),'KMBDI ')
        CALL ZDISPQ(0,TRIA(OLMBD),'KMBD  ')
        CALL ZDISPQ(0,TRIA(OLVBDI),'KVBDI ')
        CALL ZDISPQ(0,TRIA(OLXNPR),'KXNPR ')
        CALL ZDISPQ(0,TRIA(OLEAN),'KEAN  ')
        
      END IF

C Structure is now standard FEAT style

      TRIA(OTRIFL) = 0
      
      END

***********************************************************************
* Generate/release extended fictitious boundary information arrays
*
* This routine builds up the extended fictitious boudary information
* arrays KFBVT/KFBMT from the array KXNPR. TRIA must be an extended
* triangulation structure. If the arrays KFBVT/KFBMT exist, they are
* overwritten.
*
* In:
*   TRIA   - array [1..SZTRIA] of integer
*            Source triangulation structure, which should be extended.
*   IFLAG  - integer
*            Specifies which array to generate.
*            = 0: Release KFBVT/KFBMT if it's allocated
*            = 1: Generate NFBVT/KFBVT
*            = 2: Generate NFBMT/KFBMT
*            = 3: Generate both
*
* Out:
*   TRIA.NFBVT - Number of fictitious boundary vertices
*   TRIA.NFBMT - Number of fictitious boundary edges
*   TRIA.KFBVT - array wíth all fictitious boundary vertices
*   TRIA.KFBMT - array wíth all fictitious boundary edges
*
* When KFBVT/KFBMT do not exist, arrays for all vertices/edges are
* allocated on the heap. KFBVT is allocated as array [1..NVT] and
* KMBMT as array [1..NMT], as the fictitious boundary infomation
* might change later. although allocated as 1..NVT/NMT, only the
* first NFBVT/NFBMT entries of the arrays are used.
*
* If the arrays already exist, they are overwritten without being
* reallocated. If one of these arrays is not large enough,
* the corresponding NFBxT value is set to -1 to indicate that
* the routine broke down.
***********************************************************************

      SUBROUTINE GENFBA (TRIA,IFLAG)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cout.inc'
      
      INTEGER TRIA(SZTRIA),IFLAG
      
      INTEGER NFBVT,NFBMT,NVT,NMT,KXNPR
      INTEGER FBVTL, FBMTL, KFBVT,KFBMT,I,J
      
C     TRIA must be an extended structure

      IF (TRIA(OTRIFL).EQ.0) RETURN
      
C     Should we generate or release the arrays?

      IF (IFLAG.EQ.0) THEN
        IF (TRIA(OLFBVT).NE.0) CALL ZDISP(0,TRIA(OLFBVT),'KFBVT ')
        IF (TRIA(OLFBMT).NE.0) CALL ZDISP(0,TRIA(OLFBMT),'KFBMT ')
        RETURN
      END IF
      
C     How many nodes/midpoints do we have altogether?

      NVT = TRIA(ONVT)
      NMT = TRIA(ONMT)
      
C     Loop through all elements in KXNPR to find the fictitious
C     boundary nodes/midpoints

      KXNPR = L(TRIA(OLXNPR))
      NFBVT = 0
      NFBMT = 0
      
C     Handle vertices:
      
      IF (IAND(IFLAG,1).NE.0) THEN
      
C       Is the arrays allocated? If not, allocate!
C       Save the number of free entries in FBxTL

        IF (TRIA(OLFBVT).EQ.0) THEN
          CALL ZNEW (NVT,3,TRIA(OLFBVT),'KFBVT ')
          FBVTL = NVT
        ELSE
          CALL ZLEN (TRIA(OLFBVT),FBVTL)
        END IF
        KFBVT = L(TRIA(OLFBVT))  

C       Loop through the vertices

        DO I=0,NVT-1
        
C         Do we have a fictitious boundary element?
        
          J = KWORK(KXNPR+2*I)
          IF (IAND(J,2**13).NE.0) THEN
          
            NFBVT = NFBVT + 1
            
C           Cancel if the memory is not large enough
            
            IF (NFBVT.GT.FBVTL) THEN
              NFBVT = -1
              GOTO 10
            END IF
            
C           Add the vertex to the end of the list

            KWORK(KFBVT+NFBVT) = I+1
          
          END IF
        
        END DO ! I

10      CONTINUE

C       Save the number of vertices to the structure

        TRIA(ONFBVT) = NFBVT
        
      END IF
      
C     Handle edges:

      IF (IAND(IFLAG,2).NE.0) THEN

C       Is the arrays allocated? If not, allocate!
C       Save the number of free entries in FBxTL
      
        IF (TRIA(OLFBMT).EQ.0) THEN
          CALL ZNEW (NMT,3,TRIA(OLFBMT),'KFBMT ')
          FBMTL = NMT
        ELSE
          CALL ZLEN (TRIA(OLFBMT),FBMTL)
        END IF
        KFBMT = L(TRIA(OLFBMT))  

C       Loop through the edges

        DO I=NVT,NVT+NMT-1
        
C         Do we have a fictitious boundary element?
        
          J = KWORK(KXNPR+2*I)
          IF (IAND(J,2**13).NE.0) THEN
          
            NFBMT = NFBMT + 1
            
C           Cancel if the memory is not large enough
            
            IF (NFBMT.GT.FBMTL) THEN
              NFBMT = -1
              GOTO 20
            END IF
            
C           Add the edge to the end of the list
            
            KWORK(KFBMT+NFBMT) = I+1
          
          END IF
        
        END DO ! I
        
20      CONTINUE
      
C       Save the number of edges to the structure

        TRIA(ONFBMT) = NFBMT

      END IF
      
      END

***********************************************************************
* Backup/duplicate triangulation structure
*
* This routine makes a copy of a triangulation structure in memory.
* The variable IDPFLG decides on which arrays are copied in memory
* and which not.
*
* By setting the corresponding bit in IDPFLG to 0, the array is
* duplicated in memory, and any change to the new array will not harm
* the original one.
* By setting the corresponding bit in IDPFLG to 1, the array is not
* duplicated. The handle of the original structure is simply put into
* the new structure to make the information accessable. Then this array
* is shared between two triangulation structures!
*
* In:
*   TRIA   - array [1..SZTRIA] of integer
*            Original triangulation structure
*   IUPD   - integer
*            Defines how to create the backup.
*            = 0: Treat TRIADP as empty destination structure.
*                 Clear TRIADP and build it according to TRIA
*                 and IDPFLG.
*            = 1: Treat TRIADP as existing copy of TRIA which has to
*                 be updated by the help of TRIA. Recover all
*                 vectors of TRIADP by those of TRIA.
*                 (I.e. those arrays which were duplicated by a previous
*                 call to TRIDUP with IUPD=0.)
*                 IDPFLG can used to specify which data do copy
*                 from TRIA to TRIADP:
*                 =0 : Copy all data that was copied previously
*                 >0:  Copy all arrays where the corresponding bit
*                      in IDPFLG is =0 (see below) and which exist as
*                      duplicates. Arrays corresponding
*                      to bits =1 or where the handles in TRIA and
*                      TRIDP coincide are not touched.
*   IDPFLG - integer
*            Bitfield that decides which handles are a copy of another
*            structure, thus which arrays are shared between the new
*            and the old structure. Is only used if IUPD=0.
*             Bit  0: Don't copy DCORVG , simply copy its handle
*             Bit  1: Don't copy DCORMG , simply copy its handle
*             Bit  2: Don't copy KVERT  , simply copy its handle
*             Bit  3: Don't copy KMID   , simply copy its handle
*             Bit  4: Don't copy KADJ   , simply copy its handle
*             Bit  5: Don't copy KVEL   , simply copy its handle
*             Bit  6: Don't copy KMEL   , simply copy its handle
*             Bit  7: Don't copy KNPR   , simply copy its handle
*             Bit  8: Don't copy KMM    , simply copy its handle
*             Bit  9: Don't copy KVBD   , simply copy its handle
*             Bit 10: Don't copy KEBD   , simply copy its handle
*             Bit 11: Don't copy KBCT   , simply copy its handle
*             Bit 12: Don't copy DVBDP  , simply copy its handle
*             Bit 13: Don't copy DMBDP  , simply copy its handle
*             Bit 14: Don't copy KMBD   , simply copy its handle
*             Bit 15: Don't copy KXNPR  , simply copy its handle
*             Bit 16: Don't copy KEAN   , simply copy its handle
*             Bit 17: Don't copy KVBDI  , simply copy its handle
*             Bit 18: Don't copy KMBDI  , simply copy its handle
*             Bit 19: Don't copy DAREA  , simply copy its handle
*             Bit 20: Don't copy KFBVT  , simply copy its handle
*             Bit 21: Don't copy KFBMT  , simply copy its handle
*
* Out:
*   TRIADP - array [1..SZTRIA] of integer
*            Will receive a copy of TRIA. Depending in IDPFLG, some
*            of the arrays of TRIA are duplicated in memory.
***********************************************************************

      SUBROUTINE TRIDUP (TRIA,TRIADP,IUPD,IDPFLG)
      
      IMPLICIT NONE
      
      INCLUDE 'stria.inc'
      
      INTEGER TRIA(SZTRIA),TRIADP(SZTRIA),IDPFLG
      
      INTEGER I,IDP,IUPD
      
      IF (IUPD.EQ.0) THEN
      
C       TRIADP is an empty structure which has to be initialized.
C      
C       Clear the new structure

        CALL LCL3(TRIADP,SZTRIA)
        
C       Decide on IDPFLG which arrays to copy

        TRIADP(OIDPFLG) = IDPFLG
        IDP = IDPFLG
        
      ELSE
      
C       Create a bitfield what to copy by ORing IDPFLG with what
C       we have in TRIADP. That way, only arrays that exist as
C       real duplicates are copied from TRIA to TRIADP.
      
        IDP = IOR(IDPFLG,TRIADP(OIDPFLG))
        
      END IF
      
C     As TRIADP(.)=0, a new handle is allocated and the old array is
C     copied into that:
      
      IF ((IAND(IDP,2** 0).EQ.0).AND.(TRIA(OLCORVG).NE.0))
     *  CALL ZCPY(TRIA(OLCORVG),'TRISRC',TRIADP(OLCORVG),'TRIDST')
      
      IF ((IAND(IDP,2** 1).EQ.0).AND.(TRIA(OLCORMG).NE.0))
     *  CALL ZCPY(TRIA(OLCORMG),'TRISRC',TRIADP(OLCORMG),'TRIDST')
      
      IF ((IAND(IDP,2** 2).EQ.0).AND.(TRIA(OLVERT).NE.0))
     *  CALL ZCPY(TRIA(OLVERT ),'TRISRC',TRIADP(OLVERT ),'TRIDST')
      
      IF ((IAND(IDP,2** 3).EQ.0).AND.(TRIA(OLMID).NE.0))
     *  CALL ZCPY(TRIA(OLMID  ),'TRISRC',TRIADP(OLMID  ),'TRIDST')
      
      IF ((IAND(IDP,2** 4).EQ.0).AND.(TRIA(OLADJ).NE.0))
     *  CALL ZCPY(TRIA(OLADJ  ),'TRISRC',TRIADP(OLADJ  ),'TRIDST')
      
      IF ((IAND(IDP,2** 5).EQ.0).AND.(TRIA(OLVEL).NE.0))
     *  CALL ZCPY(TRIA(OLVEL  ),'TRISRC',TRIADP(OLVEL  ),'TRIDST')
      
      IF ((IAND(IDP,2** 6).EQ.0).AND.(TRIA(OLMEL).NE.0))
     *  CALL ZCPY(TRIA(OLMEL  ),'TRISRC',TRIADP(OLMEL  ),'TRIDST')
      
      IF ((IAND(IDP,2** 7).EQ.0).AND.(TRIA(OLNPR).NE.0))
     *  CALL ZCPY(TRIA(OLNPR  ),'TRISRC',TRIADP(OLNPR  ),'TRIDST')
      
      IF ((IAND(IDP,2** 8).EQ.0).AND.(TRIA(OLMM).NE.0))
     *  CALL ZCPY(TRIA(OLMM   ),'TRISRC',TRIADP(OLMM   ),'TRIDST')
      
      IF ((IAND(IDP,2** 9).EQ.0).AND.(TRIA(OLVBD).NE.0))
     *  CALL ZCPY(TRIA(OLVBD  ),'TRISRC',TRIADP(OLVBD  ),'TRIDST')
      
      IF ((IAND(IDP,2**10).EQ.0).AND.(TRIA(OLEBD).NE.0))
     *  CALL ZCPY(TRIA(OLEBD  ),'TRISRC',TRIADP(OLEBD  ),'TRIDST')
      
      IF ((IAND(IDP,2**11).EQ.0).AND.(TRIA(OLBCT).NE.0))
     *  CALL ZCPY(TRIA(OLBCT  ),'TRISRC',TRIADP(OLBCT  ),'TRIDST')
      
      IF ((IAND(IDP,2**12).EQ.0).AND.(TRIA(OLVBDP).NE.0))
     *  CALL ZCPY(TRIA(OLVBDP ),'TRISRC',TRIADP(OLVBDP ),'TRIDST')
      
      IF ((IAND(IDP,2**13).EQ.0).AND.(TRIA(OLMBDP).NE.0))
     *  CALL ZCPY(TRIA(OLMBDP ),'TRISRC',TRIADP(OLMBDP ),'TRIDST')
      
      IF ((IAND(IDP,2**14).EQ.0).AND.(TRIA(OLMBD).NE.0))
     *  CALL ZCPY(TRIA(OLMBD  ),'TRISRC',TRIADP(OLMBD  ),'TRIDST')
      
      IF ((IAND(IDP,2**15).EQ.0).AND.(TRIA(OLXNPR).NE.0))
     *  CALL ZCPY(TRIA(OLXNPR ),'TRISRC',TRIADP(OLXNPR ),'TRIDST')
      
      IF ((IAND(IDP,2**16).EQ.0).AND.(TRIA(OLEAN).NE.0))
     *  CALL ZCPY(TRIA(OLEAN  ),'TRISRC',TRIADP(OLEAN  ),'TRIDST')
      
      IF ((IAND(IDP,2**17).EQ.0).AND.(TRIA(OLVBDI).NE.0))
     *  CALL ZCPY(TRIA(OLVBDI ),'TRISRC',TRIADP(OLVBDI ),'TRIDST')
      
      IF ((IAND(IDP,2**18).EQ.0).AND.(TRIA(OLMBDI).NE.0))
     *  CALL ZCPY(TRIA(OLMBDI ),'TRISRC',TRIADP(OLMBDI ),'TRIDST')
      
      IF ((IAND(IDP,2**19).EQ.0).AND.(TRIA(OLAREA).NE.0))
     *  CALL ZCPY(TRIA(OLAREA ),'TRISRC',TRIADP(OLAREA ),'TRIDST')
      
      IF ((IAND(IDP,2**20).EQ.0).AND.(TRIA(OLFBVT).NE.0))
     *  CALL ZCPY(TRIA(OLFBVT ),'TRISRC',TRIADP(OLFBVT ),'TRIDST')
      
      IF ((IAND(IDP,2**21).EQ.0).AND.(TRIA(OLFBMT).NE.0))
     *  CALL ZCPY(TRIA(OLFBMT ),'TRISRC',TRIADP(OLFBMT ),'TRIDST')
      
C     Copy the missing information from the old structure to the new one
      
      DO I=1,SZTRIA
        IF (TRIADP(I).EQ.0) TRIADP(I)=TRIA(I)
      END DO
      
C     That's it
      
      END
      
***********************************************************************
* Correct handles of shared information
*
* With this routine, destroyed handles of shared information between
* two triangulation structures can be recovered. TRIA is assumed to
* be a valid triangulation structure. The handles of all arrays
* of information that TRIA shares with TRIADS are copied to TRIADS.
*
* This routine is typically used in refinement processes when e.g.
* coordinates of grid points are shared between two levels of
* refinement. Upon refining the coarser grid, the handle in the
* source grid gets invalid, while the handle of the fine grid is
* the new valid one. TRICSH can now be used to copy the handle
* back from the fine grid triangulation structure to the coarse
* grid triangulation structure.
*
*  Example: DCORVG is shared between all refinement levels.
*  Then the caller can create the grid by:
*
*    CALL GENTRI (IMETH=2, filename='anything.tri')  -> read coarse gr.
*    DO I=NLMIN+1,NLMAX                              -> refine to NLMAX
*      CALL GENTRI (TRIA(I-1),TRIA(I),IMETH=1,IDPFLG=1)
*      --> generates TRIA(I), whereas destroys LCORVG on level I-1!
*    END DO
*    DO I=NLMAX-1,NLMIN,-1             -> write LCORVG of lv. n+1 into
*      CALL TRICSH (TRIA(I),TRIA(I+1))    LCORVG of level n
*    END DO
*
*  Afterwards, the information is again shared correctly between all
*  levels.
*
* In:
*   TRIA   - array [1..SZTRIA] of integer
*            Valid triangulation structure, which is a modified
*            copy of TRIADS.
*
* Out:
*   TRIADS - array [1..SZTRIA] of integer
*            The handles of all arrays that are shared with TRIA
*            are copied to TRIADS. Old handles in TRIADS are assumed
*            to be invalid and will be overwritten by those of TRIA.
***********************************************************************

      SUBROUTINE TRICSH (TRIA,TRIADS)
      
      IMPLICIT NONE
      
      INCLUDE 'stria.inc'
      
      INTEGER TRIA(SZTRIA),TRIADS(SZTRIA)
      
      INTEGER IDPFLG
      
C     Get the bitfield that tells us about which arrays are shared

      IDPFLG = TRIA(OIDPFLG)
      
C     Restore the handles
      
      IF (IAND(IDPFLG,2** 0).NE.0) TRIADS(OLCORVG) = TRIA(OLCORVG)
      IF (IAND(IDPFLG,2** 1).NE.0) TRIADS(OLCORMG) = TRIA(OLCORMG)
      IF (IAND(IDPFLG,2** 2).NE.0) TRIADS(OLVERT ) = TRIA(OLVERT )
      IF (IAND(IDPFLG,2** 3).NE.0) TRIADS(OLMID  ) = TRIA(OLMID  )
      IF (IAND(IDPFLG,2** 4).NE.0) TRIADS(OLADJ  ) = TRIA(OLADJ  )
      IF (IAND(IDPFLG,2** 5).NE.0) TRIADS(OLVEL  ) = TRIA(OLVEL  )
      IF (IAND(IDPFLG,2** 6).NE.0) TRIADS(OLMEL  ) = TRIA(OLMEL  )
      IF (IAND(IDPFLG,2** 7).NE.0) TRIADS(OLNPR  ) = TRIA(OLNPR  )
      IF (IAND(IDPFLG,2** 8).NE.0) TRIADS(OLMM   ) = TRIA(OLMM   )
      IF (IAND(IDPFLG,2** 9).NE.0) TRIADS(OLVBD  ) = TRIA(OLVBD  )
      IF (IAND(IDPFLG,2**10).NE.0) TRIADS(OLEBD  ) = TRIA(OLEBD  )
      IF (IAND(IDPFLG,2**11).NE.0) TRIADS(OLBCT  ) = TRIA(OLBCT  )
      IF (IAND(IDPFLG,2**12).NE.0) TRIADS(OLVBDP ) = TRIA(OLVBDP )
      IF (IAND(IDPFLG,2**13).NE.0) TRIADS(OLMBDP ) = TRIA(OLMBDP )
      IF (IAND(IDPFLG,2**14).NE.0) TRIADS(OLMBD  ) = TRIA(OLMBD  )
      IF (IAND(IDPFLG,2**15).NE.0) TRIADS(OLXNPR ) = TRIA(OLXNPR )
      IF (IAND(IDPFLG,2**16).NE.0) TRIADS(OLEAN  ) = TRIA(OLEAN  )
      IF (IAND(IDPFLG,2**17).NE.0) TRIADS(OLVBDI ) = TRIA(OLVBDI )
      IF (IAND(IDPFLG,2**18).NE.0) TRIADS(OLMBDI ) = TRIA(OLMBDI )
      IF (IAND(IDPFLG,2**19).NE.0) TRIADS(OLAREA ) = TRIA(OLAREA )
      IF (IAND(IDPFLG,2**20).NE.0) TRIADS(OLFBVT ) = TRIA(OLFBVT )
      IF (IAND(IDPFLG,2**21).NE.0) TRIADS(OLFBMT ) = TRIA(OLFBMT )
      
      END

***********************************************************************
* Restore triangulation structure
*
* This routine restores data of a triangulation structure. All
* information arrays not shared between TRIA and another triangulation
* structure are copied (back) into the TRIADS. 
*
* In:
*   TRIA   - array [1..SZTRIA] of integer
*            Backup of a triangulation structure
*
* Out:
*   TRIADS - array [1..SZTRIA] of integer
*            All arrays where a duplicates exist in TRIA are copied
*            to TRIA, overwriting the old information arrays.
***********************************************************************

      SUBROUTINE TRIRST (TRIA,TRIADS)
      
      IMPLICIT NONE
      
      INCLUDE 'stria.inc'
      
      INTEGER TRIA(SZTRIA),TRIADS(SZTRIA)
      
      INTEGER IDPFLG
      
C     Get the bitfield that tells us about which arrays we have a
C     duplicate from

      IDPFLG = TRIA(OIDPFLG)
      
C     Restore the arrays
      
      IF ( (IAND(IDPFLG,2** 0).EQ.0) .AND. (TRIADS(OLCORVG).NE.0) ) 
     *  CALL ZCPY(TRIADS(OLCORVG),'TRISRC',TRIA(OLCORVG),'TRIDST')
     
      IF ( (IAND(IDPFLG,2** 1).EQ.0) .AND. (TRIADS(OLCORMG).NE.0) )
     *  CALL ZCPY(TRIADS(OLCORMG),'TRISRC',TRIA(OLCORMG),'TRIDST')
     
      IF ( (IAND(IDPFLG,2** 2).EQ.0) .AND. (TRIADS(OLVERT ).NE.0) )
     *  CALL ZCPY(TRIADS(OLVERT ),'TRISRC',TRIA(OLVERT ),'TRIDST')
     
      IF ( (IAND(IDPFLG,2** 3).EQ.0) .AND. (TRIADS(OLMID  ).NE.0) )
     *  CALL ZCPY(TRIADS(OLMID  ),'TRISRC',TRIA(OLMID  ),'TRIDST')
     
      IF ( (IAND(IDPFLG,2** 4).EQ.0) .AND. (TRIADS(OLADJ  ).NE.0) )
     *  CALL ZCPY(TRIADS(OLADJ  ),'TRISRC',TRIA(OLADJ  ),'TRIDST')
     
      IF ( (IAND(IDPFLG,2** 5).EQ.0) .AND. (TRIADS(OLVEL  ).NE.0) )
     *  CALL ZCPY(TRIADS(OLVEL  ),'TRISRC',TRIA(OLVEL  ),'TRIDST')
     
      IF ( (IAND(IDPFLG,2** 6).EQ.0) .AND. (TRIADS(OLMEL  ).NE.0) )
     *  CALL ZCPY(TRIADS(OLMEL  ),'TRISRC',TRIA(OLMEL  ),'TRIDST')
     
      IF ( (IAND(IDPFLG,2** 7).EQ.0) .AND. (TRIADS(OLNPR  ).NE.0) )
     *  CALL ZCPY(TRIADS(OLNPR  ),'TRISRC',TRIA(OLNPR  ),'TRIDST')
     
      IF ( (IAND(IDPFLG,2** 8).EQ.0) .AND. (TRIADS(OLMM   ).NE.0) )
     *  CALL ZCPY(TRIADS(OLMM   ),'TRISRC',TRIA(OLMM   ),'TRIDST')
     
      IF ( (IAND(IDPFLG,2** 9).EQ.0) .AND. (TRIADS(OLVBD  ).NE.0) )
     *  CALL ZCPY(TRIADS(OLVBD  ),'TRISRC',TRIA(OLVBD  ),'TRIDST')
     
      IF ( (IAND(IDPFLG,2**10).EQ.0) .AND. (TRIADS(OLEBD  ).NE.0) )
     *  CALL ZCPY(TRIADS(OLEBD  ),'TRISRC',TRIA(OLEBD  ),'TRIDST')
     
      IF ( (IAND(IDPFLG,2**11).EQ.0) .AND. (TRIADS(OLBCT  ).NE.0) )
     *  CALL ZCPY(TRIADS(OLBCT  ),'TRISRC',TRIA(OLBCT  ),'TRIDST')
     
      IF ( (IAND(IDPFLG,2**12).EQ.0) .AND. (TRIADS(OLVBDP ).NE.0) )
     *  CALL ZCPY(TRIADS(OLVBDP ),'TRISRC',TRIA(OLVBDP ),'TRIDST')
     
      IF ( (IAND(IDPFLG,2**13).EQ.0) .AND. (TRIADS(OLMBDP ).NE.0) )
     *  CALL ZCPY(TRIADS(OLMBDP ),'TRISRC',TRIA(OLMBDP ),'TRIDST')
     
      IF ( (IAND(IDPFLG,2**14).EQ.0) .AND. (TRIADS(OLMBD  ).NE.0) )
     *  CALL ZCPY(TRIADS(OLMBD  ),'TRISRC',TRIA(OLMBD  ),'TRIDST')
     
      IF ( (IAND(IDPFLG,2**15).EQ.0) .AND. (TRIADS(OLXNPR ).NE.0) )
     *  CALL ZCPY(TRIADS(OLXNPR ),'TRISRC',TRIA(OLXNPR ),'TRIDST')
     
      IF ( (IAND(IDPFLG,2**16).EQ.0) .AND. (TRIADS(OLEAN  ).NE.0) )
     *  CALL ZCPY(TRIADS(OLEAN  ),'TRISRC',TRIA(OLEAN  ),'TRIDST')
     
      IF ( (IAND(IDPFLG,2**17).EQ.0) .AND. (TRIADS(OLVBDI ).NE.0) )
     *  CALL ZCPY(TRIADS(OLVBDI ),'TRISRC',TRIA(OLVBDI ),'TRIDST')
     
      IF ( (IAND(IDPFLG,2**18).EQ.0) .AND. (TRIADS(OLMBDI ).NE.0) )
     *  CALL ZCPY(TRIADS(OLMBDI ),'TRISRC',TRIA(OLMBDI ),'TRIDST')
     
      IF ( (IAND(IDPFLG,2**19).EQ.0) .AND. (TRIADS(OLAREA ).NE.0) )
     *  CALL ZCPY(TRIADS(OLAREA ),'TRISRC',TRIA(OLAREA ),'TRIDST')
     
      IF ( (IAND(IDPFLG,2**20).EQ.0) .AND. (TRIADS(OLFBVT ).NE.0) )
     *  CALL ZCPY(TRIADS(OLFBVT ),'TRISRC',TRIA(OLFBVT ),'TRIDST')
     
      IF ( (IAND(IDPFLG,2**21).EQ.0) .AND. (TRIADS(OLFBMT ).NE.0) )
     *  CALL ZCPY(TRIADS(OLFBMT ),'TRISRC',TRIA(OLFBMT ),'TRIDST')
      
      END

***********************************************************************
* Delete triangulation structure
*
* This routine deletes a triangulation structure. All arrays that
* are a real duplicates or that exist only once are deleted from memory. 
* All arrays that are maintained by another triangulation structure
* are not touched.
*
* Note that if some arrays arised from the duplication of a structure,
* deleting the "root" structure will destroy the arrays in all "child"
* structures as well without notification!!!
*
* In:
*   TRIA   - array [1..SZTRIA] of integer
*            Backup of a triangulation structure
*
* Out:
*   All arrays in TRIA wich are duplicates of another triangulation
*   structure are deleted from memory. TRIA is filled with 0.
***********************************************************************

      SUBROUTINE TRIDEL (TRIA)
      
      IMPLICIT NONE
      
      INCLUDE 'stria.inc'
      
      INTEGER TRIA(SZTRIA)
      
      INTEGER IDPFLG
      
C     Get the bitfield that tells us about which arrays we have a
C     duplicate from

      IDPFLG = TRIA(OIDPFLG)
      
C     Delete the arrays
      
      IF (IAND(IDPFLG,2** 0).EQ.0) CALL ZDISPQ(0,TRIA(OLCORVG),'TRIDUP')
      IF (IAND(IDPFLG,2** 1).EQ.0) CALL ZDISPQ(0,TRIA(OLCORMG),'TRIDUP')
      IF (IAND(IDPFLG,2** 2).EQ.0) CALL ZDISPQ(0,TRIA(OLVERT ),'TRIDUP')
      IF (IAND(IDPFLG,2** 3).EQ.0) CALL ZDISPQ(0,TRIA(OLMID  ),'TRIDUP')
      IF (IAND(IDPFLG,2** 4).EQ.0) CALL ZDISPQ(0,TRIA(OLADJ  ),'TRIDUP')
      IF (IAND(IDPFLG,2** 5).EQ.0) CALL ZDISPQ(0,TRIA(OLVEL  ),'TRIDUP')
      IF (IAND(IDPFLG,2** 6).EQ.0) CALL ZDISPQ(0,TRIA(OLMEL  ),'TRIDUP')
      IF (IAND(IDPFLG,2** 7).EQ.0) CALL ZDISPQ(0,TRIA(OLNPR  ),'TRIDUP')
      IF (IAND(IDPFLG,2** 8).EQ.0) CALL ZDISPQ(0,TRIA(OLMM   ),'TRIDUP')
      IF (IAND(IDPFLG,2** 9).EQ.0) CALL ZDISPQ(0,TRIA(OLVBD  ),'TRIDUP')
      IF (IAND(IDPFLG,2**10).EQ.0) CALL ZDISPQ(0,TRIA(OLEBD  ),'TRIDUP')
      IF (IAND(IDPFLG,2**11).EQ.0) CALL ZDISPQ(0,TRIA(OLBCT  ),'TRIDUP')
      IF (IAND(IDPFLG,2**12).EQ.0) CALL ZDISPQ(0,TRIA(OLVBDP ),'TRIDUP')
      IF (IAND(IDPFLG,2**13).EQ.0) CALL ZDISPQ(0,TRIA(OLMBDP ),'TRIDUP')
      IF (IAND(IDPFLG,2**14).EQ.0) CALL ZDISPQ(0,TRIA(OLMBD  ),'TRIDUP')
      IF (IAND(IDPFLG,2**15).EQ.0) CALL ZDISPQ(0,TRIA(OLXNPR ),'TRIDUP')
      IF (IAND(IDPFLG,2**16).EQ.0) CALL ZDISPQ(0,TRIA(OLEAN  ),'TRIDUP')
      IF (IAND(IDPFLG,2**17).EQ.0) CALL ZDISPQ(0,TRIA(OLVBDI ),'TRIDUP')
      IF (IAND(IDPFLG,2**18).EQ.0) CALL ZDISPQ(0,TRIA(OLMBDI ),'TRIDUP')
      IF (IAND(IDPFLG,2**19).EQ.0) CALL ZDISPQ(0,TRIA(OLAREA ),'TRIDUP')
      IF (IAND(IDPFLG,2**20).EQ.0) CALL ZDISPQ(0,TRIA(OLFBVT ),'TRIDUP')
      IF (IAND(IDPFLG,2**21).EQ.0) CALL ZDISPQ(0,TRIA(OLFBMT ),'TRIDUP')
      
C     Clear TRIA

      CALL LCL3(TRIA,SZTRIA)
      
      END

***********************************************************************
* Edge-to-nodes
*
* This routine transforms the number of an edge inside of the geometry
* into the both adjacent nodes on the beginning and the end of that
* edge.
* This routine does not necessarily use the KEAN-array. If KEAN does
* not exist, the vertex information is calculated directly. Thus it 
* ED2NDS is a general alternative for using KEAN directly.
*
* In:
*   TRIA   - array [1..SZTRIA] of integer
*            Triangulation structure
*   IEDG   - number of the edge; 1..NMT
*
* Out:
*   IVT1   - number of the first vertex adjacent to that edge
*   IVT2   - number of the second vertex adjacent to that edge
***********************************************************************

      SUBROUTINE ED2NDS (TRIA,IEDG,IVT1,IVT2)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      INCLUDE 'cmem.inc'
      
      INTEGER TRIA(SZTRIA),IEDG,IVT1,IVT2
      
C local variables

      INTEGER IEL,KMEL,KVERT,KMID,ILED,NVT,KEAN
      
C Do we have extended information to speed up our work?

      IF ((IAND(TRIA(OTRIFL),1).NE.0).AND.(TRIA(OLEAN).NE.0)) THEN

C       We are lucky, now that's easy. Simply use KEAN to get the
C       necessary information

        KEAN = L(TRIA(OLEAN))
        IVT1 = KWORK(KEAN+2*(IEDG-1))
        IVT2 = KWORK(KEAN+2*(IEDG-1)+1)
        RETURN

      END IF
      
C We don't have extended information, so we have to go an alternative
C way :(
      
      KMEL = L(TRIA(OLMEL))
      KVERT = L(TRIA(OLVERT))
      KMID = L(TRIA(OLMID))
      NVT = TRIA(ONVT)
      
C Use KMEL to find out the first element adjacent to that edge

      IEL = KWORK(KMEL+2*(IEDG-1))-1
      
C Loop through the edges of that element to find the local
C edge number 1..4 of out edge

      DO ILED = 0,NNVE-1
        IF ((KWORK(KMID+IEL*NNVE+ILED)-NVT).EQ.IEDG) THEN
          
C Ok, edge was found. Now take the adjacent vertices from KVERT

          IVT1 = KWORK(KVERT+IEL*NNVE+ILED)
          IVT2 = KWORK(KVERT+IEL*NNVE+MOD(ILED+1,NNVE))
          RETURN
          
        END IF
      END DO
      
C Oops, something went wrong. Set IVT1/IVT2 to 0 to indicate

      IVT1 = 0
      IVT2 = 0
      
      END
      
***********************************************************************
* Search boundary node index
*
* This routine accepts a vertex or an edge number of a boundary vertex/
* edge and determines the appropriate index of that vertex in the 
* KVBD/KMBD-array.
* If the enhanced TRIA-structure is set up, this is done by a quick
* bisection search, otherwise a linear search through the KVBD-array
* is performed.
*
* In:
*   TRIA   - array [1..SZTRIA] of integer
*            Triangulation structure
*   IVT    - 1..NVT:     vertex number of a vertex on the boundary
*            NVT+1..NMT: number of an edge on the boundary
* Out:
*   IDX    - index of the vertex/edge in the KVBD/KMBD-array,
*            =-1, if IVT is not a boundary vertex/edge or if IVT was
*                 not found
*   IBC    - number of the boundary component that contains IVT
*            =-1, if IVT is not a boundary vertex/edge or if IVT was
*                 not found
***********************************************************************

      SUBROUTINE SRBDNI (TRIA, IVT, IDX, IBC)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      INCLUDE 'cmem.inc'
      
      INTEGER TRIA(SZTRIA),IVT,IDX,IBC,IRVT,IRVI
      
C local variables

      INTEGER I,J,K,KBCT,KVBD,KEAN,KMBD,KVBDI,KMID
      INTEGER IEL,KEBD
      
      KBCT = L(TRIA(OLBCT))
      KVBD = L(TRIA(OLVBD))
      KEBD = L(TRIA(OLEBD))
      KMBD = L(TRIA(OLMBD))
      KMID = L(TRIA(OLMID))
      
      IDX = -1
      IBC = -1
      
C Is the TRIA-structure enhanced?

      IF ((TRIA(OTRIFL).EQ.0).OR.(TRIA(OLEAN).EQ.0).OR.
     *    (TRIA(OLVBDI).EQ.0)) THEN

C Without the enhanced structures we have no chance of being fast;
C a linear search is to do. Do we have a vertex or an edge?

        IF (IVT.LE.TRIA(ONVT)) THEN
        
C Search the KVBD-array for the node

          DO IBC=1,TRIA(ONBCT)
            DO IDX=KWORK(KBCT+IBC-1),KWORK(KBCT+IBC)-1
              IF (KWORK(KVBD+IDX-1).EQ.IVT) THEN
C Got the vertex
                RETURN
              END IF
            END DO
          END DO
          
C Vertex not found
          IBC = -1
          IDX = -1
          
        ELSE
        
C Search all boundary elements to find the element that contains
C our boundary edge.

          DO IBC=1,TRIA(ONBCT)
            DO IDX=KWORK(KBCT+IBC-1),KWORK(KBCT+IBC)-1

C             The corresponding adjacent element has number (- 1)

              IEL = KWORK(KEBD+IDX-1)-1
              
C             Loop through the edges of that element to find our
C             edge:

              DO J=0,TRIA(ONVE)-1
                IF (KWORK(KMID+IEL*NNVE+J).EQ.(IVT-TRIA(ONVT))) THEN

C                 Hooray, we found the element and the edge. IDX and IBC 
C                 point to where we found the edge - so we can
C                 return now. 
                
                  RETURN
                
                END IF
              END DO
 
C             No, that was the wrong element. Next one...
              
            END DO
          END DO
          
C Edge not found

          IBC = -1
          IDX = -1
          
        END IF      
      
      ELSE
      
C If we should search for an edge, we instead search for the vertex
C preceding this edge. Use KEAN to find out the preceding vertex

        KEAN = L(TRIA(OLEAN))

        IF (IVT.LE.TRIA(ONVT)) THEN
          IRVT = IVT
        ELSE
          IRVT = KWORK(KEAN+2*(IVT-TRIA(ONVT)))
        END IF
        
C IRVT contains the number of the "real" vertex we are searching for.
C In the enhanced structure there's the sorted array KVBDI...

        KVBDI = L(TRIA(OLVBDI))
        
C this array is sorted for vertex numbers. Use bisection search to
C find our node in that array

        I=0
        J=TRIA(ONVBD)-1
        IRVI = 0
        
        IF (KWORK(KVBDI+2*I).EQ.IRVT) THEN
          IRVI = I
        ELSE IF (KWORK(KVBDI+2*J).EQ.IRVT) THEN
          IRVI = J
        ELSE
C Crippled Do-While-Loop
10        CONTINUE
            K = (I+J)/2
            IF (KWORK(KVBDI+2*K).GT.IRVT) THEN
              J = K
            ELSE
              I = K
            END IF
          IF ((I.NE.J).AND.(KWORK(KVBDI+2*K).NE.IRVT)) GOTO 10
          
          IF (KWORK(KVBDI+2*K).EQ.IRVT) THEN
            IRVI = K
          ELSE
            IRVI = -1
          END IF
        END IF
      
        IF (IRVI.NE.-1) THEN

C Vertex found; determine index by the second entry of the structure
          
          IDX = KWORK(KVBDI+2*IRVI+1)
          
C Finally determine the boundary component. This is a linear
C search, but as the number of "real" boundary components is
C typically low, this is quick!

          DO IBC=1,TRIA(ONBCT)
            IF ( (IDX.GE.KWORK(KBCT+IBC-1)).AND.
     *           (IDX.LT.KWORK(KBCT+IBC)) ) RETURN
          END DO
          
C The following does not happen when all structures are ok:

          WRITE (*,*) 'SRBDNI: KBCT-structure destroyed!'
          STOP
        
        ELSE
C not found :(
          IBC = -1
          IDX = -1
        END IF
      
      END IF
      
      END
      
***********************************************************************
* Generate regularly distributed vertices on edges
*
* The routine generates for every edge in the triangulation the
* coordinates of a number of points on the edge. For every edge,
* NVPED vertices will be created in DCORMG on that edge.
* The position of each vertex on each edge will be decided by
* the PAR-array: Let P,Q be the endpoint of edge I an arbitrary edge 
* and s=PAI(I), then GXEDMP creates a vertex at coordinates
* s*Q + (1-s)*P, so PAR describes the "parameter" value of each
* of the NVPED vertices to create on each edge.
*
* In:
*   TRIA   - array [1..SZTRIA] of integer
*            Triangulation structure
*   NVPED  - Number of points to create/calculate on each edge
*   PAR    - array [1..NVPED] of double
*            Parameter value (0..1) of each point to create on
*            each edge
*
* Out:
*  If DCORMG does not exist, it will be cerated.
*  If DCORMG is too large, is will be shrinked.
*  If DCORMG has appropriate size, it will be overwritten.
*  If DCORMG is too small, it will be reallocated.
*  In any case DCORMG is updated to hold on every edge the coordinates
*  of the points specified by PAR(.).
*
* Remark: If the routine is called with NVPED=0, all data of
*  vertices on edges is deleted from DCORMG!
***********************************************************************
      
      SUBROUTINE GXVTED (TRIA,NVPED,PAR)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      INCLUDE 'cmem.inc'
      
C parameters
      
      INTEGER TRIA(SZTRIA),NVPED
      DOUBLE PRECISION PAR(NVPED)
      
C local variables

      INTEGER I,J,KCORMG,IVT1,IVT2
      INTEGER CNVPED,CNVEDT,CNIELV,CNIEVT,CNANT,CNMT
      DOUBLE PRECISION X1,Y1,X2,Y2
      
      CNVPED = TRIA(ONVPED)
      CNVEDT = TRIA(ONVEDT)
      CNIELV = TRIA(ONIELV)
      CNIEVT = TRIA(ONIEVT)
      CNANT  = TRIA(ONANT )
      CNMT   = TRIA(ONMT  )
      
C Is DCORMG existing?

      IF (TRIA(OLCORMG).EQ.0) THEN
      
C       No, create it with the correct size:
        
        IF (CNMT*NVPED.NE.0) THEN
          CALL ZNEW (2*CNMT*NVPED,-1,TRIA(OLCORMG),'DCORMG ')
        END IF
        TRIA(ONVPED) = NVPED
        TRIA(ONVEDT) = CNMT*NVPED
      
      ELSE IF (TRIA(ONVPED).LT.NVPED) THEN
      
C Is DCORMG large enough?
      
C       Reallocate DCORMG and transfer the old data - if DCORMG 
C       should exist at all:

        IF (CNMT*NVPED+CNIEVT+CNANT.GT.0) THEN

          CALL ZNEW (2*(CNMT*NVPED+CNIEVT+CNANT),
     *               -1,J,'DCORMG ')
        
          KCORMG = L(TRIA(OLCORMG))
        
          CALL LCP1 (DWORK(KCORMG+2*CNVEDT),
     *             DWORK(L(J)+2*CNMT*NVPED),
     *             2*(CNIEVT+CNANT))

          TRIA(OLCORMG) = J
     
        ELSE
        
          J = 0
          
        END IF

        CALL ZDISP (0,TRIA(OLCORMG),'DCORMG ')
          
        TRIA(ONVPED)  = NVPED
        TRIA(ONVEDT)  = CNMT*NVPED
      
      ELSE IF (TRIA(ONVPED).LT.NVPED) THEN
      
C Is DCORMG too large?

C       Copy the existing data to the correct position:

        KCORMG = L(TRIA(OLCORMG))

        CALL LCP1 (DWORK(KCORMG+2*CNVEDT),
     *             DWORK(KCORMG+2*CNMT*NVPED),
     *             2*(CNIEVT+CNANT))
     
        TRIA(ONVPED)  = NVPED
        TRIA(ONVEDT)  = CNMT*NVPED

C       Release unnecessary memory

        CALL ZDISP (2*(CNMT*NVPED+CNIEVT+CNANT),
     *              TRIA(OLCORMG),'DCORMG ')

      END IF
      
C Ok, DCORMG is prepared properly.
C Now we can overwrite some parts of DCORMG with the coordinates
C of the vertices on the edges. Loop over all edges:

      KCORMG = L(TRIA(OLCORMG))
      
      DO I=0,CNMT-1
      
C       On each edge loop over the parameter values of the points,
C       given by PAR:

        DO J=0,NVPED-1
        
C         Get the coordinates of the vertices adjacent to that edge:

          CALL ED2NDS (TRIA,I+1,IVT1,IVT2)
          CALL NDE2XY (IVT1,TRIA,X1,Y1)
          CALL NDE2XY (IVT2,TRIA,X2,Y2)
          
C         Write the coordinates of the point into the DCORMG-array
        
          DWORK(KCORMG+2*(NVPED*I+J))   = PAR(J+1)*X2+(1D0-PAR(J+1))*X1
          DWORK(KCORMG+2*(NVPED*I+J)+1) = PAR(J+1)*Y2+(1D0-PAR(J+1))*Y1
        
        END DO
      
      END DO
      
      END 
      
***********************************************************************
* Generate regularly distributed inner-element vertices
*
* The routine generates for every element in the triangulation the
* coordinates of a number of points in the inner of that element,
* e.g. element midpoints,.... For every element,
* NIELV vertices will be created in DCORMG.
* The position of each vertex in each element will be decided by
* the PAR-array. Every element in the PAR-array is a coordinate pair
* of a point on the reference quadrilateral [-1,1]^2, which will
* be transferred via a bilinear transformation onto the actual
* element.
*
* In:
*   TRIA   - array [1..SZTRIA] of integer
*            Triangulation structure
*   NVPED  - Number of points to create/calculate on each element
*   PAR    - array [1..2,1..NIELV] of double
*            (XI1,XI2)=(PAR(1,I),PAR(2,I)) defines coordinates
*            of the I'th vertex in coordinates on the reference 
*            element.
*
* Out:
*  If DCORMG does not exist, it will be cerated.
*  If DCORMG is too large, is will be shrinked.
*  If DCORMG has appropriate size, it will be overwritten.
*  If DCORMG is too small, it will be reallocated.
*  In any case DCORMG is updated to hold on every element the 
*  coordinates of the points specified by PAR(.,.).
*
* Remark: If the routine is called with NIELV=0, all data of
*  vertices on edges is deleted from DCORMG!
***********************************************************************
      
      SUBROUTINE GXIEVT (TRIA,NIELV,PAR)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      INCLUDE 'cmem.inc'
      
C parameters
      
      INTEGER TRIA(SZTRIA),NIELV
      DOUBLE PRECISION PAR(2,NIELV)
      
C local variables

      INTEGER I,J,KCORMG,KVERT
      INTEGER CNVPED,CNVEDT,CNIELV,CNIEVT,CNANT,CNEL
      DOUBLE PRECISION DCOORD(2,4),DJF(2,2),X,Y,DJAC(2,2),DETJ
      
      CNVPED = TRIA(ONVPED)
      CNVEDT = TRIA(ONVEDT)
      CNIELV = TRIA(ONIELV)
      CNIEVT = TRIA(ONIEVT)
      CNANT  = TRIA(ONANT )
      CNEL   = TRIA(ONEL  )

C Is DCORMG existing?

      IF (TRIA(OLCORMG).EQ.0) THEN
      
C       No, create it with the correct size:

        IF (CNVEDT+NIELV*CNEL+CNANT.NE.0) THEN
          CALL ZNEW (2*(CNVEDT+NIELV*CNEL+CNANT),
     *               -1,TRIA(OLCORMG),'DCORMG ')
        END IF
        
        TRIA(ONIELV) = NIELV
        TRIA(ONIEVT) = CNEL*NIELV
      
      ELSE IF (TRIA(ONIELV).LT.NIELV) THEN
      
C Is DCORMG large enough?
      
C       Reallocate DCORMG and transfer the old data - if there should
C       exist that array at all

        IF (CNVEDT+NIELV*CNEL+CNANT.GT.0) THEN
          CALL ZNEW (2*(CNVEDT+NIELV*CNEL+CNANT),
     *               -1,J,'DCORMG ')

C         Copy vertices on edges

          KCORMG =L(TRIA(OLCORMG))

          CALL LCP1 (DWORK(KCORMG),DWORK(L(J)),
     *               2*CNVEDT)
        
C         Copy additional vertices
        
          CALL LCP1 (DWORK(KCORMG+2*(CNVEDT+CNIEVT)),
     *               DWORK(L(J)+2*(CNVEDT+NIELV*CNEL)),
     *               2*CNANT)

          TRIA(OLCORMG) = J
      
        END IF
     
        CALL ZDISP (0,TRIA(OLCORMG),'DCORMG ')
        
        TRIA(ONIELV) = NIELV
        TRIA(ONIEVT) = CNEL*NIELV
      
      ELSE IF (TRIA(ONIELV).LT.NIELV) THEN
      
C Is DCORMG too large?

C       Copy the existing additional vertex data to the correct position:

        KCORMG =L(TRIA(OLCORMG))

        CALL LCP1 (DWORK(KCORMG+2*(CNVEDT+CNIEVT)),
     *             DWORK(KCORMG+2*(CNVEDT+NIELV*CNEL)),
     *             2*CNANT)
     
        TRIA(ONIELV) = NIELV
        TRIA(ONIEVT) = CNEL*NIELV

C       Release unnecessary memory

        CALL ZDISP (2*(CNVEDT+NIELV*CNEL+CNANT),
     *              TRIA(OLCORMG),'DCORMG ')

      END IF
      
C Ok, DCORMG is prepared properly.
C Now we can overwrite some parts of DCORMG with the coordinates
C of the vertices on the elements. Loop over all elements:

      KCORMG = L(TRIA(OLCORMG))
      KVERT =  L(TRIA(OLVERT))
      
      DO I=0,CNEL-1
      
C       On each element grab the four corver vertices for the transformation

        DO J=1,4
          CALL NDE2XY (KWORK(KVERT+I*NNVE)+J-1,
     *                 TRIA,DCOORD(1,J),DCOORD(2,J))
        END DO
        
C       Prepare auxiliary Jacobian factors

        CALL QINIJF (DCOORD,DJF)
        
C       Loop over all "reference coordinates":

        DO J=0,NIELV
        
C         Calculate the transformation if the coordinates given in 
C         PAR(J+1) and store the calculated coordinates in DCORMG.
        
          CALL QTRAF (DCOORD,DJF,DJAC,DETJ,PAR(1,J+1),PAR(2,J+1),X,Y)
        
          DWORK(KCORMG+2*(CNVEDT+NIELV*I+J))   = X
          DWORK(KCORMG+2*(CNVEDT+NIELV*I+J)+1) = Y

        END DO
      
      END DO
      
      END 
      