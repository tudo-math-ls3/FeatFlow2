************************************************************************
* This file contains a couple of routines to handle boundary
* conditions: Implementation of boundary conditions into matrix, RHS-
* vectors, Dirichlet, Neumann, Pressure-drop boundaries,...
*
* To perform its task, user defined callback routines in INDAT2D.F
* are called which have to provide necessary information to the routines
* here.
*
* For a proper handling of geometry information, all routines here
* assume some special data in the user defined parts of the
* triangulation structures. In detail, the user defined block of
* a TRIA structure is assumed to have the following form:
*
*   TRIA.TRIUD(1) = NBDMT
*                 = Number of boundary vertices
*   TRIA.TRIUD(2) = LSCNPR
*                 = Handle to array [1..NMT] of integer = KSCNPR
*     Shortcut nodal property for all edges. This contains a list of
*     all boundary nodes. The first NVBD entries contain the numbers
*     of all boundary edges. Then all fictitious boundary edges
*     follow, attached to that list. NBDMT is the actual size of
*     this array, as it may vary over time if the fictitious boundary
*     changes. As the maximum number of edges in the geometry is NMT,
*     NBDMT<=NMT. The entries follow the following rule:
*      |Entry|   : Number of the edge on the boundary
*      Entry > 0 : The edge is a Dirichlet edge
*      Entry < 0 : The edge is a Neumann edge
*   TRIA.TRIUD(3) = LIARR
*                   Handle to precalculated integer information about
*                   the current geometry or 0 if not used
*   TRIA.TRIUD(4) = LDARR
*                   Handle to precalculated double information about
*                   the current geometry or 0 if not used
************************************************************************

************************************************************************
* Update boundary information
*
* Extended calling convention
*
* This routine updates/marks DIRICHLET and NEUMANN boundary components
* in a triangulation structure. There is a loop through all vertices
* and edges that tests which are DIRICHLET elements, NEUMANN elements
* and fictitious boundary elements.
*
* In:
*   TRIA   - Triangulation structure which KXNPR element is to
*            be build. KXNPR must be in a "clean" state, i.e. without
*            any marks about fictitious boundary or Neumann components.
*            Precalculated information about the geometry is assumed
*            to be stored as follows:
*              TRIA.TRIUD[3] = LIARR
*              TRIA.TRIUD[4] = LDARR
*   TIMENS - Current simulation time; 
*            =0 for stationary simulation
*   RE     - Reynolds-number of the problem
*   IPARAM - array [1..*] of integer 
*   DPARAM - array [1..*] of integer 
*            TIntAssembly/TDoubleAssembly assembly structures; gives
*            additional information about the discretization.
*            This is passed to user defined callback routines so that 
*            they can access more detailed information about the
*            problem. Not used in this routine.
*   IGEOM  - array [1..*] of integer 
*   DGEOM  - array [1..*] of double 
*            Integer- and double-precision parameter blocks with
*            geometry information. Passed to fictitious boundary
*            routines. Not used in this routine.
*
* If either LIARR or LDARR is 0, precalculated information is not used.
*
* Out:
*   INEUM  - >0, if there are Neumann boundaries in the geometry
* 
*   KXNPR  - Is build according to the geometry information.
*            Fictitious boundary Dirichlet vertices/edges and 
*            Neumann edges are marked according to the
*            definition is STRIA.INC.
*   NBDMT  - Returns the number of elements in KSCNPR.
*   KSCNPR - array [1..NMT]
*            The first NBDMT elements are filled with the shortcut
*            nodal property of all boundary edges; see above.
*
* KXNPR may be identical to TRIA.KXNPR on entry of this routine.
* It is modified according to the boundary information. If KXNPR and
* TRIA.KXNPR do not coincide, only the vertices and edges corresponding
* to Neumann boundary parts and Dirichlet fictitious boundary parts
* are modified in KXNPR. TRIA.KXNPR is not touched.
************************************************************************

      SUBROUTINE UPDBDX (TRIA,TIMENS,INEUM,KXNPR,
     *                   NBDMT,KSCNPR, IPARAM,DPARAM,IGEOM,DGEOM)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INTEGER INEUM,IPARAM(*),IGEOM(*)
      INTEGER TRIA(SZTRIA),KXNPR(2,*),NBDMT,KSCNPR(*)
      DOUBLE PRECISION TIMENS,DPARAM(*),DGEOM(*)

C     externals

      INTEGER ISFBDY,ISFBDM,NFBDYC
      EXTERNAL ISFBDY,ISFBDM,NFBDYC

C     local variables

      DOUBLE PRECISION DPARV,DPARM,DPARN1,DPARN2
      INTEGER KMBD,IMBD,KMBDP,KVBD,KVBDP,INPART,IMT,IVT
      INTEGER KADJ,KVERT,KCORVG,KMID,KXNPRO
      INTEGER NPART,INPRN,INPRV,INPRM,IEL,IVE,NEL,IVT1,IVE1,IVT2,IMBCP
      INTEGER IM1,NVT
      DOUBLE PRECISION XX,YY,XX1,YY1,XX2,YY2

      IF (TRIA(OTRIFL).EQ.0) THEN
        WRITE (*,'(A)') 'UPDBDR: Triangulation structure invalid!'
        STOP
      END IF
      
C     At first we handle the standard boundary nodes and mark,
C     which nodes are Neumann boundary nodes:

      INEUM=0
      
      NEL = TRIA(ONEL)
      NVT = TRIA(ONVT)

      KMBD  = L(TRIA(OLMBD))
      KMBDP = L(TRIA(OLMBDP))
      KVBD  = L(TRIA(OLVBD))
      KVBDP = L(TRIA(OLVBDP))
      
      KADJ  = L(TRIA(OLADJ))
      KVERT = L(TRIA(OLVERT))
      KCORVG = L(TRIA(OLCORVG))
      KMID   = L(TRIA(OLMID))
      
C     Get the original KXNPR-array from the TRIA structure

      KXNPRO = L(TRIA(OLXNPR))

C     At first ask the user-defined routines how many
C     Neumann boundary components exist:
      
      NPART = 0
      INPRN = 0
      CALL NEUDAT(NPART,INPRN,DPARN1,DPARN2,TIMENS,
     *            IPARAM,DPARAM)
     
C     NPART is now the number of Neumann components.
C
C     Clear the number of boundary edges

      NBDMT = 0

C     Loop through the vertices(edges on the boundary:

      DO IMBD=0,TRIA(ONVBD)-1
      
C       What's the parameter value of the midpoint of that boundary?

        DPARM = DWORK(KMBDP+IMBD)
        IMT   = KWORK(KMBD+IMBD)

        DPARV = DWORK(KVBDP+IMBD)
        IVT   = KWORK(KVBD+IMBD)
        
C       And what is the boundary component?

        INPRV = KWORK(KXNPRO+2*(IVT-1)+1)
        INPRM = KWORK(KXNPRO+2*(IMT-1)+1)
        
C       Originally, each boundary point/edge is a Dirichlet node.
C       Set the Dirichlet-flag:

        KXNPR(1,IVT) = IOR(KXNPR(1,IVT),2**12)
        KXNPR(1,IMT) = IOR(KXNPR(1,IMT),2**12)

C       Loop through the NPARTS boundary parts and check, if the
C       vertex/edge belongs to a Neumann segment:
        
        DO INPART=1,NPART
        
C         Get the parameter values DPARN1,DPARN2 of the boundary segment
C         from the user-defined routines
        
          CALL NEUDAT(INPART,INPRN,DPARN1,DPARN2,TIMENS,
     *                IPARAM,DPARAM)
     
C         Check if DPARV is inside of that segment. 
     
          IF ((DPARV.GT.DPARN1).AND.(DPARV.LT.DPARN2)
     *                         .AND.(INPRV.EQ.INPRN)) THEN
     
C           Yes, we have at least one Neumann boundary segment:

            INEUM=1
            
C           Set KXNPR of that vertex to indicate the Neumann vertex:
            
            KXNPR(1,IVT) = IAND(KXNPR(1,IVT),NOT(2**12))
            
          ENDIF
          
C         Check if DPARM is inside of that segment. 
     
          IF ((DPARM.GT.DPARN1).AND.(DPARM.LT.DPARN2)
     *                         .AND.(INPRM.EQ.INPRN)) THEN
     
C           Yes, we have at least one Neumann boundary segment:

            INEUM=1
            
C           Set KXNPR of that edge to indicate the Neumann segment:
            
            KXNPR(1,IMT) = IAND(KXNPR(1,IMT),NOT(2**12))
            
          ENDIF
          
        END DO ! INPART
        
C       Now check the KXNPR flag of the edge. Add the edge to the shortcut 
C       nodal property array - either as Neumann or as Dirichlet node:

        NBDMT = NBDMT+1

        IF (IAND(KXNPR(1,IMT),2**12).EQ.0) THEN

C         Add the edge to the list of boundary edges as Neumann boundary

          KSCNPR(NBDMT) = -(IMT-NVT)
          
        ELSE            
          
C         Add the edge to the list of boundary edges as Dirichlet boundary

          KSCNPR(NBDMT) = IMT-NVT
          
        END IF
        
      END DO ! IMBD
      
C     The standard boundary edges are completed. Now handle
C     Fictitious boundary vertices/edges.
C     The array KFBMT is build manually - this is much faster.

C     If there are no fict. boundary components, we skip this...

      IF (NFBDYC(IGEOM,DGEOM).GT.0) THEN
        
C     Loop through all elements and all edges to find the nodes 
C     belonging to elements where a moving boundary goes through...
     
        DO IEL=0,NEL-1
        
          DO IVE=0,NNVE-1
          
C           Only handle edges (with adjacent vertices) adjacent to
C           elements with larger number than the current element.
C           This ensures that all edges are only handles once.
          
            IF (KWORK(KADJ+NNVE*IEL+IVE).LT.(IEL+1)) THEN
      
C             take the X/Y-coordinate of the current vertex
      
              IVT1 = KWORK(KVERT+NNVE*IEL+IVE)
              XX   = DWORK(KCORVG+2*(IVT1-1))
              YY   = DWORK(KCORVG+2*(IVT1-1)+1)
        
C             Test, if the point is in a fictitious boundary component. If yes,
C             set KXNPR(1,.) = 1 and KXNPR(2,.) to the number of the fictitious
C             boundary component.
C
C             Remark: We first handle the (corner) vertices, then the midpoints.
C             Only (corner) vertices are exported later to GMV files.
C             Using the KXNPR-array he GMV output routine can decide which points
C             belong to fictitious boundary components to "stamp out" this region.
C             This has nothing to do with any calculation! As we are using
C             nonconforming finite elements that are midpoint oriented, only
C             the KXNPR-information of the midpoints are concerned during
C             the modification the matrix!

              IF ((TRIA(OTRIUD+2).NE.0).AND.(TRIA(OTRIUD+3).NE.0)) THEN

C               When using precalculated information, LIARR/LDARR is stored
C               in the TRIA-structure at TRIA.TRIUD[3..4]:
              
                IMBCP = ISFBDM (IVT1, TRIA, 
     *                          TRIA(OTRIUD+2), TRIA(OTRIUD+3), 0,
     *                          IGEOM,DGEOM)
     
              ELSE
              
C               No precalculated information
              
                IMBCP = ISFBDY (XX,YY,0,IGEOM,DGEOM)
                
              END IF
              
C             Originally, each vertex does not belong to a fictitious
C             boundary component...

              KXNPR(1,IVT1) = IAND(KXNPR(1,IVT1),NOT(2**13))
              
C             And if the vertex is an inner vertex, it's therefore also
C             not of Dirichlet type:
              
              IF (IAND(KWORK(KXNPRO+2*(IVT1-1)),2**8).EQ.0) 
     *          KXNPR(1,IVT1) = IAND(KXNPR(1,IVT1),NOT(2**12))

C             IMCP <> 0 indicates a fictitious boundary node.
C             IMCP > 0 indicates a Dirichlet fictitious boundary node.
C             IMCP < 0 indicates a Neumann fictitious boundary node.
              
              IF (IMBCP.NE.0) THEN
              
C               We have a fictitious boundary node
              
                KXNPR(1,IVT1) = IOR(KXNPR(1,IVT1),2**13)
                
C               Don't change the KXNPR(2,.) if the vertex is a boundary vertex!
C               Real boundary information has higher priority than fictitious
C               boundary information!

                IF (IAND(KWORK(KXNPRO+2*(IVT1-1)),2**8).EQ.0) 
     *            KXNPR(2,IVT1) = ABS(IMBCP)
     
C               Is that even a Dirichlet node?

                IF (IMBCP.GT.0) THEN
                  KXNPR(1,IVT1) = IOR(KXNPR(1,IVT1),2**12)
                END IF
     
              ENDIF

C             Take a look to the midpoint of the current element, following the
C             current vertex. 
C             Warnung: For efficienty:
C               IVT2, IVE1, IVE2 - 0-based
C               IVT1             - 1-based

              IVE1=IVE+1
              IF (IVE1.EQ.4) IVE1=0
              IVT2=KWORK(KVERT+NNVE*IEL+IVE1)-1
              XX1=DWORK(KCORVG+2*IVT2)
              YY1=DWORK(KCORVG+2*IVT2+1)
              XX2=0.5D0*(XX+XX1)
              YY2=0.5D0*(YY+YY1)

C             Get the number of that edge

              IM1=KWORK(KMID+NNVE*IEL+IVE)    
      
C             Test also for this edge if it's in a fictitious boundary
              
              IF ((TRIA(OTRIUD+2).NE.0).AND.(TRIA(OTRIUD+3).NE.0)) THEN
                IMBCP = ISFBDM (IM1, TRIA, 
     *                          TRIA(OTRIUD+2), TRIA(OTRIUD+3), 0,
     *                          IGEOM,DGEOM)
              ELSE
                IMBCP = ISFBDY (XX2,YY2,0,IGEOM,DGEOM)
              END IF

C             Originally, the edge does not belong to a fictitious
C             boundary component...

              KXNPR(1,IM1) = IAND(KXNPR(1,IM1),NOT(2**13))
              
C             And if the edge is an inner edge, it's therefore also
C             not of Dirichlet type:
              
              IF (IAND(KWORK(KXNPRO+2*(IM1-1)),2**8).EQ.0) 
     *          KXNPR(1,IM1) = IAND(KXNPR(1,IM1),NOT(2**12))

C             If this edge is in a fictitious boundary component, set
C             KXNPR appropriately.

              IF (IMBCP.NE.0) THEN
              
                KXNPR(1,IM1) = IOR(KXNPR(1,IM1),2**13)
                
C               Is that even a Dirichlet edge?

                IF (IMBCP.GT.0) THEN
                  KXNPR(1,IM1) = IOR(KXNPR(1,IM1),2**12)
                END IF
                
C               If the edge is an inner edge, we do some more work:

                IF (IAND(KWORK(KXNPRO+2*(IM1-1)),2**9).EQ.0) THEN
                
C                 Save the fict. boundary component                
      
                  KXNPR(2,IM1) = ABS(IMBCP)

C                 Add the edge to the list of edges - as Dirichlet edge 
C                 (KBDMT>0) if it is one. If the edge os of Neumann-type,
C                 or if the edge is part of the real boundary, 
C                 don't add it to the list, so the edge is treated as
C                 standard inner edge. (The last case prevents
C                 we add more than NMT edges altogether to KSCNPR
C                 by adding some edges twice - above and here!)
C
C                 Note that if the edge is of Neumann type, we mustn't
C                 add it to KSCNPR with a value < 0! This might lead to
C                 unpredictable (or more concrete: slightly wrong)
C                 results in the solution, as there might be a special
C                 handling of boundary edges of Neumann type in
C                 modifying matrices etc. - this does not hold for
C                 inner edges!

                  IF ((IMBCP.GT.0).AND.
     *            (IAND(KXNPR(1,IVT1),2**8).NE.2**8)) THEN
                    NBDMT = NBDMT+1
                    KSCNPR(NBDMT) = IM1-NVT
                  END IF
               
                END IF ! ( and(KXNPRO(1,IM1),2**8) = 0 )

              END IF ! (IMBCP.GT.0)

            END IF ! (KWORK(KADJ+NNVE*IEL+IVE).LT.IEL)
            
          END DO ! IVE
          
        END DO ! IEL
        
      END IF ! (NFBDYC().GT.0)
        
      END

************************************************************************
* Implement PRESSURE DROP values into the velocity RHS-vectors DF1/DF2.
*
* Extended calling convention.
*
* In:
*   DF1,
*   DF2    - array [1..*] of double
*            X- any Y-RHS vectors which are to be modified. These
*            are assumed to be assembled by E030/EM30/E031/E031 finite
*            element.
*   TIMENS - Current simulation time
*   TSTEPB - Length of current time step in time stepping scheme.
*   RE     - Reynolds number; from the DAT file
*   TRIA   - Triangulation structure which KXNPR element is to
*            be build. KXNPR must be in a "clean" state, i.e. without
*            any marks about fictitious boundary or Neumann components
*   IPARAM - array [1..*] of integer 
*   DPARAM - array [1..*] of integer 
*            TIntAssembly/TDoubleAssembly assembly structures; gives
*            additional information about the discretization.
*            This is passed to user defined callback routines so that 
*            they can access more detailed information about the
*            problem. Not used in this routine.
*   IGEOM  - array [1..*] of integer 
*   DGEOM  - array [1..*] of double 
*            Integer- and double-precision parameter blocks with
*            geometry information. Passed to boundary
*            routines. Not used in this routine.
*
* Out:
*   DF1,
*   DF2    - Modified RHS vectors
************************************************************************

      SUBROUTINE PDSETX (DF1,DF2,TIMENS,TSTEPB,RE,TRIA,IPARAM,DPARAM,
     *                   IGEOM,DGEOM)

C     ??? Remark: This should be rewritten using real line integration
C     instead of taking the midpoint rule... Would be much more exact.
C     A cubature formula should then be given by a parameter.
C     The boundary values can then be interpreted as "normal stress
C     in a point" instead of "mean pressure along an edge".

      IMPLICIT NONE
      
C main COMMON blocks

      INCLUDE 'cmem.inc'

      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'

C parameters 

      DOUBLE PRECISION DF1(*),DF2(*)
      INTEGER IPARAM(*),TRIA(SZTRIA),IGEOM(*)
      DOUBLE PRECISION DPARAM(*),TSTEPB,TIMENS,RE,DGEOM(*)

C local variables

      INTEGER KEBD,KVERT,KMID,KXNPR

      INTEGER IMID, IV1, IV2, INPR, INPART, INPRN, NPART
      INTEGER KCORVG,KMBDP,IMBD,IMT,KMBD,INPRM,KEAN,NVT
      DOUBLE PRECISION DPAR, DPARN1, DPARN2, PX1, PY1, PX2, PY2
      DOUBLE PRECISION DN1, DN2, PMEAN, DPARM,FDATIN
      EXTERNAL FDATIN
      
C     Dereference some handles

      KMBD   = L(TRIA(OLMBD))
      KEAN   = L(TRIA(OLEAN))
      KEBD   = L(TRIA(OLEBD))
      KVERT  = L(TRIA(OLVERT))
      KMID   = L(TRIA(OLMID))
      KXNPR  = L(TRIA(OLXNPR))
      KCORVG = L(TRIA(OLCORVG))
      KMBDP  = L(TRIA(OLMBDP))
      NVT    = TRIA(ONVT)
      
C     At first ask the user-defined routines how many
C     Neumann boundary components exist:
      
      NPART = 0
      INPRN = 0
      CALL NEUDAT(NPART,INPRN,DPARN1,DPARN2,TIMENS,
     *            IPARAM,DPARAM)

C     NPART is now the number of Neumann components.
C     
C     Loop through the vertices(edges) on the boundary:

      DO IMBD=0,TRIA(ONVBD)-1

C       What's the parameter value of the edge of that boundary
C       (i.e. the midpoint of the edge)?

        DPARM = DWORK(KMBDP+IMBD)
        IMT   = KWORK(KMBD+IMBD)
        
C       And what is the boundary component?

        INPRM = KWORK(KXNPR+2*(IMT-1)+1)
          
C       Loop through the NPARTS boundary parts and check, if the
C       edge belongs to a Neumann segment:
        
        DO INPART=1,NPART
        
C         Get the parameter values DPARN1,DPARN2 of the boundary segment
C         from the user-defined routines

          CALL NEUDAT(INPART,INPRN,DPARN1,DPARN2,TIMENS,IPARAM,DPARAM)

C         Check if DPARV is inside of that segment. 
          
          IF ((DPAR.GT.DPARN1).AND.(DPAR.LT.DPARN2)
     *                        .AND.(INPR.EQ.INPRN)) THEN
     
C           Get the vertices adjacent to the edge:

            IV1  = DWORK(KEAN+2*(IMT-NVT-1))
            IV2  = DWORK(KEAN+2*(IMT-NVT-1)+1)
     
C           Calculate the outer normal vector to the edge (IV1,IV2):
     
            PX1  = DWORK(KCORVG+2*IV1)
            PY1  = DWORK(KCORVG+2*IV1+1)
            PX2  = DWORK(KCORVG+2*IV2)
            PY2  = DWORK(KCORVG+2*IV2+1)
            DN1  =-PY2+PY1
            DN2  = PX2-PX1
            
C           Call the user defined callback routine to get the
C           mean pressure value:
            
            PMEAN = FDATIN(7,INPR,DPAR,DPAR,TIMENS,RE,
     *              IMT,TRIA,IPARAM,DPARAM,IGEOM,DGEOM)
     
C           Include that into the RHS vectors - cf. p. 257 (235) in 
C           Turek's book.
C           The pressure drop condition is a time independent
C           right hand side condition! Therefore it's not dependent
C           on the THETA-value of any THETA-Scheme, but only
C           dependent on the length of the current time step.
     
            DF1(IMID) = DF1(IMID)+PMEAN*DN1*TSTEPB
            DF2(IMID) = DF2(IMID)+PMEAN*DN2*TSTEPB
            
          ENDIF
          
        END DO ! INPART

      END DO ! IMBD

      END

************************************************************************
* Implement DIRICHLET values into the vector (DX1,DX2).
*
* Extended calling convention
*
* This routine can be used for updating the velocity solution
* vectors DU1/DU2 and velocity RHS vectors DF1/DF2 with the fixed
* values on the Dirichlet boundary.
*
* It's assumed that DF1/DF2 is assembled with E030/E031/EM30/EM31.
*
* In:
*   DX1,
*   DX2    - array [1..*] of double
*            X- any Y-velocity vectors which are to be modified.
*            Both vectors are assumed to correspond to values in the
*            midpoints of edges!
*   TRIA   - Triangulation structure which KXNPR element is to
*            be build. KXNPR must be in a "clean" state, i.e. without
*            any marks about fictitious boundary or Neumann components
*   PARX,
*   PARY   - Subroutines of the used parametrization to retrieve
*            X/Y-coordinate of a parameter value.
*   UE     - Subroutine that must return the Dirichlet function value
*            in a point X/Y.
*   TIMENS - Current time in instationary Navier-Stokes calculation
*   RE     - Reynolds number; from the DAT file
*   IPARAM - array [1..*] of integer 
*   DPARAM - array [1..*] of integer 
*            TIntAssembly/TDoubleAssembly assembly structures; gives
*            additional information about the discretization.
*            This is passed to user defined callback routines so that 
*            they can access more detailed information about the
*            problem. Not used in this routine.
*   IGEOM  - array [1..*] of integer 
*   DGEOM  - array [1..*] of double 
*            Integer- and double-precision parameter blocks with
*            geometry information. Passed to boundary
*            routines. Not used in this routine.
*
* Return:
*   DX1, 
*   DX2    - Modified vectors. The Dirichlet values are implemented.
************************************************************************

      SUBROUTINE BDRSTX (DX1,DX2,KXNPR,
     *                   PARX,PARY,UE,TIMENS,RE,
     *                   TRIA,IPARAM,DPARAM,IGEOM,DGEOM)

      IMPLICIT NONE
      
C main COMMON blocks

      INCLUDE 'cmem.inc'

      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'

C parameters      
      
      DOUBLE PRECISION DX1(*),DX2(*),DGEOM(*)
      DOUBLE PRECISION DPARAM(*),TIMENS,RE
      INTEGER TRIA(*),IPARAM(*),KXNPR(2,*),IGEOM(*)
      
C externals

      DOUBLE PRECISION PARX, PARY, UE, FBINDT,FBMIND
      EXTERNAL PARX,PARY,UE,FBINDT,FBMIND 

C local variables

      INTEGER NVT, NEL,IADJEL
      INTEGER KCORVG,KADJ,KMID,KVERT
      INTEGER INPR, IEL, IVE, IVT1, IVT2, IVE1, IM1
      DOUBLE PRECISION U1, U2, XX, YY, XX1, YY1, XX2, YY2

      NVT = TRIA(ONVT)
      NEL = TRIA(ONEL)
      
C     Dereference some handles
      
      KCORVG = L(TRIA(OLCORVG))
      KADJ   = L(TRIA(OLADJ))
      KVERT  = L(TRIA(OLVERT))
      KMID   = L(TRIA(OLMID))

C Process all vertices -> loop through all elements and on
C each element through all vertices.

      DO IEL = 0,NEL-1
      
        DO IVE = 0,NNVE-1
        
          IADJEL = KWORK(KADJ+IEL*NNVE+IVE)
          
          IF (IADJEL.LT.(IEL+1)) THEN
            
C           Get the corner vertex
            
            IVT1 = KWORK(KVERT+IEL*NNVE+IVE)-1
            XX   = DWORK(KCORVG+2*IVT1)
            YY   = DWORK(KCORVG+2*IVT1+1)
            
            IVE1=IVE+1
            IF(IVE1.EQ.4) IVE1=0
            
C           Calculate the midpoint of the edge
            
            IVT2 = KWORK(KVERT+IEL*NNVE+IVE1)-1
            XX1  = DWORK(KCORVG+2*IVT2)
            YY1  = DWORK(KCORVG+2*IVT2+1)
            XX2  = 0.5D0*(XX+XX1)
            YY2  = 0.5D0*(YY+YY1)
            
            IM1  = KWORK(KMID+IEL*NNVE+IVE)
            
C           Get the nodal property identifier:
            
            INPR = KXNPR(1,IM1)
            
C           Should we handle it as Dirichlet?

            IF (IAND(INPR,2**12).NE.0) THEN
            
              IF (IAND(INPR,2**13).NE.0) THEN
              
C               Fictitious boundary Dirichlet node.
C               Check if we have precalculated data:

                IF ((TRIA(OTRIUD+2).NE.0).AND.
     *              (TRIA(OTRIUD+3).NE.0)) THEN
     
                  U1 = FBMIND (1, IM1, TIMENS,RE, TRIA, 
     *                         TRIA(OTRIUD+2), TRIA(OTRIUD+3),
     *                         IPARAM,DPARAM,IGEOM,DGEOM)
     
                  U2 = FBMIND (2, IM1, TIMENS,RE, TRIA, 
     *                         TRIA(OTRIUD+2), TRIA(OTRIUD+3),
     *                         IPARAM,DPARAM,IGEOM,DGEOM)

                ELSE

                  U1 = FBINDT (1,XX2,YY2,TIMENS,RE,IM1,TRIA,
     *                         IPARAM,DPARAM,IGEOM,DGEOM)
                  U2 = FBINDT (2,XX2,YY2,TIMENS,RE,IM1,TRIA,
     *                         IPARAM,DPARAM,IGEOM,DGEOM)    
     
                END IF
                
                DX1(IM1-NVT)=U1
                DX2(IM1-NVT)=U2
              
              ELSE 
            
C               General Dirichlet-Node; use analytic function

                U1 = UE(XX2,YY2,1,TIMENS,RE,IM1,TRIA,IPARAM,DPARAM,
     *                  IGEOM,DGEOM)
                U2 = UE(XX2,YY2,2,TIMENS,RE,IM1,TRIA,IPARAM,DPARAM,
     *                  IGEOM,DGEOM)
                DX1(IM1-NVT) = U1
                DX2(IM1-NVT) = U2
              
              END IF
              
            END IF
            
          END IF ! (IADJEL.LT.IEL)
          
        END DO ! IVE
        
      END DO ! IEL

      END

************************************************************************
* Multiply the NEUMANN-components of the vector DX with A1
*
* This will search the geometry for all nodes that belong to Neumann
* boundary of real boundary components. Fictitious boundary Neumann
* nodes are ignored. The corresponding entries in the vector
* DX (which is assumed to be assembled with E030/E031/EM30/EM31) are
* multiplied by A1.
*
* In:
*   DX     - source vector
*   NBDMT  - Number of edges on the real boundary
*   KSCNPR - array [1..NBDMT] of integer
*            The first NBDMT elements are filled with the shortcut
*            nodal property of all boundary edges.
*   A1     - double
*
* Out:
*   DX     - modified vector
************************************************************************

      SUBROUTINE BDRDEF (DX,KSCNPR,NBDMT,A1)

      IMPLICIT NONE
      
C parameter

      INTEGER NBDMT
      DOUBLE PRECISION DX(*),A1
      INTEGER KSCNPR(*)
      
C     local variables

      INTEGER IMBD, IMID     

C     loop over all boundary vertices

      DO IMBD=1,NBDMT
      
C       Get the shortcut nodal property
      
        IMID=KSCNPR(IMBD)
        
C       Values > 0 are Dirichlet, values < 0 Neumann boundary edges
        
        IF (IMID.LT.0) THEN
          DX(-IMID) = A1*DX(-IMID)
        END IF
        
      END DO ! IMBD
      
      END

************************************************************************
* Multiply rows corresponding to the NEUMANN-components of 
* the matrix DA with A1
*
* This will search the geometry for all nodes that belong to Neumann
* boundary of real boundary components. Fictitious boundary Neumann
* nodes are ignored. The corresponding entries in the matrix
* DA (which is expected in matrix structure 7 or 9) are
* multiplied by A1.
*
* In:
*   DA,
*   KLD    - matrix to be modified
*   NBDMT  - Number of edges on the real boundary
*   KSCNPR - array [1..NBDMT] of integer
*            The first NBDMT elements are filled with the shortcut
*            nodal property of all boundary edges.
*   A1     - double
*
* Out:
*   DX     - modified vector
************************************************************************

      SUBROUTINE BDRDFA (DA,KLD,KSCNPR,NBDMT,A1)

      IMPLICIT NONE
      
C parameter

      INTEGER NBDMT
      DOUBLE PRECISION DA(*),A1
      INTEGER KSCNPR(*),KLD(*)
      
C     local variables

      INTEGER IMBD, IMID, I

C     loop over all boundary vertices

      DO IMBD=1,NBDMT
      
C       Get the shortcut nodal property
      
        IMID=KSCNPR(IMBD)
        
C       Values > 0 are Dirichlet, values < 0 Neumann boundary edges
        
        IF (IMID.LT.0) THEN
        
C         Where is the corresponding line?
C         Multiply the entries...
      
          DO I=KLD(-IMID),KLD(-IMID+1)-1
            DA(I) = A1*DA(I)
          END DO
          
        END IF
        
      END DO ! IMBD
      
      END

************************************************************************
* Updates the matrix entries for all DIRICHLET boundary nodes.
*
* Replaces all matrix lines corresponding to DIRICHLET nodes by
* unit vectors.
*
* In:
*   DA,
*   KCOL,
*   KLD    - System matrix, to be modified
*   NBDMT  - Number of edges on the real boundary
*   KSCNPR - array [1..NMT] of integer
*            The first NBDMT elements are filled with the shortcut
*            nodal property of all boundary edges.
*
* Out:
*   DA     - modified matrix
************************************************************************

      SUBROUTINE BDRYA (DA,KCOL,KLD,KSCNPR,NBDMT)

      IMPLICIT NONE
      
C main COMMON blocks

      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'

C parameters 

      DOUBLE PRECISION DA(*)
      INTEGER KCOL(*),KLD(*),NBDMT,KSCNPR(NBDMT)

C local variables

      INTEGER IMBD, IMID, ICOL

      DO IMBD=1,NBDMT
      
C       Get the shortcut nodal property

        IMID = KSCNPR(IMBD)

C       Values > 0 are Dirichlet, values < 0 Neumann boundary edges

        IF (IMID.GT.0) THEN

C         The number IMID is the number of the edge, similar to KMBD.
C
C         The diagonal element is set to 1, the other elements
C         are set to 0.
C         We don't use LCL1 here, since the number of elements
C         to set to 0 is usually so small, that even an LCL1 would
C         take too much time.

          DA(KLD(IMID)) = 1D0

          DO ICOL=KLD(IMID)+1,KLD(IMID+1)-1
            DA(ICOL) = 0D0
          END DO
          
        END IF

      END DO
      
      END

************************************************************************
* Updates the matrix entries of a global matrix for all 
* DIRICHLET boundary nodes.
*
* Replaces all matrix lines corresponding to DIRICHLET nodes by
* unit vectors. The matrix is assumed to be a global matrix, given in
* matrix structure 9. The first NEQU lines correspond to the X-velocity,
* the second NEQU to the Y-velocity and the rest for pressure.
*
* In:
*   NEQU   - Length of each velocity component in the global matrix
*   DA,
*   KCOL,
*   KLD    - Global system matrix, to be modified
*   NBDMT  - Number of edges on the real boundary
*   KSCNPR - array [1..NMT] of integer
*            The first NBDMT elements are filled with the shortcut
*            nodal property of all boundary edges.
*
* Out:
*   DA     - modified matrix
************************************************************************

      SUBROUTINE BDRYGA (NEQU,DA,KCOL,KLD,KSCNPR,NBDMT)

      IMPLICIT NONE
      
C main COMMON blocks

      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'

C parameters 

      DOUBLE PRECISION DA(*)
      INTEGER KCOL(*),KLD(*),NBDMT,KSCNPR(NBDMT),NEQU

C local variables

      INTEGER IMBD, IMID, ICOL

      DO IMBD=1,NBDMT
      
C       Get the shortcut nodal property

        IMID = KSCNPR(IMBD)

C       Values > 0 are Dirichlet, values < 0 Neumann boundary edges

        IF (IMID.GT.0) THEN

C         The number IMID is the number of the edge, similar to KMBD.
C
C         The diagonal element is set to 1, the other elements
C         are set to 0.

C         Modify first velocity part:

          DO ICOL=KLD(IMID),KLD(IMID+1)-1
            IF (KCOL(ICOL).EQ.IMID) THEN
              DA(ICOL) = 1D0
            ELSE
              DA(ICOL) = 0D0
            END IF
          END DO
          
C         Modify second velocity part:

          DO ICOL=KLD(IMID+NEQU),KLD(IMID+NEQU+1)-1
            IF (KCOL(ICOL).EQ.(IMID+NEQU)) THEN
              DA(ICOL) = 1D0
            ELSE
              DA(ICOL) = 0D0
            END IF
          END DO
          
        END IF

      END DO
      
      END

************************************************************************
* Set the DIRICHLET-components of the vector (D1,D2) to zero
*
* This is typically used to force entries in the defect vector 
* corresponding to Dirichlet nodes to zero.
*
* In:
*   D1     - array [1..*] of double
*   D2     - array [1..*] of double
*   NBDMT  - Number of edges on the real boundary
*   KSCNPR - array [1..NMT] of integer
*            The first NBDMT elements are filled with the shortcut
*            nodal property of all boundary edges.
*
* Out:
*   D1,
*   D2     - modified vectors
************************************************************************

      SUBROUTINE BDRY0 (D1,D2,KSCNPR,NBDMT)

      IMPLICIT NONE
      
C main COMMON blocks

      INCLUDE 'cbasictria.inc'

C parameters

      DOUBLE PRECISION D1(*),D2(*)
      INTEGER NBDMT,KSCNPR(NBDMT)

C local variables

      INTEGER IMBD, IMID

      DO IMBD = 1,NBDMT

C       Get the shortcut nodal property

        IMID = KSCNPR(IMBD)

C       Values > 0 are Dirichlet, values < 0 Neumann boundary edges

        IF (IMID.GT.0) THEN
          D1(IMID) = 0D0
          D2(IMID) = 0D0
        END IF
        
      END DO

      END
