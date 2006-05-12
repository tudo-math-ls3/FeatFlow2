************************************************************************
* This file contains routines to support tracers, which can be
* written to a GMV file in the postprocessing.
*
* "Tracers" are points inside of the domain with no mass that move
* with the fluid. They are characterised by different quantities as
* described in the "TRACERS.INC" include file.
************************************************************************

************************************************************************
* Initialise tracer structure
*
* This routine allocates a tracer structure for NTRCS tracers on the
* heap. All tracers are initialised in position (0,0) with no velocity.
*
* In:
*   NTRCS  : Number of tracers
*
* Out:
*   LTRCS  : Handle to TDTracers structure for all the tracers on the
*            heap.
*
* To clean up with tracers, the caller can simply delete the allocated
* structure from the heap!
************************************************************************

      SUBROUTINE INTRCS (NTRCS,LTRCS)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'stracers.inc'
      
      INTEGER NTRCS,LTRCS
      
C     Cancel if < 1 tracer

      IF (NTRCS.LT.1) RETURN
      
C     Allocate the structure on the heap:

      CALL ZNEW(SZTRCS + SZTRCR*NTRCS,1,LTRCS,'TRCS  ')
      
C     Initialise number of tracers

      DWORK(L(LTRCS)+OTRCNT-1) = DBLE(NTRCS)
      
      END
      
************************************************************************
* Initialise tracer structure for regularly distributed tracers
*
* This routine allocates a tracer structure for NTRCS tracers on the
* heap. Each element receives NTRCPE x NTRCPE tracers.
*
* In:
*   NTRCPE : SQRT(Number of tracers per element).
*   DCORVG,
*   KVERT,
*   NEL    : Geometry information; used to position the tracers
*            
*
* Out:
*   LTRCS  : Handle to TDTracers structure for all the tracers on the
*            heap.
*
* To clean up with tracers, the caller can simply delete the allocated
* structure from the heap!
************************************************************************

      SUBROUTINE INRTRC (NTRCPE,LTRCS,DCORVG,KVERT,NEL)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'stracers.inc'
      INCLUDE 'cbasictria.inc'
      
      INTEGER NTRCPE,LTRCS,KVERT(NNVE,*),NEL
      DOUBLE PRECISION DCORVG(2,*)
      
C     local variables

      INTEGER IEL,I,J,ITRC,KTRC
      DOUBLE PRECISION DCOORD (2,4),DJF(2,2),DLEN,X,Y,DJAC(2,2),DETJ
      DOUBLE PRECISION DOFS
      
C     Cancel if < 1 tracer

      IF (NTRCPE.LT.1) RETURN
      
C     Allocate the structure on the heap:

      CALL ZNEW(SZTRCS + SZTRCR*NTRCPE*NTRCPE*NEL,1,LTRCS,'TRCS  ')
      
C     Initialise number of tracers

      DWORK(L(LTRCS)+OTRCNT-1) = DBLE(NTRCPE*NTRCPE*NEL)
      DWORK(L(LTRCS)+OTRPEL-1) = DBLE(NTRCPE)
      
C     We place NTRCPE horizontal and vertical tracers in the
C     reference element and map them to the real element.
C     The reference element [-1,1^]^2 can therefore be divided into a
C     "grid" with edge length...

      DLEN = 2D0/DBLE(NTRCPE)
      
C     DOFS is added as a "boarder" to all X- and Y-coordinates on
C     the reference element:
      
      DOFS = -1D0-0.5D0*DLEN
      
C     ITRC counts the number of initialised tracers

      ITRC = 0
      
C     Loop through the elements

      DO IEL=1,NEL
      
C       Calculate auxiliary Jacobian factors from the transformation
C       of the reference element to the real element

        DO I=1,4
          J = KVERT(I,IEL)
          DCOORD(1,I) = DCORVG(1,J)
          DCOORD(2,I) = DCORVG(2,J)
        END DO
        CALL QINIJF (DCOORD,DJF)
        
C       On each element, calculate the tracer coordinates

        DO I=1,NTRCPE
          DO J=1,NTRCPE

            CALL QTRAF (DCOORD,DJF,DJAC,DETJ,
     *                  DOFS+DBLE(I)*DLEN,DOFS+DBLE(J)*DLEN,X,Y)

C           Where is the tracer structure?

            KTRC = L(LTRCS) + OTRDEF-1 + ITRC*SZTRCR - 1
            
C           Initialise the tracer structure

            DWORK(KTRC+OTRXPOS) = X
            DWORK(KTRC+OTRYPOS) = Y
            DWORK(KTRC+OTRELEM) = DBLE(IEL)
            
C           Next tracer

            ITRC = ITRC+1
            
          END DO
        END DO
      
      END DO
      
      END
      
************************************************************************
* Move tracers
*
* This routine moves all the tracers in a tracer structure according
* to their current velocity and the given time step. Tracers outside 
* of the domain are not touched. Tracers inside of the domain are moved 
* in a way to prevent leaving the domain through Dirichlet boundaries.
* If a tracer leaves the Domain through a Neumann boudary segment, 
* the element status of the tracer is set to -1.
*
* In:
*   LTRCS  : Handle to TDTracers structure on the heap
*   TRIAS  : Triangulation structure on all levels
*   NLMIN  : Minimum level in TRIAS
*   NLMAX  : Maximum level in TRIAS
*   NEQ    : Number of equations in the X- and Y-velocity field
*   TSTEP  : Time step
*
* Out:
*   The tracers in the TSTracers structure are moved by an Explicit
*   Euler approach according to the velocity.
************************************************************************
      
      SUBROUTINE TRCMOV(LTRCS,TRIAS,NLMIN,NLMAX,TSTEP)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'stracers.inc'
      INCLUDE 'stria.inc'
      
C     parameters
      
      INTEGER LTRCS,NLMIN,NLMAX,TRIAS(SZTRIA,*)
      DOUBLE PRECISION TSTEP
      
C     local variables

      INTEGER I,KTRC,IEL,IELOLD,IEL2
      DOUBLE PRECISION X,Y,XV,YV,XNEW,YNEW,X0,Y0
      INTEGER IRES,IEDG,NVT,KXNPR
      DOUBLE PRECISION X1,Y1,X2,Y2
      
C     Where are the tracers?

      KTRC = L(LTRCS)+OTRDEF-1
      
      NVT = TRIAS(ONVT,NLMAX)
      KXNPR = L(TRIAS(OLXNPR,NLMAX))
      
C     Loop through all the tracers

      DO I=0,FLOOR(DWORK(L(LTRCS)+OTRCNT-1))-1
      
C       Where is the tracer structure?

        KTRC = L(LTRCS) + OTRDEF-1 + I*SZTRCR - 1
      
C       What's the element of that tracer? Where is it?

        IELOLD = FLOOR(DWORK(KTRC+OTRELEM))
        X   = DWORK(KTRC+OTRXPOS)
        Y   = DWORK(KTRC+OTRYPOS)
        XV  = DWORK(KTRC+OTRXVEL)
        YV  = DWORK(KTRC+OTRYVEL)
        
        IF (IELOLD.GE.0) THEN
        
C         The tracer is in the domain. Make sure we have the correct
C         element that contains the tracer:

          IF (IELOLD.EQ.0) THEN
          
C           Element completely unknown - hierarchical search

            CALL PSRCH5(TRIAS,NLMIN,NLMAX,X,Y,IELOLD,IEL)
            
          ELSE
          
C           Raytracing Search, beginning at element IELOLD

            CALL PSRCH4(X,Y,
     *                  DWORK(L(TRIAS(OLCORVG,NLMAX))),
     *                  KWORK(L(TRIAS(OLVERT,NLMAX))),
     *                  KWORK(L(TRIAS(OLMID,NLMAX))),
     *                  KWORK(L(TRIAS(OLADJ,NLMAX))),
     *                  NVT,
     *                  IELOLD,IEL,
     *                  IRES,IEDG,X1,Y1,X2,Y2)
     
            IF (IRES.NE.0) IEL = 0
          
          END IF
          
C         If the element was not found, the point is not in the
C         domain.

          IF (IEL.EQ.0) THEN
            DWORK(KTRC+OTRXVEL) = 0D0
            DWORK(KTRC+OTRYVEL) = 0D0
            DWORK(KTRC+OTRELEM) = DBLE(-1)
          ELSE
          
C           Ok, we have the element. Now move the point
C           with explicit Euler:

            XNEW = X + XV*TSTEP
            YNEW = Y + YV*TSTEP

C           Find the new element of the point. Have we left the 
C           domain?

            IELOLD = IEL          
            CALL PSRCH4(XNEW,YNEW,
     *                  DWORK(L(TRIAS(OLCORVG,NLMAX))),
     *                  KWORK(L(TRIAS(OLVERT,NLMAX))),
     *                  KWORK(L(TRIAS(OLMID,NLMAX))),
     *                  KWORK(L(TRIAS(OLADJ,NLMAX))),
     *                  NVT,
     *                  IELOLD,IEL2,
     *                  IRES,IEDG,X1,Y1,X2,Y2)
            
            IF (IRES.EQ.0) THEN

C             We known the new element.
C             Update the new coordinates.

              DWORK(KTRC+OTRXPOS) = XNEW
              DWORK(KTRC+OTRYPOS) = YNEW
              DWORK(KTRC+OTRELEM) = DBLE(IEL2)
              
            ELSE IF (IRES.NE.2) THEN

C             Raytracing failed, remove the tracer from the domain.

              DWORK(KTRC+OTRXPOS) = XNEW
              DWORK(KTRC+OTRYPOS) = YNEW
              DWORK(KTRC+OTRXVEL) = 0D0
              DWORK(KTRC+OTRYVEL) = 0D0
              DWORK(KTRC+OTRELEM) = DBLE(-1)

            ELSE
            
C             The domain was left on element IEL through edge IEDG.
C             X1,Y1,X2,Y2 gives the coordinates of the corners of
C             that edge. 
C             Is that edge a Neumann edge?
              
              IF (IAND(KWORK(KXNPR+2*(NVT+IEDG-1)),2**12).EQ.0) THEN
              
C               Neumann edge - the point leaves the domain

                DWORK(KTRC+OTRXPOS) = XNEW
                DWORK(KTRC+OTRYPOS) = YNEW
                DWORK(KTRC+OTRXVEL) = 0D0
                DWORK(KTRC+OTRYVEL) = 0D0
                DWORK(KTRC+OTRELEM) = DBLE(-1)
              
              ELSE
              
C               Bad luck, Dirichlet edge. 
C               Calculate the intersection point on the boudary
C               into (X0,Y0).

                CALL LINSCT (X1,Y1,X2,Y2,XNEW,YNEW,X,Y,X0,Y0,IRES)
                
C               Put the tracer inbetween (XNEW,YNEW) and (X0,Y0)

                XNEW = 0.75D0*X0 + 0.25D0*XNEW
                YNEW = 0.75D0*Y0 + 0.25D0*YNEW
                
C               Hopefully we find the corresponding element of the
C               new position. We start at IEL.

                CALL PSRCH4(XNEW,YNEW,
     *                      DWORK(L(TRIAS(OLCORVG,NLMAX))),
     *                      KWORK(L(TRIAS(OLVERT,NLMAX))),
     *                      KWORK(L(TRIAS(OLMID,NLMAX))),
     *                      KWORK(L(TRIAS(OLADJ,NLMAX))),
     *                      NVT,
     *                      IEL,IEL2,
     *                      IRES,IEDG,X1,Y1,X2,Y2)
     
                IF (IRES.NE.0) THEN
                
C                 Hmmm, perhaps the hierarchical search works...

                  CALL PSRCH5(TRIAS,NLMIN,NLMAX,XNEW,YNEW,IEL,IEL)
                
                END IF
                
                IF (IEL.NE.0) THEN
                
C                 Uff, we made it.
                
                  DWORK(KTRC+OTRXPOS) = XNEW
                  DWORK(KTRC+OTRYPOS) = YNEW
                  DWORK(KTRC+OTRELEM) = IEL
                  
                ELSE
                
C                 No chance. Something baaaaad happened here.
C                 Remove the tracer.
                
                  DWORK(KTRC+OTRXPOS) = XNEW
                  DWORK(KTRC+OTRYPOS) = YNEW
                  DWORK(KTRC+OTRXVEL) = 0D0
                  DWORK(KTRC+OTRYVEL) = 0D0
                  DWORK(KTRC+OTRELEM) = DBLE(-1)

                END IF
              
              END IF
            
            END IF 
          
          END IF
        
        END IF
      
      END DO
      
      END

************************************************************************
* Recycle tracers
*
* This routine searches in the tracer structure for "unused" tracers
* outside of the domain. These tracers are reintroduced in boundary
* elements not containing tracers or containing less tracers than
* defined in the initialisation phase.
*
* The routine only works properly if the tracer structure was
* initialised by INRTRC. In any other case, nothing will happen.
************************************************************************
      
      SUBROUTINE TRCREC(LTRCS,NEL,DCORVG,KVERT,NVBD,KEBD)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'stracers.inc'
      INCLUDE 'cbasictria.inc'
      
C     parameters

      INTEGER LTRCS,NVBD,KEBD(*),NEL,KVERT(NNVE,*)
      DOUBLE PRECISION DCORVG(2,*)
      
C     local variables

      INTEGER LTRCNT,KTRC,IEL,NEBD,LBDCNT,KBDCNT,IFBDE,ITRMAX,ITRIDX
      INTEGER PEL,CNT,I,J,K,KTRCNT
      DOUBLE PRECISION DLEN,DOFS,X,Y,XACT,YACT

      DOUBLE PRECISION DCOORD (2,4),DJF(2,2),DJAC(2,2),DETJ
      
      INTEGER IR
      DOUBLE PRECISION USRAN
      EXTERNAL USRAN
      
C     Check if we have regularly distribuited tracers. If not: cancel.

      PEL = FLOOR(DWORK(L(LTRCS)+OTRPEL-1))
      IF (PEL.EQ.0) RETURN
      
C     We count how many tracers are present in each element.
C     Allocate a working array:

      CALL ZNEW(NEL,3,LTRCNT,'KTRCNT')
      KTRCNT=L(LTRCNT)
      
      DO I=1,FLOOR(DWORK(L(LTRCS)+OTRCNT-1))
      
C       Where is the tracer structure?

        KTRC = L(LTRCS) + OTRDEF-1 + I*SZTRCR - 1
      
C       Only count those tracers in the domain where we know where 
C       they are:

        IEL = FLOOR(DWORK(KTRC+OTRELEM))
        
        IF (IEL.GE.1) THEN
          KWORK(KTRCNT+IEL-1) = KWORK(KTRCNT+IEL-1) + 1
        END IF
      
      END DO
      
C     In the next step we extract the elements adjacent to the 
C     boundary. Allocate a 2nd working array, 2-dimensional,
C     which will be build as:
C       KBDCNT(1,.) = Tracers in the element
C       KBDCNT(2,.) = Element number
C
C     Allocate that working array:

      CALL ZNEW (2*NVBD,3,LBDCNT,'KBDCNT')
      KBDCNT = L(LBDCNT)

C     Loop through the elements on the boundary, don't handle an
C     element twice. (At least not an element which has more than
C     one edge on the same boundary component; elements that
C     hit more than one boundary component might be hit twice,
C     but we handle that later).
C     Build KBDCNT.

      NEBD = 0
      DO I=1,NVBD
        IF ((I.EQ.1).OR.(KEBD(I).NE.KEBD(I-1))) THEN
          KWORK(KBDCNT+2*NEBD+1) = KEBD(I)
          KWORK(KBDCNT+2*NEBD)   = KWORK(KTRCNT+KEBD(I)-1)
          NEBD = NEBD+1
        END IF
      END DO
      
C     NEBD is now the number of elements in KBDCNT.
C      
C     Sort that array with Heapsort for the element number.

      CALL HSRTIN (NEBD,2,KWORK(KBDCNT),2)
      
C     Write "-1" into the tracer counter of each element that exists
C     more than once. This happens rather seldom but can happen if
C      - one cell touches two boundary components (in a coarse grid)
C      - the last cell of a boundary component is the first cell
C        (what it usually is because the domain is connected)

      DO I=2,NEBD
        IF (KWORK(KBDCNT+2*(I-1)+1).EQ.
     *      KWORK(KBDCNT+2*(I-2)+1)) THEN
          KWORK(KBDCNT+2*(I-1)) = -1
        END IF
      END DO

C     Sort that array with Heapsort for the number of tracers in
C     each element.

      CALL HSRTIN (NEBD,2,KWORK(KBDCNT),1)
      
C     Search for the first non-duplicate element

      DO IFBDE=1,NEBD
        IF (KWORK(KBDCNT+2*(IFBDE-1)).GE.0) GOTO 10
      END DO
10    CONTINUE

C     Might be that IFBDE is = NEBD+1, but then, nothing will happen
C     in the following... highly unlikely...
C
C     We now run for every tracer outside of the domain through all
C     boundary elements and "feed" the tracer into the element with
C     the least number of tracers in it. ITRMAX describes the
C     current maximum number of tracers in an element that are allowed,
C     and ITRIDX the next tracer index.
      
      ITRIDX = IFBDE
      ITRMAX = KWORK(KBDCNT+2*(ITRIDX-1))
      
      IR = 1000
      
      DO I=1,FLOOR(DWORK(L(LTRCS)+OTRCNT-1))
      
C       Where is the tracer structure?

        KTRC = L(LTRCS) + OTRDEF-1 + I*SZTRCR - 1
        
C       Is this an unassigned tracer which must be reintroduced? 
C       (element=-1)
C       An element number <=-2 indicates a user defined tracer
C       status; these tracers are ignored completely.

        IF (FLOOR(DWORK(KTRC+OTRELEM)).EQ.-1) THEN
          
C         Insert the tracer into the ITRIDX'th element in the
C         list of boundary elements.

          IEL = KWORK(KBDCNT+2*(ITRIDX-1)+1)
          DWORK(KTRC+OTRELEM) = DBLE(IEL)
          
C         How many tracers are in that element?
C         I there are less than original, equally distribute the
C         tracers there. If there are more, take a random position
C         in that element.
          
          DLEN = 2D0/DBLE(PEL)
          DOFS = -1D0+0.5D0*DLEN
          CNT = KWORK(KBDCNT+2*(ITRIDX-1))
          
          IF (CNT.LT.PEL*PEL) THEN
          
C           Calculate the point's coordinates on the reference element

            X = DBLE(CNT/PEL)*DLEN+DOFS
            Y = DBLE(MOD(CNT,PEL))*DLEN+DOFS

          ELSE
          
C           Put it anywhere into the reference element

            X = USRAN(IR)
            Y = USRAN(IR)
          
          END IF
          
C         Use the transformation from the reference element
C         to calculate the real position
          
          DO K=1,4
            J = KVERT(K,IEL)
            DCOORD(1,K) = DCORVG(1,J)
            DCOORD(2,K) = DCORVG(2,J)
          END DO
          CALL QINIJF (DCOORD,DJF)
          CALL QTRAF (DCOORD,DJF,DJAC,DETJ,
     *                X,Y,XACT,YACT)
          
C         Store the position, increase the number of tracers in that
C         element
          
          KWORK(KBDCNT+2*(ITRIDX-1)) = KWORK(KBDCNT+2*(ITRIDX-1))+1

          DWORK(KTRC+OTRXPOS) = XACT
          DWORK(KTRC+OTRYPOS) = YACT

C         Increment ITRIDX to feed the next tracer into the
C         following element in the list -- as long as this
C         element does not have more tracers than ITRMAX!

          IF (ITRIDX.LT.NEBD) THEN
          
            IF (KWORK(KBDCNT+2*(ITRIDX+1-1)).LE.ITRMAX) THEN
              ITRIDX = ITRIDX + 1
            ELSE
            
C             If the next element contains more tracers than the
C             current one, start again with the first element.
C             This ensures a "homogenuous" distribution of the
C             tracers in the boundary elements.

              ITRIDX = IFBDE
              ITRMAX = KWORK(KBDCNT+2*(ITRIDX-1))
            
            END IF
          
          ELSE
          
C           All elements full (unlikely to happen). Start at 
C           the  beginning.
          
            ITRIDX = IFBDE
            ITRMAX = KWORK(KBDCNT+2*(ITRIDX-1))
            
          END IF
          
        END IF
      
      END DO
      
C     Release auxiliary memory

      CALL ZDISP (0,LBDCNT,'KBDCNT')
      CALL ZDISP (0,LTRCNT,'KTRCNT')
      
      END 

      
************************************************************************
* Calculate tracer velocity
*
* This routine calculates the velocity of all tracers inside of
* the domain from a given finite element solution.
*
* In:
*   LTRCS  : Handle to TDTracers structure on the heap
*   TRIAS  : Triangulation structure on all levels
*   NLMIN  : Minimum level in TRIAS
*   NLMAX  : Maximum level in TRIAS
*   NEQ    : Length of the X- and Y-solution vectors
*   DU     : array [1..NEQ] of double
*            X-velocity field
*   DV     : array [1..NEQ] of double
*            Y-velocity field
*   ELE    : Element that was used for the discretisation of DU and DV
*   BNONPR : =true, if ELE is a nonparametric element, otherwise =false
* 
* Out:
*   The velocity vectors of all tracers in the domain are updated
*   according to DU and DV.
************************************************************************

      SUBROUTINE TRCVUP (LTRCS,TRIAS,NLMIN,NLMAX,NEQ,DU,DV,ELE,BNONPR)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'stria.inc'
      INCLUDE 'stracers.inc'
      
      INTEGER LTRCS,TRIAS(SZTRIA,*),NEQ,NLMIN,NLMAX
      DOUBLE PRECISION DU(NEQ),DV(NEQ)
      EXTERNAL ELE
      LOGICAL BNONPR
      
C     local variables

      INTEGER I,IRES,IEL,IELOLD,KTRC,NVT,IEDG
      DOUBLE PRECISION X,Y,XV,YV,DUMMY1,DUMMY2,X1,Y1,X2,Y2

      NVT = TRIAS(ONVT,NLMAX)

C     Loop through all the tracers

      DO I=0,FLOOR(DWORK(L(LTRCS)+OTRCNT-1))-1
      
C       Where is the tracer structure?

        KTRC = L(LTRCS) + OTRDEF-1 + I*SZTRCR - 1
      
C       What's the element of that tracer? Where is it?

        IELOLD = FLOOR(DWORK(KTRC+OTRELEM))
        X   = DWORK(KTRC+OTRXPOS)
        Y   = DWORK(KTRC+OTRYPOS)
      
C       Ignore tracers with element number <=-1
C       (-> tracers outside of the domain or tracers with
C       user defined tracer status)
      
        IF (IELOLD.GE.0) THEN
        
C         The tracer is in the domain. Make sure we have the correct
C         element that contains the tracer:

          IF (IELOLD.EQ.0) THEN
          
C           Element completely unknown - hierarchical search

            CALL PSRCH5(TRIAS,NLMIN,NLMAX,X,Y,IELOLD,IEL)
            
          ELSE
          
C           Raytracing Search, beginning at element IELOLD

            CALL PSRCH4(X,Y,
     *                  DWORK(L(TRIAS(OLCORVG,NLMAX))),
     *                  KWORK(L(TRIAS(OLVERT,NLMAX))),
     *                  KWORK(L(TRIAS(OLMID,NLMAX))),
     *                  KWORK(L(TRIAS(OLADJ,NLMAX))),
     *                  NVT,
     *                  IELOLD,IEL,
     *                  IRES,IEDG,X1,Y1,X2,Y2)
     
            IF (IRES.NE.0) IEL = 0
          
          END IF
          
C         If the element was not found, the point is not in the
C         domain. Ignore it. Otherwise...

          IF (IEL.NE.0) THEN
          
C           ...evaluate the FE function to get the X- and Y-velocity

            CALL SCEVLQ (NEQ,DU,ELE,BNONPR,X,Y,IEL,
     *                   TRIAS(1,NLMAX),DWORK(L(TRIAS(OLCORVG,NLMAX))),
     *                   KWORK(L(TRIAS(OLVERT,NLMAX))),
     *                   KWORK(L(TRIAS(OLMID,NLMAX))),
     *                   XV, DUMMY1, DUMMY2)

            CALL SCEVLQ (NEQ,DV,ELE,BNONPR,X,Y,IEL,
     *                   TRIAS(1,NLMAX),DWORK(L(TRIAS(OLCORVG,NLMAX))),
     *                   KWORK(L(TRIAS(OLVERT,NLMAX))),
     *                   KWORK(L(TRIAS(OLMID,NLMAX))),
     *                   YV, DUMMY1, DUMMY2)
     
C           and save the velocity.
             
            DWORK(KTRC+OTRXVEL) = XV
            DWORK(KTRC+OTRYVEL) = YV
          
          END IF
          
        END IF ! IEL <> -1
      
      END DO
      
      END
      
************************************************************************
* Write simple tracers
*
* This writes tracers to the GMV file.
*
* Remark: This realises only limited support for tracers, as no
* function values are assigned to the tracers. The only "data" that is
* assigned to each tracer is the TAG element in the tracer structure.
*
* In:
*   MFILE : The file handle of the GMV file; must be open.
*   LTRCS : Handle to TDTracers structure on the heap
*
* Remark: The routine will write all tracers regardless of their
* status, position or tag.
************************************************************************
      
      SUBROUTINE GMVTRS (MFILE,LTRCS)

      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'stracers.inc'
      
      INTEGER MFILE,LTRCS
      
C     local variables

      DOUBLE PRECISION X,Y,XV,YV,TAG
      INTEGER I,NTRC,KTRC
      
C     Introduce tracers in the GMV file

      NTRC = FLOOR(DWORK(L(LTRCS)+OTRCNT-1))
      
      WRITE (MFILE,'(A,I10)') 'tracers ',NTRC
      
C     Loop through all the tracers, write out X-position

      DO I=0,NTRC-1
      
C       Where is the tracer structure?

        KTRC = L(LTRCS) + OTRDEF-1 + I*SZTRCR - 1
      
C       Get information about the tracer
      
        X   = DWORK(KTRC+OTRXPOS)
        WRITE (MFILE,'(F20.10)') X
      
      END DO

C     Write Y-position

      DO I=0,NTRC-1
      
C       Where is the tracer structure?

        KTRC = L(LTRCS) + OTRDEF-1 + I*SZTRCR - 1
      
C       Get information about the tracer
      
        Y   = DWORK(KTRC+OTRYPOS)
        WRITE (MFILE,'(F20.10)') Y
      
      END DO
      
C     Write a dummy as Z-position

      DO I=0,NTRC-1
        WRITE (MFILE,'(F20.10)') 0D0
      END DO
      
C     Write X-velocity

      WRITE (MFILE,'(A)') 'X_velocity'
      DO I=0,NTRC-1
      
C       Where is the tracer structure?

        KTRC = L(LTRCS) + OTRDEF-1 + I*SZTRCR - 1
      
C       Get information about the tracer
      
        XV  = DWORK(KTRC+OTRXVEL)
        WRITE (MFILE,'(F20.10)') XV
      
      END DO
      
C     Write Y-velocity

      WRITE (MFILE,'(A)') 'Y_velocity'
      
      DO I=0,NTRC-1
      
C       Where is the tracer structure?

        KTRC = L(LTRCS) + OTRDEF-1 + I*SZTRCR - 1
      
C       Get information about the tracer
      
        YV  = DWORK(KTRC+OTRYVEL)
        WRITE (MFILE,'(F20.10)') YV
      
      END DO

C     Write user defined tag

      WRITE (MFILE,'(A)') 'Tag'

      DO I=0,NTRC-1
      
C       Where is the tracer structure?

        KTRC = L(LTRCS) + OTRDEF-1 + I*SZTRCR - 1
      
C       Get information about the tracer
      
        TAG = DWORK(KTRC+OTRTAG)
        WRITE (MFILE,'(F20.10)') TAG
      
      END DO
      
      WRITE (MFILE,'(A)') 'endtrace'
      
      END
