***********************************************************************
* This file provides drag-/lift calculation routines like the file
* draglift.f. In contrast to draglift.f the routines here
* build up upon the new calculation method with the general
* volume integration method VOLINT in bdvol.f.
* Also the line integration method is improved to work independent
* of any COMMON block.
***********************************************************************

************************************************************************
* Calculate Drag and Lift on the real boundary.
*
* Extended calling convention
*
* This routine performs a line integration on parts of the real boundaty
* to calculate the body forces that arise there. 
*
* The body forces are defined as the integrals
*
*     dfw=2 int_s [dpf1 dut/dn n_y - p n_x] ds / dpf2
*     daw=2 int_s [dpf1 dut/dn n_x + p n_y] ds / dpf2
*
* with dpf1 and dpf2 being user-defined coefficients according to the
* geometry. Usually in benchmark-like geometries there is:
*   dpf1 = RHO*NU = (density of the fluid)*viscosity
*   dpf2 = RHO*DIST*UMEAN**2
*        = (density of the fluid)*(length of the obstacle facing the flow)
*          *(mean velocity of the fluid)^2
*  
* In:
*   DU1,
*   DU2,
*   DP     - velocity- and pressure vectors of the solution,
*            corresponding to ELE
*   KVERT,
*   KMID,
*   KVBD,
*   KEBD,
*   KBCT,
*   DVBDP,
*   DCORVG - usual geometry information
*   ELE    - Element that is used for the discretization of the velocity
*            vectors; corresponds to DU1/DU2.
*   BNONPR - Decides on parametric/nonparametric element.
*            =false, if ELE is a parametric element
*            =true, if ELE is a nonparametric element
*   TRIA   - array [1..SZTRIA] of integer
*            Triangulation structure of the underlying mesh
*
*   IBDC   - Number of the boundary component where the force should
*            be calculated
*   DMIN,
*   DMAX   - minimum and maximum parameter value on boundary component
*            IBDC between where the forces should be calculated.
*
*   DPF1,
*   DPF2   - coefficients in the integral for the computation
*   ICUB   - Number of the line cubature formula to use
* 
* Out:
*   DFW    - Horizontal force
*   DAW    - Vertical force
*   DLEN   - Length of the boundary component between parameter
*            values DMIN and DMAX
************************************************************************
      
      SUBROUTINE BDFORX(DU1,DU2,DP,
     *                  TRIA,KVERT,KMID,KVBD,KEBD,KBCT,DCORVG,
     *                  DVBDP,ELE,BNONPR,
     *                  IBDC,DMIN,DMAX,DPF1,DPF2,ICUB,
     *                  DFW,DAW,DLEN)
     
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'ccub.inc'
      
      INCLUDE 'cbasictria.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      
      INCLUDE 'stria.inc'

C parameters

      DOUBLE PRECISION DU1(*),DU2(*),DP(*),DCORVG(2,*),DFW,DAW,DLEN
      INTEGER KVERT(NNVE,*),KMID(NNVE,*),KVBD(*),KEBD(*)
      INTEGER TRIA(SZTRIA),IBDC,ICUB,KBCT(*)
      DOUBLE PRECISION DMIN,DMAX,DPF1,DPF2,DVBDP(*)
      EXTERNAL ELE
      LOGICAL BNONPR

C externals

      INTEGER NDFL
      EXTERNAL NDFL
      DOUBLE PRECISION TMAX
      EXTERNAL TMAX
      
C local variables

      INTEGER IELTYP,IVT1,IVT2,II,IEBD
      INTEGER IVE,JP,JDFL
      INTEGER KDFG(NNBAS),KDFL(NNBAS),IDFL
      DOUBLE PRECISION PX1,PX2,PY1,PY2,DLH,XI1,XI2
      DOUBLE PRECISION DTX,DTY,DNX,DNY,DPCONT,DJ1,DJ2,DJ3,DJ4,XX,YY,DUT
      DOUBLE PRECISION DXILIN(NNCUBP),DPAR1,DPAR2,OM

C     Initialize the output variables to 0:

      DFW  = 0D0
      DAW  = 0D0
      DLEN = 0D0

C     Return immediately if the number of the boundary component
C     is out of bounds

      IF ((IBDC.LT.1).OR.(IBDC.GT.TRIA(ONBCT))) RETURN

C     Don't calculate anything if one of the coefficients is 0:

      IF ((DPF1.EQ.0D0).OR.(DPF2.EQ.0D0)) RETURN

C     In the cubature we need function values and derivatives:

      BDER(1) = .TRUE.
      BDER(2) = .TRUE.
      BDER(3) = .TRUE.

C     Get the element type

      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      
C     Number of DOF's per element?
      
      IDFL = NDFL(IELTYP)

C     Initialize the cubature formula. This sets the weights as well 
C     as NCUBP according to the cubature rule on the unit interval
C     [-1,1].

      CALL CB1 (ICUB)
      
C     Transfer the coordinates of the cubature points on [-1,1] to
C     DXILIN. We'll have to transfer them later on to the cubature
C     points on the line of the real element:

      DO II=1,NCUBP
        DXILIN(II) = DXI(II,1)
      END DO
      
C     Loop over all edges belonging to our current boundary component.
C     IEBD counts the index of the edges inside of the KEBD array
C     and the corresponding adjacend vertices in KVBD:

      DO IEBD = KBCT(IBDC),KBCT(IBDC+1)-1
      
C       Get the element adjacent to the current boundary edge:
      
        IEL = KEBD(IEBD)
        
C       Get the two vertices adjacent to the edge. The first vertex
C       can be obtained by KVBD(IEBD). The second vertex is either
C       KVBD(IEBD+1) or (if we are at the last edge) the first
C       vertex of the whole boundary component:
        
        IVT1 = KVBD(IEBD)
        
        IF (IEBD.NE.KBCT(IBDC+1)-1) THEN
          IVT2 = KVBD(IEBD+1)
        ELSE
          IVT2 = KVBD(KBCT(IBDC))
        END IF
        
C       Get the parameter values of the two vertices. If one of
C       the vertices is in the desired range [DMIN,DMAX], do the
C       calculation. If both are outside, we mustn't calculate the
C       integral contribution on that edge!

        DPAR1 = DVBDP(IEBD)
        IF (IEBD.NE.KBCT(IBDC+1)-1) THEN
          DPAR2 = DVBDP(IEBD+1)
        ELSE
          DPAR2 = TMAX(IBDC)
        END IF
        
        IF ( ((DMIN.LE.DPAR1).AND.(DPAR1.LE.DMAX)) .OR.
     *       ((DMIN.LE.DPAR2).AND.(DPAR2.LE.DMAX)) ) THEN
        
C         Now quickly determine which edge 1..4 of the current element
C         is the one we are treating at the moment; the edge starts with
C         vertex IVT1!

          DO II=1,4
            IF (KVERT(II,IEL).EQ.IVT1) GOTO 102
          END DO
          WRITE(MTERM,*) 
     *         'BDFORX: 1st vertex of edge not found. KVERT destroyed?'
          RETURN
102       CONTINUE

C         The number of the edge is now in II. We have to transfer
C         the coordinates of the cubature points from 1D to 2D depending
C         on this edge.

          IF (II.EQ.1) THEN
C           Edge 1 is on the bottom of the reference element
            DO II = 1,NCUBP
              DXI(II,1) = DXILIN(II)
              DXI(II,2) = -1D0
            END DO
          ELSE IF (II.EQ.2) THEN
C           Edge 2 is on the right of the reference element
            DO II = 1,NCUBP
              DXI(II,1) = 1D0
              DXI(II,2) = DXILIN(II)
            END DO
          ELSE IF (II.EQ.3) THEN
C           Edge 3 is on the top of the reference element
            DO II = 1,NCUBP
              DXI(II,1) = -DXILIN(II)
              DXI(II,2) = 1D0
            END DO
          ELSE IF (II.EQ.4) THEN
C           Edge 4 is on the left of the reference element
            DO II = 1,NCUBP
              DXI(II,1) = -1D0
              DXI(II,2) = -DXILIN(II)
            END DO
          END IF
          
C         Fetch the start- and end-coordinates of the edge as well as
C         its length:

          PX1 = DCORVG(1,IVT1)
          PX2 = DCORVG(1,IVT2)
          PY1 = DCORVG(2,IVT1)
          PY2 = DCORVG(2,IVT2)

          DLH = SQRT((PX2-PX1)**2+(PY2-PY1)**2)

C         Sum up the length of all edges to DLEN:

          DLEN = DLEN+DLH
          
C         Calculate the tangential as well as the normal vector
C         of the edge:

          DTX = (PX2-PX1)/DLH
          DTY = (PY2-PY1)/DLH
          DNX = -DTY
          DNY =  DTX

C         Fetch the pressure of the current element:

          DPCONT = DP(IEL)

C         Initialize the local and global DOF's on our current element
C         into KDFG and KDFL:

          CALL NDFGLX(TRIA,IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)

C         Initialize the variables in the element structure
C         according to our current element

          DO IVE = 1,4
            JP=KVERT(IVE,IEL)
            KVE(IVE)=JP
            DX(IVE)=DCORVG(1,JP)
            DY(IVE)=DCORVG(2,JP)
          END DO

C         Calculate auxiliary Jacobian factors

          DJ1 = 0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
          DJ2 = 0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
          DJ3 = 0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
          DJ4 = 0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))

C         Call the element subroutine to adapt to our current
C         configuration of cubature points:

          CALL ELE(0D0,0D0,-2)
      
C         Loop through the cubature points on the current line on
C         the current element to calculate the contributions to the
C         integral:

          DO ICUBP = 1,ICUB
          
C           Get the coordinates of the cubature point on the reference
C           element
          
            XI1=DXI(ICUBP,1)
            XI2=DXI(ICUBP,2)

C           Calculate the Jacobian of the bilinear mapping onto the
C           reference element and its determinant.
C           This is not used here, but maybe inside of the element
C           itself for the calculation of derivatives!

            DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
            DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
            DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
            DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
            
            DETJ = DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
            
C           Calculate the OMEGA for the integration by multiplication
C           of the integration coefficient by the Jacobian of the
C           mapping.
C           The determinant of the mapping of the unit interval [-1,1]
C           to the real line is 0.5*length of the line!
            
            OM   = DOMEGA(ICUBP)*0.5D0*DLH

C           Call the element to calculate the values in the current
C           cubature point:

            IF (BNONPR) THEN
            
C             Calculate the real coordinates of the cubature point
C             with the bilinear transformation

              XX=0.5D0*(DX(1)+DX(2)+DJ1)+0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *          +0.5D0*DJ1*XI2+0.5D0*DJ2*XI1*XI2
              YY=0.5D0*(DY(1)+DY(3)+DJ3)+0.5D0*DJ4*XI1
     *          +0.5D0*(DY(3)-DY(1)-DJ4)*XI2-0.5D0*DJ3*XI1*XI2
     
              CALL ELE(XX,YY,-3)
              
            ELSE
            
              CALL ELE(XI1,XI2,-3)
              
            END IF ! BNONPR

C           Loop through the DOF's on our element and calculate
C           "dut" on our current element:

            DUT=0D0
            DO JDFL=1,IDFL
              DUT = DUT + DU1(KDFG(JDFL))*DBAS(KDFL(JDFL),2)*DTX*DNX
     *                  + DU2(KDFG(JDFL))*DBAS(KDFL(JDFL),2)*DTY*DNX
     *                  + DU1(KDFG(JDFL))*DBAS(KDFL(JDFL),3)*DTX*DNY
     *                  + DU2(KDFG(JDFL))*DBAS(KDFL(JDFL),3)*DTY*DNY
            END DO

C           Sum this up to the two integral contributions in DFW and DAW 

            DFW = DFW + OM*(DPF1*DUT*DNY-DPCONT*DNX)
            DAW = DAW - OM*(DPF1*DUT*DNX+DPCONT*DNY)

          END DO ! ICUBP
          
        END IF ! DMIN<DPAR1<DMAX or DMIN<DPAR2<DMAX

      END DO ! IEBD

C     Final modification to DFW and DAW gives the drag and lift.

      DFW = 2D0*DFW/DPF2      
      DAW = 2D0*DAW/DPF2

99999 END

************************************************************************
* Calculate Integral of function on the real bondary.
*
* Extended calling convention
*
* This routine performs a line integration on parts of the real boundaty
* to calculate the Iintegral of a function that arise there. 
*
* In:
*   DX     - Solution vectors of the solution,
*            corresponding to ELE
*   KVERT,
*   KMID,
*   KVBD,
*   KEBD,
*   KBCT,
*   DVBDP,
*   DCORVG - usual geometry information
*   ELE    - Element that is used for the discretization of the velocity
*            vectors; corresponds to DP
*   BNONPR - Decides on parametric/nonparametric element.
*            =false, if ELE is a parametric element
*            =true, if ELE is a nonparametric element
*   TRIA   - array [1..SZTRIA] of integer
*            Triangulation structure of the underlying mesh
*
*   IBDC   - Number of the boundary component where the force should
*            be calculated
*   DMIN,
*   DMAX   - minimum and maximum parameter value on boundary component
*            IBDC between where the forces should be calculated.
*
*   ICUB   - Number of the line cubature formula to use
* 
* Out:
*   DVAL   - Integral of the FEM-function
*   DLEN   - Length of the boundary component between parameter
*            values DMIN and DMAX
************************************************************************
      
      SUBROUTINE BDINTX(DF,TRIA,KVERT,KMID,KVBD,KEBD,KBCT,DCORVG,
     *                  DVBDP,ELE,BNONPR,
     *                  IBDC,DMIN,DMAX,ICUB,
     *                  DVAL,DLEN)
     
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'ccub.inc'
      
      INCLUDE 'cbasictria.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      
      INCLUDE 'stria.inc'

C parameters

      DOUBLE PRECISION DF(*),DCORVG(2,*),DLEN
      INTEGER KVERT(NNVE,*),KMID(NNVE,*),KVBD(*),KEBD(*)
      INTEGER TRIA(SZTRIA),IBDC,ICUB,KBCT(*)
      DOUBLE PRECISION DMIN,DMAX,DVBDP(*),DVAL
      EXTERNAL ELE
      LOGICAL BNONPR

C externals

      INTEGER NDFL
      EXTERNAL NDFL
      DOUBLE PRECISION TMAX
      EXTERNAL TMAX
      
C local variables

      INTEGER IELTYP,IVT1,IVT2,II,IEBD
      INTEGER IVE,JP,JDFL
      INTEGER KDFG(NNBAS),KDFL(NNBAS),IDFL
      DOUBLE PRECISION PX1,PX2,PY1,PY2,DLH,XI1,XI2
      DOUBLE PRECISION DJ1,DJ2,DJ3,DJ4,XX,YY,DPT
      DOUBLE PRECISION DXILIN(NNCUBP),DPAR1,DPAR2,OM

C     Initialize the output variables to 0:

      DVAL   = 0D0
      DLEN   = 0D0
      
C     Return immediately if the number of the boundary component
C     is out of bounds

      IF ((IBDC.LT.1).OR.(IBDC.GT.TRIA(ONBCT))) RETURN

C     In the cubature we need only function values:

      BDER(1) = .TRUE.
      BDER(2) = .FALSE.
      BDER(3) = .FALSE.

C     Get the element type

      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      
C     Number of DOF's per element?
      
      IDFL = NDFL(IELTYP)

C     Initialize the cubature formula. This sets the weights as well 
C     as NCUBP according to the cubature rule on the unit interval
C     [-1,1].

      CALL CB1 (ICUB)
      
C     Transfer the coordinates of the cubature points on [-1,1] to
C     DXILIN. We'll have to transfer them later on to the cubature
C     points on the line of the real element:

      DO II=1,NCUBP
        DXILIN(II) = DXI(II,1)
      END DO
      
C     Loop over all edges belonging to our current boundary component.
C     IEBD counts the index of the edges inside of the KEBD array
C     and the corresponding adjacend vertices in KVBD:

      DO IEBD = KBCT(IBDC),KBCT(IBDC+1)-1
      
C       Get the element adjacent to the current boundary edge:
      
        IEL = KEBD(IEBD)
        
C       Get the two vertices adjacent to the edge. The first vertex
C       can be obtained by KVBD(IEBD). The second vertex is either
C       KVBD(IEBD+1) or (if we are at the last edge) the first
C       vertex of the whole boundary component:
        
        IVT1 = KVBD(IEBD)
        
        IF (IEBD.NE.KBCT(IBDC+1)-1) THEN
          IVT2 = KVBD(IEBD+1)
        ELSE
          IVT2 = KVBD(KBCT(IBDC))
        END IF
        
C       Get the parameter values of the two vertices. If both of
C       the vertices are in the desired range [DMIN,DMAX], do the
C       calculation. If ine is both outside, we mustn't calculate the
C       integral contribution on that edge!

        DPAR1 = DVBDP(IEBD)
        IF (IEBD.NE.KBCT(IBDC+1)-1) THEN
          DPAR2 = DVBDP(IEBD+1)
        ELSE
          DPAR2 = TMAX(IBDC)
        END IF
        
        IF ( ((DMIN.LE.DPAR1).AND.(DPAR1.LE.DMAX)) .OR.
     *       ((DMIN.LE.DPAR2).AND.(DPAR2.LE.DMAX)) ) THEN
        
C         Now quickly determine which edge 1..4 of the current element
C         is the one we are treating at the moment; the edge starts with
C         vertex IVT1!

          DO II=1,4
            IF (KVERT(II,IEL).EQ.IVT1) GOTO 102
          END DO
          WRITE(MTERM,*) 
     *         'BDIPRX: 1st vertex of edge not found. KVERT destroyed?'
          RETURN
102       CONTINUE

C         The number of the edge is now in II. We have to transfer
C         the coordinates of the cubature points from 1D to 2D depending
C         on this edge.

          IF (II.EQ.1) THEN
C           Edge 1 is on the bottom of the reference element
            DO II = 1,NCUBP
              DXI(II,1) = DXILIN(II)
              DXI(II,2) = -1D0
            END DO
          ELSE IF (II.EQ.2) THEN
C           Edge 2 is on the right of the reference element
            DO II = 1,NCUBP
              DXI(II,1) = 1D0
              DXI(II,2) = DXILIN(II)
            END DO
          ELSE IF (II.EQ.3) THEN
C           Edge 3 is on the top of the reference element
            DO II = 1,NCUBP
              DXI(II,1) = -DXILIN(II)
              DXI(II,2) = 1D0
            END DO
          ELSE IF (II.EQ.4) THEN
C           Edge 4 is on the left of the reference element
            DO II = 1,NCUBP
              DXI(II,1) = -1D0
              DXI(II,2) = -DXILIN(II)
            END DO
          END IF
          
C         Fetch the start- and end-coordinates of the edge as well as
C         its length:

          PX1 = DCORVG(1,IVT1)
          PX2 = DCORVG(1,IVT2)
          PY1 = DCORVG(2,IVT1)
          PY2 = DCORVG(2,IVT2)

          DLH = SQRT((PX2-PX1)**2+(PY2-PY1)**2)

C         Sum up the length of all edges to DLEN:

          DLEN = DLEN+DLH
          
C         Initialize the local and global DOF's on our current element
C         into KDFG and KDFL:

          CALL NDFGLX(TRIA,IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)

C         Initialize the variables in the element structure
C         according to our current element

          DO IVE = 1,4
            JP=KVERT(IVE,IEL)
            KVE(IVE)=JP
            DX(IVE)=DCORVG(1,JP)
            DY(IVE)=DCORVG(2,JP)
          END DO

C         Calculate auxiliary Jacobian factors

          DJ1 = 0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
          DJ2 = 0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
          DJ3 = 0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
          DJ4 = 0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))

C         Call the element subroutine to adapt to our current
C         configuration of cubature points:

          CALL ELE(0D0,0D0,-2)
      
C         Loop through the cubature points on the current line on
C         the current element to calculate the contributions to the
C         integral:

          DO ICUBP = 1,ICUB
          
C           Get the coordinates of the cubature point on the reference
C           element
          
            XI1=DXI(ICUBP,1)
            XI2=DXI(ICUBP,2)

C           Calculate the Jacobian of the bilinear mapping onto the
C           reference element and its determinant.
C           This is not used here, but maybe inside of the element
C           itself for the calculation of derivatives!

            DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
            DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
            DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
            DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
            
            DETJ = DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
            
C           Calculate the OMEGA for the integration by multiplication
C           of the integration coefficient by the Jacobian of the
C           mapping.
C           The determinant of the mapping of the unit interval [-1,1]
C           to the real line is 0.5*length of the line!
            
            OM   = DOMEGA(ICUBP)*0.5D0*DLH

C           Call the element to calculate the values in the current
C           cubature point:

            IF (BNONPR) THEN
            
C             Calculate the real coordinates of the cubature point
C             with the bilinear transformation

              XX=0.5D0*(DX(1)+DX(2)+DJ1)+0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *          +0.5D0*DJ1*XI2+0.5D0*DJ2*XI1*XI2
              YY=0.5D0*(DY(1)+DY(3)+DJ3)+0.5D0*DJ4*XI1
     *          +0.5D0*(DY(3)-DY(1)-DJ4)*XI2-0.5D0*DJ3*XI1*XI2
     
              CALL ELE(XX,YY,-3)
              
            ELSE
            
              CALL ELE(XI1,XI2,-3)
              
            END IF ! BNONPR

C           Loop through the DOF's on our element and calculate
C           "dut" on our current element:

            DPT=0D0
            DO JDFL=1,IDFL
              DPT = DPT + DF(KDFG(JDFL))*DBAS(KDFL(JDFL),1)
            END DO

C           Sum this up to the integral pressure

            DVAL = DVAL + OM*DPT

          END DO ! ICUBP
          
        END IF ! DMIN<DPAR1<DMAX or DMIN<DPAR2<DMAX

      END DO ! IEBD

99999 END

***********************************************************************
* Initialise ALPHA-vector
*
* This routine prepares the ALPHA-vector DALPHA such that it is 1
* on all entries belonging to nodes in the (fictitious) boundary
* component IFBC and 0 everywhere else. If IFBC=0, ALPHA is prepared
* such that it is 1 on all nodes belonging to any fictitious boundary
* component.
* This routine is designed to work with midpoints, which are
* nodes if EM30/EM31 is used. In this case there must be NU=NMT.
*
* In:
*  DALPHA  - array [1..NU] of double; ALPHA-vector to be filled
*  KXNPR,
*  KMID    - information from the triangulation
*  IFBC    - Number of boundary component to prepare the ALPHA
*            vector for:
*            0 = all fictitious boundary components
*            NBCT+1..NBCT+NFBDYC: prepare for fictitious boundary
*                component IFBC-NBCT (=1..NFBDYC)
*
* Out:
*  DALPHA  - the ALPHA-vector filled with 1/0
***********************************************************************

      SUBROUTINE INALPH (DALPHA,KXNPR,KMID,NEL,NVT,NBCT,NU,IFBC)

      IMPLICIT NONE

      INCLUDE 'cbasictria.inc'
      
C parameters

      INTEGER KXNPR(2,*),KMID(NNVE,*),IFBC,NU,NEL,NVT,NBCT
      DOUBLE PRECISION DALPHA(NU)

C local variables

      INTEGER IEL0,IVE,IM1,IVT
      
      CALL  LCL1(DALPHA, NU) 
      
      DO IEL0=1,NEL
        DO IVE=1,4
          
          IM1=KMID(IVE,IEL0)
          
          IF (IFBC.EQ.0) THEN
          
C           Handle the case of "all fictitious boundaries":

            IF (IAND(KXNPR(1,IM1),2**13).NE.0) THEN
              DALPHA(IM1-NVT) = 1.0D0
            END IF
          
          ELSE

C           Case: A specific fictitious boundary component:
          
            IF ( (IFBC.GT.NBCT).AND.
     *           (IAND(KXNPR(1,IM1),2**13).NE.0).AND.
     *           (IAND(KXNPR(1,IM1),2**8).EQ.0).AND.
     *           (KXNPR(2,IM1).EQ.IFBC)) THEN
     
     
              IVT=KXNPR(1,IM1)
              DALPHA(IM1-NVT) = 1.0D0
              
            ENDIF
            
          ENDIF
       
        ENDDO
        
      ENDDO

      END      

***********************************************************************
* Relaxate ALPHA-vector
*
* This routine is an additional relaxation routine for setting up
* the ALPHA vector. Normally ALPHA is set to 1 for those points
* that belong to a fictitious boundary component, and is set to 0
* for all other points.
* The problem is that this perhaps sets up ALPHA the wrong way
* if points are not inside the boundary component, but very near to
* it. This routine searches for edge midpoints (these are the points
* which correspond to entries in ALPHA) and analyses the distance
* between the edge midpoint and the corners of the element.
* If the (reconstructed) intersection point between the boundary
* component and the edge midpoint is less that TOL percent of
* the length of the line to the corner, ALPHA is set to 1 for
* that edge midpoint as well.
*
*            |-----o-----|
* interface  |           |
* with -> ---:-----      |     
* inters.  ..O.....\     M     <- point M is included as well if
* points   ..|......-----x-- }        |M-x| / |M-C| <= TOL
*          ..|...........|..  
*          ..|-----O-----C..
*                  ^ point is included anyway
*
* In:
*  DALPHA  - array [1..NU] of double; ALPHA-vector to be filled
*  DCORVG, DCORMG, KVERT, KMID, NVT,
*          - information from the triangulation
*  IFBC    - Number of boundary component to prepare the ALPHA
*            vector for:
*            0 = all fictitious boundary components
*            NBCT+1..NBCT+NFBDYC: prepare for fictitious boundary
*                component IFBC-NBCT (=1..NFBDYC)
*  TOL     - 0 <= TOL < 1. Tolerance factor. If TOL > 0, the setup
*            of ALPHA is slightly modified such that ALPHA=1
*            also for points that are outside, but near the
*            reconstructed boundary.
*  IGEOM  - array [1..*] of integer 
*  DGEOM  - array [1..*] of double 
*           Integer- and double-precision parameter blocks with
*           geometry information. Passed to fictitious boundary
*           routines. Not used in this routine.
*
* Out:
*  DALPHA  - the ALPHA-vector filled with 1/0
***********************************************************************

      SUBROUTINE RXALPH (DALPHA,DCORVG,KVERT,DCORMG,KMID,NEL,NVT,NU,
     *                   IFBC,TOL,IGEOM,DGEOM)

      IMPLICIT NONE

      INCLUDE 'cbasictria.inc'
      
C parameters

      INTEGER KVERT(NNVE,*),KMID(NNVE,*),IFBC,NU,NVT,NEL,IGEOM(*)
      DOUBLE PRECISION DCORVG(2,*),DCORMG(2,*),DALPHA(NU),TOL,DGEOM(*)

C externals

      INTEGER ISFBDY
      EXTERNAL ISFBDY

C local variables

      INTEGER IEL,KV1,KV2
      INTEGER IEDGE1, IEDGE2
      DOUBLE PRECISION X1,Y1,X2,Y2,MX1,MY1,MX2,MY2
      DOUBLE PRECISION DIST1,DIST2, DEDGE1,DEDGE2
      LOGICAL BFOUND

C loop through all elements

      DO IEL=1,NEL

C Reconstruct the interface on the current element 
C with 8 approximation steps

        CALL RCFBLI (DCORVG, KVERT, IEL, 8, IFBC, .FALSE.,
     *               X1,Y1,IEDGE1, X2,Y2,IEDGE2, BFOUND,
     *               IGEOM,DGEOM)
     
        IF (BFOUND) THEN
        
C Ok, interface is found and given by 2 points, together with the numbers
C of the edges containing these points. 
C We are interested in the case that 2 edge midpoints are outside
C of the domain and one midpoint is near the domain. As we have the number
C of the edges where the fictitious boundary intersects, we know
C on which edges we have to search for these points...

C Determine if both edge midpoints are outside:
          
          MX1 = DCORMG(1,KMID(IEDGE1,IEL)-NVT)
          MY1 = DCORMG(2,KMID(IEDGE1,IEL)-NVT)
          MX2 = DCORMG(1,KMID(IEDGE2,IEL)-NVT)
          MY2 = DCORMG(2,KMID(IEDGE2,IEL)-NVT)
          
          IF ((ISFBDY(MX1,MY1,IFBC,IGEOM,DGEOM).EQ.0).AND.
     *        (ISFBDY(MX2,MY2,IFBC,IGEOM,DGEOM).EQ.0)) THEN

C Determine both distances - of the midpoint to the reconstructed
C intersection point

            DIST1 = ( (MX1-X1)**2 + (MY1-Y1)**2 )
            DIST2 = ( (MX2-X2)**2 + (MY2-Y2)**2 )
            
C and the length of the edges - to be more exact: the half length of
C the edges, sinde we are interested at the relative distance to the
C *midpoint* of the edge.
            
            KV1 = KVERT(IEDGE1,IEL)
            KV2 = KVERT(MOD(IEDGE1,NNVE)+1,IEL)
            DEDGE1 = DSQRT( (DCORVG(1,KV2)-DCORVG(1,KV1))**2 +
     *                      (DCORVG(2,KV2)-DCORVG(2,KV1))**2 ) * 0.5D0

            KV1 = KVERT(IEDGE2,IEL)
            KV2 = KVERT(MOD(IEDGE2,NNVE)+1,IEL)
            DEDGE2 = DSQRT( (DCORVG(1,KV2)-DCORVG(1,KV1))**2 +
     *                      (DCORVG(2,KV2)-DCORVG(2,KV1))**2 ) * 0.5D0
            
C If the distance is <= TOL, treat the midpoint like it's contained 
C in the boundary component, so set the ALPHA-value to 1.

            IF ( DIST1.LE.(TOL*DEDGE1) ) THEN
              DALPHA (KMID(IEDGE1,IEL)-NVT) = 1D0
            END IF

            IF ( DIST2.LE.(TOL*DEDGE2) ) THEN
              DALPHA (KMID(IEDGE2,IEL)-NVT) = 1D0
            END IF
            
          END IF

        END IF
        
      END DO

      END      

************************************************************************
* Calculates lift (dfwy) and drag (dfwx)
*
* Calculates the X- and Y- body-forces for the fictitious boundary
* component IFBC. If IFBC=0, the sum of the body forces of all
* fictitious boundary components will be calculated.
*
* Volume integration, constant pressure
*
* In:
*  DU1/DU2/DP - Solution vector(s) for velocity and pressure
*  DPF1/DPF2 - are correction factors for the drag/lift coefficient.
*              For standard drag/lift-calculation these are provided by
*              the DPF(1)/DPF(2) parameters that are initialised in 
*              INDAT2D.F by the user initialisation routine.
*  TRIA      - array [1..SZTRIA] of integer
*              Triangulation structure of the underlying mesh
*  DCORVG,
*  KVERT,
*  DCORMG,
*  KMID,
*  KXNPR     - Trangulation arrays; must correspond to TRIA!
*  IFBC      - Number of fictitious boundary component to calculate the
*              body forces for, or =0 to calculate the sum of all body 
*              forces.
*  ELE       - Element function
*  BNONPR    - true, if ELE is a nonparametric element
*  LAMBDA    - Relaxation parameter in [0..1] for relaxation 
*              of ALPHA vector by reconstructed interface, see RXALPH.
*              0D0=standard=no relaxation
*              1D0=treat diagonal intersection lines like horizontal
*
* Out:
*  DFWX      - Drag coefficient
*  DFWY      - Lift coefficient
************************************************************************

      SUBROUTINE BDFVLG (DU1,DU2,DP,NEQU,DCORVG,KVERT,DCORMG,
     *                   KMID,KXNPR,TRIA,
     *                   ELE,BNONPR,DFWX,DFWY,DPF1,DPF2,
     *                   IFBC,LAMBDA,IGEOM,DGEOM)
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      
      INCLUDE 'stria.inc'
      
C parameters

      DOUBLE PRECISION DU1(*),DU2(*),DP(*),DCORVG(2,*),DCORMG(2,*)
      INTEGER KXNPR(2,*),KVERT(NNVE,*),KMID(NNVE,*),IFBC
      INTEGER TRIA(SZTRIA),IGEOM(*)
      INTEGER NEQU
      DOUBLE PRECISION DFWX,DFWY,DPF1,DPF2,LAMBDA,DGEOM(*)
      LOGICAL BNONPR
      
C For calculation of number of relaxed midpoiints:
C     DOUBLE PRECISION AL

C externals
      EXTERNAL ELE
      EXTERNAL CINCP,CVDCP,CISCP,CDOCP,CVDCE
      
C local variables

      INTEGER IINFO(6+SZTRIA)
      
C Handle to DINFO-Array for integration routines.
C
C Structure:
C   double DRAG-coefficient
C   double LIFT-coefficient
C   double temp [3]
C   double DPF1 (from parameter)
C   double DPF2 (from parameter)
C   double ALPHA [NU]

      INTEGER LDINFO
      
C no calculation if one of the DPF-parameters that are provided
C in PTSDAT (file indat2d.f) is 0.

      IF ((DPF1.EQ.0D0).OR.(DPF2.EQ.0D0)) RETURN
      
C Allocate memory for DINFO-arrays. 

      CALL ZNEW (NEQU+7,1,LDINFO,'LDINFO')
      
c Use the vector dalpha as a test function in the weak formulation
c (only the difusion part is taken, should take even the
c convection but in featflow it seems to be ok) ... 
c this needs to be checked in timedependent case

C Prepare the ALPHA-vector for all fict. boundary components:

      CALL INALPH (DWORK(L(LDINFO)+7),KXNPR,KMID,
     *             TRIA(ONEL),TRIA(ONVT),TRIA(ONBCT),NEQU,IFBC)
      
C      AL=0D0
C      DO I=1,NU
C        AL=AL+DWORK(L(LDINFO)+7+I-1)
C      END DO
      
C Relaxate the alpha for nearly-horizontal/vertical lines

      CALL RXALPH (DWORK(L(LDINFO)+7),DCORVG,KVERT,DCORMG,KMID,
     *             TRIA(ONEL),TRIA(ONVT),
     *             NEQU,IFBC,LAMBDA,IGEOM,DGEOM)

C      DO I=1,NU
C        AL=AL-DWORK(L(LDINFO)+7+I-1)
C      END DO
C      PRINT *,'#Relaxated nodes = ',INT(-AL+0.5D0)

C provide the parameters DPF1, DPF2

      DWORK(L(LDINFO)+5) = DPF1
      DWORK(L(LDINFO)+6) = DPF2
      
C IINFO(1..5): reserved for future use.
C In IINFO(6) we provide the desired fictitious boundary number.
C Actually we don't use it (because ALPHA contains information about
C th boundary component), but...

      IINFO(6) = IFBC
      
C In IINFO(7..7+SZTRIA-1] we save the triangulation structure.
C This allowes the callback routines to access the triangulation

      CALL LCP3 (TRIA,IINFO(7),SZTRIA)
      
C Call the calculation routine; 
C IINFO(1..5) saves the number of elements with 0..4 edge 
C midpoints in a fictitious boundary object
      
      CALL VOLINT (DU1,DU2,DP,TRIA,KVERT,KMID,
     *    DCORVG,ELE,BNONPR,8,
     *    IINFO, 6, DWORK(L(LDINFO)), NEQU+7,
     *    CINCP,CVDCE,CVDCP,CISCP,CDOCP,.FALSE.,IGEOM,DGEOM)
      
C Obtain the values of the drag- and lift-coefficient:

      DFWX = DWORK(L(LDINFO))
      DFWY = DWORK(L(LDINFO)+1)
      
C Release memory, finish

      CALL ZDISP (0,LDINFO,'LDINFO')
      
99999 END


************************************************************************
* Volume-integration callback-routines for constant pressure:
************************************************************************

************************************************************************
* Initialisation
************************************************************************
      SUBROUTINE CINCP (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
     *                  TRIA,ELE,BNONPR)
      
      IMPLICIT NONE
      
C parameters
      
      INTEGER NIINFO,IINFO(NIINFO),NDINFO,TRIA(*)
      LOGICAL BNONPR
      DOUBLE PRECISION DINFO(NDINFO),DU1(*),DU2(*),DP(*)
      EXTERNAL ELE
      
C initialise drag/lift-values with 0
      
      DINFO(1) = 0D0
      DINFO(2) = 0D0
      
      IINFO(1) = 0
      IINFO(2) = 0
      IINFO(3) = 0
      IINFO(4) = 0
      IINFO(5) = 0
      
C IINFO(6) holds the number of the fictitious bondary component and
C must not be initialised with 0!
      
      END

************************************************************************
* Postprocessing
************************************************************************
      SUBROUTINE CDOCP (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
     *                  TRIA,ELE,BNONPR)
      
      IMPLICIT NONE
      
C parameters
      
      INTEGER NIINFO,IINFO(NIINFO),NDINFO,TRIA(*)
      DOUBLE PRECISION DINFO(NDINFO),DU1(*),DU2(*),DP(*)
      LOGICAL BNONPR
      EXTERNAL ELE
      
C local variables 

      DOUBLE PRECISION DPF1,DPF2

C obtain drag-/lift correction parameters from the info field:

      DPF1 = DINFO(6)
      DPF2 = DINFO(7)

C Correct the values according to the formulas of the
C drag/lift coefficients with the help of the information
C provided by PTSDAT:
      
      IF (DPF2.NE.0D0) THEN
        DINFO(1)=2D0*DINFO(1)/DPF2
        DINFO(2)=2D0*DINFO(2)/DPF2
      END IF      

      END
      
************************************************************************
* Preparation of calculation on an element
************************************************************************

      SUBROUTINE CVDCE (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
     *                  IELEM,TRIA,ELE,BNONPR,IDFL,KDFG,KDFL,
     *                  IGEOM,DGEOM)
      IMPLICIT NONE
      
!      INCLUDE 'cmem.inc'
!      INCLUDE 'cbasictria.inc'
!      INCLUDE 'cbasicelem.inc'
!      INCLUDE 'celem.inc'
!      
!      INCLUDE 'stria.inc'
!      
!C parameters
!
      INTEGER NIINFO,IINFO(NIINFO),NDINFO,IELEM,TRIA(*),IGEOM(*)
      INTEGER IDFL,KDFG(*),KDFL(*),DGEOM(*)
      DOUBLE PRECISION DINFO(NDINFO),DU1(*),DU2(*),DP(*)
      LOGICAL BNONPR
      EXTERNAL ELE
!
!C local variables
!
!      INTEGER I,J,IM
!      DOUBLE PRECISION X,Y
!      INTEGER KMID,KCORMG,NVT
!      
!      INTEGER ISFBDY
!      EXTERNAL ISFBDY
!
!C Fetch some triangulation information to local variables.
!
!      KMID   = L(TRIA(OLMID))
!      KCORMG = L(TRIA(OLCORMG))
!      NVT    = TRIA(ONVT)
!
!C Calculate how many edge midpoints of the current element
!C are in a boundary object
!
!      J=0
!      
!      DO I=1,NNVE
!C Get element midpoint number
!        
!        IM = KWORK(KMID+(IELEM-1)*NNVE+I-1) - NVT - 1
!        
!        IF (IM.NE.0) THEN
!        
!C is it in in a boundary object?
!
!          X = DWORK(KCORMG+2*IM)
!          Y = DWORK(KCORMG+2*IM+1)
!
!          IF (ISFBDY (X,Y,IINFO(6),IGEOM,DGEOM).NE.0) THEN
!            J=J+1
!          END IF
!          
!        END IF
!      
!      END DO
!      
!C Increment IINFO according to J
!
!      IINFO(J+1) = IINFO(J+1)+1
      
      END
      
************************************************************************
* Calculation of (partial) values in a cubature point
************************************************************************
      SUBROUTINE CVDCP (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
     *                  XX,YY,XX1,YY1,ICUBP,IELEM,TRIA,ELE,BNONPR,
     *                  IDFL,KDFG,KDFL)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      
C parameters
      
      INTEGER NIINFO,IINFO(NIINFO),NDINFO,IELEM,ICUBP,TRIA(*)
      DOUBLE PRECISION DINFO(NDINFO),DU1(*),DU2(*),DP(*)
      DOUBLE PRECISION XX,YY,XX1,YY1
      INTEGER IDFL,KDFG(*),KDFL(*)
      LOGICAL BNONPR
      EXTERNAL ELE

C local variables

      DOUBLE PRECISION DAV,DAX,DAY,ALPH
      INTEGER I,IG

C The ALPHA-vector is stored in DINFO(8..NU+8)

      DAV=0D0 ! VALUE
      DAX=0D0 ! X-der.
      DAY=0D0 ! Y-der.
      DO I=1,IDFL
        IG=KDFG(I)
        ALPH=DINFO(7+IG)
        DAV=DAV+ALPH*DBAS(KDFL(I),1)
        DAX=DAX+ALPH*DBAS(KDFL(I),2)
        DAY=DAY+ALPH*DBAS(KDFL(I),3)
      END DO

C save these values to DINFO for later use

      DINFO(3) = DAV
      DINFO(4) = DAX
      DINFO(5) = DAY

      END
      
************************************************************************
* Summing up the values from the cubature point to the values of
* interest.
************************************************************************
      SUBROUTINE CISCP (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
     *                  DETJ,ICUBP,IELEM,OMEGA,
     *                  DU1V,DU1X,DU1Y,DU2V,DU2X,DU2Y,DPV)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cbasicelem.inc'
      
C parameters
      
      INTEGER NIINFO,IINFO(NIINFO),NDINFO,ICUBP,IELEM
      DOUBLE PRECISION DINFO(NDINFO),DU1(*),DU2(*),DP(*)
      DOUBLE PRECISION DETJ,OMEGA
      DOUBLE PRECISION DU1V,DU1X,DU1Y,DU2V,DU2X,DU2Y,DPV
      
      DOUBLE PRECISION DUT,DT1,DT2
      
C local variables

      DOUBLE PRECISION DAV,DAX,DAY,DN1,DN2,AH1,AH2,DPF1,DPF2

C get value/derivative

      DAV = DINFO(3)
      DAX = DINFO(4)
      DAY = DINFO(5)
      
      DPF1 = DINFO(6)
      DPF2 = DINFO(7)
      
C Form the integrand. Don't normalise the DN1/DN2, since it is not
C really a normal vector but a scaled one, calculated with the
C help of ALPHA!

      DN1=-DAX
      DN2=-DAY
      
C Calculate the force vector. This is defined as:
C
C (FD) = int ( sigma * grad(alpha) ) dx
C (FL)    V
C
C There are now different possibilities about that sigma.
C The original and most correct formulation is the following, which
C is also known as "deformation formulation":
C
C     sigma = -p*I + dpf1*[ grad(u) + grad(u)^T ]
C
C Unfortunately this gives not the best results in the benchmark
C tests. Another formulation is the following, which is also called
C "gradient formulation":
C
C     sigma = -p*I + dpf1*[ grad(u) ]
C
C This can be used in the case that div(u)=0 - in this case
C the strong formulation of the Navier Stokes Equation, which
C can be derived using this sigma, are the same because of
C                 0 = div(u_h) = grad(u_h)^T
C for the discretised matrices!
C 
C Unfortunately using nonconforming finite elements the deformation
C formulation does not give very accurate results, the results with
C the gradient formulation are much better. This is not true anymore
C if other finite elements are used: Using Q2 the deformation
C formulation gives far better results!
C
C We give all implementations here - comment the formulation in
C that you want to use/test!
C
C For volume integration we use:
C
C alpha = c*normal vector 
C -> up to multiplicative constant, so not normalised
C -> must not be normalised because multiplicative constant
C    is necessary for the integral because of Green's formula!
C    (see construction of this)
      
C Implementation of the deformation formulation of the stress tensor:
      
      AH1 = -DPV*DN1 + DPF1*(2D0*DU1X*DN1+(DU2X*DN2+DU1Y*DN2))
      AH2 = -DPV*DN2 + DPF1*((DU2X*DN1+DU1Y*DN1)+2D0*DU2Y*DN2)
      
C Implementation of the gradient formulation of the stress tensor:

C      AH1=-DPV*DN1+DPF1*(DU1X*DN1+DU1Y*DN2)
C      AH2=-DPV*DN2+DPF1*(DU2X*DN1+DU2Y*DN2)
     
C Sum up to the integral

      DINFO(1)=DINFO(1)+AH1*OMEGA
      DINFO(2)=DINFO(2)+AH2*OMEGA
      
      END
      
************************************************************************
* Calculates lift (dfwy), drag (dfwx) 
*
* Calculates the X- and Y- body-forces for the fictitious boundary
* component IFBC. If IFBC=0, the sum of the body forces of all
* fictitious boundary components will be calculated.
*
* Reconstructed interface integration, constant pressure
*
* In:
*  DU1/DU2/DP - Solution vector(s) for velocity and pressure
*  DPF1/DPF2 - are correction factors for the drag/lift coefficient.
*              For standard drag/lift-calculation these are provided by
*              the DPF(1)/DPF(2) parameters that are initialised in 
*              INDAT2D.F by the user initialisation routine.
*  IFBC      - Number of fictitious boundary component to calculate the
*              body forces for, or =0 to calculate the sum of all body 
*              forces.
*  ELE       - Element function
*  BNONPR    - true, if the element is a nonparametric element
*  IGEOM  - array [1..*] of integer 
*  DGEOM  - array [1..*] of double 
*           Integer- and double-precision parameter blocks with
*           geometry information. Passed to boundary
*           routines. Not used in this routine.
*
* Out:
*  DFWX      - Drag coefficient
*  DFWY      - Lift coefficient
************************************************************************

      SUBROUTINE BDFRIG (DU1,DU2,DP,DCORVG,KVERT,DCORMG,KMID,
     *                   TRIA,ELE,BNONPR,DFWX,DFWY,DPF1,DPF2,IFBC,
     *                   IGEOM,DGEOM)


      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasictria.inc'

      INCLUDE 'ccub.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      
      INCLUDE 'stria.inc'
      
C parameters
      DOUBLE PRECISION DU1(*),DU2(*),DP(*),DCORVG(2,*),DCORMG(2,*)
      DOUBLE PRECISION DGEOM(*)
      INTEGER KVERT(NNVE,*),KMID(NNVE,*),IFBC,TRIA(*),IGEOM(*)
      DOUBLE PRECISION DFWX,DFWY,DPF1,DPF2
      LOGICAL BNONPR

C externals

      EXTERNAL ELE
      INTEGER NDFL
      EXTERNAL NDFL
      
C local variables

      INTEGER IELTYP,I,IVE,JP,IG,IFBC1,IEDGE1,IEDGE2
      DOUBLE PRECISION DU1V,DU1X,DU1Y,DU2V,DU2X,DU2Y,DPV
      DOUBLE PRECISION DCOORD(2,4),OM
      DOUBLE PRECISION DCUBP(2,NNCUBP)
      DOUBLE PRECISION XX,YY,XI1,XI2
      DOUBLE PRECISION X1,Y1,X2,Y2
      LOGICAL BDOINT
      INTEGER IDFL,KDFG(NNCUBP),KDFL(NNCUBP)
      
      DOUBLE PRECISION AH1,AH2,DN,DN1,DN2,DT1,DT2
      
C     No calculation if one of the DPF-parameters that are provided
C     in PTSDAT (file indat2d.f) is 0.

      IF ((DPF1.EQ.0D0).OR.(DPF2.EQ.0D0)) RETURN
      
c     preparation - evaluation of parameters
      IER=0

c     which derivatives of basis functions are needed?
      DO  I = 1,NNDER
        BDER(I)=.FALSE.
      ENDDO
      
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.
      
c     Dummy call of ELE, sets number of element

      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)

      IDFL=NDFL(IELTYP)

C Prepare conforming element for cubature. Perform a dummy call 
C to the element in the conforming case.
C This is done for saving arithmetic operations in later calls.
C
C In the nonconforming case the element has to be initialised
C separately on every element, so there's no use in initialising it
C globally.
C
C Normally we have to set ICUBP before calling the element.
C But because our later line cubature formula uses quadrature points
C that are more or less arbitrary placed in the element,
C we only support elements here not depending on the type
C of the cubature formula!
C Therefore we don't set ICUBP here before we call the element.
C Elements that depend on this information will either give
C false values or halt the program!

      ICUBP=0
      
      IF (.NOT.BNONPR) THEN
        CALL ELE(0D0,0D0,-2)
      END IF

C User defined initialisation

      DFWX = 0D0
      DFWY = 0D0

C Loop over all elements.
C
C The variable IEL is stored in the COMMON block /ELEM/ in case that
C any procedure has to know in which element we are...

      DO IEL=1,TRIA(ONEL)
      
C Get the degrees of freedom:

        CALL NDFGLX(TRIA,IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
        IF (IER.LT.0) RETURN

C Store the coordinates of the corners of the current element in 
C the COMMON block variables DX/DY (necessary for the element) 
C as well as in the DCOORD array for the coordinate transformation.

        DO IVE = 1, TRIA(ONVE)
          JP=KVERT(IVE,IEL)
          KVE(IVE)=JP
          DX(IVE)=DCORVG(1,JP)
          DY(IVE)=DCORVG(2,JP)
          DCOORD (1,IVE) = DCORVG(1,JP)
          DCOORD (2,IVE) = DCORVG(2,JP)
        END DO
        
C Is there at all a line on the current element?
C Reconstruct the interface using 8 approximation steps

        BDOINT = .FALSE.
        CALL RCFBLI (DCORVG, KVERT, IEL, 
     *               8, IFBC, .FALSE.,
     *               X1,Y1, IEDGE1, X2,Y2, IEDGE2, BDOINT,
     *               IGEOM,DGEOM)
     
C If the two points are identical, we have the case that exactly the
C corner of one element is hit. Of course there is nothing to integrate
C then! Btw.: There must not be done anything as the length of the
C line is 0, which would produce NAN when dividing by that!

        IF ((X1.EQ.X2).AND.(Y1.EQ.Y2)) BDOINT=.FALSE.

        IF (BDOINT) THEN

C Initialise the line-cubature formula für the current element.
C This will set up the
C quadrature points on the real element as well as on the reference
C element. The cubature points on the reference element are
C stored in DXI/DYI, which can be found in the COMMON block.
C The coordinates on the cubature points on the real element
C will be stored in DCUBP.

C The CB1 routine will initialise the weighting factors, the NCUBP
C information and the distribution of the cubature points on the
C reference interval [-1,1]:

          CALL CB1 (4)
      
C Transform the parametrisation into real coordinates.
C Call the fictitious boundary projection routine to calculate the point.
C
C DETJ will receive the length of the line and must be used
C as a weighting factor in addition to the weighting factors
C of the quadrature points.

          CALL FBCBPR (DCOORD,IEL,IFBC,X1,Y1,X2,Y2,DCUBP,DETJ,
     *                 IGEOM,DGEOM)

C Transform the real coordinates back into coordinates on 
C the reference element - to support conforming
C FE-approaches. This requires the inverse transformation...
C This of course is not necessary in the case of non-parametric elements,
C but we do this in all cases here for sure. Maybe that another
C user-provided routine depend on this.

          DO I=1,NCUBP
            CALL QBTRAF (DCOORD, DXI(I,1), DXI(I,2),
     *                   DCUBP(1,I), DCUBP(2,I))
          END DO

C Prepare nonparametric elements for cubature on current element.
C The COMMON block variable ICUBP was initialised earlier...
C
C We set ICUBP to 0 for the same reason as explained above.
C The later loop will change ICUBP as it is used for looping through
C the cubature points on the element, so we have to reset it here.

          ICUBP=0
          IF (BNONPR) THEN
            CALL ELE(0D0,0D0,-2)
          END IF

c Loop over all cubature points

          DO ICUBP = 1, NCUBP
          
c Cubature point on the reference element

            XI1 = DXI(ICUBP,1)
            XI2 = DXI(ICUBP,2)
          
C Cubature point on the real element:

            XX = DCUBP(1,ICUBP)
            YY = DCUBP(2,ICUBP)
            
C Calculate the weighting factor for the current cubature point
C with the help of the Jacobian determinant = (length of the line)/2

            OM = DOMEGA(ICUBP)*DETJ
          
C Evaluate the basis functions in the cubature point
C for the velocities

            IF(BNONPR) THEN
              CALL ELE(XX,YY,-3)
            ELSE
              CALL ELE(XI1,XI2,-3)
            ENDIF
            IF (IER.LT.0) RETURN
      
C At first build the normal vector from the start-/endpoint of the line.
C Remember that the coordinates of the endpoints of the line are oriented 
C to direct in the tangential direction of the boundary component, see
C definition of RCFBLI!

            DN1=-(Y2-Y1)
            DN2=(X2-X1)

C           normalise it
      
            DN = 1D0/DSQRT(DN1*DN1+DN2*DN2)
            DN1 = DN*DN1
            DN2 = DN*DN2
      
C           Turn it back to get the (normalised) tangential vector:

            DT1 = DN2
            DT2 = -DN1
            
C           Calculate the values of the function in the cubature 
C           points:
C            
C           X-Component:
            
            DU1V=0D0 ! value
            DU1X=0D0 ! x dreiv.
            DU1Y=0D0 ! y deriv
            DO I=1,IDFL
              IG=KDFG(I)
              DU1V=DU1V+DU1(IG)*DBAS(KDFL(I),1)
              DU1X=DU1X+DU1(IG)*DBAS(KDFL(I),2)
              DU1Y=DU1Y+DU1(IG)*DBAS(KDFL(I),3)
            END DO
          
C           Y-Component:
          
            DU2V=0D0 ! value
            DU2X=0D0 ! x dreiv.
            DU2Y=0D0 ! y deriv
            DO I=1,IDFL
              IG=KDFG(I)
              DU2V=DU2V+DU2(IG)*DBAS(KDFL(I),1)
              DU2X=DU2X+DU2(IG)*DBAS(KDFL(I),2)
              DU2Y=DU2Y+DU2(IG)*DBAS(KDFL(I),3)
            END DO
          
C           Pressure:
          
            DPV=DP(IEL)
      
C           Calculate the force vector. This is defined as:
C
C           (FD) = int ( sigma * alpha ) dx
C           (FL)    V
C
C           There are now different possibilities about that sigma.
C           The original and most correct formulation is the following, which
C           is also known as "deformation formulation":
C
C               sigma = -p*I + dpf1*[ grad(u) + grad(u)^T ]
C
C           Unfortunately this gives not the best results in the 
C           benchmark tests. Another formulation is the following, 
C           which is also called "gradient formulation":
C
C               sigma = -p*I + dpf1*[ grad(u) ]
C
C           This can be used in the case that div(u)=0 - in this case
C           the strong formulation of the Navier Stokes Equation, which
C           can be derived using this sigma, are the same because of
C                           0 = div_h(u) = grad_h(u)^T 
C           for the discretised matrices!
C
C           In case of line integrals there's a third possibility how to
C           define that integral. This possibility was that one which  
C           was proposed in the original paper where the DFG-benchmark 
C           was proposed.
C           It should be equivalent with the avove formulation, but to be
C           honest, we haven't found out why...:
C
C               sigma = -p*I + dpf1 * [ Du_T/Dn_S ] * n
C                     = -p*I + dpf1 * <Du*n,t> * t
C
C           with n the normal vector and t the tangential vector to 
C           the surface.
C           
C           Unfortunately using nonconforming Finite Elements, the 
C           deformation formulation does not give very accurate results, 
C           the results with the gradient formulation are much better. 
C           This is not true anymore if other finite elements are used:
C           Using Q2 the deformation formulation gives far better 
C           results!
C
C           The third method was the original method used by FeatFlow 
C           for line integrals around a non-fictitious boundary 
C           components. For fictitious boundary components it should
C           give roughly thesame result as the deformation
C           formulation...

C           Implementation of the deformation formulation of the
C           stress tensor:
      
            AH1 = -DPV*DN1 + DPF1*(2*DU1X*DN1+(DU2X*DN2+DU1Y*DN2))
            AH2 = -DPV*DN2 + DPF1*((DU2X*DN1+DU1Y*DN1)+2*DU2Y*DN2)
  
C           Implementation of the gradient formulation of the
C           stress tensor:

C            AH1=-DPV*DN1 + DPF1*(DU1X*DN1+DU1Y*DN2)
C            AH2=-DPV*DN2 + DPF1*(DU2X*DN1+DU2Y*DN2)
  
C           Implementation of the surface integral with the help of 
C           tangential vectors:
  
C            DUT = (DU1X*DN1+DU1Y*DN2)*DT1 + (DU2X*DN1+DU2Y*DN2)*DT2
C            AH1 = -DPV*DN1 + DPF1*(DUT*DT1)
C            AH2 = -DPV*DN2 + DPF1*(DUT*DT2)
      
C           Sum up to the integral with the weighting factor OMEGA

            DFWX=DFWX+AH1*OM
            DFWY=DFWY+AH2*OM
        
          END DO
          
        END IF
        
      END DO

C Correct the values according to the formulas of the
C drag/lift coefficients with the help of the information
C provided by PTSDAT:

      IF (DPF2.EQ.0D0) THEN
        DFWX=0D0
        DFWY=0D0
      ELSE
        DFWX=2D0*DFWX/DPF2
        DFWY=2D0*DFWY/DPF2
      END IF

      END
