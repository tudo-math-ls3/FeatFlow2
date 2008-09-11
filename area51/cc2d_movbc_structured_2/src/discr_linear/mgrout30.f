* **********************************************************************
* This file contains extended prolongation/restriction routines for
* the E030/EM30 element. They support adaptive switching to constant 
* pronongation/restriction depending on the aspect ratio of the element.
*
* The file also contains examples of the routines YPRxxx that can
* be used with the M011 solver for prolongation/restriction, but
* these routines are commented out. In the present implementation
* the YPROLU/YPROLP routines in MGROUT.F are the actual routines
* which call the correct prol/rest.
* 
* **********************************************************************

************************************************************************
* MP030  - Standard prolongation with averaging of the values of a node
*          on the edge of two elements  
*
* Parameters:
* DU1    - coarse grid vector       
* KVERT1 - array with node numbers of the element-nodes on the 
*          coarse grid
* KVERT1 - array with node numbers of the element-nodes on the 
*          fine grid
* KMID1  - array with node numbers of the element-midpoints on the 
*          coarse grid
* KMID2  - array with node numbers of the element-midpoints on the 
*          fibne grid
* KADJ1  - array describing the neighbourhood of an element on the
*          coarse grid
* KADJ2  - array describing the neighbourhood of an element on the
*          fine rid
* NVT1   - number of nodes on the coarse grid
* NVT2   - number of nodes on the fine grid
* NEL1   - number of elements on the coarse grid
* NEL2   - number of elements on the fine grid
* NMT2   - number of midpoints on the fine grid
*
* Output:
* DU2    - fine grid vector
************************************************************************

      SUBROUTINE MP030(DU1,DU2,KVERT1,KVERT2,KMID1,KMID2,
     *                 KADJ1,KADJ2,NVT1,NVT2,NEL1,NEL2,NMT2)
     
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
     
C parameters
     
      DOUBLE PRECISION DU1(*),DU2(*)
      INTEGER KVERT1(NNVE,*),KVERT2(NNVE,*),KMID1(NNVE,*),KMID2(NNVE,*)
      INTEGER KADJ1(NNVE,*),KADJ2(NNVE,*)
      INTEGER NVT1, NVT2, NEL1, NEL2, NMT2
      
C local variables

      DOUBLE PRECISION A1, A2, A3, A4, A5, A6, A7, A8
      INTEGER IEL1
      INTEGER IM1,IM2,IM3,IM4
      INTEGER IELH1,IELH2,IELH3,IELH4
      INTEGER IA, IB, IC
      DOUBLE PRECISION DUH1,DUH2,DUH3,DUH4

C Die folgenden Parameter unterscheiden die Prolongation/Restriktion
C von EM30 zu EM31 !
C Alle Koeffizienten hier sind bereits halbiert, so dass bei der 
C Mittelung der Funktionswerte in den neuen Knoten nicht mehr halbiert 
C werden muss:
      PARAMETER (A1=0.5D0,A2=-0.0625D0,A3=0D0,A4=0.0625D0)
C      PARAMETER (A1=0.5D0,A2=-0.125D0,A3=0D0,A4=0.125D0)
      PARAMETER (A5=0.625D0,A6=0.125D0,A7=0.125D0,A8=0.125D0)

C *** Zero initialization of (DU2)
      CALL  LCL1 (DU2,NMT2)
C
      DO 10 IEL1=1,NEL1
C
      IM1=KMID1(1,IEL1)-NVT1
      IM2=KMID1(2,IEL1)-NVT1
      IM3=KMID1(3,IEL1)-NVT1
      IM4=KMID1(4,IEL1)-NVT1
C
      DUH1=DU1(IM1)
      DUH2=DU1(IM2)
      DUH3=DU1(IM3)
      DUH4=DU1(IM4)
C
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)
C
C *** The edge IM1 and the corresponding fine inner node
C
      IF (KADJ1(1,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH1)-NVT2
       IB=KMID2(4,IELH2)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH1+A2*DUH2+A3*DUH3+A4*DUH4
       DU2(IB)=DU2(IB)+   A1*DUH1+A4*DUH2+A3*DUH3+A2*DUH4
      ELSE
       IA=KMID2(1,IELH1)-NVT2
       IB=KMID2(4,IELH2)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH1+A2*DUH2+A3*DUH3+A4*DUH4)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH1+A4*DUH2+A3*DUH3+A2*DUH4)
      ENDIF
      IC=KMID2(2,IELH1)-NVT2
      DU2(IC)=A5*DUH1+A6*(DUH2+DUH4)+A7*DUH3
C
C *** The edge IM2 and the corresponding fine inner node
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH2)-NVT2
       IB=KMID2(4,IELH3)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH2+A2*DUH3+A3*DUH4+A4*DUH1
       DU2(IB)=DU2(IB)+   A1*DUH2+A4*DUH3+A3*DUH4+A2*DUH1
      ELSE
       IA=KMID2(1,IELH2)-NVT2
       IB=KMID2(4,IELH3)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH2+A2*DUH3+A3*DUH4+A4*DUH1)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH2+A4*DUH3+A3*DUH4+A2*DUH1)
      ENDIF
      IC=KMID2(2,IELH2)-NVT2
      DU2(IC)=A5*DUH2+A6*(DUH3+DUH1)+A7*DUH4
C
C *** The edge IM3 and the corresponding fine inner node
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH3)-NVT2
       IB=KMID2(4,IELH4)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH3+A2*DUH4+A3*DUH1+A4*DUH2
       DU2(IB)=DU2(IB)+   A1*DUH3+A4*DUH4+A3*DUH1+A2*DUH2
      ELSE
       IA=KMID2(1,IELH3)-NVT2
       IB=KMID2(4,IELH4)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH3+A2*DUH4+A3*DUH1+A4*DUH2)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH3+A4*DUH4+A3*DUH1+A2*DUH2)
      ENDIF
      IC=KMID2(2,IELH3)-NVT2
      DU2(IC)=A5*DUH3+A6*(DUH4+DUH2)+A7*DUH1
C
C *** The edge IM4 and the corresponding fine inner node
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH4)-NVT2
       IB=KMID2(4,IELH1)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH4+A2*DUH1+A3*DUH2+A4*DUH3
       DU2(IB)=DU2(IB)+   A1*DUH4+A4*DUH1+A3*DUH2+A2*DUH3
      ELSE
       IA=KMID2(1,IELH4)-NVT2
       IB=KMID2(4,IELH1)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH4+A2*DUH1+A3*DUH2+A4*DUH3)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH4+A4*DUH1+A3*DUH2+A2*DUH3)
      ENDIF
      IC=KMID2(2,IELH4)-NVT2
      DU2(IC)=A5*DUH4+A6*(DUH1+DUH3)+A7*DUH2
C
C
10    CONTINUE
C
C
      END



************************************************************************
* MP030X - Extended prolongartion. Takes care of the aspect ratio
*          and allows different weights for the two contributions 
*          of the elements adjacent to an edge for the value of a node.
*
* Parameters:
* DU1    - coarse grid vector       
* KVERT1 - array with node numbers of the element-nodes on the 
*          coarse grid
* KVERT1 - array with node numbers of the element-nodes on the 
*          fine grid
* KMID1  - array with node numbers of the element-midpoints on the 
*          coarse grid
* KMID2  - array with node numbers of the element-midpoints on the 
*          fibne grid
* KADJ1  - array describing the neighbourhood of an element on the
*          coarse grid
* KADJ2  - array describing the neighbourhood of an element on the
*          fine rid
* NVT1   - number of nodes on the coarse grid
* NVT2   - number of nodes on the fine grid
* NEL1   - number of elements on the coarse grid
* NEL2   - number of elements on the fine grid
* NMT1   - number of midpoints on the coarse grid
* NMT2   - number of midpoints on the fine grid
* DCORVG - array with coordinates of the nodes in the coarse grid
* AREA1  - array containing the areas of the elements on the coarse 
*          grid
* IINT1  - type of the averaging on the element edges:             
*     <=0: standard averaging of both contributions by 1/2          
*      =1: weighted averaging of the interpolated function values:
*          The area of the current coarse grid element determines 
*          the weight. (L2-projection, standard)
*      =2: weighted averaging of the interpolated function values:
*          The area of the neightbour element of the coarse grid 
*          the weight. 
* PARP   - upper bound aspect ratio; for all elements with higher AR
*          the prolongation is switched to constant prolongation 
* IADPR1 - controls switching to constant prolongation:     
*      =1: switch depending on aspect ratio of current element
*      =2: switch depending on aspect ratio of current element and
*          neighbour element
*
* Output:
* DU2    - fine grid vector
************************************************************************
 
      SUBROUTINE MP030X(DU1,DU2,KVERT1,KVERT2,KMID1,KMID2,KADJ1,KADJ2,
     *               NVT1,NVT2,NEL1,NEL2,NMT1,NMT2,DCORVG,
     *               AREA1,IINT1,PARP,IADPR1)

      IMPLICIT NONE

      INTEGER NNVE
      PARAMETER (NNVE=4)

C Common blocks for output (in case of errors, debug,...)

      INTEGER M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      
C local variables
      
      DOUBLE PRECISION DU1(*),DU2(*),DCORVG(2,*)
      DOUBLE PRECISION AREA1(*)
      INTEGER KVERT1(NNVE,*),KVERT2(NNVE,*),
     *        KMID1(NNVE,*),KMID2(NNVE,*),KADJ1(NNVE,*),KADJ2(NNVE,*)
      INTEGER NVT1,NVT2,NEL1,NEL2,NMT1,NMT2
      
      INTEGER IEL1, IA, IB, IC, NCEL
      INTEGER IM1, IM2, IM3, IM4
      INTEGER IELA1, IELA2, IELA3, IELA4, IELH1, IELH2, IELH3, IELH4
      INTEGER IGRID1, IGRID2, IGRID3, IGRID4
      DOUBLE PRECISION DAREA, DAREA1, DAREA2, DAREA3, DAREA4
      DOUBLE PRECISION ARIEL, ARIEL1, ARIEL2, ARIEL3, ARIEL4
      DOUBLE PRECISION DUH1, DUH2, DUH3, DUH4
      DOUBLE PRECISION WEIGH1, WEIGH2, WEIGH3, WEIGH4

      INTEGER IINT1, IADPR1
      DOUBLE PRECISION PARP

C Die folgenden Parameter unterscheiden die Prolongation/Restriktion
C von EM30 zu EM31 !
C PRWEIG (.,1) gibt die Konstanten fuer die Standard-Prolongation an,
C PRWEIG (.,2) die Konstanten fuer die konstante Prolongation.

      DOUBLE PRECISION PRWEIG (8,2)
      DATA PRWEIG    /1D0,-0.25D0,0D0,0.25D0,
     *                0.625D0,0.125D0,0.125D0,0.125D0,
     *                1D0,0D0,0D0,0D0,
     *                1D0,0D0,0D0,0D0/ 
      SAVE PRWEIG

C Zielvektor mit 0 vorbelegen:

      CALL  LCL1 (DU2,NMT2)

C Anzahl der Elemente mit zu grossem Aspect-Ratio zaehlen

      NCEL = 0

C Alle Elemente mit IEL1 durchlaufen

      DO IEL1=1,NEL1

C Nummern der Seitenmittelpunkte IM1, IM2, IM3, IM4 auf dem Element
C im Grobgitter besorgen. NVT1 davon abziehen, damit global erster 
C Seitenmittelpunkt Nummer 1 hat (dieser steht im Vektor an Position 1!).

      IM1=KMID1(1,IEL1)-NVT1
      IM2=KMID1(2,IEL1)-NVT1
      IM3=KMID1(3,IEL1)-NVT1
      IM4=KMID1(4,IEL1)-NVT1

C Funktionswerte auf diesen Mittelpunkten ermitteln.
C Der Vektor enthaelt die Funktionswerte in den Seitenmitten
C da das Element keine Funktionswerte in den Ecken besitzt:
C DU1(IMj) ist der Funktionswert im j'ten Seitenmittelpunkt:

      DUH1=DU1(IM1)
      DUH2=DU1(IM2)
      DUH3=DU1(IM3)
      DUH4=DU1(IM4)

C Nummern der benachbarten Elemente im Grobgitter ermitteln:
C
C           +--------+
C           |        |
C           | IELA3  |
C           |        |
C  +--------4--------3--------+
C  |        |        |        |
C  | IELA4  |  IEL   | IELA2  |
C  |        |        |        |
C  +--------1--------2--------+
C           |        |
C           | IELA1  |
C           |        |
C           +--------+

      IELA1=KADJ1(1,IEL1)
      IELA2=KADJ1(2,IEL1)
      IELA3=KADJ1(3,IEL1)
      IELA4=KADJ1(4,IEL1)

C Berechnung des Verzerrungsfaktors des aktuellen Elementes IEL1
C mit Hilfe der Funktion ARCALC in der Datei ARCALC.F.
C Ergebnis wird in ARIEL zurueckgeliefert.

      CALL ARCALC(ARIEL ,IEL1 ,KVERT1,DCORVG)
C Evtl. Kehrwert bilden, falls Element in falsche Koordinaten-
C richtung verzerrt ist, damit ARIEL >=1 ist:
      IF (ARIEL.LT.1D0) ARIEL=1D0/ARIEL

C Flaeche des aktuellen Elementes auf dem Grobgitter bestimmen fuer
C die spaetere Berechnung der Gewichtungsfaktoren:

      DAREA=AREA1(IEL1)

C Betrachten wir das erste Nachbarelement IELA1.
C Falls es so ein Element gibt, dann berechnen wir ebenfalls genauso
C wie oben den Verzerrungsfaktor und den Flaecheninhalt:

      IF (IELA1.NE.0) THEN
       CALL ARCALC(ARIEL1,IELA1,KVERT1,DCORVG)
       IF (ARIEL1.LT.1D0) ARIEL1=1D0/ARIEL1
       DAREA1=AREA1(IELA1)
      ELSE
       ARIEL1=0D0
       DAREA1=0D0
      ENDIF

C Das gleiche dann natuerlich auch fuer das Nachbarelement IELA2:

      IF (IELA2.NE.0) THEN
       CALL ARCALC(ARIEL2,IELA2,KVERT1,DCORVG)
       IF (ARIEL2.LT.1D0) ARIEL2=1D0/ARIEL2
       DAREA2=AREA1(IELA2)
      ELSE
       ARIEL2=0D0
       DAREA2=0D0
      ENDIF

C fuer das Nachbarelement IELA3:

      IF (IELA3.NE.0) THEN
       CALL ARCALC(ARIEL3,IELA3,KVERT1,DCORVG)
       IF (ARIEL3.LT.1D0) ARIEL3=1D0/ARIEL3
       DAREA3=AREA1(IELA3)
      ELSE
       ARIEL3=0D0
       DAREA3=0D0
      ENDIF

C und fuer das Nachbarelement IELA4:

      IF (IELA4.NE.0) THEN
       CALL ARCALC(ARIEL4,IELA4,KVERT1,DCORVG)
       IF (ARIEL4.LT.1D0) ARIEL4=1D0/ARIEL4
       DAREA4=AREA1(IELA4)
      ELSE
       ARIEL4=0D0
       DAREA4=0D0
      ENDIF

C Nummern der Elemente auf dem Feingitter bestimmen, die im aktuellen
C Element auf dem Grobgitter enthalten sind.
C Jedes Element auf dem Grobgitter enthaelt (wg. der regulaeren 
C Verfeinerung) 4 kleine Elemente auf dem Feingitter. 
C Eines der Elemente hat die gleiche Nummer wie das Element auf dem
C Grobgitter. An die anderen Nummern kommen wir ueber das Array
C KADJ2 dran:
C
C IELH1 := Nummer des aktuellen Elementes 
C        = Nummer des Feingitterelementes rechts unten im 
C          Grobgitterelement
C IELH2 := Nummer des an der 2. Kante von IELH1 anliegenden Elementes
C        = Nummer des Feingitterelementes rechts unten
C IELH3 := Nummer des an der 2. Kante von IELH2 anliegenden Elementes
C        = Nummer des Feingitterelementes rechts oben
C IELH4 := Nummer des an der 2. Kante von IELH3 anliegenden Elementes
C        = Nummer des Feingitterelementes links oben

      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)

C Nummern der Seitenmittelpunkte im Feingitter bestimmen, die bei
C der Prolongation benutzt werden.
C
C NVT2 abziehen, damit global erster Seitenmittelpunkt Nummer 1 hat
C (dieser steht im Vektor an Position 1!).

      IA=KMID2(1,IELH1)-NVT2
      IB=KMID2(4,IELH2)-NVT2
      IC=KMID2(2,IELH1)-NVT2

C Es liegt nun folgende Situation vor:
C
C   4               IM3                3
C     ===============X================
C     |              |               |
C     |              |               |
C     |    IELH4     |     IELH3     |
C     |              |               |
C     |                              |
C IM4 X----------- IEL1 -------------X IM2
C     |                              |
C     |              |               |
C     |    IELH1     o IC  IELH2     |
C     |              |               |
C     |              |               |
C   1 =======o=======X=======o======== 2
C     |     IA      IM1      IB      |
C     |                              |
C     |                              |
C     |                              |
C     |                              |
C     |            IELA1             |
C     |                              |
C     |                              |
C     |                              |
C     |                              |
C     |                              |
C     ================================

C Gewichtungsfaktoren fuer die Interpolation berechnen:

C IINT1<=0: Ungewichtete Mittelung der interpolierten Funktionswerte 
C           bei der Berechnung des neuen Funktionswertes im
C           Seitenmittelpunkt auf dem Feingitter
C IINT1=1: Gewichtete Mittelung der interpolierten Funktionswerte:
C          Flaecheninhalt des aktuellen Grobgitterelementes bestimmt
C          das Gewicht.
C IINT1=2: Gewichtete Mittelung der interpolierten Funktionswerte:
C          Flaecheninhalt des benachbarten Grobgitterelementes bestimmt
C          das Gewicht.

      IF (IINT1.LE.0) THEN 
        WEIGH1=0.5D0
        WEIGH2=0.5D0
        WEIGH3=0.5D0
        WEIGH4=0.5D0
      ELSE IF (IINT1.EQ.1) THEN
        WEIGH1=DAREA /(DAREA+DAREA1)
        WEIGH2=DAREA /(DAREA+DAREA2)
        WEIGH3=DAREA /(DAREA+DAREA3)
        WEIGH4=DAREA /(DAREA+DAREA4)
      ELSE IF (IINT1.GE.2) THEN
        WEIGH1=DAREA1/(DAREA+DAREA1)
        WEIGH2=DAREA2/(DAREA+DAREA2)
        WEIGH3=DAREA3/(DAREA+DAREA3)
        WEIGH4=DAREA4/(DAREA+DAREA4)
      END IF

C Wenn an einer Kante kein Nachbarelement existiert, so muss der
C Gewichtungsfaktor auf 1 gesetzt werden, da dort nichts zusammengezaehlt
C werden muss/kann/darf:

      IF (IELA1.EQ.0) WEIGH1=1D0
      IF (IELA2.EQ.0) WEIGH2=1D0
      IF (IELA3.EQ.0) WEIGH3=1D0
      IF (IELA4.EQ.0) WEIGH4=1D0

C Umschalter zwischen Standard- und konstanter Prolongation:
C IGRIDi=1=Standard-Prolongation an der Elementseite i
C IGRIDi=2=konstante Prolongation an der Elementseite i

      IGRID1=1
      IGRID2=1
      IGRID3=1
      IGRID4=1
      
C Nun vergleichen: falls der Verzerrungsfaktor des aktuellen Elementes 
C (oder im Fall von IADPR1>=2 das auch des Nachbarelementes)
C zu gross ist, so wird die Restriktion an dieser
C Seite des Elementes auf Konstant umgeschaltet:

      IF ((IADPR1.GE.1).AND.(PARP.NE.-1D0)) THEN
        IF ((ARIEL.GT.PARP)) IGRID1=2
        IF ((ARIEL.GT.PARP)) IGRID2=2
        IF ((ARIEL.GT.PARP)) IGRID3=2
        IF ((ARIEL.GT.PARP)) IGRID4=2
        IF (IADPR1.GE.2) THEN
          IF (ARIEL1.GT.PARP) IGRID1=2
          IF (ARIEL2.GT.PARP) IGRID2=2
          IF (ARIEL3.GT.PARP) IGRID3=2
          IF (ARIEL4.GT.PARP) IGRID4=2
        END IF
        IF ((IGRID1.EQ.2).OR.(IGRID2.EQ.2).OR.(IGRID3.EQ.2).OR.
     *      (IGRID4.EQ.2)) NCEL = NCEL+1
      END IF  

C Prolongation durchfuehren: Beitraege des aktuellen Elementes zu
C den Funktionswerten im Feingitter in IA, IB und IC berechnen und
C auf den prolongierten Vektor draufaddieren:

      DU2(IA)=DU2(IA)+WEIGH1*
     *        (PRWEIG(1,IGRID1)*DUH1+PRWEIG(2,IGRID1)*DUH2+
     *         PRWEIG(3,IGRID1)*DUH3+PRWEIG(4,IGRID1)*DUH4)
      DU2(IB)=DU2(IB)+WEIGH1*
     *        (PRWEIG(1,IGRID1)*DUH1+PRWEIG(4,IGRID1)*DUH2+
     *         PRWEIG(3,IGRID1)*DUH3+PRWEIG(2,IGRID1)*DUH4)
      DU2(IC)=PRWEIG(5,IGRID1)*DUH1+PRWEIG(6,IGRID1)*(DUH2+DUH4)+
     *        PRWEIG(7,IGRID1)*DUH3

C Fertig. Dies war die Prolongation um den Seitenmittelpunkt IM1.
C
C Als naechstes folgen gleichen Prolongationsbefehle nochmals fuer 
C die neuen Seitenmittelpunkte um die Punkte IM2, IM3 und IM4 herum:

C Neue Seitenmitten um IM2:

      IA=KMID2(1,IELH2)-NVT2
      IB=KMID2(4,IELH3)-NVT2
      IC=KMID2(2,IELH2)-NVT2

      DU2(IA)=DU2(IA)+WEIGH2*
     *        (PRWEIG(1,IGRID2)*DUH2+PRWEIG(2,IGRID2)*DUH3+
     *         PRWEIG(3,IGRID2)*DUH4+PRWEIG(4,IGRID2)*DUH1)
      DU2(IB)=DU2(IB)+WEIGH2*
     *        (PRWEIG(1,IGRID2)*DUH2+PRWEIG(4,IGRID2)*DUH3+
     *         PRWEIG(3,IGRID2)*DUH4+PRWEIG(2,IGRID2)*DUH1)
      DU2(IC)=PRWEIG(5,IGRID2)*DUH2+PRWEIG(6,IGRID2)*(DUH3+DUH1)+
     *        PRWEIG(7,IGRID2)*DUH4

C Neue Seitenmitten um IM3:

      IA=KMID2(1,IELH3)-NVT2
      IB=KMID2(4,IELH4)-NVT2
      IC=KMID2(2,IELH3)-NVT2

      DU2(IA)=DU2(IA)+WEIGH3*
     *        (PRWEIG(1,IGRID3)*DUH3+PRWEIG(2,IGRID3)*DUH4+
     *         PRWEIG(3,IGRID3)*DUH1+PRWEIG(4,IGRID3)*DUH2)
      DU2(IB)=DU2(IB)+WEIGH3*
     *        (PRWEIG(1,IGRID3)*DUH3+PRWEIG(4,IGRID3)*DUH4+
     *         PRWEIG(3,IGRID3)*DUH1+PRWEIG(2,IGRID3)*DUH2)
      DU2(IC)=PRWEIG(5,IGRID3)*DUH3+PRWEIG(6,IGRID3)*(DUH4+DUH2)+
     *        PRWEIG(7,IGRID3)*DUH1

C Neue Seitenmitten um IM4:

      IA=KMID2(1,IELH4)-NVT2
      IB=KMID2(4,IELH1)-NVT2
      IC=KMID2(2,IELH4)-NVT2

      DU2(IA)=DU2(IA)+WEIGH4*
     *        (PRWEIG(1,IGRID4)*DUH4+PRWEIG(2,IGRID4)*DUH1+
     *         PRWEIG(3,IGRID4)*DUH2+PRWEIG(4,IGRID4)*DUH3)
      DU2(IB)=DU2(IB)+WEIGH4*
     *        (PRWEIG(1,IGRID4)*DUH4+PRWEIG(4,IGRID4)*DUH1+
     *         PRWEIG(3,IGRID4)*DUH2+PRWEIG(2,IGRID4)*DUH3)
      DU2(IC)=PRWEIG(5,IGRID4)*DUH4+PRWEIG(6,IGRID4)*(DUH1+DUH3)+
     *        PRWEIG(7,IGRID4)*DUH2

      END DO

      END


************************************************************************
* MR030  - Standard restriction with averaging of the values of a node
*          on the edge of two elements  
*
* Parameters:
* DU2    - fine grid vector
* KVERT1 - array with node numbers of the element-nodes on the 
*          coarse grid
* KVERT1 - array with node numbers of the element-nodes on the 
*          fine grid
* KMID1  - array with node numbers of the element-midpoints on the 
*          coarse grid
* KMID2  - array with node numbers of the element-midpoints on the 
*          fibne grid
* KADJ1  - array describing the neighbourhood of an element on the
*          coarse grid
* KADJ2  - array describing the neighbourhood of an element on the
*          fine rid
* NVT1   - number of nodes on the coarse grid
* NVT2   - number of nodes on the fine grid
* NEL1   - number of elements on the coarse grid
* NEL2   - number of elements on the fine grid
* NMT1   - number of midpoints on the coarse grid
*
* Output:
* DU1    - coarse grid vector       
************************************************************************

      SUBROUTINE MR030(DU2,DU1,KVERT2,KVERT1,KMID2,KMID1,
     *                 KADJ2,KADJ1,NVT2,NVT1,NEL2,NEL1)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
     
C parameters
     
      DOUBLE PRECISION DU1(*),DU2(*)
      INTEGER KVERT1(NNVE,*),KVERT2(NNVE,*),KMID1(NNVE,*),KMID2(NNVE,*)
      INTEGER KADJ1(NNVE,*),KADJ2(NNVE,*)
      INTEGER NVT1, NVT2, NEL1, NEL2, NMT2
      
C local variables

      DOUBLE PRECISION A1, A2, A3, A4, A5, A6, A7, A8
      INTEGER IEL1
      INTEGER IM1,IM2,IM3,IM4
      INTEGER IELH1,IELH2,IELH3,IELH4
      INTEGER IA, IB, IC
      INTEGER I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12

      DOUBLE PRECISION DUH1,DUH2,DUH3,DUH4,DUH5,DUH6,DUH7,DUH8
      DOUBLE PRECISION DUH9,DUH10,DUH11,DUH12

C Die folgenden Parameter unterscheiden die Prolongation/Restriktion
C von EM30 zu EM31 !:
      PARAMETER (A1=1D0,A2=-0.0625D0,A3=0D0,A4=0.0625D0)

C      PARAMETER (A1=1D0,A2=-0.125D0,A3=0D0,A4=0.125D0)
      PARAMETER (A5=0.625D0,A6=0.125D0,A7=0.125D0,A8=0.125D0)

      DO 10 IEL1=1,NEL1
C
      IM1=KMID1(1,IEL1)-NVT1
      IM2=KMID1(2,IEL1)-NVT1
      IM3=KMID1(3,IEL1)-NVT1
      IM4=KMID1(4,IEL1)-NVT1
C
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)
C
      I1=KMID2(1,IELH1)-NVT2
      I2=KMID2(4,IELH2)-NVT2
      I3=KMID2(1,IELH2)-NVT2
      I4=KMID2(4,IELH3)-NVT2
      I5=KMID2(1,IELH3)-NVT2
      I6=KMID2(4,IELH4)-NVT2
      I7=KMID2(1,IELH4)-NVT2
      I8=KMID2(4,IELH1)-NVT2
      I9=KMID2(2,IELH1)-NVT2
      I10=KMID2(2,IELH2)-NVT2
      I11=KMID2(2,IELH3)-NVT2
      I12=KMID2(2,IELH4)-NVT2
C
      DUH1= DU2(I1)
      DUH2= DU2(I2)
      DUH3= DU2(I3)
      DUH4= DU2(I4)
      DUH5= DU2(I5)
      DUH6= DU2(I6)
      DUH7= DU2(I7)
      DUH8= DU2(I8)
      DUH9= DU2(I9)
      DUH10=DU2(I10)
      DUH11=DU2(I11)
      DUH12=DU2(I12)
C
C
C *** The edge IM1
C
      IF (KADJ1(1,IEL1).NE.0) THEN
C     case of an inner edge
        IF (KADJ1(1,IEL1).GT.IEL1) THEN
           DU1(IM1)= A1*(DUH1+DUH2)+A2*(DUH4+DUH7)
     *                 +A3*(DUH5+DUH6)+A4*(DUH3+DUH8)
     *                   +A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
        ELSE
           DU1(IM1)=DU1(IM1)+A2*(DUH4+DUH7)
     *                        +A3*(DUH5+DUH6)+A4*(DUH3+DUH8)
     *                          +A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
        ENDIF
      ELSE
       DU1(IM1)=     A1*(DUH1+DUH2)+2D0*A2*(DUH4+DUH7)
     *          +2D0*A3*(DUH5+DUH6)+2D0*A4*(DUH3+DUH8)
     *          +    A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
      ENDIF
C
C *** The edge IM2
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
C     case of an inner edge
        IF (KADJ1(2,IEL1).GT.IEL1) THEN
           DU1(IM2)= A1*(DUH3+DUH4)+A2*(DUH6+DUH1)
     *                 +A3*(DUH7+DUH8)+A4*(DUH5+DUH2)
     *                   +A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
        ELSE
           DU1(IM2)=DU1(IM2)+A2*(DUH6+DUH1)
     *                        +A3*(DUH7+DUH8)+A4*(DUH5+DUH2)
     *                          +A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
        ENDIF
      ELSE
       DU1(IM2)=     A1*(DUH3+DUH4)+2D0*A2*(DUH6+DUH1)
     *          +2D0*A3*(DUH7+DUH8)+2D0*A4*(DUH5+DUH2)
     *          +    A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
      ENDIF
C
C *** The edge IM3
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
C     case of an inner edge
        IF (KADJ1(3,IEL1).GT.IEL1) THEN
           DU1(IM3)= A1*(DUH5+DUH6)+A2*(DUH8+DUH3)
     *                 +A3*(DUH1+DUH2)+A4*(DUH7+DUH4)
     *                   +A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
        ELSE
           DU1(IM3)=DU1(IM3)+A2*(DUH8+DUH3)
     *                        +A3*(DUH1+DUH2)+A4*(DUH7+DUH4)
     *                          +A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
        ENDIF
      ELSE
       DU1(IM3)=     A1*(DUH5+DUH6)+2D0*A2*(DUH8+DUH3)
     *          +2D0*A3*(DUH1+DUH2)+2D0*A4*(DUH7+DUH4)
     *          +    A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
      ENDIF
C
C *** The edge IM4
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
C     case of an inner edge
        IF (KADJ1(4,IEL1).GT.IEL1) THEN
           DU1(IM4)= A1*(DUH7+DUH8)+A2*(DUH2+DUH5)
     *                 +A3*(DUH3+DUH4)+A4*(DUH1+DUH6)
     *                   +A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
        ELSE
           DU1(IM4)=DU1(IM4)+A2*(DUH2+DUH5)
     *                        +A3*(DUH3+DUH4)+A4*(DUH1+DUH6)
     *                          +A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
        ENDIF
      ELSE
       DU1(IM4)=     A1*(DUH7+DUH8)+2D0*A2*(DUH2+DUH5)
     *          +2D0*A3*(DUH3+DUH4)+2D0*A4*(DUH1+DUH6)
     *          +    A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
      ENDIF

10    CONTINUE

      END

************************************************************************
* MR030X - Extended restriction. Takes care of the aspect ratio
*          and allows different weights for the two contributions 
*          of the elements adjacent to an edge for the value of a node.
*
* Parameters:
* DU2    - fine grid vector       
* KVERT1 - array with node numbers of the element-nodes on the 
*          coarse grid
* KVERT1 - array with node numbers of the element-nodes on the 
*          fine grid
* KMID1  - array with node numbers of the element-midpoints on the 
*          coarse grid
* KMID2  - array with node numbers of the element-midpoints on the 
*          fibne grid
* KADJ1  - array describing the neighbourhood of an element on the
*          coarse grid
* KADJ2  - array describing the neighbourhood of an element on the
*          fine rid
* NVT1   - number of nodes on the coarse grid
* NVT2   - number of nodes on the fine grid
* NEL1   - number of elements on the coarse grid
* NEL2   - number of elements on the fine grid
* NMT1   - number of midpoints on the coarse grid
* NMT2   - number of midpoints on the fine grid
* DCORVG - array with coordinates of the nodes in the coarse grid
* AREA1  - array containing the areas of the elements on the coarse 
*          grid
* IINT1  - type of the averaging on the element edges:             
*     <=0: standard averaging of both contributions by 1/2          
*      =1: weighted averaging of the interpolated function values:
*          The area of the current coarse grid element determines 
*          the weight. (L2-projection, standard)
*      =2: weighted averaging of the interpolated function values:
*          The area of the neightbour element of the coarse grid 
*          the weight. 
* PARP   - upper bound aspect ratio; for all elements with higher AR
*          the prolongation is switched to constant prolongation 
* IADPR1 - controls switching to constant prolongation:     
*      =1: switch depending on aspect ratio of current element
*      =2: switch depending on aspect ratio of current element and
*          neighbour element
*
* Output:
* DU1    - coarse grid vector
************************************************************************

      SUBROUTINE MR030X(DU2,DU1,KVERT2,KVERT1,KMID2,KMID1,KADJ2,KADJ1,
     *               NVT2,NVT1,NEL2,NEL1,NMT2,NMT1,DCORVG,
     *               AREA1,IINT1,PARR,IADPR1)

      IMPLICIT NONE

      INTEGER NNVE
      PARAMETER (NNVE=4)

C Common blocks for output (in case of errors, debug,...)

      INTEGER M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8

C local variables
      
      DOUBLE PRECISION DU1(*),DU2(*),DCORVG(2,*)
      DOUBLE PRECISION AREA1(*)
      INTEGER KVERT1(NNVE,*),KVERT2(NNVE,*),
     *        KMID1(NNVE,*),KMID2(NNVE,*),KADJ1(NNVE,*),KADJ2(NNVE,*)
      INTEGER NVT1,NVT2,NEL1,NEL2,NMT1,NMT2
      
      INTEGER IEL1, IA, IB, IC, NCEL
      INTEGER IM1, IM2, IM3, IM4
      INTEGER IELA1, IELA2, IELA3, IELA4, IELH1, IELH2, IELH3, IELH4
      INTEGER IGRID1, IGRID2, IGRID3, IGRID4
      INTEGER I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, I11, I12
      DOUBLE PRECISION DAREA, DAREA1, DAREA2, DAREA3, DAREA4
      DOUBLE PRECISION ARIEL, ARIEL1, ARIEL2, ARIEL3, ARIEL4
      DOUBLE PRECISION DUH1, DUH2, DUH3, DUH4, DUH5, DUH6
      DOUBLE PRECISION DUH7, DUH8, DUH9, DUH10, DUH11, DUH12
      DOUBLE PRECISION WEIGH1, WEIGH2, WEIGH3, WEIGH4

      INTEGER IINT1, IADPR1
      DOUBLE PRECISION PARR

C Die folgenden Parameter unterscheiden die Prolongation/Restriktion
C von EM30 zu EM31 !
C PRWEIG (.,1) gibt die Konstanten fuer die Standard-Prolongation an,
C PRWEIG (.,2) die Konstanten fuer die konstante Prolongation.

      DOUBLE PRECISION PRWEIG (8,2)
      DATA PRWEIG    /1D0,-0.25D0,0D0,0.25D0,
     *                0.625D0,0.125D0,0.125D0,0.125D0,
     *                1D0,0D0,0D0,0D0,
     *                1D0,0D0,0D0,0D0/ 
      SAVE PRWEIG

      CALL  LCL1 (DU1,NMT1)

C Alle Elemente mit IEL1 durchlaufen

      DO IEL1=1,NEL1

C Nummern der Seitenmittelpunkte IM1, IM2, IM3, IM4 auf dem Element
C im Grobgitter besorgen. NVT1 davon abziehen, damit global erster 
C Seitenmittelpunkt Nummer 1 hat (dieser steht im Vektor an Position 1!).

      IM1=KMID1(1,IEL1)-NVT1
      IM2=KMID1(2,IEL1)-NVT1
      IM3=KMID1(3,IEL1)-NVT1
      IM4=KMID1(4,IEL1)-NVT1

C Nummern der benachbarten Elemente im Grobgitter ermitteln:
C
C           +--------+
C           |        |
C           | IELA3  |
C           |        |
C  +--------4--------3--------+
C  |        |        |        |
C  | IELA4  |  IEL   | IELA2  |
C  |        |        |        |
C  +--------1--------2--------+
C           |        |
C           | IELA1  |
C           |        |
C           +--------+

      IELA1=KADJ1(1,IEL1)
      IELA2=KADJ1(2,IEL1)
      IELA3=KADJ1(3,IEL1)
      IELA4=KADJ1(4,IEL1)

C Berechnung des Verzerrungsfaktors des aktuellen Elementes IEL1
C mit Hilfe der Funktion ARCALC in der Datei ARCALC.F.
C Ergebnis wird in ARIEL zurueckgeliefert.

      CALL ARCALC(ARIEL ,IEL1 ,KVERT1,DCORVG)
C Evtl. Kehrwert bilden, falls Element in falsche Koordinaten-
C richtung verzerrt ist, damit ARIEL >=1 ist:
      IF (ARIEL.LT.1D0) ARIEL=1D0/ARIEL

C Flaeche des aktuellen Elementes auf dem Grobgitter bestimmen fuer
C die spaetere Berechnung der Gewichtungsfaktoren:

      DAREA=AREA1(IEL1)

C Betrachten wir das erste Nachbarelement IELA1.
C Falls es so ein Element gibt, dann berechnen wir ebenfalls genauso
C wie oben den Verzerrungsfaktor und den Flaecheninhalt:

      IF (IELA1.NE.0) THEN
       CALL ARCALC(ARIEL1,IELA1,KVERT1,DCORVG)
       IF (ARIEL1.LT.1D0) ARIEL1=1D0/ARIEL1
       DAREA1=AREA1(IELA1)
      ELSE
       ARIEL1=0D0
       DAREA1=0D0
      ENDIF

C Das gleiche dann natuerlich auch fuer das Nachbarelement IELA2:

      IF (IELA2.NE.0) THEN
       CALL ARCALC(ARIEL2,IELA2,KVERT1,DCORVG)
       IF (ARIEL2.LT.1D0) ARIEL2=1D0/ARIEL2
       DAREA2=AREA1(IELA2)
      ELSE
       ARIEL2=0D0
       DAREA2=0D0
      ENDIF

C fuer das Nachbarelement IELA3:

      IF (IELA3.NE.0) THEN
       CALL ARCALC(ARIEL3,IELA3,KVERT1,DCORVG)
       IF (ARIEL3.LT.1D0) ARIEL3=1D0/ARIEL3
       DAREA3=AREA1(IELA3)
      ELSE
       ARIEL3=0D0
       DAREA3=0D0
      ENDIF

C und fuer das Nachbarelement IELA4:

      IF (IELA4.NE.0) THEN
       CALL ARCALC(ARIEL4,IELA4,KVERT1,DCORVG)
       IF (ARIEL4.LT.1D0) ARIEL4=1D0/ARIEL4
       DAREA4=AREA1(IELA4)
      ELSE
       ARIEL4=0D0
       DAREA4=0D0
      ENDIF

C Nummern der Elemente auf dem Feingitter bestimmen, die im aktuellen
C Element auf dem Grobgitter enthalten sind.
C Jedes Element auf dem Grobgitter enthaelt (wg. der regulaeren 
C Verfeinerung) 4 kleine Elemente auf dem Feingitter. 
C Eines der Elemente hat die gleiche Nummer wie das Element auf dem
C Grobgitter. An die anderen Nummern kommen wir ueber das Array
C KADJ2 dran:
C
C IELH1 := Nummer des aktuellen Elementes 
C        = Nummer des Feingitterelementes rechts unten im 
C          Grobgitterelement
C IELH2 := Nummer des an der 2. Kante von IELH1 anliegenden Elementes
C        = Nummer des Feingitterelementes rechts unten
C IELH3 := Nummer des an der 2. Kante von IELH2 anliegenden Elementes
C        = Nummer des Feingitterelementes rechts oben
C IELH4 := Nummer des an der 2. Kante von IELH3 anliegenden Elementes
C        = Nummer des Feingitterelementes links oben

      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)

C Nummern der ganzen Seitenmittelpunkte auf dem Feingitter ermitteln,
C an denen die Funktionswerte gegeben sind. I1-I8 sind die Nummern der
C Seitenmittelpunkte auf dem Rand des Vierecks (gegen den Uhrzeitersinn)
C und I9-I12 die der Seitenmitten innerhalb des Vierecks (ebenfalls
C gegen den Uhzeigersinn)

      I1=KMID2(1,IELH1)-NVT2
      I2=KMID2(4,IELH2)-NVT2
      I3=KMID2(1,IELH2)-NVT2
      I4=KMID2(4,IELH3)-NVT2
      I5=KMID2(1,IELH3)-NVT2
      I6=KMID2(4,IELH4)-NVT2
      I7=KMID2(1,IELH4)-NVT2
      I8=KMID2(4,IELH1)-NVT2
      I9=KMID2(2,IELH1)-NVT2
      I10=KMID2(2,IELH2)-NVT2
      I11=KMID2(2,IELH3)-NVT2
      I12=KMID2(2,IELH4)-NVT2

C Funktionswerte auf diesen Mittelpunkten ermitteln.
C Der Vektor enthaelt die Funktionswerte in den Seitenmitten
C da das Element keine Funktionswerte in den Ecken besitzt:
C DU2(IMj) ist der Funktionswert im j'ten Seitenmittelpunkt:

      DUH1= DU2(I1)
      DUH2= DU2(I2)
      DUH3= DU2(I3)
      DUH4= DU2(I4)
      DUH5= DU2(I5)
      DUH6= DU2(I6)
      DUH7= DU2(I7)
      DUH8= DU2(I8)
      DUH9= DU2(I9)
      DUH10=DU2(I10)
      DUH11=DU2(I11)
      DUH12=DU2(I12)

C Es liegt nun folgende Situation vor:
C
C   4       I6      IM3      I5      3
C     =======o=======X=======o========
C     |              |               |
C     |              |               |
C  I7 o    IELH4     o I11 IELH3     o I4
C     |              |               |
C     |                              |
C IM4 X------o---- IEL1 -----o-------X IM2
C     |    I12               I10     |
C     |              |               |
C  I8 o    IELH1     o I9  IELH2     o I3
C     |              |               |
C     |              |               |
C   1 =======o=======X=======o======== 2
C     |     I1      IM1      I2      |
C     |                              |
C     |                              |
C     |                              |
C     |                              |
C     |            IELA1             |
C     |                              |
C     |                              |
C     |                              |
C     |                              |
C     |                              |
C     ================================

C Gewichtungsfaktoren fuer die Interpolation berechnen:

C IINT1<=0: Ungewichtete Mittelung der interpolierten Funktionswerte 
C           bei der Berechnung des neuen Funktionswertes im
C           Seitenmittelpunkt auf dem Feingitter
C IINT1=1: Gewichtete Mittelung der interpolierten Funktionswerte:
C          Flaecheninhalt des aktuellen Grobgitterelementes bestimmt
C          das Gewicht.
C IINT1=2: Gewichtete Mittelung der interpolierten Funktionswerte:
C          Flaecheninhalt des benachbarten Grobgitterelementes bestimmt
C          das Gewicht.

      IF (IINT1.LE.0) THEN 
        WEIGH1=0.5D0
        WEIGH2=0.5D0
        WEIGH3=0.5D0
        WEIGH4=0.5D0
      ELSE IF (IINT1.EQ.1) THEN
        WEIGH1=DAREA /(DAREA+DAREA1)
        WEIGH2=DAREA /(DAREA+DAREA2)
        WEIGH3=DAREA /(DAREA+DAREA3)
        WEIGH4=DAREA /(DAREA+DAREA4)
      ELSE IF (IINT1.GE.2) THEN
        WEIGH1=DAREA1/(DAREA+DAREA1)
        WEIGH2=DAREA2/(DAREA+DAREA2)
        WEIGH3=DAREA3/(DAREA+DAREA3)
        WEIGH4=DAREA4/(DAREA+DAREA4)
      END IF

C Wenn an einer Kante kein Nachbarelement existiert, so muss der
C Gewichtungsfaktor auf 1 gesetzt werden, da dort nichts zusammengezaehlt
C werden muss/kann/darf:

      IF (IELA1.EQ.0) WEIGH1=1D0
      IF (IELA2.EQ.0) WEIGH2=1D0
      IF (IELA3.EQ.0) WEIGH3=1D0
      IF (IELA4.EQ.0) WEIGH4=1D0

C Umschalter zwischen Standard- und konstanter Restriktion:
C IGRIDi=1=Standard-Restriktion an der Elementseite i
C IGRIDi=2=konstante Restriktionan der Elementseite i

      IGRID1=1
      IGRID2=1
      IGRID3=1
      IGRID4=1
      
C Nun vergleichen: falls der Verzerrungsfaktor des aktuellen Elementes 
C (oder im Fall von IADPR1>=2 das auch des Nachbarelementes)
C zu gross ist, so wird die Restriktion an dieser
C Seite des Elementes auf Konstant umgeschaltet:

      IF ((IADPR1.GE.1).AND.(PARR.NE.-1D0)) THEN
        IF ((ARIEL.GT.PARR)) IGRID1=2
        IF ((ARIEL.GT.PARR)) IGRID2=2
        IF ((ARIEL.GT.PARR)) IGRID3=2
        IF ((ARIEL.GT.PARR)) IGRID4=2
        IF (IADPR1.GE.2) THEN
          IF (ARIEL1.GT.PARR) IGRID1=2
          IF (ARIEL2.GT.PARR) IGRID2=2
          IF (ARIEL3.GT.PARR) IGRID3=2
          IF (ARIEL4.GT.PARR) IGRID4=2
        END IF
      END IF  

C Restriktion durchfuehren: Beitraege des aktuellen Elementes zu
C den Funktionswerten im Grobgitter in IM1 berechnen und
C auf den restringierten Vektor draufaddieren:

      DU1(IM1)=DU1(IM1)+WEIGH1*
     *      (PRWEIG(1,IGRID1)*(DUH1+DUH2)+PRWEIG(2,IGRID1)*(DUH4+DUH7)
     *      +PRWEIG(3,IGRID1)*(DUH5+DUH6)+PRWEIG(4,IGRID1)*(DUH3+DUH8))
     *      +PRWEIG(5,IGRID1)*DUH9+PRWEIG(6,IGRID1)*(DUH10+DUH12)
     *      +PRWEIG(7,IGRID1)*DUH11

C Fertig. Dies war die Restriktion auf den Seitenmittelpunkt IM1
C des aktuellen Elementes IEL1.
C
C Als naechstes folgen gleichen Restriktionsbefehle nochmals fuer 
C die anderen Seitenmitten IM2, IM3 und IM4:

C Restriktion auf IM2:

      DU1(IM2)=DU1(IM2)+WEIGH2*
     *      (PRWEIG(1,IGRID2)*(DUH3+DUH4)+PRWEIG(2,IGRID2)*(DUH6+DUH1)
     *      +PRWEIG(3,IGRID2)*(DUH7+DUH8)+PRWEIG(4,IGRID2)*(DUH5+DUH2))
     *      +PRWEIG(5,IGRID2)*DUH10+PRWEIG(6,IGRID2)*(DUH11+DUH9)
     *      +PRWEIG(7,IGRID2)*DUH12

C Restriktion auf IM3:

      DU1(IM3)=DU1(IM3)+WEIGH3*
     *      (PRWEIG(1,IGRID3)*(DUH5+DUH6)+PRWEIG(2,IGRID3)*(DUH8+DUH3)
     *      +PRWEIG(3,IGRID3)*(DUH1+DUH2)+PRWEIG(4,IGRID3)*(DUH7+DUH4))
     *      +PRWEIG(5,IGRID3)*DUH11+PRWEIG(6,IGRID3)*(DUH12+DUH10)
     *      +PRWEIG(7,IGRID3)*DUH9

C Restriktion auf IM4:

      DU1(IM4)=DU1(IM4)+WEIGH4*
     *      (PRWEIG(1,IGRID4)*(DUH7+DUH8)+PRWEIG(2,IGRID4)*(DUH2+DUH5)
     *      +PRWEIG(3,IGRID4)*(DUH3+DUH4)+PRWEIG(4,IGRID4)*(DUH1+DUH6))
     *      +PRWEIG(5,IGRID4)*DUH12+PRWEIG(6,IGRID4)*(DUH9+DUH11)
     *      +PRWEIG(7,IGRID4)*DUH10

      END DO

      END

************************************************************************
*
* YPRL30
*
*   Purpose: - performs the prolongation   DUF:=p(DUC)
*              with
*                  DUF   - fine correction vector on level ILEV
*                  DUC   - coarse correction vector on level ILEV-1
*  
*            - DUF and DUC have the structure  DU=(DU1,DU2,DP)
*  
************************************************************************
!      SUBROUTINE YPRL30(DUC,DUF)  
!      INCLUDE 'commons.for'
!      INCLUDE 'commonsmg.for'
!      INCLUDE 'commonsini.for'
!      DIMENSION DUF(*),DUC(*)
!
!C DEBUG-Code zum Prolongations-/Restriktionstest
!C      DO I=1,KNEQ(ILEV-1)
!C        DUC(I) = SIN(DBLE(I))*I-EXP(-LOG(DBLE(I)))
!C      END DO
!C      print *
!C      print *,'Prolongation vector input:'
!C      DO I=1,KNEQ(ILEV-1)
!C        WRITE (*,'(D24.14)') DUC(I)
!C      END DO
!
!C *** addresses for the fine level ILEV;
!C Don't change the current level ILEV but update the address fields
!      ISETLV=2
!      CALL  SETLEV (ILEV,ISETLV)
!C
!      KVERTF=L(LVERT)
!      KMIDF =L(LMID )
!      KADJF =L(LADJ )
!C
!C *** addresses for the coarse level ILEV-1
!      I1=ILEV-1
!      NXC=KNEQ(I1)
!      NELC=KNEL(I1)
!      NVTC=KNVT(I1)
!      NMTC=KNMT(I1)
!
!      KVERTC=L(KLVERT(I1))
!      KMIDC =L(KLMID (I1))
!      KADJC =L(KLADJ (I1))
!      KAREA =L(KLAREA (I1))
!      KCORVG=L(KLCVG(I1))
!      
!C evtl. ist IPR < 0, falls diese Prol. fuer EM31 verwendet wird !      
!      
!      IF (ABS(IPR).EQ.1) THEN
!C Standard-Prolongation ohne jegliche Gewichtung
!        CALL MP030(DUC,DUF,KWORK(KVERTC),KWORK(KVERTF),
!     *             KWORK(KMIDC),KWORK(KMIDF),KWORK(KADJC),KWORK(KADJF),
!     *             NVTC,NVT,NELC,NEL,NMT)
!      ELSE IF (ABS(IPR).EQ.2) THEN
!        CALL MP030C(DUC,DUF,KWORK(KVERTC),KWORK(KVERTF),
!     *              KWORK(KMIDC),KWORK(KMIDF),KWORK(KADJC),KWORK(KADJF),
!     *              NVTC,NVT,NELC,NEL,NMT,DWORK(KAREA),IAVRTP)
!      ELSE IF ((ABS(IPR).EQ.3).OR.(ABS(IPR).EQ.0)) THEN
!C Bitfeld für Umschaltung auf konstante P/R auswerten
!        IADPR1=0
!        IF(IAND(IADPRM,1).EQ.1) IADPR1=1
!        IF(IAND(IADPRM,9).EQ.9) IADPR1=2
!C Bitfeld für Gewichtung auswerten
!        IAVRT2=0
!        IF(IAND(IAVRTP,1).EQ.1) IAVRT2=1
!        IF(IAND(IAVRTP,5).EQ.5) IAVRT2=2
!        CALL MP030X(DUC,DUF,KWORK(KVERTC),KWORK(KVERTF),
!     *              KWORK(KMIDC),KWORK(KMIDF),KWORK(KADJC),KWORK(KADJF),
!     *              NVTC,NVT,NELC,NEL,NMTC,NMT,DWORK(KCORVG),
!     *              DWORK(KAREA),IAVRT2,DPREPS,IADPR1)
!      END IF
!
!C DEBUG-Code zum Prolongations-/Restriktionstest
!C      print *
!C      print *,'Prolongated to:'
!C      DO I=1,KNEQ(ILEV)
!C        WRITE (*,'(D24.14)') DUF(I)
!C      END DO
!C      CALL EXIT (255)
!
!      END
!
!************************************************************************
!*
!*   Purpose: - performs the defect restriction   DDC:=r(DDF)
!*              with
!*                  DDF - fine defect vector on level ILEV+1
!*                  DDC - coarse defect vector on level ILEV
!*  
!*            - DDF and DDC have the structure  DD=(DD1,DD2,DDP)
!*  
!************************************************************************
!      SUBROUTINE YRST30(DDF,DDC)  
!      INCLUDE 'commons.for'
!      INCLUDE 'commonsmg.for'
!      INCLUDE 'commonsini.for'
!      DIMENSION DDF(*),DDC(*)
!
!C DEBUG-Code zum Prolongations-/Restriktionstest
!C ----------------------------------------------
!C      DO I=1,KNEQ(ILEV+1)
!C        DDF(I) = SIN(DBLE(I))*I-EXP(-LOG(DBLE(I)))
!C      END DO
!C      print *
!C      print *,'Restriction vector input:'
!C      DO I=1,KNEQ(ILEV+1)
!C        WRITE (*,'(D24.14)') DDF(I)
!C      END DO
!C ----------------------------------------------
!C DEBUG Ende
!
!C Don't change the current level ILEV but update the address fields
!      ISETLV=2
!      CALL  SETLEV (ILEV,ISETLV)
!C
!      KVERTC=L(LVERT)
!      KMIDC =L(LMID )
!      KADJC =L(LADJ )
!C
!C *** addresses for the fine level ILEV+1
!      I1=ILEV+1
!      NXF=KNEQ(I1)
!      NELF=KNEL(I1)
!      NVTF=KNVT(I1)
!      NMTF=KNMT(I1)
!
!      KVERTF=L(KLVERT(I1))
!      KMIDF =L(KLMID (I1))
!      KADJF =L(KLADJ (I1))
!      KCORVG=L(KLCVG(I1))
!
!C Hier die Aspect Ratios des groeberen Gitters uebergeben !
!      KAREA =L(KLAREA (ILEV))
!  
!C evtl. ist IPR < 0, falls diese Restriktion fuer EM31 verwendet wird !
!
!      IF (ABS(IPR).EQ.1) THEN    
!        CALL MR030(DDF,DDC,KWORK(KVERTF),KWORK(KVERTC),
!     *             KWORK(KMIDF),KWORK(KMIDC),KWORK(KADJF),KWORK(KADJC),
!     *             NVTF,NVT,NELF,NEL)
!      ELSE IF (ABS(IPR).EQ.2) THEN
!        CALL MR030C(DDF,DDC,KWORK(KVERTF),KWORK(KVERTC),
!     *             KWORK(KMIDF),KWORK(KMIDC),KWORK(KADJF),KWORK(KADJC),
!     *             NVTF,NVT,NELF,NEL,NMT,DWORK(KAREA),IAVRTP)
!      ELSE IF ((ABS(IPR).EQ.3).OR.(ABS(IPR).EQ.0)) THEN
!C Bitfeld für Umschaltung auf konstante P/R auswerten
!        IADPR1=0
!        IF(IAND(IADPRM,2).EQ.2) IADPR1=1
!        IF(IAND(IADPRM,18).EQ.18) IADPR1=2
!C Bitfeld für Gewichtung auswerten
!        IAVRT2=0
!        IF(IAND(IAVRTP,2).EQ.2) IAVRT2=1
!        IF(IAND(IAVRTP,10).EQ.10) IAVRT2=2
!        CALL MR030X(DDF,DDC,KWORK(KVERTF),KWORK(KVERTC),
!     *             KWORK(KMIDF),KWORK(KMIDC),KWORK(KADJF),KWORK(KADJC),
!     *             NVTF,NVT,NELF,NEL,NMTF,NMT,DWORK(KCORVG),
!     *             DWORK(KAREA),IAVRT2,DPREPS,IADPR1)
!      END IF
!
!C DEBUG-Code zum Prolongations-/Restriktionstest
!C ----------------------------------------------
!C      print *
!C      print *,'Restricted to:'
!C      DO I=1,KNEQ(ILEV)
!C        WRITE (*,'(D24.14)') DDC(I)
!C      END DO
!C ----------------------------------------------
!C DEBUG Ende
!
!99999 END
!C
!************************************************************************
!C
!************************************************************************
!      SUBROUTINE YSTP30(DX,DD,DB,NEQ,ALPHA)  
!************************************************************************
!*
!*   PURPOSE: - PERFORMS STEP SIZE CONTROL FOR PROLONGATION WITH
!*                  DX    - OLD FINE SOLUTION VECTOR ON LEVEL ILEV
!*                  DD    - FINE CORRECTION VECTOR ON LEVEL ILEV
!*                  DB    - FINE RIGHT HAND SIDE VECTOR ON LEVEL ILEV
!*                  ALPHA - RELAXATION PARAMETER ACCORDING TO SOME
!*                          OPTIMIZATION CRITERION BUT WITHIN THE LIMITS
!*                          AMINU AND AMAXU (COMMON /RPARM/)
!*  
!*            - ADAPTIVE STRATEGIE NOT YET IMPLEMENTED
!*  
!************************************************************************
!      INCLUDE 'commons.for'
!      INCLUDE 'commonsmg.for'
!      INCLUDE 'commonsini.for'
!      DIMENSION DX(*),DD(*),DB(*)
!      
!C      print *
!C      print *,'YSTEP: Vektor'
!C      do i=1,neq
!C        write (*,'(F12.4)') dx (i)
!C      end do
!C      print *,'-----'      
!      
!C Zaehler/Nenner des Step-length-control-Parameters vorbelegen
!        DBX = 1D0
!        DBY = 1D0
!        
!        IF (ISLCTL.EQ.1) THEN
!C Step-length-control mit Energieminimierung
!
!C Zuerst ins entsprechende Level wechseln
!          ISETLV=2
!          ILOLD=ILEV
!          CALL SETLEV (ILEV,ISETLV)
!          LSLC1 = KLSLC1 (ILEV)
!          LSLC2 = KLSLC2 (ILEV)
!
!C Zaehler/Nenner fuer den Step-length-control-Parameter berechnen
!          CALL LCP1 (DB,DWORK(L(LSLC1)),NEQ)
!          CALL YAX(DX,DWORK(L(LSLC1)),NEQ,-1D0,1D0)
!
!          CALL LSP1(DD,DWORK(L(LSLC1)),NEQ,DBY)
!
!          CALL YAX(DD,DWORK(L(LSLC1)),NEQ,1D0,0D0)
!          CALL LSP1(DD,DWORK(L(LSLC1)),NEQ,DBX)
!
!C Level zurueckschalten
!          CALL SETLEV (ILOLD,ISETLV)
!
!        END IF
!        
!        IF (ISLCTL.EQ.2) THEN
!C Step-length-control mit Defektnorm-Minimierung
!
!C Zuerst ins entsprechende Level wechseln
!          ISETLV=2
!          ILOLD=ILEV
!          CALL SETLEV (ILEV,ISETLV)
!          LSLC1 = KLSLC1 (ILEV)
!          LSLC2 = KLSLC2 (ILEV)
!
!C Zaehler/Nenner fuer den Step-length-control-Parameter berechnen
!          CALL YAX(DD,DWORK(L(LSLC2)),NEQ,1D0,0D0)
!          CALL LSP1(DWORK(L(LSLC2)),DWORK(L(LSLC2)),NEQ,DBX)
!
!          CALL LCP1 (DB,DWORK(L(LSLC1)),NEQ)
!          CALL YAX(DX,DWORK(L(LSLC1)),NEQ,-1D0,1D0)
!
!          CALL LSP1(DWORK(L(LSLC1)),DWORK(L(LSLC2)),NEQ,DBY)
!
!C Level zurueckschalten
!          CALL SETLEV (ILOLD,ISETLV)
!
!        END IF
!        
!C Spezialfaelle abfangen und Schranken einhalten.
!C Falls ISLCTL=0 ist, so erzwingt die Vorbelegung von DBX/DBY den Parameter
!C ALPHA=1D0
!        IF (DBX.NE.0D0) THEN
!          ALPHA=DBY/DBX
!        ELSE
!C DBX=0 bei 0 als Startvektor der Iteration !
!          ALPHA=1D0
!        END IF
!        IF (ALPHA.GT.DSLMAX) ALPHA=DSLMAX
!        IF (ALPHA.LT.DSLMIN) ALPHA=DSLMIN
!
!        IF (ISLCTL.LE.-1) THEN
!C Fester step-length-control-Parameter; ueberschreibt alles
!          ALPHA = DSLALP
!        END IF
!
!        IF (IAND(MDBGT,1).NE.0) THEN
!          WRITE (MTERM,'(A,F21.14)') 'Step-length parameter ALPHA=',
!     *          ALPHA       
!        END IF
!        
!99999 END
!
