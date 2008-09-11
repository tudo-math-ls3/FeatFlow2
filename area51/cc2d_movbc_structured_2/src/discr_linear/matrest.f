
C **********************************************************************
C IGXYIX - Get X/Y-Index
C Performs a search for the index ILD such that matrix(IX,IY)=KLA(ILD) 
C holds. I.e. searches in the matrix array for the index belonging
C to the position IX/IY, so that the caller can modify this matrix
C element directly.
C **********************************************************************

      INTEGER FUNCTION IGXYIX (IX, IY, KCOLA1, KLDA1)

        IMPLICIT NONE
        
        INTEGER M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      
        COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
        
        INTEGER KLDA1, KCOLA1, ICOL, ILD, IX, IY

        DIMENSION KLDA1(*),KCOLA1(*)
        
C Look into row IY:

        DO ILD=KLDA1(IY)+1,KLDA1(IY+1)-1
          ICOL=KCOLA1(ILD)
C If there's column IX in this row we can stop here
          IF (ICOL.EQ.IX) GOTO 1
        END DO
        
C Otherwise: error - this element does not exist in our matrix
        
        WRITE(MTERM,'(A,I6,I6,A))') 'Error in IGXYIX: Entry (',IX,',',
     *                              IY,') not in matrix !'
        
1       IGXYIX = ILD
        
      END 



************************************************************************
*                                                                      
* MAREST                                                               
*                                                                      
* Purpose: Restricrion of the fine grid matrix to the coarse grid
*          matrix. Searches in the matrix for rows that belong to
*          elements with an aspect ratio that is too large. These
*          rows are rebuild by constant restriction of the fine grid
*          matrix to the coarse grid matrix.
*
* Important Parameters:
* 
*  IADM1 - configuration parameter for switching to constant 
*          restriction:
*          =1: switch depending on aspect ratio of current element
*          =2: switch depending on aspect ratio of current element
*              and aspect ratio of neighbour element
*  PARM  - maximum aspect ratio. Rows corresponding to elements
*          with AR>PARM (and neighbour elements with AR>PARM resp.)
*          are rebuild by restriction
*
************************************************************************

      SUBROUTINE MAREST(KVERT1,KVERT2,KMID1,KMID2,KADJ1,KADJ2,
     *                  KLDA1,KLDA2,KCOLA1,KCOLA2,VA1,VA2,DCORV1,AREA1,
     *                  NEL1,NEL2,NVT1,NVT2,PARM,IADM)
************************************************************************
C
      IMPLICIT NONE

      INTEGER NNVE
      PARAMETER (NNVE=4)

      DOUBLE PRECISION PARM
      INTEGER NEL1,NEL2,NVT1,NVT2,IADM
      INTEGER KVERT1(NNVE,*),KVERT2(NNVE,*),KMID1(NNVE,*),
     *        KMID2(NNVE,*),KADJ1(NNVE,*),KADJ2(NNVE,*)
      INTEGER KLDA1(*),KLDA2(*),KCOLA1(*),KCOLA2(*)
      DOUBLE PRECISION DCORV1(2,*)
      DOUBLE PRECISION VA1(*),VA2(*)
      DOUBLE PRECISION AREA1(*)
      INTEGER KIEL(NNVE)
      
      INTEGER IEL, IVE, IVE1, IVE2, IVE3, IVE4, ILD
      INTEGER JVE, JVE1, JVE2, JVE3, JVE4
      INTEGER IMID1, IMID2, IMID3, IMID4, IMID5 
      INTEGER IADJ1, IADJ2, IADJ3, IADJ4
      INTEGER JADJ1, JADJ2, JADJ3, JADJ4
      INTEGER IM1, IM2, IM3, IM4, IM5, IM6, IM7, IM8, IM9,
     *        IM10, IM11, IM12
      INTEGER IEL1, IEL2, JEL1, JEL2
      DOUBLE PRECISION DVAL1, DVAL2, DVAL3, DVAL4, DVAL5 
      DOUBLE PRECISION PV1, PV2, PV3, PV4, PV5, PV6, PV7, PV8,
     *                 PV9, PV10, PV11, PV12
      DOUBLE PRECISION ARIEL, ARIELA
      INTEGER IGXYIX
      EXTERNAL IGXYIX

C Alle Elemente IEL=1..NEL1 des Grobgitters durchlaufen:

      DO 10 IEL=1,NEL1
C
C      BPRINT=((iel.eq.9).and.(ilev.eq.1))
c

C Verzerrungsfaktor des aktuellen Elementes mit ARCALC berechnen.
C Ergebnis wird in ARIEL abgelegt.

      CALL ARCALC(ARIEL,IEL,KVERT1,DCORV1)

C Evtl. Kehrwert bilden, da wir an dem Verhaeltnis
C "laengere/kuerzere" Seite des Elementes intreressiert sind

      IF (ARIEL.LT.1D0) ARIEL=1D0/ARIEL
      
C Die Elemente auf dem Feingitter bestimmen, die in dem aktuellen
C Element auf dem Grobgitter enthalten sind

      KIEL(1)=IEL
      KIEL(2)=KADJ2(2,KIEL(1))
      KIEL(3)=KADJ2(2,KIEL(2))
      KIEL(4)=KADJ2(2,KIEL(3))

C Zur Zeit liegt damit folgende Variablenbelegung vor:
C                 
C                                    
C     4==============X===============3
C     |              |               |
C     |              |               |
C     |    KIEL(4)   |     KIEL(3)   |
C     |              |               |
C     |                              |
C     X----------- IEL1 -------------X 
C     |                              |
C     |              |               |
C     |    KIEL(1)   |     KIEL(2)   |
C     |              |               |
C     |              |               |
C     1==============X===============2
C                       
C
C mit Verzerrungsfaktor ARIEL.

C Nun werden die vier Knoten des Elementes durchlaufen:
C IVE=1..4:

      DO 20 IVE=1,4

C IVE1 erhaelt die Nummer des aktuellen Knotens (z.B. 1)
C auf der Seitenmitte, bzw. Nummer der aktuell betrachteten 
C Kante (in einer lokalen Nummerierung 1..4)
C IVE2 erhaelt die Nummer des rechten Nachbarknotens/
C rechte Nachbarkante auf dem Element (im Beispiel: 2)
C IVE4 erhaelt die Nummer des linken Nachbarknotens/
C rechte Nachbarkante auf dem Element (im Beispiel: 4)

      IVE1=IVE
      IVE2=MOD(IVE1,4)+1
      IVE4=MOD(IVE1+2,4)+1

C Also mit IVE=1:
C
C IVE4 O----------- IEL1 -------------O IVE2=2
C  =4  |                              |
C      |              |               |
C      |    KIEL(1)   |     KIEL(2)   |
C      |              |               |
C      |              |               |
C      ===============O================  
C                   IVE1=1
C
C Nun bekommen IADJ1, IADJ2 und IADJ4 die Nummern der
C Nachbarelemente zugewiesen, die den Knoten IVEx
C mit unserem aktuellen Element gemeinsam haben:
C
      IADJ1=KADJ1(IVE1,IEL)
      IADJ2=KADJ1(IVE2,IEL)
      IADJ4=KADJ1(IVE4,IEL)

C Also:
C
C  +-------+-------+-------+
C  |       |       |       |
C  | IADJ4 4  IEL  2 IADJ2 |
C  |       |       |       |
C  +-------+---1---+-------+
C          |       |
C          | IADJ1 |
C          |       |
C          +-------+


C IEL1 und IEL2 erhalten nun die Nummern der beiden
C Elemente im Feingitter, die den aktuellen Seitenmittelpunkt
C IVE=IVE1 gemeinsam haben:

      IEL1=KIEL(IVE1)
      IEL2=KIEL(IVE2)
      
C Also:
C
C      O------------ IEL ------------O     
C      |                             |
C      |              |              |
C      |     IEL1     |     IEL2     |
C      |              |              |
C      |              |              |
C      1==============O==============2
C                    IVE1

C Noch ein paar Variablen fuer die folgenden Berechnungen 
C initialisieren:

      IMID1=0
      IMID2=0
      IMID3=0
      IMID4=0
      IMID5=0

      DVAL1=0D0
      DVAL2=0D0
      DVAL3=0D0
      DVAL4=0D0
      DVAL5=0D0
C
C Wir betrachten das Nachbarelement ADJ1 unseres aktuellen Knotens
C IEL1: 
C  +-------+-------+-------+
C  |       |       |       |
C  |       4 IEL1  2       |
C  |       |       |       |
C  +-------+---1---+-------+
C          |       |
C          | IADJ1 |
C          |       |
C          +-------+
C Falls dieses Element existiert und eine kleinere Nummer hat,
C als unser aktuelles Element, so haben wir die Matrix an der jeweiligen
C Stelle schon modifiziert (da die Matrix ja die Ueberlappungen der
C Basisfunktionen der beiden Elemente beschreibt). Daher brauchen
C wir die Matrix hier nicht nochmal zu modifizieren und koennen
C die Modifikationsroutine ueberspringen:

      IF ((IADJ1.LT.IEL).AND.(IADJ1.NE.0)) GOTO 20

C Falls das Nachbarelement auf dem Grobgitter existiert und eine 
C groessere Nummer hat, als unser aktuelles Element, so pruefen wir 
C den Verzerrungsfaktor - sowohl den des aktuellen Elementes wie
C auch den des Nachbarelementes:

      IF (IADJ1.NE.0) THEN
C Verzerrungsfaktor berechnen, bei Bedarf Kehrwert:
         CALL ARCALC(ARIELA,IADJ1,KVERT1,DCORV1)
         IF (ARIELA.LT.1D0) ARIELA=1D0/ARIELA

C Falls sowohl das aktuelle Element (bzw. bei IADM1 auch das 
C Nachbarelement)
C keine grosse Anisotropie zeigen (Verzerrungsfaktor <= PARM),
C so brauchen wir die Matrix nicht zu modifizieren, d.h.
C in diesem Fall ueberspringen wir die Modifikationsroutine:

        IF (IADM.LE.0) GOTO 20

C Prüfen, ob eines der beiden Elemente ein schlechtes AR hat

        IF ((ARIEL.LT.PARM).AND.(ARIELA.LT.PARM)) GOTO 20      
        
C evtl. prüfen, ob es das aktuelle Element ist und sein soll

        IF ((IADM.EQ.1).AND.(ARIEL.LT.PARM)) GOTO 20

      ELSE
      
C Wenn kein Nachbarelement vorhanden ist, so entscheidet
C der Verzerrungsfaktor des aktuellen Elementes, ob 
C die Matrix modifiziert wird:

        IF (IADM.LE.0) GOTO 20

        IF (ARIEL.LT.PARM) GOTO 20

      ENDIF

C Wir sind in dem Fall, dass wir die Matrix modifizieren muessen
C aufgrund zu grosser Anisotropie des aktuellen oder des Nachbar-
C elementes.
C
C Zuerst werden die globalen Nummern der Seitenmittelpunkte IMID1..3
C zu den lokalen Nummeern IVE1,2,4 auf dem aktuellen Element IEL1 
C des Grobgitters ermittelt.
C Davon NVT1 abziehen, damit der global erste Seitenmittelpunkt die
C Nummer 1 bekommt. Die Nummer korrespondiert mit der Nummer
C des Freiheitsgrades im Loesungsvektor.

      IMID1 =KMID1(IVE1,IEL)-NVT1
      IMID2 =KMID1(IVE2,IEL)-NVT1
      IMID3 =KMID1(IVE4,IEL)-NVT1

C Also:
C
C IMID3 O------------ IEL ------------O IMID2=IVE2
C =IVE3 |                             |
C       |              |              |
C       |     IEL1     |     IEL2     |
C       |              |              |
C       |              |              |
C       1==============O==============2
C                 IMID1=IVE1

C Nun bestimmen wir die globalen Nummern der Seitenmittelpunkte 
C auf dem Feingitter (um -NVT korrigiert, wie oben):

      IM1 =KMID2(1,IEL1)-NVT2
      IM2 =KMID2(4,IEL2)-NVT2
      IM3 =KMID2(2,IEL1)-NVT2
      IM4 =KMID2(4,IEL1)-NVT2
      IM5 =KMID2(1,IEL2)-NVT2
      IM6 =KMID2(3,IEL1)-NVT2
      IM7 =KMID2(2,IEL2)-NVT2
      
C Also:
C
C             IM6            IM7 
C       +------O----- IEL ----O-------+
C       |4            3|3            2|
C       |     IEL1     |     IEL2     |
C   IM4 O              O IM3          O IM5
C       |              |              |
C       |1            2|4            1|
C       1=======O======+=======O======2
C              IM1            IM2      

C Falls es nun ein unteres Nachbarelement gibt, so brauchen
C wir auch die Knotennummern der Knoten des Nachbarelementes
C auf dem Fein- und Grobgitter:

      IF (IADJ1.NE.0) THEN
      
C Zuerst mal eben die 4 Nachbarelementes durchlaufen, um unser
C aktuelles Element wiederzufinden:

        DO JVE=1,4
          IF (KADJ1(JVE,IADJ1).EQ.IEL) GOTO 32
        END DO  

C           4--------3
C           |        |
C           |  IEL   |
C           |   ^    |
C  +--------1---^----2--------+
C  |        |   ^    |        |
C  |   ?    | IADJ1  |   ?    |
C  |        |        |        |
C  +--------+--------+--------+
C           |        |
C           |   ?    |
C           |        |
C           +--------+

C Danach enthaelt JVE die lokale Nummer (1..4) der Kante, welche
C von IADJ1 aus an unser aktuelles Element IEL stoesst. Nun suchen wir
C fuer IADJ1 genauso wie fuer IEL die ganzen lokalen und globalen
C Knotennummern raus:

32      JVE1=JVE
        JVE2=MOD(JVE1,4)+1
        JVE4=MOD(JVE1+2,4)+1

C In unserem Beispiel waere JVE=JVE1=3!
C Dazu kommt der rechte Nachbar JVE2 und der linke Nachbar JVE4:
C Dies bedeutet:
C
C IVE4 O----         IEL1         ----O IVE2  
C      |                              |
C      |                              |
C      |                              |
C      |              |               |
C      |              |               |
C      4==============O===============3  
C      |         JVE=JVE1=3           |
C      |                              |
C      |                              |
C      |                              |
C      |                              |
C JVE2 O                              O JVE4=2
C =4   |                              |
C      |                              |
C      |                              |
C      |                              |
C      |                              |
C      1==============O===============2  
C
C Nummern der an IADJ1 angrenzenden Elemente auf dem  Grobgitter
C ermitteln:

        JADJ2=KADJ1(JVE2,IADJ1)
        JADJ4=KADJ1(JVE4,IADJ1)

C Also:
C
C           +--------+
C           |   |    |
C           |- IEL --|
C           |   |    |
C  +--------+--------+--------+
C  |        |        |        |
C  | JADJ2  | IADJ1  | JADJ4  |
C  |        |        |        |
C  +--------+--------+--------+
C           |        |
C           |   ?    |
C           |        |
C           +--------+

C Nun auf das Feingitter gehen. Dort betrachten wir die Elemente
C IEL1 und IEL2 aus dem Grobgitterelement IEL. Wir schauen nun
C nach, welche Feingitterelemente aus dem Grobgitterelement
C IADJ1 an diese beiden Elemente angrenzen:

        JEL1=KADJ2(4,IEL2)
        JEL2=KADJ2(1,IEL1)

C Also:
C
C      O------------ IEL -------------O       
C      |                              |
C      |              |               |
C      |     IEL1     |     IEL2      |
C      |      v       |      v        |
C      |1     v       |      v       2|
C      +======v=======O======v========+  
C      |4     v       |      v       3|
C      |      v       |      v        |
C      |     JEL2     |     JEL1      |
C      |              |               |
C      |                              |
C      O----------- IADJ1 ------------O       
C      |                              |
C      |              |               |
C      |              |               |
C      |              |               |
C      |1             |              2|
C      +==============O===============+  C

C Nun ermitteln wir noch die globalen Nummern der Seitenmittelpunkte
C des linken und rechten Nachbars von JVE1 auf dem Grobgitter-
C element IADJ1 (um -NVT1 korrigiert, damit der global erste Seiten-
C mittelpunkt Nummer 1 hat):
 
        IMID4 =KMID1(JVE2,IADJ1)-NVT1
        IMID5 =KMID1(JVE4,IADJ1)-NVT1
        
C Also:
C
C       O----         IEL          ----O       
C       |                              |
C       |                              |
C       |                              |
C       |              |               |
C       |              |               |
C       4==============O===============3  
C       |         JVE=JVE1=3           |
C       |                              |
C       |                              |
C       |                              |
C IMID4 |                              | IMID5
C =JVE2 O            IADJ1             O =JVE4
C       |                              |
C       |                              |
C       |                              |
C       |                              |
C       |                              |
C       1==============O===============2  

C Schliesslich ermitteln wir noch die globalen Nummern der 
C Seitenmitten von JEL1 und JEL2 auf dem Feingitter:
 
        IM8 =KMID2(2,JEL1)-NVT2
        IM9 =KMID2(4,JEL1)-NVT2
        IM10=KMID2(1,JEL2)-NVT2
        IM11=KMID2(3,JEL1)-NVT2
        IM12=KMID2(2,JEL2)-NVT2

C Also:
C
C       O----         IEL          ----O       
C       |                              |
C       |                              |
C       |                              |
C       |              |               |
C       |              |               |
C       4==============+===============3  
C       |1            4|2             1|
C       |     JEL2     |     JEL1      |
C  IM10 O              O IM8           O IM9
C       |              |               |
C       |2            3|3             4|      
C       +------O---- IADJ1 ----O-------+      
C       |4   IM12     3|3    IM11     2|
C       |              |               |
C       |              |               |
C       |              |               |
C       |1            2|4             1|
C       1==============+===============2  


      ELSE

C Falls es kein Nachbarelement gibt, so setzen wir als Zeichen
C dafuer IM8 auf 0

        IM8 =0
      ENDIF

C Die IMx-Variablen enthalten somit nun zusammenfassend die 
C globalen Nummern (um -NVT korrigiert) der folgenden Seitenmitten
C auf dem Feingitter:
C
C       |     IM6             IM7      |
C       +------O----- IEL ----O--------+
C       |                              |
C       |              |               |
C   IM4 O              O IM3           O IM5
C       |              |               |
C       |              |               |
C       +=======O======+=======O=======+
C       |      IM1     |      IM2      |
C       |              |               |
C  IM10 O              O IM8           O IM9
C       |              |               |
C       |                              |      
C       +------O---- IADJ1 ----O-------+      
C       |    IM12            IM11      |


C Nun benutzen wir die Werte der Feingittermatrix an den Stellen IMx,
C um per Galerkin-Zugang die Eintraege auf der Grobgittermatrxix VA1
C an den Stellen IMID1,...,IMID5 zu berechnen.
C Dies entspricht der Anwendung eines Progongations und eines
C Restriktionsoperators auf die Feingittermatrix VA2.

      ILD = IGXYIX (IM1,IM3,KCOLA2,KLDA2)
      PV1=VA2(KLDA2(IM1))+VA2(ILD)
 
      IF (IM8.NE.0) THEN
        ILD = IGXYIX (IM1,IM8,KCOLA2,KLDA2)
        PV1=PV1+VA2(ILD)
      END IF

      ILD = IGXYIX (IM2,IM3,KCOLA2,KLDA2)
      PV2=VA2(KLDA2(IM2))+VA2(ILD)

      IF (IM8.NE.0) THEN
        ILD = IGXYIX (IM2,IM8,KCOLA2,KLDA2)
        PV2=PV2+VA2(ILD)
      END IF

      ILD = IGXYIX (IM3,IM1,KCOLA2,KLDA2)
      PV3=VA2(KLDA2(IM3))+VA2(ILD)

      ILD = IGXYIX (IM3,IM2,KCOLA2,KLDA2)
      PV3=PV3+VA2(ILD)

      ILD = IGXYIX (IM4,IM1,KCOLA2,KLDA2)
      PV4=VA2(ILD)

      ILD = IGXYIX (IM4,IM3,KCOLA2,KLDA2)
      PV4=PV4+VA2(ILD)

      ILD = IGXYIX (IM5,IM2,KCOLA2,KLDA2)
      PV5=VA2(ILD)

      ILD = IGXYIX (IM5,IM3,KCOLA2,KLDA2)
      PV5=PV5+VA2(ILD)

      ILD = IGXYIX (IM6,IM1,KCOLA2,KLDA2)
      PV6=VA2(ILD)

      ILD = IGXYIX (IM6,IM3,KCOLA2,KLDA2)
      PV6=PV6+VA2(ILD)

      ILD = IGXYIX (IM7,IM2,KCOLA2,KLDA2)
      PV7=VA2(ILD)

      ILD = IGXYIX (IM7,IM3,KCOLA2,KLDA2)
      PV7=PV7+VA2(ILD)

      IF (IM8.NE.0) THEN
        ILD = IGXYIX (IM8,IM1,KCOLA2,KLDA2)
        PV8=VA2(KLDA2(IM8))+VA2(ILD)

        ILD = IGXYIX (IM8,IM2,KCOLA2,KLDA2)
        PV8=PV8+VA2(ILD)

        ILD = IGXYIX (IM9,IM8,KCOLA2,KLDA2)
        PV9=VA2(ILD)

        ILD = IGXYIX (IM9,IM2,KCOLA2,KLDA2)
        PV9=PV9+VA2(ILD)

        ILD = IGXYIX (IM10,IM8,KCOLA2,KLDA2)
        PV10=VA2(ILD)

        ILD = IGXYIX (IM10,IM1,KCOLA2,KLDA2)
        PV10=PV10+VA2(ILD)

        ILD = IGXYIX (IM11,IM8,KCOLA2,KLDA2)
        PV11=VA2(ILD)

        ILD = IGXYIX (IM11,IM2,KCOLA2,KLDA2)
        PV11=PV11+VA2(ILD)

        ILD = IGXYIX (IM12,IM8,KCOLA2,KLDA2)
        PV12=VA2(ILD)

        ILD = IGXYIX (IM12,IM1,KCOLA2,KLDA2)
        PV12=PV12+VA2(ILD)
      ENDIF

      IF (IADJ1.EQ.0) THEN
        DVAL1=1D0/1D0*(PV1+PV2+PV3)
      ELSE
        DVAL1=1D0/1D0*(PV1+PV2+PV3+PV8)
      ENDIF

      IF (IADJ2.EQ.0) THEN
        DVAL2=1D0/1D0*(PV5+PV7)
      ELSE
        DVAL2=1D0/1D0*(PV5+PV7)
      ENDIF

      IF (IADJ4.EQ.0) THEN
        DVAL3=1D0/1D0*(PV4+PV6)
      ELSE
        DVAL3=1D0/1D0*(PV4+PV6)
      ENDIF

      IF (IADJ1.NE.0) THEN

       IF (JADJ2.EQ.0) THEN
         DVAL4=1D0/1D0*(PV10+PV12)
       ELSE
         DVAL4=1D0/1D0*(PV10+PV12)
       ENDIF

       IF (JADJ4.EQ.0) THEN
         DVAL5=1D0/1D0*(PV9+PV11)
       ELSE
         DVAL5=1D0/1D0*(PV9+PV11)
       ENDIF

      ENDIF

C Nun werden die eigentlichen Matrixeintraege fuer VA1 berechnet.
C Zuerst muss die zu IMID1 gehoerende Zeile in der Matrix geloescht
C werden:

      DO ILD=KLDA1(IMID1),KLDA1(IMID1+1)-1
        VA1(ILD)=0E0
      END DO  

C Dann werden die Eintraege in der betroffenen Zeile geaendert:

      VA1(KLDA1(IMID1))=DVAL1

      ILD = IGXYIX (IMID2,IMID1,KCOLA1,KLDA1)
      VA1(ILD)=DVAL2
    
      ILD = IGXYIX (IMID3,IMID1,KCOLA1,KLDA1)
      VA1(ILD)=DVAL3

      IF (IADJ1.NE.0) THEN
        ILD = IGXYIX (IMID4,IMID1,KCOLA1,KLDA1)
        VA1(ILD)=DVAL4

        ILD = IGXYIX (IMID5,IMID1,KCOLA1,KLDA1)
        VA1(ILD)=DVAL5
      ENDIF

C Aktueller Seitenmittelpunkt IMID1 abgehandelt. Wir schalten nun
C (im Beispiel gegen den Uhrzeigersinn) auf den naechsten 
C Seitenmittelpunkt um und machen dort genauso weiter

20    CONTINUE

C Aktuelles Element fertig. Auf naechstes Element schalten

10    CONTINUE

9999  END
      
