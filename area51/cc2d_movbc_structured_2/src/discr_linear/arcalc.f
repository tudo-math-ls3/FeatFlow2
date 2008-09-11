*+**********************************************************************
* This file contains code to calculate the aspect ratio of the elements
* in the grid:
*
* WHCALC - Calculates width and height
* ARCALC - Calculates the aspect ratios
* SETARE - Calculates the area of all elements
************************************************************************

      SUBROUTINE WHCALC (DM1M3,DM2M4,IELM,KVERT,DCORVG)
C
C returns:
C   DM1M3 - distance between midpoints M1, M3
C   DM2M4 - distance between midpoints M2, M4
C 
      IMPLICIT NONE

      DOUBLE PRECISION DCORVG(2,*)
      INTEGER KVERT(4,*)
      
      DOUBLE PRECISION DM1M3,DM2M4
      INTEGER IELM
      
      DOUBLE PRECISION XE1, XE2, XE3, XE4, XM1, XM2, XM3, XM4
      DOUBLE PRECISION YE1, YE2, YE3, YE4, YM1, YM2, YM3, YM4
      INTEGER IECK1, IECK2, IECK3, IECK4
     
      IECK1=KVERT(1,IELM)    
      IECK2=KVERT(2,IELM)    
      IECK3=KVERT(3,IELM)    
      IECK4=KVERT(4,IELM)    

      XE1=DCORVG(1,IECK1)
      YE1=DCORVG(2,IECK1)
      XE2=DCORVG(1,IECK2)
      YE2=DCORVG(2,IECK2)
      XE3=DCORVG(1,IECK3)
      YE3=DCORVG(2,IECK3)
      XE4=DCORVG(1,IECK4)
      YE4=DCORVG(2,IECK4)

      XM1=(XE2+XE1)/2D0
      YM1=(YE2+YE1)/2D0
      XM2=(XE3+XE2)/2D0
      YM2=(YE3+YE2)/2D0
      XM3=(XE4+XE3)/2D0
      YM3=(YE4+YE3)/2D0
      XM4=(XE1+XE4)/2D0
      YM4=(YE1+YE4)/2D0

      DM1M3 = DSQRT((XM3-XM1)**2+(YM3-YM1)**2)
      DM2M4 = DSQRT((XM4-XM2)**2+(YM4-YM2)**2)

      END


      SUBROUTINE ARCALC (RATIO,IELM,KVERT,DCORVG)
****************************************************************
c    by Rainer Schmachtel
C=======================================================================
C    Calculates the aspect ratio (AR) of element IELM,
C    and returns it on RATIO.
C=======================================================================
C-----------------------------------------------------------------------
      IMPLICIT NONE

      DOUBLE PRECISION RATIO, DCORVG(2,*)
      INTEGER KVERT(4,*)

      INTEGER IELM
      
      DOUBLE PRECISION NENNER, ZAEHLER
      DOUBLE PRECISION XE1, XE2, XE3, XE4, XM1, XM2, XM3, XM4
      DOUBLE PRECISION YE1, YE2, YE3, YE4, YM1, YM2, YM3, YM4
      INTEGER IECK1, IECK2, IECK3, IECK4

      IECK1=KVERT(1,IELM)    
      IECK2=KVERT(2,IELM)    
      IECK3=KVERT(3,IELM)    
      IECK4=KVERT(4,IELM)    

      XE1=DCORVG(1,IECK1)
      YE1=DCORVG(2,IECK1)
      XE2=DCORVG(1,IECK2)
      YE2=DCORVG(2,IECK2)
      XE3=DCORVG(1,IECK3)
      YE3=DCORVG(2,IECK3)
      XE4=DCORVG(1,IECK4)
      YE4=DCORVG(2,IECK4)

      XM1=(XE2+XE1)/2D0
      YM1=(YE2+YE1)/2D0
      XM2=(XE3+XE2)/2D0
      YM2=(YE3+YE2)/2D0
      XM3=(XE4+XE3)/2D0
      YM3=(YE4+YE3)/2D0
      XM4=(XE1+XE4)/2D0
      YM4=(YE1+YE4)/2D0

      ZAEHLER=(XM3-XM1)**2+(YM3-YM1)**2
      NENNER =(XM4-XM2)**2+(YM4-YM2)**2

      RATIO=DSQRT(ZAEHLER/MAX(NENNER,1D-10))
C      RATIO=DSQRT(ZAEHLER/NENNER)

      END

************************************************************************
*
* SETARE
*
* Flaecheninhalt aller Elemente bestimmen.
*
*   Purpose: - writes on  AREA(IEL)  the area of the element IEL,
*              IEL=1,...,NEL
*            - writes on  AREA(NEL+1) the sum of all  AREA(IEL)
*            - KVERT,DCORVG are the usual FEAT arrays
*
************************************************************************
C      SUBROUTINE   SETARE  (AREA,NEL,KVERT,DCORVG)
C=======================================================================
C     Declarations
C=======================================================================
C      IMPLICIT NONE

C      INTEGER NNVE
C      PARAMETER (NNVE=4)
C      DOUBLE PRECISION DCORVG(2,*)
C      INTEGER NEL,KVERT(NNVE,*)

C      DOUBLE PRECISION AREA(*)
      
C      DOUBLE PRECISION SUM,X1,X2,X3,X4,Y1,Y2,Y3,Y4,AAA
C      INTEGER IEL,I1,I2,I3,I4
C=======================================================================
C      SUM=0.D0
C      DO  11  IEL=1,NEL
C
C      I1=KVERT(1,IEL)
C      I2=KVERT(2,IEL)
C      I3=KVERT(3,IEL)
C      I4=KVERT(4,IEL)
C
C      X1=DCORVG(1,I1)
C      X2=DCORVG(1,I2)
C      X3=DCORVG(1,I3)
C      X4=DCORVG(1,I4)
C
C      Y1=DCORVG(2,I1)
C      Y2=DCORVG(2,I2)
C      Y3=DCORVG(2,I3)
C      Y4=DCORVG(2,I4)
C
C Berechnung des Flaecheninhaltes des Rechtecks.
C Dazu wird der Flaecheninhalt der beiden Dreiecke
C (X1,Y1)->(X2,Y2)->(X3,Y3) und (X1,Y1)->(X3,Y3)->(X4,Y4)
C berechnet, indem die Determinante der Matrix aus den aufspannenden
C Vektoren berechnet und halbiert wird.:
C det ( X1-X2 X3-X2 ) + det (X1-X4 X3-X4 )
C       Y1-Y2 Y1-Y3          Y1-Y4 Y3-Y4
C
C      AAA=0.5D0*(  DABS((X1-X2)*(Y3-Y2)-(Y1-Y2)*(X3-X2))
C     *            +DABS((X1-X4)*(Y3-Y4)-(Y1-Y4)*(X3-X4)) )
C      AREA(IEL)=AAA
C      SUM=SUM+AAA
C  11  CONTINUE
C
C      AREA(NEL+1)=SUM
C
C      END
