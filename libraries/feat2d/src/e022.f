************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.3)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek                     *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* E022                                                                 *
*                                                                      *
* Purpose  Calculation of values and derivatives of element  22        *
*          Piecewice linear element (4P1/P1 element)                   *
*          (same #dof as quadratic element)                            *
*          Numbering convention:                                       *
*                                  3                                   *
*                                  +                                   *
*                                  *E4*                                *
*                                6 +*****+ 5                           *
*                                  *E2*E1*E3*                          *
*                                  +*****+*****+                       *
*                                  1     4     2                       *
*                                                                      *
* Caution: To evaluate function values use piecewise 1x1 Gauss cubature*
*          formula or piecewise trapezoidal rule,                      *
*          to evaluate derivatives use piecewise 1x1 Gauss cubature    *
*          formula in subroutine CB2T only                             *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  03/21/89                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* XI1      R*8                                                         *
* XI2      R*8    Evaluation point in barycentric coordinates          *
* XI3      R*8                                                         *
* IPAR     I*4    Switch for desired action                            *
*                  0  Calculation of the desired values of DBAS        *
*                 -1  Set number of element                            *
*                 -2  Calculate values on the reference element        *
*                     for all cubature points and save them            *
*                 -3  same as 0, but use the values saved before       *
*                                                                      *
* BDER     LOG    Derivative J is calculated if BDER(J).EQ..TRUE.      *
*                 Multiindices are enumbered from 1 to 6               *
* DJAC     R*8    Jacobian                                             *
* DETJ     R*8    Determinant of the jacobian                          *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DBAS     R*8    Array of function values and derivatives for         *
*                 all basis functions (block /ELEM/)                   *
* IPAR     I*4    Set to  22  if called with IPAR = -1                 *
* IER      I*4    Error indicator                                      *
*                 -131  Desired derivative not available               *
*                 -132  Triangle has vanishing area                    *
*                 -130  Clockwise ordering of the corner points        *
*                                                                      *
************************************************************************
C
      SUBROUTINE E022 (XI1,XI2,XI3,IPAR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NDOFL=6)
      PARAMETER (Q2=.5D0,Q3=1D0/3D0)
      DIMENSION KCASE(NNCUBP)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      SAVE /OUTPUT/,/ERRCTL/,/ELEM/,/CHAR/,/CUB/
      SAVE KCASE
C
      SUB='E022'
      IF (ICHECK.GE.998) CALL OTRC('E022  ','03/21/89')
C
      IER=0
      GOTO (10,1,2,3) -IPAR+1
C
10    ICUBP0=1
      CALL E022A(XI1,XI2,XI3,KCASE,ICUBP0)
      GOTO 2200
C
1     IPAR=22
      GOTO 99999
C
2     IF (ICUBP.EQ.12.OR.ICUBP.EQ.13) THEN
C *** piecewise trapezoidal or piecewise midpoint rule
       DO 21 ICUBP1=1,4
       IAUX=(ICUBP1-1)*NCUBP/4
       DO 21 ICUBP0=IAUX+1,IAUX+NCUBP/4
21     KCASE(ICUBP0)=ICUBP1
      ELSE
C *** Other integration rules not defined piecewise
C *** or with interior cubature points only
       DO 22 ICUBP0=1,NCUBP
22     CALL E022A (DXI(ICUBP0,1),DXI(ICUBP0,2),DXI(ICUBP0,3),
     *             KCASE,ICUBP0)
      ENDIF
      GOTO 99999
C
3     ICUBP0=ICUBP
C
2200  IF (ICHECK.EQ.0) GOTO 2201
      IF (DETJ.LE.0D0) CALL WERR(-130,'E022  ')
      IF (DETJ.GT.0D0.AND.DETJ.LT.1D-70) CALL WERR(-132,'E022  ')
C *** No second order derivatives available
C *** Used for second order problems only
      IF (BDER(4).OR.BDER(5).OR.BDER(6)) CALL WERR(-131,'E022  ')
      IF (IER.NE.0) GOTO 99999
C
2201  XJ2=2D0/DETJ
      XJ22=XJ2*Q2
      XJ23=XJ2*Q3
C
      GOTO (2211,2212,2213,2214,2215,2216,2217,2218,2219,2220)
     *      KCASE(ICUBP0)
C
2211  IF (BDER(1)) THEN
       DBAS(1,1)=0D0
       DBAS(2,1)=0D0
       DBAS(3,1)=0D0
       DBAS(4,1)= XI1+XI2-XI3
       DBAS(5,1)=-XI1+XI2+XI3
       DBAS(6,1)= XI1-XI2+XI3
      ENDIF
      IF (BDER(2)) THEN
       DBAS(1,2)=0D0
       DBAS(2,2)=0D0
       DBAS(3,2)=0D0
       DBAS(4,2)= DJAC(2,1)*XJ2
       DBAS(5,2)= (DJAC(2,2)-DJAC(2,1))*XJ2
       DBAS(6,2)=-DJAC(2,2)*XJ2
      ENDIF
      IF (BDER(3)) THEN
       DBAS(1,3)=0D0
       DBAS(2,3)=0D0
       DBAS(3,3)=0D0
       DBAS(4,3)=-DJAC(1,1)*XJ2
       DBAS(5,3)=-(DJAC(1,2)-DJAC(1,1))*XJ2
       DBAS(6,3)= DJAC(1,2)*XJ2
      ENDIF
      GOTO 99999
C
2212  IF (BDER(1)) THEN
       DBAS(1,1)= XI1-XI2-XI3
       DBAS(2,1)=0D0
       DBAS(3,1)=0D0
       DBAS(4,1)= 2D0*XI2
       DBAS(5,1)=0D0
       DBAS(6,1)= 2D0*XI3
      ENDIF
      IF (BDER(2)) THEN
       DBAS(1,2)=-(DJAC(2,2)-DJAC(2,1))*XJ2
       DBAS(2,2)=0D0
       DBAS(3,2)=0D0
       DBAS(4,2)= DJAC(2,2)*XJ2
       DBAS(5,2)=0D0
       DBAS(6,2)=-DJAC(2,1)*XJ2
      ENDIF
      IF (BDER(3)) THEN
       DBAS(1,3)= (DJAC(1,2)-DJAC(1,1))*XJ2
       DBAS(2,3)=0D0
       DBAS(3,3)=0D0
       DBAS(4,3)=-DJAC(1,2)*XJ2
       DBAS(5,3)=0D0
       DBAS(6,3)= DJAC(1,1)*XJ2
      ENDIF
      GOTO 99999
C
2213  IF (BDER(1)) THEN
       DBAS(1,1)=0D0
       DBAS(2,1)=-XI1+XI2-XI3
       DBAS(3,1)=0D0
       DBAS(4,1)= 2D0*XI1
       DBAS(5,1)= 2D0*XI3
       DBAS(6,1)=0D0
      ENDIF
      IF (BDER(2)) THEN
       DBAS(1,2)=0D0
       DBAS(2,2)= DJAC(2,2)*XJ2
       DBAS(3,2)=0D0
       DBAS(4,2)=-(DJAC(2,2)-DJAC(2,1))*XJ2
       DBAS(5,2)=-DJAC(2,1)*XJ2
       DBAS(6,2)=0D0
      ENDIF
      IF (BDER(3)) THEN
       DBAS(1,3)=0D0
       DBAS(2,3)=-DJAC(1,2)*XJ2
       DBAS(3,3)=0D0
       DBAS(4,3)= (DJAC(1,2)-DJAC(1,1))*XJ2
       DBAS(5,3)= DJAC(1,1)*XJ2
       DBAS(6,3)=0D0
      ENDIF
      GOTO 99999
C
2214  IF (BDER(1)) THEN
       DBAS(1,1)=0D0
       DBAS(2,1)=0D0
       DBAS(3,1)=-XI1-XI2+XI3
       DBAS(4,1)=0D0
       DBAS(5,1)= 2D0*XI2
       DBAS(6,1)= 2D0*XI1
      ENDIF
      IF (BDER(2)) THEN
       DBAS(1,2)=0D0
       DBAS(2,2)=0D0
       DBAS(3,2)=-DJAC(2,1)*XJ2
       DBAS(4,2)=0D0
       DBAS(5,2)= DJAC(2,2)*XJ2
       DBAS(6,2)=-(DJAC(2,2)-DJAC(2,1))*XJ2
      ENDIF
      IF (BDER(3)) THEN
       DBAS(1,3)=0D0
       DBAS(2,3)=0D0
       DBAS(3,3)= DJAC(1,1)*XJ2
       DBAS(4,3)=0D0
       DBAS(5,3)=-DJAC(1,2)*XJ2
       DBAS(6,3)= (DJAC(1,2)-DJAC(1,1))*XJ2
      ENDIF
      GOTO 99999

2215  IF (BDER(1)) THEN
       DBAS(1,1)=0D0
       DBAS(2,1)=0D0
       DBAS(3,1)=0D0
       DBAS(4,1)=1D0
       DBAS(5,1)=0D0
       DBAS(6,1)=0D0
      ENDIF
      IF (BDER(2)) THEN
       DBAS(1,2)=-(DJAC(2,2)-DJAC(2,1))*XJ23
       DBAS(2,2)= DJAC(2,2)*XJ23
       DBAS(3,2)=0D0
       DBAS(4,2)= 2D0*DJAC(2,1)*XJ23
       DBAS(5,2)= (DJAC(2,2)-2D0*DJAC(2,1))*XJ23
       DBAS(6,2)=-(DJAC(2,2)+DJAC(2,1))*XJ23
      ENDIF
      IF (BDER(3)) THEN
       DBAS(1,3)= (DJAC(1,2)-DJAC(1,1))*XJ23
       DBAS(2,3)=-DJAC(1,2)*XJ23
       DBAS(3,3)=0D0
       DBAS(4,3)=-2D0*DJAC(1,1)*XJ23
       DBAS(5,3)=-(DJAC(1,2)-2D0*DJAC(1,1))*XJ23
       DBAS(6,3)= (DJAC(1,2)+DJAC(1,1))*XJ23
      ENDIF
      GOTO 99999
C
2216  IF (BDER(1)) THEN
       DBAS(1,1)=0D0
       DBAS(2,1)=0D0
       DBAS(3,1)=0D0
       DBAS(4,1)=0D0
       DBAS(5,1)=1D0
       DBAS(6,1)=0D0
      ENDIF
      IF (BDER(2)) THEN
       DBAS(1,2)=0D0
       DBAS(2,2)= DJAC(2,2)*XJ23
       DBAS(3,2)=-DJAC(2,1)*XJ23
       DBAS(4,2)= (2D0*DJAC(2,1)-DJAC(2,2))*XJ23
       DBAS(5,2)= 2D0*(DJAC(2,2)-DJAC(2,1))*XJ23
       DBAS(6,2)=-(2D0*DJAC(2,2)-DJAC(2,1))*XJ23
      ENDIF
      IF (BDER(3)) THEN
       DBAS(1,3)=0D0
       DBAS(2,3)=-DJAC(1,2)*XJ23
       DBAS(3,3)= DJAC(1,1)*XJ23
       DBAS(4,3)= (DJAC(1,2)-2D0*DJAC(1,1))*XJ23
       DBAS(5,3)=-2D0*(DJAC(1,2)-DJAC(1,1))*XJ23
       DBAS(6,3)= (2D0*DJAC(1,2)-DJAC(1,1))*XJ23
      ENDIF
      GOTO 99999
C
2217  IF (BDER(1)) THEN
       DBAS(1,1)=0D0
       DBAS(2,1)=0D0
       DBAS(3,1)=0D0
       DBAS(4,1)=0D0
       DBAS(5,1)=0D0
       DBAS(6,1)=1D0
      ENDIF
      IF (BDER(2)) THEN
       DBAS(1,2)=-(DJAC(2,2)-DJAC(2,1))*XJ23
       DBAS(2,2)=0D0
       DBAS(3,2)=-DJAC(2,1)*XJ23
       DBAS(4,2)= (DJAC(2,1)+DJAC(2,2))*XJ23
       DBAS(5,2)= (2D0*DJAC(2,2)-DJAC(2,1))*XJ23
       DBAS(6,2)=-2D0*DJAC(2,2)*XJ23
      ENDIF
      IF (BDER(3)) THEN
       DBAS(1,3)= (DJAC(1,2)-DJAC(1,1))*XJ23
       DBAS(2,3)=0D0
       DBAS(3,3)= DJAC(1,1)*XJ23
       DBAS(4,3)=-(DJAC(1,1)+DJAC(1,2))*XJ23
       DBAS(5,3)=-(2D0*DJAC(1,2)-DJAC(1,1))*XJ23
       DBAS(6,3)= 2D0*DJAC(1,2)*XJ23
      ENDIF
      GOTO 99999
C
2218  IF (BDER(1)) THEN
       DBAS(1,1)=0D0
       DBAS(2,1)=0D0
       DBAS(3,1)=0D0
       DBAS(4,1)= .5D0+XI2-XI3
       DBAS(5,1)=0D0
       DBAS(6,1)= .5D0-XI2+XI3
      ENDIF
      IF (BDER(2)) THEN
       DBAS(1,2)=-(DJAC(2,2)-DJAC(2,1))*XJ22
       DBAS(2,2)=0D0
       DBAS(3,2)=0D0
       DBAS(4,2)= (DJAC(2,1)+DJAC(2,2))*XJ22
       DBAS(5,2)= (DJAC(2,2)-DJAC(2,1))*XJ22
       DBAS(6,2)=-(DJAC(2,2)+DJAC(2,1))*XJ22
      ENDIF
      IF (BDER(3)) THEN
       DBAS(1,3)= (DJAC(1,2)-DJAC(1,1))*XJ22
       DBAS(2,3)=0D0
       DBAS(3,3)=0D0
       DBAS(4,3)=-(DJAC(1,1)+DJAC(1,2))*XJ22
       DBAS(5,3)=-(DJAC(1,2)-DJAC(1,1))*XJ22
       DBAS(6,3)= (DJAC(1,2)+DJAC(1,1))*XJ22
      ENDIF
      GOTO 99999
C
2219  IF (BDER(1)) THEN
       DBAS(1,1)=0D0
       DBAS(2,1)=0D0
       DBAS(3,1)=0D0
       DBAS(4,1)= XI1+.5D0-XI3
       DBAS(5,1)=-XI1+.5D0+XI3
       DBAS(6,1)=0D0
      ENDIF
      IF (BDER(2)) THEN
       DBAS(1,2)=0D0
       DBAS(2,2)= DJAC(2,2)*XJ22
       DBAS(3,2)=0D0
       DBAS(4,2)=-(DJAC(2,2)-2D0*DJAC(2,1))*XJ22
       DBAS(5,2)= (DJAC(2,2)-2D0*DJAC(2,1))*XJ22
       DBAS(6,2)=-DJAC(2,2)*XJ22
      ENDIF
      IF (BDER(3)) THEN
       DBAS(1,3)=0D0
       DBAS(2,3)=-DJAC(1,2)*XJ22
       DBAS(3,3)=0D0
       DBAS(4,3)= (DJAC(1,2)-2D0*DJAC(1,1))*XJ22
       DBAS(5,3)=-(DJAC(1,2)-2D0*DJAC(1,1))*XJ22
       DBAS(6,3)= DJAC(1,2)*XJ22
      ENDIF
      GOTO 99999
C
2220  IF (BDER(1)) THEN
       DBAS(1,1)=0D0
       DBAS(2,1)=0D0
       DBAS(3,1)=0D0
       DBAS(4,1)=0D0
       DBAS(5,1)=-XI1+XI2+.5D0
       DBAS(6,1)= XI1-XI2+.5D0
      ENDIF
      IF (BDER(2)) THEN
       DBAS(1,2)=0D0
       DBAS(2,2)=0D0
       DBAS(3,2)=-DJAC(2,1)*XJ22
       DBAS(4,2)= DJAC(2,1)*XJ22
       DBAS(5,2)= (2D0*DJAC(2,2)-DJAC(2,1))*XJ22
       DBAS(6,2)=-(2D0*DJAC(2,2)-DJAC(2,1))*XJ22
      ENDIF
      IF (BDER(3)) THEN
       DBAS(1,3)=0D0
       DBAS(2,3)=0D0
       DBAS(3,3)= DJAC(1,1)*XJ22
       DBAS(4,3)=-DJAC(1,1)*XJ22
       DBAS(5,3)=-(2D0*DJAC(1,2)-DJAC(1,1))*XJ22
       DBAS(6,3)= (2D0*DJAC(1,2)-DJAC(1,1))*XJ22
      ENDIF
C
99999 END
C
C
C
      SUBROUTINE E022A(XI1,XI2,XI3,KCASE,ICUBP0)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (Q2=.5D0)
      DIMENSION KCASE(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /OUTPUT/,/ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('E022A ','03/21/89')
C
C *** Determine position if cubature point
C *** Determine local element
      BXI1=XI1.EQ.Q2
      BXI2=XI2.EQ.Q2
      BXI3=XI3.EQ.Q2
      IF (.NOT.(BXI1.OR.BXI2.OR.BXI3)) THEN
       KCASE(ICUBP0)=1
       IF (XI1.GT.Q2) KCASE(ICUBP0)=2
       IF (XI2.GT.Q2) KCASE(ICUBP0)=3
       IF (XI3.GT.Q2) KCASE(ICUBP0)=4
      ELSE IF (BXI1) THEN
       KCASE(ICUBP0)=8
       IF (BXI2) KCASE(ICUBP0)=5
       IF (BXI3) KCASE(ICUBP0)=7
      ELSE IF (BXI2) THEN
       KCASE(ICUBP0)=9
       IF (BXI3) KCASE(ICUBP0)=6
      ELSE
       KCASE(ICUBP0)=10
      ENDIF
      END
