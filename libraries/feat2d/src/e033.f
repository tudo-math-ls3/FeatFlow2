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
* E033                                                                 *
*                                                                      *
* Purpose  Calculation of values and derivatives of element  33        *
*          Piecewice bilinear element (4Q1/Q1 element)                 *
*          (same #dof as biquadratic element)                          *
*          Numbering convention:                                       *
*                                  4    7    3                         *
*                                  +****+****+                         *
*                                  * E4 * E3 *                         *
*                                 8+****9****+6                        *
*                                  * E1 * E2 *                         *
*                                  +****+****+                         *
*                                  1    5    2                         *
*                                                                      *
* Caution: To evaluate function values use piecewise 1x1 Gauss cubature*
*          formula or piecewise trapezoidal rule,                      *
*          to evaluate derivatives use piecewise 1x1 Gauss cubature    *
*          formula in subroutine CB2Q only                             *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  01/03/90                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* XI1      R*8                                                         *
* XI2      R*8    Evaluation point in cartesian coordinates            *
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
* IPAR     I*4    Set to  33  if called with IPAR = -1                 *
* IER      I*4    Error indicator                                      *
*                 -131  Desired derivative not available               *
*                 -132  Triangle has vanishing area                    *
*                 -130  Clockwise ordering of the corner points        *
*                                                                      *
************************************************************************
C
      SUBROUTINE E033 (XI1,XI2,IPAR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NDOFL=9)
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
      SUB='E033'
      IF (ICHECK.GE.998) CALL OTRC('E033  ','03/21/89')
C
      IER=0
      GOTO (10,1,2,3) -IPAR+1
C
10    ICUBP0=1
      CALL E033A(XI1,XI2,KCASE,ICUBP0)
      GOTO 3300
C
1     IPAR=33
      GOTO 99999
C
2     IF (ICUBP.EQ.13) THEN
C *** piecewise trapezoidal rule
       DO 21 ICUBP1=1,4
       IAUX=(ICUBP1-1)*NCUBP/4
       DO 21 ICUBP0=IAUX+1,IAUX+NCUBP/4
21     KCASE(ICUBP0)=ICUBP1
      ELSE
C *** Other integration rules not defined piecewise
C *** or with interior cubature points only
       DO 22 ICUBP0=1,NCUBP
22     CALL E033A (DXI(ICUBP0,1),DXI(ICUBP0,2),KCASE,ICUBP0)
      ENDIF
      GOTO 99999
C
3     ICUBP0=ICUBP
C
3300  IF (ICHECK.EQ.0) GOTO 3301
      IF (DETJ.LT.0D0) CALL WERR(-130,'E033  ')
      IF (DETJ.GT.0D0.AND.DETJ.LT.1D-70) CALL WERR(-132,'E033  ')
C *** No second order derivatives available
C *** Used for second order problems only
      IF (BDER(4).OR.BDER(5).OR.BDER(6)) CALL WERR(-131,'E033  ')
      IF (IER.NE.0) GOTO 99999
C
3301  XJ1=1D0/DETJ
C
      GOTO (3311,3312,3313,3314,3315,3316,3317,3318,3319) KCASE(ICUBP0)
C
3311  IF (BDER(1)) THEN
       DBAS(1,1)= XI1*XI2
       DBAS(2,1)= 0D0
       DBAS(3,1)= 0D0
       DBAS(4,1)= 0D0
       DBAS(5,1)=-XI1*XI2-XI2
       DBAS(6,1)= 0D0
       DBAS(7,1)= 0D0
       DBAS(8,1)=-XI1*XI2-XI1
       DBAS(9,1)= XI1*XI2+XI1+XI2+1D0
      ENDIF
      IF (BDER(2)) THEN
       DBAS(1,2)= (XI2*DJAC(2,2)-XI1*DJAC(2,1))*XJ1
       DBAS(2,2)= 0D0
       DBAS(3,2)= 0D0
       DBAS(4,2)= 0D0
       DBAS(5,2)=(-XI2*DJAC(2,2)+(XI1+1D0)*DJAC(2,1))*XJ1
       DBAS(6,2)= 0D0
       DBAS(7,2)= 0D0
       DBAS(8,2)=((-XI2-1D0)*DJAC(2,2)+XI1*DJAC(2,1))*XJ1
       DBAS(9,2)=((XI2+1D0)*DJAC(2,2)-(XI1+1D0)*DJAC(2,1))*XJ1
      ENDIF
      IF (BDER(3)) THEN
       DBAS(1,3)= (-XI2*DJAC(1,2)+XI1*DJAC(1,1))*XJ1
       DBAS(2,3)= 0D0
       DBAS(3,3)= 0D0
       DBAS(4,3)= 0D0
       DBAS(5,3)=(XI2*DJAC(1,2)-(XI1+1D0)*DJAC(1,1))*XJ1
       DBAS(6,3)= 0D0
       DBAS(7,3)= 0D0
       DBAS(8,3)=((XI2+1D0)*DJAC(1,2)-XI1*DJAC(1,1))*XJ1
       DBAS(9,3)=(-(XI2+1D0)*DJAC(1,2)+(XI1+1D0)*DJAC(1,1))*XJ1
      ENDIF
      GOTO 99999
C
3312  IF (BDER(1)) THEN
       DBAS(1,1)= 0D0
       DBAS(2,1)=-XI1*XI2
       DBAS(3,1)= 0D0
       DBAS(4,1)= 0D0
       DBAS(5,1)= XI1*XI2-XI2
       DBAS(6,1)= XI1*XI2+XI1
       DBAS(7,1)= 0D0
       DBAS(8,1)= 0D0
       DBAS(9,1)=-XI1*XI2-XI1+XI2+1D0
      ENDIF
      IF (BDER(2)) THEN
       DBAS(1,2)= 0D0
       DBAS(2,2)= (-XI2*DJAC(2,2)+XI1*DJAC(2,1))*XJ1
       DBAS(3,2)= 0D0
       DBAS(4,2)= 0D0
       DBAS(5,2)= (XI2*DJAC(2,2)-(XI1-1D0)*DJAC(2,1))*XJ1
       DBAS(6,2)= ((XI2+1D0)*DJAC(2,2)-XI1*DJAC(2,1))*XJ1
       DBAS(7,2)= 0D0
       DBAS(8,2)= 0D0
       DBAS(9,2)= ((-XI2-1D0)*DJAC(2,2)+(XI1-1D0)*DJAC(2,1))*XJ1
      ENDIF
      IF (BDER(3)) THEN
       DBAS(1,3)= 0D0
       DBAS(2,3)= (XI2*DJAC(1,2)-XI1*DJAC(1,1))*XJ1
       DBAS(3,3)= 0D0
       DBAS(4,3)= 0D0
       DBAS(5,3)= (-XI2*DJAC(1,2)+(XI1-1D0)*DJAC(1,1))*XJ1
       DBAS(6,3)= ((-XI2-1D0)*DJAC(1,2)+XI1*DJAC(1,1))*XJ1
       DBAS(7,3)= 0D0
       DBAS(8,3)= 0D0
       DBAS(9,3)= ((XI2+1D0)*DJAC(1,2)-(XI1-1D0)*DJAC(1,1))*XJ1
      ENDIF
      GOTO 99999
C
3313  IF (BDER(1)) THEN
       DBAS(1,1)= 0D0
       DBAS(2,1)= 0D0
       DBAS(3,1)= XI1*XI2
       DBAS(4,1)= 0D0
       DBAS(5,1)= 0D0
       DBAS(6,1)=-XI1*XI2+XI1
       DBAS(7,1)=-XI1*XI2+XI2
       DBAS(8,1)= 0D0
       DBAS(9,1)= XI1*XI2-XI1-XI2+1D0
      ENDIF
      IF (BDER(2)) THEN
       DBAS(1,2)= 0D0
       DBAS(2,2)= 0D0
       DBAS(3,2)= (XI2*DJAC(2,2)-XI1*DJAC(2,1))*XJ1
       DBAS(4,2)= 0D0
       DBAS(5,2)= 0D0
       DBAS(6,2)= ((-XI2+1D0)*DJAC(2,2)+XI1*DJAC(2,1))*XJ1
       DBAS(7,2)= (-XI2*DJAC(2,2)+(XI1-1D0)*DJAC(2,1))*XJ1
       DBAS(8,2)= 0D0
       DBAS(9,2)= ((XI2-1D0)*DJAC(2,2)-(XI1-1D0)*DJAC(2,1))*XJ1
      ENDIF
      IF (BDER(3)) THEN
       DBAS(1,3)= 0D0
       DBAS(2,3)= 0D0
       DBAS(3,3)= (-XI2*DJAC(1,2)+XI1*DJAC(1,1))*XJ1
       DBAS(4,3)= 0D0
       DBAS(5,3)= 0D0
       DBAS(6,3)= ((XI2-1D0)*DJAC(1,2)-XI1*DJAC(1,1))*XJ1
       DBAS(7,3)= (XI2*DJAC(1,2)-(XI1-1D0)*DJAC(1,1))*XJ1
       DBAS(8,3)= 0D0
       DBAS(9,3)= ((-XI2+1D0)*DJAC(1,2)+(XI1-1D0)*DJAC(1,1))*XJ1
      ENDIF
      GOTO 99999
C
3314  IF (BDER(1)) THEN
       DBAS(1,1)= 0D0
       DBAS(2,1)= 0D0
       DBAS(3,1)= 0D0
       DBAS(4,1)=-XI1*XI2
       DBAS(5,1)= 0D0
       DBAS(6,1)= 0D0
       DBAS(7,1)= XI1*XI2+XI2
       DBAS(8,1)= XI1*XI2-XI1
       DBAS(9,1)=-XI1*XI2+XI1-XI2+1D0
      ENDIF
      IF (BDER(2)) THEN
       DBAS(1,2)= 0D0
       DBAS(2,2)= 0D0
       DBAS(3,2)= 0D0
       DBAS(4,2)= (-XI2*DJAC(2,2)+XI1*DJAC(2,1))*XJ1
       DBAS(5,2)= 0D0
       DBAS(6,2)= 0D0
       DBAS(7,2)= (XI2*DJAC(2,2)-(XI1+1D0)*DJAC(2,1))*XJ1
       DBAS(8,2)= ((XI2-1D0)*DJAC(2,2)-XI1*DJAC(2,1))*XJ1
       DBAS(9,2)= ((-XI2+1D0)*DJAC(2,2)+(XI1+1D0)*DJAC(2,1))*XJ1
      ENDIF
      IF (BDER(3)) THEN
       DBAS(1,3)= 0D0
       DBAS(2,3)= 0D0
       DBAS(3,3)= 0D0
       DBAS(4,3)= (XI2*DJAC(1,2)-XI1*DJAC(1,1))*XJ1
       DBAS(5,3)= 0D0
       DBAS(6,3)= 0D0
       DBAS(7,3)= (-XI2*DJAC(1,2)+(XI1+1D0)*DJAC(1,1))*XJ1
       DBAS(8,3)= ((-XI2+1D0)*DJAC(1,2)+XI1*DJAC(1,1))*XJ1
       DBAS(9,3)= ((XI2-1D0)*DJAC(1,2)-(XI1+1D0)*DJAC(1,1))*XJ1
      ENDIF
      GOTO 99999
C
3315  IF (BDER(1)) THEN
       DBAS(1,1)= 0D0
       DBAS(2,1)= 0D0
       DBAS(3,1)= 0D0
       DBAS(4,1)= 0D0
       DBAS(5,1)=-XI2
       DBAS(6,1)= 0D0
       DBAS(7,1)= 0D0
       DBAS(8,1)= 0D0
       DBAS(9,1)= XI2+1D0
      ENDIF
      IF (BDER(2)) THEN
       DBAS(1,2)= XI2*DJAC(2,2)*XJ1*.5D0
       DBAS(2,2)=-XI2*DJAC(2,2)*XJ1*.5D0
       DBAS(3,2)= 0D0
       DBAS(4,2)= 0D0
       DBAS(5,2)= DJAC(2,1)*XJ1
       DBAS(6,2)= (XI2+1D0)*DJAC(2,2)*XJ1*.5D0
       DBAS(7,2)= 0D0
       DBAS(8,2)= (-XI2-1D0)*DJAC(2,2)*XJ1*.5D0
       DBAS(9,2)=-DJAC(2,1)*XJ1
      ENDIF
      IF (BDER(3)) THEN
       DBAS(1,3)=-XI2*DJAC(1,2)*XJ1*.5D0
       DBAS(2,3)= XI2*DJAC(1,2)*XJ1*.5D0
       DBAS(3,3)= 0D0
       DBAS(4,3)= 0D0
       DBAS(5,3)=-DJAC(1,1)*XJ1
       DBAS(6,3)= (-XI2-1D0)*DJAC(1,2)*XJ1*.5D0
       DBAS(7,3)= 0D0
       DBAS(8,3)= (XI2+1D0)*DJAC(1,2)*XJ1*.5D0
       DBAS(9,3)= DJAC(1,1)*XJ1
      ENDIF
      GOTO 99999
C
3316  IF (BDER(1)) THEN
       DBAS(1,1)= 0D0
       DBAS(2,1)= 0D0
       DBAS(3,1)= 0D0
       DBAS(4,1)= 0D0
       DBAS(5,1)= 0D0
       DBAS(6,1)= XI1
       DBAS(7,1)= 0D0
       DBAS(8,1)= 0D0
       DBAS(9,1)=-XI1+1D0
      ENDIF
      IF (BDER(2)) THEN
       DBAS(1,2)= 0D0
       DBAS(2,2)= XI1*DJAC(2,1)*XJ1*.5D0
       DBAS(3,2)=-XI1*DJAC(2,1)*XJ1*.5D0
       DBAS(4,2)= 0D0
       DBAS(5,2)=-(XI1-1D0)*DJAC(2,1)*XJ1*.5D0
       DBAS(6,2)= DJAC(2,2)*XJ1
       DBAS(7,2)= (XI1-1D0)*DJAC(2,1)*XJ1*.5D0
       DBAS(8,2)= 0D0
       DBAS(9,2)= -DJAC(2,2)*XJ1
      ENDIF
      IF (BDER(3)) THEN
       DBAS(1,3)= 0D0
       DBAS(2,3)=-XI1*DJAC(1,1)*XJ1*.5D0
       DBAS(3,3)= XI1*DJAC(1,1)*XJ1*.5D0
       DBAS(4,3)= 0D0
       DBAS(5,3)= (XI1-1D0)*DJAC(1,1)*XJ1*.5D0
       DBAS(6,3)= -DJAC(1,2)*XJ1
       DBAS(7,3)= -(XI1-1D0)*DJAC(1,1)*XJ1*.5D0
       DBAS(8,3)= XI1*DJAC(1,1)*XJ1*.5D0
       DBAS(9,3)= DJAC(1,2)*XJ1
      ENDIF
      GOTO 99999
C
3317  IF (BDER(1)) THEN
       DBAS(1,1)= 0D0
       DBAS(2,1)= 0D0
       DBAS(3,1)= 0D0
       DBAS(4,1)= 0D0
       DBAS(5,1)= 0D0
       DBAS(6,1)= 0D0
       DBAS(7,1)= XI2
       DBAS(8,1)= 0D0
       DBAS(9,1)=-XI2+1D0
      ENDIF
      IF (BDER(2)) THEN
       DBAS(1,2)= 0D0
       DBAS(2,2)= 0D0
       DBAS(3,2)= XI2*DJAC(2,2)*XJ1*.5D0
       DBAS(4,2)=-XI2*DJAC(2,2)*XJ1*.5D0
       DBAS(5,2)= 0D0
       DBAS(6,2)= (-XI2+1D0)*DJAC(2,2)*XJ1*.5D0
       DBAS(7,2)= -DJAC(2,1)*XJ1
       DBAS(8,2)= (XI2-1D0)*DJAC(2,2)*XJ1*.5D0
       DBAS(9,2)= DJAC(2,1)*XJ1
      ENDIF
      IF (BDER(3)) THEN
       DBAS(1,3)= 0D0
       DBAS(2,3)= 0D0
       DBAS(3,3)=-XI2*DJAC(1,2)*XJ1*.5D0
       DBAS(4,3)= XI2*DJAC(1,2)*XJ1*.5D0
       DBAS(5,3)= 0D0
       DBAS(6,3)= (XI2-1D0)*DJAC(1,2)*XJ1*.5D0
       DBAS(7,3)= DJAC(1,1)*XJ1
       DBAS(8,3)= (-XI2+1D0)*DJAC(1,2)*XJ1*.5D0
       DBAS(9,3)= -DJAC(1,1)*XJ1
      ENDIF
      GOTO 99999
C
3318  IF (BDER(1)) THEN
       DBAS(1,1)= 0D0
       DBAS(2,1)= 0D0
       DBAS(3,1)= 0D0
       DBAS(4,1)= 0D0
       DBAS(5,1)= 0D0
       DBAS(6,1)= 0D0
       DBAS(7,1)= 0D0
       DBAS(8,1)=-XI1
       DBAS(9,1)= XI1+1D0
      ENDIF
      IF (BDER(2)) THEN
       DBAS(1,2)=-XI1*DJAC(2,1)*XJ1*.5D0
       DBAS(2,2)= 0D0
       DBAS(3,2)= 0D0
       DBAS(4,2)= XI1*DJAC(2,1)*XJ1*.5D0
       DBAS(5,2)= (XI1+1D0)*DJAC(2,1)*XJ1*.5D0
       DBAS(6,2)= 0D0
       DBAS(7,2)=-(XI1+1D0)*DJAC(2,1)*XJ1*.5D0
       DBAS(8,2)=-DJAC(2,2)*XJ1
       DBAS(9,2)= DJAC(2,2)*XJ1
      ENDIF
      IF (BDER(3)) THEN
       DBAS(1,3)= XI1*DJAC(1,1)*XJ1*.5D0
       DBAS(2,3)= 0D0
       DBAS(3,3)= 0D0
       DBAS(4,3)=-XI1*DJAC(1,1)*XJ1*.5D0
       DBAS(5,3)=-(XI1+1D0)*DJAC(1,1)*XJ1*.5D0
       DBAS(6,3)= 0D0
       DBAS(7,3)= (XI1+1D0)*DJAC(1,1)*XJ1*.5D0
       DBAS(8,3)= DJAC(1,2)*XJ1
       DBAS(9,3)=-DJAC(1,2)*XJ1
      ENDIF
      GOTO 99999
C
3319  IF (BDER(1)) THEN
       DBAS(1,1)= 0D0
       DBAS(2,1)= 0D0
       DBAS(3,1)= 0D0
       DBAS(4,1)= 0D0
       DBAS(5,1)= 0D0
       DBAS(6,1)= 0D0
       DBAS(7,1)= 0D0
       DBAS(8,1)= 0D0
       DBAS(9,1)= 1D0
      ENDIF
      IF (BDER(2)) THEN
       DBAS(1,2)= 0D0
       DBAS(2,2)= 0D0
       DBAS(3,2)= 0D0
       DBAS(4,2)= 0D0
       DBAS(5,2)= DJAC(2,1)*XJ1*.5D0
       DBAS(6,2)= DJAC(2,2)*XJ1*.5D0
       DBAS(7,2)=-DJAC(2,1)*XJ1*.5D0
       DBAS(8,2)=-DJAC(2,2)*XJ1*.5D0
       DBAS(9,2)= 0D0
      ENDIF
      IF (BDER(3)) THEN
       DBAS(1,3)= 0D0
       DBAS(2,3)= 0D0
       DBAS(3,3)= 0D0
       DBAS(4,3)= 0D0
       DBAS(5,3)=-DJAC(1,1)*XJ1*.5D0
       DBAS(6,3)=-DJAC(1,2)*XJ1*.5D0
       DBAS(7,3)= DJAC(1,1)*XJ1*.5D0
       DBAS(8,3)= DJAC(1,2)*XJ1*.5D0
       DBAS(9,3)= 0D0
      ENDIF
C
99999 END
C
C
C
      SUBROUTINE E033A(XI1,XI2,KCASE,ICUBP0)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KCASE(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /OUTPUT/,/ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('E033A ','01/03/90')
C
C *** Determine position of cubature point
C *** BC means center of gravity
      BC=XI1.EQ.0D0.AND.XI2.EQ.0D0
      IF (.NOT.BC) THEN
C *** BINT means that cubature point is not on the edges of subdomains
       BINT=XI1.NE.0D0.AND.XI2.NE.0D0
       BXI1L=XI1.LT.0D0
       BXI2L=XI2.LT.0D0
       BXI1G=XI1.GT.0D0
       BXI2G=XI2.GT.0D0
       IF (BINT) THEN
        IF (BXI1L.AND.BXI2L) KCASE(ICUBP0)=1
        IF (BXI1G.AND.BXI2L) KCASE(ICUBP0)=2
        IF (BXI1G.AND.BXI2G) KCASE(ICUBP0)=3
        IF (BXI1L.AND.BXI2G) KCASE(ICUBP0)=4
       ELSE
        IF (BXI2L) KCASE(ICUBP0)=5
        IF (BXI1G) KCASE(ICUBP0)=6
        IF (BXI2G) KCASE(ICUBP0)=7
        IF (BXI1L) KCASE(ICUBP0)=8
       ENDIF
      ELSE
       KCASE(ICUBP0)=9
      ENDIF
C
      END
