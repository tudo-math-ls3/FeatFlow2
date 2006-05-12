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
* E034                                                                 *
*                                                                      *
* Purpose  Calculation of values and derivatives of element  34        *
*          Piecewice constant element (4Q0 element)                    *
*          Numbering convention:                                       *
*                                  ***********                         *
*                                  * E4 * E3 *                         *
*                                  ***********                         *
*                                  * E1 * E2 *                         *
*                                  ***********                         *
*                                                                      *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  03/21/89                                               *
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
      SUBROUTINE E034 (XI1,XI2,IPAR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NDOFL=4)
      PARAMETER (Q2=.5D0,Q4=.25D0)
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
      SUB='E034'
      IF (ICHECK.GE.998) CALL OTRC('E034  ','03/21/89')
C
      IER=0
      GOTO (10,1,2,3) -IPAR+1
C
10    ICUBP0=1
      CALL E034A(XI1,XI2,KCASE,ICUBP0)
      GOTO 3400
C
1     IPAR=34
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
22     CALL E034A (DXI(ICUBP0,1),DXI(ICUBP0,2),KCASE,ICUBP0)
      ENDIF
      GOTO 99999
C
3     ICUBP0=ICUBP
C
3400  IF (ICHECK.EQ.0) GOTO 3401
      IF (DETJ.LT.0D0) CALL WERR(-130,'E034  ')
      IF (DETJ.GT.0D0.AND.DETJ.LT.1D-70) CALL WERR(-132,'E034  ')
C *** No first and second order derivatives available
      IF (BDER(2).OR.BDER(3).OR.BDER(4).OR.BDER(5).OR.BDER(6))
     *    CALL WERR(-131,'E034  ')
      IF (IER.NE.0) GOTO 99999
C
3401  GOTO (3411,3412,3413,3414,3415,3416,3417,3418,3419) KCASE(ICUBP0)
C
3411  DBAS(1,1)=1D0
      DBAS(2,1)=0D0
      DBAS(3,1)=0D0
      DBAS(4,1)=0D0
      GOTO 99999
C
3412  DBAS(1,1)=0D0
      DBAS(2,1)=1D0
      DBAS(3,1)=0D0
      DBAS(4,1)=0D0
      GOTO 99999
C
3413  DBAS(1,1)=0D0
      DBAS(2,1)=0D0
      DBAS(3,1)=1D0
      DBAS(4,1)=0D0
      GOTO 99999
C
3414  DBAS(1,1)=0D0
      DBAS(2,1)=0D0
      DBAS(3,1)=0D0
      DBAS(4,1)=1D0
      GOTO 99999
C
3415  DBAS(1,1)=Q2
      DBAS(2,1)=Q2
      DBAS(3,1)=0D0
      DBAS(4,1)=0D0
      GOTO 99999
C
3416  DBAS(1,1)=0D0
      DBAS(2,1)=Q2
      DBAS(3,1)=Q2
      DBAS(4,1)=0D0
      GOTO 99999
C
3417  DBAS(1,1)=0D0
      DBAS(2,1)=0D0
      DBAS(3,1)=Q2
      DBAS(4,1)=Q2
      GOTO 99999
C
3418  DBAS(1,1)=Q2
      DBAS(2,1)=0D0
      DBAS(3,1)=0D0
      DBAS(4,1)=Q2
C
3419  DBAS(1,1)=Q4
      DBAS(2,1)=Q4
      DBAS(3,1)=Q4
      DBAS(4,1)=Q4
      GOTO 99999
C
99999 END
C
C
C
      SUBROUTINE E034A(XI1,XI2,KCASE,ICUBP0)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KCASE(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /OUTPUT/,/ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('E034A ','03/21/89')
C
C *** Determine position if cubature point
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
