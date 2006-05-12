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
* E002                                                                 *
*                                                                      *
* Purpose  Calculation of values and derivatives of element  2         *
*          Quadratic Lagrangian element                                *
*                                                                      *
* Subroutines/functions called   E002A                                 *
*                                                                      *
* Version from  01/02/89                                               *
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
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DBAS     R*8    Array of function values and derivatives for         *
*                 all basis functions (block /ELEM/)                   *
* IPAR     I*4    Set to  2  if called with IPAR = -1                  *
* IER      I*4    Error indicator                                      *
*                 -131  Desired derivative not available               *
*                 -132  Triangle has vanishing area                    *
*                 -130  Clockwise ordering of the corner points        *
*                                                                      *
************************************************************************
C
      SUBROUTINE E002(XI1,XI2,XI3,IPAR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NDOFL=6)
      DIMENSION DHELP(NDOFL,3,NNCUBP)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      SAVE /OUTPUT/,/ERRCTL/,/ELEM/,/CHAR/,/CUB/
      SAVE DHELP
C
      SUB='E002'
      IF (ICHECK.GE.998) CALL OTRC('E002  ','01/02/89')
C
      IER=0
      GOTO (10,1,2,3) -IPAR+1
C
C *** Evaluate at point (XI1,XI2,XI3)
10    ICUBP0=1
      CALL E002A(XI1,XI2,XI3,DHELP,ICUBP0)
      GOTO 4
C
C *** Set number of element
1     IPAR=2
      GOTO 99999
C
C *** Evaluate basis functions on the reference element
C *** for each cubature point
C
2     DO 11 ICUBP0=1,NCUBP
      CALL E002A(DXI(ICUBP0,1),DXI(ICUBP0,2),DXI(ICUBP0,3),
     *            DHELP,ICUBP0)
11    CONTINUE
      GOTO 99999
C
C *** Form linear combination according to actual element coordinates
C
3     ICUBP0=ICUBP
C
4     IF (ICHECK.EQ.0) GOTO 200
      IF (DETJ.LT.0D0) CALL WERR(-130,'E002  ')
      IF (DETJ.GT.0D0.AND.DETJ.LT.1D-70) CALL WERR(-132,'E002  ')
C *** No second order derivatives available
C *** Used for second order problems only
      IF (BDER(4).OR.BDER(5).OR.BDER(6)) CALL WERR(-131,'E002  ')
      IF (IER.NE.0) GOTO 99999
C
C *** Function values
200   IF (.NOT.BDER(1)) GOTO 202
      DO 201 IDFL=1,NDOFL
201   DBAS(IDFL,1)=DHELP(IDFL,1,ICUBP0)
202   IF (.NOT.(BDER(2).OR.BDER(3))) GOTO 99999
C
C***  First order derivatives
      XJ1=1D0/DETJ
      P1=(1D0-4D0*XI1)*XJ1
      P2=(4D0*XI2-1D0)*XJ1
      P3=(4D0*XI3-1D0)*XJ1
C
      IF (.NOT.BDER(2)) GOTO 203
      DBAS(1,2)= (DJAC(2,2)-DJAC(2,1))*P1
      DBAS(2,2)= DJAC(2,2)*P2
      DBAS(3,2)=-DJAC(2,1)*P3
      DBAS(4,2)= 4D0*(XI1*DJAC(2,2)-XI2*(DJAC(2,2)-DJAC(2,1)))*XJ1
      DBAS(5,2)= 4D0*(XI3*DJAC(2,2)-XI2*DJAC(2,1))*XJ1
      DBAS(6,2)= 4D0*(-XI1*DJAC(2,1)-XI3*(DJAC(2,2)-DJAC(2,1)))*XJ1
C
203   IF (.NOT.BDER(3)) GOTO 99999
      DBAS(1,3)=-(DJAC(1,2)-DJAC(1,1))*P1
      DBAS(2,3)=-DJAC(1,2)*P2
      DBAS(3,3)= DJAC(1,1)*P3
      DBAS(4,3)= 4D0*(-XI1*DJAC(1,2)+XI2*(DJAC(1,2)-DJAC(1,1)))*XJ1
      DBAS(5,3)= 4D0*(-XI3*DJAC(1,2)+XI2*DJAC(1,1))*XJ1
      DBAS(6,3)= 4D0*(XI1*DJAC(1,1)+XI3*(DJAC(1,2)-DJAC(1,1)))*XJ1
C
99999 END
C
C
C
      SUBROUTINE E002A(XI1,XI2,XI3,DHELP,ICUBP0)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNBAS=21,NNDER=6,NNVE=4,NDOFL=6)
      DIMENSION DHELP(NDOFL,3,*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      SAVE /OUTPUT/,/ERRCTL/,/ELEM/
C
      IF (ICHECK.EQ.999) CALL OTRC('E002A ','01/02/89')
C
      IF (BDER(1)) THEN
       DHELP(1,1,ICUBP0)= XI1*(XI1-XI2-XI3)
       DHELP(2,1,ICUBP0)= XI2*(XI2-XI1-XI3)
       DHELP(3,1,ICUBP0)= XI3*(XI3-XI1-XI2)
       DHELP(4,1,ICUBP0)= 4D0*XI1*XI2
       DHELP(5,1,ICUBP0)= 4D0*XI2*XI3
       DHELP(6,1,ICUBP0)= 4D0*XI1*XI3
      ENDIF
C
      END
