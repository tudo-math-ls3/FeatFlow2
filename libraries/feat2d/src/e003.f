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
* E003                                                                 *
*                                                                      *
* Purpose  Calculation of values and derivatives of element  3         *
*          Cubic Hermite element                                       *
*                                                                      *
* Subroutines/functions called   E003A                                 *
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
* IPAR     I*4    Set to  3  if called with IPAR = -1                  *
* IER      I*4    Error indicator                                      *
*                 -131  Desired derivative not available               *
*                 -132  Triangle has vanishing area                    *
*                 -130  Clockwise ordering of the corner points        *
*                                                                      *
************************************************************************
C
      SUBROUTINE E003(XI1,XI2,XI3,IPAR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NDOFL=10)
      DIMENSION DHELP(NDOFL,3,NNCUBP)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/ELEM/,/CUB/
      SAVE DHELP
C
      SUB='E003'
      IF (ICHECK.GE.998) CALL OTRC('E003  ','01/02/89')
C
      IER=0
      GOTO (10,1,2,3) -IPAR+1
C
C *** Evaluate at point (XI1,XI2,XI3)
10    ICUBP0=1
      CALL E003A(XI1,XI2,XI3,DHELP,ICUBP0)
      GOTO 4
C
C *** Set number of element
1     IPAR=3
      GOTO 99999
C
C *** Evaluate basis functions on the reference element
C *** for each cubature point
C
2     DO 11 ICUBP0=1,NCUBP
      CALL E003A(DXI(ICUBP0,1),DXI(ICUBP0,2),DXI(ICUBP0,3),
     *           DHELP,ICUBP0)
11    CONTINUE
      GOTO 99999
C
C *** Form linear combination accordiing to actual element coordinates
C
3     ICUBP0=ICUBP
C
4     IF (ICHECK.EQ.0) GOTO 300
      IF (DETJ.LT.0D0) CALL WERR(-130,'E003  ')
      IF (DETJ.GT.0D0.AND.DETJ.LT.1D-70) CALL WERR(-132,'E003  ')
C *** No second order derivatives available
C *** Used for second order problems only
      IF (BDER(4).OR.BDER(5).OR.BDER(6)) CALL WERR(-131,'E003  ')
      IF (IER.NE.0) GOTO 99999
C
C *** Function values
300   IF (.NOT.BDER(1)) GOTO 302
      DO 301 IDFL=1,NDOFL
301   DBAS(IDFL,1)=DHELP(IDFL,1,ICUBP0)
      DO 311 IDFL=2,8,3
      DBHELP=DBAS(IDFL,1)*DJAC(1,1)+DBAS(IDFL+1,1)*DJAC(1,2)
      DBAS(IDFL+1,1)=DBAS(IDFL,1)*DJAC(2,1)+DBAS(IDFL+1,1)*DJAC(2,2)
      DBAS(IDFL,1)=DBHELP
311   CONTINUE
302   IF (.NOT.(BDER(2).OR.BDER(3))) GOTO 99999
C
C***  First order derivatives
      XJ1=1D0/DETJ
      IF (.NOT.BDER(2)) GOTO 304
      DO 303 IDFL=1,NDOFL
303   DBAS(IDFL,2)=(DHELP(IDFL,2,ICUBP0)*DJAC(2,2)
     *             -DHELP(IDFL,3,ICUBP0)*DJAC(2,1))*XJ1
      DO 313 IDFL=2,8,3
      DBHELP=DBAS(IDFL,2)*DJAC(1,1)+DBAS(IDFL+1,2)*DJAC(1,2)
      DBAS(IDFL+1,2)=DBAS(IDFL,2)*DJAC(2,1)+DBAS(IDFL+1,2)*DJAC(2,2)
      DBAS(IDFL,2)=DBHELP
313   CONTINUE
304   IF (.NOT.BDER(3)) GOTO 99999
      DO 305 IDFL=1,NDOFL
305   DBAS(IDFL,3)=(-DHELP(IDFL,2,ICUBP0)*DJAC(1,2)
     *              +DHELP(IDFL,3,ICUBP0)*DJAC(1,1))*XJ1
      DO 315 IDFL=2,8,3
      DBHELP=DBAS(IDFL,3)*DJAC(1,1)+DBAS(IDFL+1,3)*DJAC(1,2)
      DBAS(IDFL+1,3)=DBAS(IDFL,3)*DJAC(2,1)+DBAS(IDFL+1,3)*DJAC(2,2)
      DBAS(IDFL,3)=DBHELP
315   CONTINUE
C
99999 END
C
C
C
      SUBROUTINE E003A(XI1,XI2,XI3,DHELP,ICUBP0)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNBAS=21,NNDER=6,NNVE=4,NDOFL=10)
      DIMENSION DHELP(NDOFL,3,*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      SAVE /OUTPUT/,/ERRCTL/,/ELEM/
C
      IF (ICHECK.EQ.999) CALL OTRC('E003A ','01/02/89')
C
      A= XI1*XI2*XI3
      D= XI3*(XI1-XI2)
      C= XI2*(XI1-XI3)
C
      IF (.NOT.BDER(1)) GOTO 1
      DHELP(1,1,ICUBP0)= XI1*XI1*(3D0-2D0*XI1)-7D0*A
      DHELP(2,1,ICUBP0)= XI2*XI1*XI1-A
      DHELP(3,1,ICUBP0)= XI3*XI1*XI1-A
      DHELP(4,1,ICUBP0)= XI2*XI2*(3D0-2D0*XI2)-7D0*A
      DHELP(5,1,ICUBP0)= XI2*XI2*(-1D0+XI2)+2D0*A
      DHELP(6,1,ICUBP0)= XI2*XI2*XI3-A
      DHELP(7,1,ICUBP0)= XI3*XI3*(3D0-2D0*XI3)-7D0*A
      DHELP(8,1,ICUBP0)= XI2*XI3*XI3-A
      DHELP(9,1,ICUBP0)= XI3*XI3*(-1D0+XI3)+2D0*A
      DHELP(10,1,ICUBP0)= 27D0*A

1     IF (.NOT.(BDER(2).OR.BDER(3))) GOTO 99999
      DHELP(1,2,ICUBP0)= 6D0*XI1*(-1D0+XI1)-7D0*D
      DHELP(2,2,ICUBP0)= XI1*(XI1-2D0*XI2)-D
      DHELP(3,2,ICUBP0)=-2D0*XI1*XI3-D
      DHELP(4,2,ICUBP0)= 6D0*XI2*(1D0-XI2)-7D0*D
      DHELP(5,2,ICUBP0)= XI2*(-2D0+3D0*XI2)+2D0*D
      DHELP(6,2,ICUBP0)= 2D0*XI2*XI3-D
      DHELP(7,2,ICUBP0)=-7D0*D
      DHELP(8,2,ICUBP0)= XI3*XI3-D
      DHELP(9,2,ICUBP0)= 2D0*D
      DHELP(10,2,ICUBP0)= 27D0*D
C
      DHELP(1,3,ICUBP0)= 6D0*XI1*(-1D0+XI1)-7D0*C
      DHELP(2,3,ICUBP0)=-2D0*XI2*XI1-C
      DHELP(3,3,ICUBP0)= XI1*(XI1-2D0*XI3)-C
      DHELP(4,3,ICUBP0)=-7D0*C
      DHELP(5,3,ICUBP0)= 2D0*C
      DHELP(6,3,ICUBP0)= XI2*XI2-C
      DHELP(7,3,ICUBP0)= 6D0*XI3*(1D0-XI3)-7D0*C
      DHELP(8,3,ICUBP0)= 2D0*XI2*XI3-C
      DHELP(9,3,ICUBP0)= XI3*(-2D0+3D0*XI3)+2D0*C
      DHELP(10,3,ICUBP0)= 27D0*C
C
99999 END
