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
* E012                                                                 *
*                                                                      *
* Purpose  Calculation of values and derivatives of element  12        *
*          Reduced biquadratic element (8 nodes)                       *
*                                                                      *
* Subroutines/functions called   E012A                                 *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* XI1      R*8    Evaluation point in cartesian coordinates            *
* XI2      R*8    with respect to the reference square                 *
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
* IPAR     I*4    Set to  12  if called with IPAR = -1                 *
* IER      I*4    Error indicator                                      *
*                 -131  Desired derivative not available               *
*                 -132  Triangle has vanishing area                    *
*                 -130  Clockwise ordering of the corner points        *
*                                                                      *
************************************************************************
C
      SUBROUTINE E012(XI1,XI2,IPAR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NDOFL=8)
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
      SUB='E012'
      IF (ICHECK.GE.998) CALL OTRC('E012  ','01/02/89')
C
      IER=0
      GOTO (10,1,2,3) -IPAR+1
C
C *** Evaluate at point (XI1,XI2)
10    ICUBP0=1
      CALL E012A(XI1,XI2,DHELP,ICUBP0)
      GOTO 4
C
C *** Set number of element
1     IPAR=12
      GOTO 99999
C
C *** Evaluate basis functions on the reference element
C *** for each cubature point
2     DO 11 ICUBP0=1,NCUBP
      CALL E012A(DXI(ICUBP0,1),DXI(ICUBP0,2),DHELP,ICUBP0)
11    CONTINUE
      GOTO 99999
C
C *** Form linear combinations of the values in DHELP
C *** for the actual element
3     ICUBP0=ICUBP
C
4     IF (ICHECK.EQ.0) GOTO 1200
      IF (DETJ.LT.0D0) CALL WERR(-130,'E012  ')
      IF (DETJ.GT.0D0.AND.DETJ.LT.1D-70) CALL WERR(-132,'E012  ')
C *** No second order derivatives available
C *** Used for second order problems only
      IF (BDER(4).OR.BDER(5).OR.BDER(6)) CALL WERR(-131,'E012  ')
      IF (IER.NE.0) GOTO 99999
C
1200  IF (.NOT.BDER(1)) GOTO 1202
      DO 1201 IDFL=1,NDOFL
1201  DBAS(IDFL,1)=DHELP(IDFL,1,ICUBP0)
1202  XJ1=1D0/DETJ
      IF (.NOT.BDER(2)) GOTO 1204
      DO 1203 IDFL=1,NDOFL
1203  DBAS(IDFL,2)=(DHELP(IDFL,2,ICUBP0)*DJAC(2,2)
     1              -DHELP(IDFL,3,ICUBP0)*DJAC(2,1))*XJ1
1204  IF (.NOT.BDER(3)) GOTO 99999
      DO 1205 IDFL=1,NDOFL
1205  DBAS(IDFL,3)=(-DHELP(IDFL,2,ICUBP0)*DJAC(1,2)
     1               +DHELP(IDFL,3,ICUBP0)*DJAC(1,1))*XJ1
C
99999 END
C
C
C
      SUBROUTINE E012A(X1,X2,DHELP,ICUBP0)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNBAS=21,NNDER=6,NNVE=4,NDOFL=8)
      PARAMETER (Q4=.25D0,Q2=.5D0)
      DIMENSION DHELP(NDOFL,3,*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      SAVE /OUTPUT/,/ERRCTL/,/ELEM/
C
      IF (ICHECK.EQ.999) CALL OTRC('E012A ','01/02/89')
C
      IF (.NOT.BDER(1)) GOTO 1
      DHELP(1,1,ICUBP0)=-Q4*(1D0-X1)*(1D0-X2)*(1D0+X1+X2)
      DHELP(2,1,ICUBP0)=-Q4*(1D0+X1)*(1D0-X2)*(1D0-X1+X2)
      DHELP(3,1,ICUBP0)=-Q4*(1D0+X1)*(1D0+X2)*(1D0-X1-X2)
      DHELP(4,1,ICUBP0)=-Q4*(1D0-X1)*(1D0+X2)*(1D0+X1-X2)
      DHELP(5,1,ICUBP0)= Q2*(1D0-X1*X1)*(1D0-X2)
      DHELP(6,1,ICUBP0)= Q2*(1D0+X1)*(1D0-X2*X2)
      DHELP(7,1,ICUBP0)= Q2*(1D0-X1*X1)*(1D0+X2)
      DHELP(8,1,ICUBP0)= Q2*(1D0-X1)*(1D0-X2*X2)
1     IF (.NOT.(BDER(2).OR.BDER(3)))  GOTO 99999
      DHELP(1,2,ICUBP0)= Q4*(1D0-X2)*(2D0*X1+X2)
      DHELP(2,2,ICUBP0)= Q4*(1D0-X2)*(2D0*X1-X2)
      DHELP(3,2,ICUBP0)= Q4*(1D0+X2)*(2D0*X1+X2)
      DHELP(4,2,ICUBP0)= Q4*(1D0+X2)*(2D0*X1-X2)
      DHELP(5,2,ICUBP0)=-(1D0-X2)*X1
      DHELP(6,2,ICUBP0)= Q2*(1D0-X2*X2)
      DHELP(7,2,ICUBP0)=-(1D0+X2)*X1
      DHELP(8,2,ICUBP0)=-Q2*(1D0-X2*X2)
      DHELP(1,3,ICUBP0)= Q4*(1D0-X1)*(X1+2D0*X2)
      DHELP(2,3,ICUBP0)=-Q4*(1D0+X1)*(X1-2D0*X2)
      DHELP(3,3,ICUBP0)= Q4*(1D0+X1)*(X1+2D0*X2)
      DHELP(4,3,ICUBP0)=-Q4*(1D0-X1)*(X1-2D0*X2)
      DHELP(5,3,ICUBP0)=-Q2*(1D0-X1*X1)
      DHELP(6,3,ICUBP0)=-(1D0+X1)*X2
      DHELP(7,3,ICUBP0)= Q2*(1D0-X1*X1)
      DHELP(8,3,ICUBP0)=-(1D0-X1)*X2
C
99999 END
