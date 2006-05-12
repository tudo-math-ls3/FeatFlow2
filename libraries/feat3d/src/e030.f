************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT3D  (Release 1.1)               *
*                                                                      *
* Authors: J. Harig, S. Turek                                          *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* E030                                                                 *
*                                                                      *
* Purpose  Calculation of values and derivatives of element  30        *
*          nonconforming element with 6 local degrees of freedom       *
*          (mean values of areas)                                      *
*                                                                      *
*          A*(X**2-Y**2)+B*(X**2-Z**2)+C*X+D*Y+E*Y+F                   *
*                                                                      *
* Subroutines/functions called   E030A                                 *
*                                                                      *
* Version from  01/17/90                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* XI1      R*8    Evaluation point in cartesian coordinates            *
* XI2      R*8    with respect to the reference cube                   *
* XI3      R*8                                                         *
* IPAR     I*4    Switch for desired action                            *
*                  0  Calculation of the desired values of DBAS        *
*                 -1  Set number of element                            *
*                 -2  Calculate values on the reference element        *
*                     for all cubature points and save them            *
*                 -3  same as 0, but use the values saved before       *
*                                                                      *
* BDER     LOG    Derivative J is calculated if BDER(J).EQ..TRUE.      *
*                 Multiindices are enumbered from 1 to 10              *
* DJAC     R*8    Jacobian                                             *
* DETJ     R*8    Determinant of the jacobian                          *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DBAS     R*8    Array of function values and derivatives for         *
*                 all basis functions (block /ELEM/)                   *
* IPAR     I*4    Set to  30  if called with IPAR = -1                 *
* IER      I*4    Error indicator                                      *
*                 -131  Desired derivative not available               *
*                 -132  Triangle has vanishing area                    *
*                 -130  Clockwise ordering of the corner points        *
*                                                                      *
************************************************************************
C
      SUBROUTINE E030(XI1,XI2,XI3,IPAR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNAE=6,NNVE=8,NNDIM=3)
      DIMENSION DHELP(NNAE,4,NNCUBP)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,
     *                DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),
     *                IEL,NDIM
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      SAVE /OUTPUT/,/ERRCTL/,/ELEM/,/CHAR/,/CUB/
      SAVE DHELP
C
      SUB='E030'
      IF (ICHECK.GE.998) CALL OTRC('E030  ','09/01/89')
C
      IER=0
      GOTO (10,1,2,3) -IPAR+1
C
C *** Evaluate at point (XI1,XI2,XI3)
10    ICUBP0=1
      CALL E030A(XI1,XI2,XI3,DHELP,ICUBP0)
      GOTO 4
C
C *** Set number of element
1     IPAR=30
      GOTO 99999
C
C *** Evaluate basis functions on the reference element
C *** for each cubature point
2     DO 30 ICUBP0=1,NCUBP
      CALL E030A(DXI(ICUBP0,1),DXI(ICUBP0,2),DXI(ICUBP0,3),DHELP,ICUBP0)
30    CONTINUE
      GOTO 99999
C
C *** Form linear combinations of the values in DHELP
C *** for the actual element
3     ICUBP0=ICUBP
C
4     IF (ICHECK.EQ.0) GOTO 3000
      IF (ABS(DETJ).LT.1D-70) CALL WERR(-132,'E030  ')
C *** No second order derivatives available
C *** Used for second order problems only
      IF (BDER(5).OR.BDER(6).OR.BDER(7).OR.BDER(8).OR.BDER(9).OR.
     *    BDER(10)) CALL WERR(-131,'E030  ')
      IF (IER.NE.0) GOTO 99999
C
3000  IF (.NOT.BDER(1)) GOTO 3002
      DO 3001 IDFL=1,NNAE
3001  DBAS(1,IDFL,1)=DHELP(IDFL,1,ICUBP0)
C
3002  XJ1=1D0/DETJ
      IF (.NOT.BDER(2)) GOTO 3004
      DO 3003 IDFL=1,NNAE
3003  DBAS(1,IDFL,2)= XJ1*(
     *   DHELP(IDFL,2,ICUBP0)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
     *  -DHELP(IDFL,3,ICUBP0)*(DJAC(2,1)*DJAC(3,3)-DJAC(3,1)*DJAC(2,3))
     *  +DHELP(IDFL,4,ICUBP0)*(DJAC(2,1)*DJAC(3,2)-DJAC(3,1)*DJAC(2,2)))
3004  IF (.NOT.BDER(3)) GOTO 3006
      DO 3005 IDFL=1,NNAE
3005  DBAS(1,IDFL,3)= XJ1*(
     *  -DHELP(IDFL,2,ICUBP0)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *  +DHELP(IDFL,3,ICUBP0)*(DJAC(1,1)*DJAC(3,3)-DJAC(3,1)*DJAC(1,3))
     *  -DHELP(IDFL,4,ICUBP0)*(DJAC(1,1)*DJAC(3,2)-DJAC(3,1)*DJAC(1,2)))
3006  IF (.NOT.BDER(4)) GOTO 99999
      DO 3007 IDFL=1,NNAE
3007  DBAS(1,IDFL,4)= XJ1*(
     *   DHELP(IDFL,2,ICUBP0)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
     *  -DHELP(IDFL,3,ICUBP0)*(DJAC(1,1)*DJAC(2,3)-DJAC(2,1)*DJAC(1,3))
     *  +DHELP(IDFL,4,ICUBP0)*(DJAC(1,1)*DJAC(2,2)-DJAC(2,1)*DJAC(1,2)))
C
99999 END
C
C
C
      SUBROUTINE E030A(X1,X2,X3,DHELP,ICUBP0)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNBAS=27,NNDER=10,NNDIM=3,NNAE=6,NNVE=8)
      PARAMETER (Q2=0.5D0,Q4=0.25D0,Q6=1D0/6D0)
      DIMENSION DHELP(NNAE,4,1)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,
     *                DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),
     *                IEL,NDIM
      SAVE /OUTPUT/,/ERRCTL/,/ELEM/
C
      IF (ICHECK.EQ.999) CALL OTRC('E030A ','01/12/90')
C
      IF (.NOT.BDER(1)) GOTO 3002
      DHELP(1,1,ICUBP0)= Q4*(X1**2-X2**2)-Q2*(X1**2-X3**2)-Q2*X3+Q6
      DHELP(2,1,ICUBP0)=-Q2*(X1**2-X2**2)+Q4*(X1**2-X3**2)-Q2*X2+Q6
      DHELP(3,1,ICUBP0)= Q4*(X1**2-X2**2)+Q4*(X1**2-X3**2)+Q2*X1+Q6
      DHELP(4,1,ICUBP0)=-Q2*(X1**2-X2**2)+Q4*(X1**2-X3**2)+Q2*X2+Q6
      DHELP(5,1,ICUBP0)= Q4*(X1**2-X2**2)+Q4*(X1**2-X3**2)-Q2*X1+Q6
      DHELP(6,1,ICUBP0)= Q4*(X1**2-X2**2)-Q2*(X1**2-X3**2)+Q2*X3+Q6
C
3002  IF (.NOT.(BDER(2).OR.BDER(3).OR.BDER(4))) GOTO 99999
      DHELP(1,2,ICUBP0)=-Q2*X1
      DHELP(2,2,ICUBP0)=-Q2*X1
      DHELP(3,2,ICUBP0)=    X1+Q2
      DHELP(4,2,ICUBP0)=-Q2*X1
      DHELP(5,2,ICUBP0)=    X1-Q2
      DHELP(6,2,ICUBP0)=-Q2*X1
C
      DHELP(1,3,ICUBP0)=-Q2*X2
      DHELP(2,3,ICUBP0)=    X2-Q2
      DHELP(3,3,ICUBP0)=-Q2*X2
      DHELP(4,3,ICUBP0)=    X2+Q2
      DHELP(5,3,ICUBP0)=-Q2*X2
      DHELP(6,3,ICUBP0)=-Q2*X2
C
      DHELP(1,4,ICUBP0)=    X3-Q2
      DHELP(2,4,ICUBP0)=-Q2*X3
      DHELP(3,4,ICUBP0)=-Q2*X3
      DHELP(4,4,ICUBP0)=-Q2*X3
      DHELP(5,4,ICUBP0)=-Q2*X3
      DHELP(6,4,ICUBP0)=    X3+Q2
C
99999 END
