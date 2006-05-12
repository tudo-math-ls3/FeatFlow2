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
* E013                                                                 *
*                                                                      *
* Purpose  Calculation of values and derivatives of element  13        *
*          Triquadratic element (27 nodes)                             *
*                                                                      *
* Subroutines/functions called   E013A                                 *
*                                                                      *
* Version from  07/22/91                                               *
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
* IPAR     I*4    Set to  13  if called with IPAR = -1                 *
* IER      I*4    Error indicator                                      *
*                 -131  Desired derivative not available               *
*                 -132  Triangle has vanishing area                    *
*                 -130  Clockwise ordering of the corner points        *
*                                                                      *
************************************************************************
C
      SUBROUTINE E013(XI1,XI2,XI3,IPAR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNDIM=3)
      DIMENSION DHELP(NNBAS,4,NNCUBP)
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
      SUB='E013'
      IF (ICHECK.GE.998) CALL OTRC('E013  ','07/22/91')
C
      IER=0
      GOTO (10,1,2,3) -IPAR+1
C
C *** Evaluate at point (XI1,XI2,XI3)
10    ICUBP0=1
      CALL E013A(XI1,XI2,XI3,DHELP,ICUBP0)
      GOTO 4
C
C *** Set number of element
1     IPAR=13
      GOTO 99999
C
C *** Evaluate basis functions on the reference element
C *** for each cubature point
2     DO 11 ICUBP0=1,NCUBP
      CALL E013A(DXI(ICUBP0,1),DXI(ICUBP0,2),DXI(ICUBP0,3),DHELP,ICUBP0)
11    CONTINUE
      GOTO 99999
C
C *** Form linear combinations of the values in DHELP
C *** for the actual element
3     ICUBP0=ICUBP
C
4     IF (ICHECK.EQ.0) GOTO 1100
      IF (ABS(DETJ).LT.1D-70) CALL WERR(-132,'E013  ')
C *** No second order derivatives available
C *** Used for second order problems only
      IF (BDER(5).OR.BDER(6).OR.BDER(7).OR.BDER(8).OR.BDER(9).OR.
     *    BDER(10)) CALL WERR(-131,'E013  ')
      IF (IER.NE.0) GOTO 99999
C
1100  IF (.NOT.BDER(1)) GOTO 1102
      DO 1101 IDFL=1,NNBAS
1101  DBAS(1,IDFL,1)=DHELP(IDFL,1,ICUBP0)
C
1102  XJ1=1D0/DETJ
      IF (.NOT.BDER(2)) GOTO 1104
      DO 1103 IDFL=1,NNBAS
1103  DBAS(1,IDFL,2)= XJ1*(
     *   DHELP(IDFL,2,ICUBP0)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
     *  -DHELP(IDFL,3,ICUBP0)*(DJAC(2,1)*DJAC(3,3)-DJAC(3,1)*DJAC(2,3))
     *  +DHELP(IDFL,4,ICUBP0)*(DJAC(2,1)*DJAC(3,2)-DJAC(3,1)*DJAC(2,2)))
1104  IF (.NOT.BDER(3)) GOTO 1106
      DO 1105 IDFL=1,NNBAS
1105  DBAS(1,IDFL,3)= XJ1*(
     *  -DHELP(IDFL,2,ICUBP0)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *  +DHELP(IDFL,3,ICUBP0)*(DJAC(1,1)*DJAC(3,3)-DJAC(3,1)*DJAC(1,3))
     *  -DHELP(IDFL,4,ICUBP0)*(DJAC(1,1)*DJAC(3,2)-DJAC(3,1)*DJAC(1,2)))
1106  IF (.NOT.BDER(4)) GOTO 99999
      DO 1107 IDFL=1,NNBAS
1107  DBAS(1,IDFL,4)= XJ1*(
     *   DHELP(IDFL,2,ICUBP0)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
     *  -DHELP(IDFL,3,ICUBP0)*(DJAC(1,1)*DJAC(2,3)-DJAC(2,1)*DJAC(1,3))
     *  +DHELP(IDFL,4,ICUBP0)*(DJAC(1,1)*DJAC(2,2)-DJAC(2,1)*DJAC(1,2)))
C
99999 END
C
C
C
      SUBROUTINE E013A(X1,X2,X3,DHELP,ICUBP0)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNBAS=27,NNDER=10,NNDIM=3,NNVE=8)
      PARAMETER (Q8=0.125D0)
      DIMENSION DHELP(NNBAS,4,*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,
     *                DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),
     *                IEL,NDIM
      SAVE /OUTPUT/,/ERRCTL/,/ELEM/
C
      IF (ICHECK.EQ.999) CALL OTRC('E013A ','07/22/91')
C
      IF (.NOT.BDER(1)) GOTO 1102
      DHELP( 1,1,ICUBP0)= Q8*X1*(1D0-X1)*X2*(1D0-X2)*X3*(1D0-X3)
      DHELP( 2,1,ICUBP0)= Q8*X1*(1D0+X1)*X2*(1D0-X2)*X3*(1D0-X3)
      DHELP( 3,1,ICUBP0)= Q8*X1*(1D0+X1)*X2*(1D0+X2)*X3*(1D0-X3)
      DHELP( 4,1,ICUBP0)= Q8*X1*(1D0-X1)*X2*(1D0+X2)*X3*(1D0-X3)
      DHELP( 5,1,ICUBP0)= Q8*X1*(1D0-X1)*X2*(1D0-X2)*X3*(1D0+X3)
      DHELP( 6,1,ICUBP0)= Q8*X1*(1D0+X1)*X2*(1D0-X2)*X3*(1D0+X3)
      DHELP( 7,1,ICUBP0)= Q8*X1*(1D0+X1)*X2*(1D0+X2)*X3*(1D0+X3)
      DHELP( 8,1,ICUBP0)= Q8*X1*(1D0-X1)*X2*(1D0+X2)*X3*(1D0+X3)
      DHELP( 9,1,ICUBP0)= (1D0-X1*X1)*(1D0-X2*X2)*(1D0-X3*X3)
C      DHELP(10,1,ICUBP0)=
C      DHELP(11,1,ICUBP0)=
C      DHELP(12,1,ICUBP0)=
C      DHELP(13,1,ICUBP0)=
C      DHELP(14,1,ICUBP0)=
C      DHELP(15,1,ICUBP0)=
C      DHELP(16,1,ICUBP0)=
C      DHELP(17,1,ICUBP0)=
C      DHELP(18,1,ICUBP0)=
C      DHELP(19,1,ICUBP0)=
C      DHELP(20,1,ICUBP0)=
C      DHELP(21,1,ICUBP0)=
C      DHELP(22,1,ICUBP0)=
C      DHELP(23,1,ICUBP0)=
C      DHELP(24,1,ICUBP0)=
C      DHELP(25,1,ICUBP0)=
C      DHELP(26,1,ICUBP0)=
C      DHELP(27,1,ICUBP0)=
C
1102  IF (.NOT.(BDER(2).OR.BDER(3).OR.BDER(4))) GOTO 99999
C
      IF (.NOT.BDER(2)) GOTO 1103
      DHELP( 1,2,ICUBP0)= Q8*(1D0-2D0*X1)*X2*(1D0-X2)*X3*(1D0-X3)
      DHELP( 2,2,ICUBP0)= Q8*(1D0+2D0*X1)*X2*(1D0-X2)*X3*(1D0-X3)
      DHELP( 3,2,ICUBP0)= Q8*(1D0+2D0*X1)*X2*(1D0+X2)*X3*(1D0-X3)
      DHELP( 4,2,ICUBP0)= Q8*(1D0-2D0*X1)*X2*(1D0+X2)*X3*(1D0-X3)
      DHELP( 5,2,ICUBP0)= Q8*(1D0-2D0*X1)*X2*(1D0-X2)*X3*(1D0+X3)
      DHELP( 6,2,ICUBP0)= Q8*(1D0+2D0*X1)*X2*(1D0-X2)*X3*(1D0+X3)
      DHELP( 7,2,ICUBP0)= Q8*(1D0+2D0*X1)*X2*(1D0+X2)*X3*(1D0+X3)
      DHELP( 8,2,ICUBP0)= Q8*(1D0-2D0*X1)*X2*(1D0+X2)*X3*(1D0+X3)
      DHELP( 9,2,ICUBP0)=-2D0*X1*(1D0-X2*X2)*(1D0-X3*X3)
C      DHELP(10,2,ICUBP0)=
C      DHELP(11,2,ICUBP0)=
C      DHELP(12,2,ICUBP0)=
C      DHELP(13,2,ICUBP0)=
C      DHELP(14,2,ICUBP0)=
C      DHELP(15,2,ICUBP0)=
C      DHELP(16,2,ICUBP0)=
C      DHELP(17,2,ICUBP0)=
C      DHELP(18,2,ICUBP0)=
C      DHELP(19,2,ICUBP0)=
C      DHELP(20,2,ICUBP0)=
C      DHELP(21,2,ICUBP0)=
C      DHELP(22,2,ICUBP0)=
C      DHELP(23,2,ICUBP0)=
C      DHELP(24,2,ICUBP0)=
C      DHELP(25,2,ICUBP0)=
C      DHELP(26,2,ICUBP0)=
C      DHELP(27,2,ICUBP0)=
C
1103  IF (.NOT.BDER(3)) GOTO 1104
      DHELP( 1,3,ICUBP0)= Q8*X1*(1D0-X1)*(1D0-2D0*X2)*X3*(1D0-X3)
      DHELP( 2,3,ICUBP0)= Q8*X1*(1D0+X1)*(1D0-2D0*X2)*X3*(1D0-X3)
      DHELP( 3,3,ICUBP0)= Q8*X1*(1D0+X1)*(1D0+2D0*X2)*X3*(1D0-X3)
      DHELP( 4,3,ICUBP0)= Q8*X1*(1D0-X1)*(1D0+2D0*X2)*X3*(1D0-X3)
      DHELP( 5,3,ICUBP0)= Q8*X1*(1D0-X1)*(1D0-2D0*X2)*X3*(1D0+X3)
      DHELP( 6,3,ICUBP0)= Q8*X1*(1D0+X1)*(1D0-2D0*X2)*X3*(1D0+X3)
      DHELP( 7,3,ICUBP0)= Q8*X1*(1D0+X1)*(1D0+2D0*X2)*X3*(1D0+X3)
      DHELP( 8,3,ICUBP0)= Q8*X1*(1D0-X1)*(1D0+2D0*X2)*X3*(1D0+X3)
      DHELP( 9,3,ICUBP0)=-(1D0-X1*X1)*2D0*X2*(1D0-X3*X3)
C      DHELP(10,3,ICUBP0)=
C      DHELP(11,3,ICUBP0)=
C      DHELP(12,3,ICUBP0)=
C      DHELP(13,3,ICUBP0)=
C      DHELP(14,3,ICUBP0)=
C      DHELP(15,3,ICUBP0)=
C      DHELP(16,3,ICUBP0)=
C      DHELP(17,3,ICUBP0)=
C      DHELP(18,3,ICUBP0)=
C      DHELP(19,3,ICUBP0)=
C      DHELP(20,3,ICUBP0)=
C      DHELP(21,3,ICUBP0)=
C      DHELP(22,3,ICUBP0)=
C      DHELP(23,3,ICUBP0)=
C      DHELP(24,3,ICUBP0)=
C      DHELP(25,3,ICUBP0)=
C      DHELP(26,3,ICUBP0)=
C      DHELP(27,3,ICUBP0)=
C
1104  IF (.NOT.BDER(4)) GOTO 99999
      DHELP( 1,4,ICUBP0)= Q8*X1*(1D0-X1)*X2*(1D0-X2)*(1D0-2D0*X3)
      DHELP( 2,4,ICUBP0)= Q8*X1*(1D0+X1)*X2*(1D0-X2)*(1D0-2D0*X3)
      DHELP( 3,4,ICUBP0)= Q8*X1*(1D0+X1)*X2*(1D0+X2)*(1D0-2D0*X3)
      DHELP( 4,4,ICUBP0)= Q8*X1*(1D0-X1)*X2*(1D0+X2)*(1D0-2D0*X3)
      DHELP( 5,4,ICUBP0)= Q8*X1*(1D0-X1)*X2*(1D0-X2)*(1D0+2D0*X3)
      DHELP( 6,4,ICUBP0)= Q8*X1*(1D0+X1)*X2*(1D0-X2)*(1D0+2D0*X3)
      DHELP( 7,4,ICUBP0)= Q8*X1*(1D0+X1)*X2*(1D0+X2)*(1D0+2D0*X3)
      DHELP( 8,4,ICUBP0)= Q8*X1*(1D0-X1)*X2*(1D0+X2)*(1D0+2D0*X3)
      DHELP( 9,4,ICUBP0)=-(1D0-X1*X1)*(1D0-X2*X2)*2D0*X3
C      DHELP(10,4,ICUBP0)=
C      DHELP(11,4,ICUBP0)=
C      DHELP(12,4,ICUBP0)=
C      DHELP(13,4,ICUBP0)=
C      DHELP(14,4,ICUBP0)=
C      DHELP(15,4,ICUBP0)=
C      DHELP(16,4,ICUBP0)=
C      DHELP(17,4,ICUBP0)=
C      DHELP(18,4,ICUBP0)=
C      DHELP(19,4,ICUBP0)=
C      DHELP(20,4,ICUBP0)=
C      DHELP(21,4,ICUBP0)=
C      DHELP(22,4,ICUBP0)=
C      DHELP(23,4,ICUBP0)=
C      DHELP(24,4,ICUBP0)=
C      DHELP(25,4,ICUBP0)=
C      DHELP(26,4,ICUBP0)=
C      DHELP(27,4,ICUBP0)=
C
99999 END
