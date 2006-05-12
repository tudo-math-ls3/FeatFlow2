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
* E023                                                                 *
*                                                                      *
* Purpose  Calculation of values and derivatives of element  23        *
*          Augmented quadratic Lagrangian element (P2+bulb)            *
*                                                                      *
* Subroutines/functions called   E023A                                 *
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
      SUBROUTINE E023(XI1,XI2,XI3,IPAR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NDOFL=7)
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
      IF (ICHECK.GE.998) CALL OTRC('E002  ','03/21/89')
C
      IER=0
      GOTO (10,1,2,3) -IPAR+1
C
C *** Evaluate at point (XI1,XI2,XI3)
10    ICUBP0=1
      CALL E023A(XI1,XI2,XI3,DHELP,ICUBP0)
      GOTO 4
C
C *** Set number of element
1     IPAR=23
      GOTO 99999
C
C *** Evaluate basis functions on the reference element
C *** for each cubature point
C
2     DO 11 ICUBP0=1,NCUBP
      CALL E023A(DXI(ICUBP0,1),DXI(ICUBP0,2),DXI(ICUBP0,3),
     *            DHELP,ICUBP0)
11    CONTINUE
      GOTO 99999
C
C *** Form linear combination according to actual element coordinates
C
3     ICUBP0=ICUBP
C
4     IF (ICHECK.EQ.0) GOTO 2300
      IF (DETJ.LT.0D0) CALL WERR(-130,'E023  ')
      IF (DETJ.GT.0D0.AND.DETJ.LT.1D-70) CALL WERR(-132,'E023  ')
C *** No second order derivatives available
C *** Used for second order problems only
      IF (BDER(4).OR.BDER(5).OR.BDER(6)) CALL WERR(-131,'E023  ')
      IF (IER.NE.0) GOTO 99999
C
C *** Function values
2300  IF (BDER(1)) THEN
       DO 2301 IDOFL=1,NDOFL
2301   DBAS(IDOFL,1)=DHELP(IDOFL,1,ICUBP0)
      ENDIF
C
C***  First order derivatives
      XJ1=1D0/DETJ
C
      IF (BDER(2)) THEN
       DO 2302 IDOFL=1,NDOFL
2302   DBAS(IDOFL,2)=(DJAC(2,2)*DHELP(IDOFL,2,ICUBP0)
     *               -DJAC(2,1)*DHELP(IDOFL,3,ICUBP0))*XJ1
      ENDIF
C
      IF (BDER(3)) THEN
       DO 2303 IDOFL=1,NDOFL
2303   DBAS(IDOFL,3)=(-DJAC(1,2)*DHELP(IDOFL,2,ICUBP0)
     *                +DJAC(1,1)*DHELP(IDOFL,3,ICUBP0))*XJ1
      ENDIF
C
99999 END
C
C
C
      SUBROUTINE E023A(XI1,XI2,XI3,DHELP,ICUBP0)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNBAS=21,NNDER=6,NNVE=4,NDOFL=7)
      DIMENSION DHELP(NDOFL,3,*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      SAVE /OUTPUT/,/ERRCTL/,/ELEM/
C
      IF (ICHECK.EQ.999) CALL OTRC('E023A ','03/21/89')
C
      IF (BDER(1)) THEN
       DHELP(1,1,ICUBP0)= XI1*(XI1-XI2-XI3)
       DHELP(2,1,ICUBP0)= XI2*(XI2-XI1-XI3)
       DHELP(3,1,ICUBP0)= XI3*(XI3-XI1-XI2)
       DHELP(4,1,ICUBP0)= 4D0*XI1*XI2
       DHELP(5,1,ICUBP0)= 4D0*XI2*XI3
       DHELP(6,1,ICUBP0)= 4D0*XI1*XI3
       DHELP(7,1,ICUBP0)= 27D0*XI1*XI2*XI3
      ENDIF
C
      IF (.NOT.BDER(1).AND..NOT.BDER(2)) GOTO 99999
C
       DHELP(1,2,ICUBP0)= 1D0-4D0*XI1
       DHELP(2,2,ICUBP0)= 4D0*XI2-1D0
       DHELP(3,2,ICUBP0)= 0D0
       DHELP(4,2,ICUBP0)= 4D0*(XI1-XI2)
       DHELP(5,2,ICUBP0)= 4D0*XI3
       DHELP(6,2,ICUBP0)=-4D0*XI3
       DHELP(7,2,ICUBP0)=27D0*XI3*(XI1-XI2)
C
       DHELP(1,3,ICUBP0)= 1D0-4D0*XI1
       DHELP(2,3,ICUBP0)= 0D0
       DHELP(3,3,ICUBP0)= 4D0*XI3-1D0
       DHELP(4,3,ICUBP0)=-4D0*XI2
       DHELP(5,3,ICUBP0)= 4D0*XI2
       DHELP(6,3,ICUBP0)= 4D0*(XI1-XI3)
       DHELP(7,3,ICUBP0)=27D0*(XI1-XI3)*XI2
C
99999 END
