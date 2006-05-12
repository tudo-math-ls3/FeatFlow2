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
* E060                                                                 *
*                                                                      *
* Purpose  Calculation of values and derivatives of element  60        *
*          Linear Hermite type element with 3 dof (triangular element) *
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
* IPAR     I*4    Set to  60  if called with IPAR = -1                 *
* IER      I*4    Error indicator                                      *
*                 -131  Desired derivative not available               *
*                 -132  Triangle has vanishing area                    *
*                 -130  Clockwise ordering of the corner points        *
*                                                                      *
************************************************************************
C
      SUBROUTINE E060 (XI1,XI2,XI3,IPAR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4)
      PARAMETER (Q3=1D0/3D0)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      SAVE /OUTPUT/,/ERRCTL/,/ELEM/,/CHAR/,/CUB/
C
      SUB='E060'
      IF (ICHECK.GE.998) CALL OTRC('E060  ','03/21/89')
C
      IER=0
      GOTO (10,1,2,3) -IPAR+1
C
1     IPAR=60
2     GOTO 99999
C
3     CONTINUE
10    IF (ICHECK.EQ.0) GOTO 6000
      IF (DETJ.LT.0D0) CALL WERR(-130,'E060  ')
      IF (DETJ.GT.0D0.AND.DETJ.LT.1D-70) CALL WERR(-132,'E060  ')
C *** No second order derivatives available
C *** Used for second order problems only
      IF (BDER(4).OR.BDER(5).OR.BDER(6)) CALL WERR(-131,'E060  ')
      IF (IER.NE.0) GOTO 99999
C
C *** Function values
6000  IF (.NOT.BDER(1)) GOTO 6002
      DBAS(1,1)= 1D0
      DBAS(2,1)= XI2-Q3
      DBAS(3,1)= XI3-Q3
C
C *** First order derivatives
6002  IF (.NOT.(BDER(2).OR.BDER(3))) GOTO 99999
      XJ1=1D0/DETJ
      IF (.NOT.BDER(2)) GOTO 6003
      DBAS(1,2)=0D0
      DBAS(2,2)= DJAC(2,2)*XJ1
      DBAS(3,2)=-DJAC(2,1)*XJ1
6003  IF (.NOT.BDER(3)) GOTO 99999
      DBAS(1,3)=0D0
      DBAS(2,3)=-DJAC(1,2)*XJ1
      DBAS(3,3)= DJAC(1,1)*XJ1
C
99999 END
