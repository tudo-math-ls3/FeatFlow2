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
* E010                                                                 *
*                                                                      *
* Purpose  Calculation of values of element 10                         *
*          Constant element                                            *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  09/01/89                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* XI1      R*8                                                         *
* XI2      R*8    Evaluation point in cartesian coordinates            *
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
* IPAR     I*4    Set to  10  if called with IPAR = -1                 *
* IER      I*4    Error indicator                                      *
*                 -131  Desired derivative not available               *
*                 -132  Triangle has vanishing area                    *
*                 -130  Clockwise ordering of the corner points        *
*                                                                      *
************************************************************************
C
      SUBROUTINE E010 (XI1,XI2,XI3,IPAR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNDIM=3)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,
     *                DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),
     *                IEL,NDIM
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      SAVE /OUTPUT/,/ERRCTL/,/ELEM/,/CHAR/,/CUB/
C
      SUB='E010'
      IF (ICHECK.GE.998) CALL OTRC('E010  ','09/01/89')
C
      IER=0
      GOTO (10,1,2,3) -IPAR+1
C
1     IPAR=10
2     GOTO 99999
C
3     CONTINUE
10    IF (ICHECK.EQ.0) GOTO 1000
      IF (ABS(DETJ).LT.1D-70) CALL WERR(-132,'E010  ')
C *** No first or second order derivatives available
      IF (BDER(2).OR.BDER(3).OR.BDER(4).OR.BDER(5).OR.BDER(6).OR.
     *    BDER(7).OR.BDER(8).OR.BDER(9).OR.BDER(10))
     *    CALL WERR(-131,'E010  ')
      IF (IER.NE.0) GOTO 99999
C
C *** Function values
1000  IF (BDER(1)) DBAS(1,1,1)=1D0
C
99999 END
