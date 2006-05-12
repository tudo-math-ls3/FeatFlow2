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
* E014                                                                 *
*                                                                      *
* Purpose  Calculation of values and derivatives of element  14        *
*          Bicubic element (16 nodes)                                  *
*                                                                      *
* Subroutines/functions called   E014A                                 *
*                                                                      *
* Version from  10/27/89                                               *
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
* IPAR     I*4    Set to  14  if called with IPAR = -1                 *
* IER      I*4    Error indicator                                      *
*                 -131  Desired derivative not available               *
*                 -132  Triangle has vanishing area                    *
*                 -130  Clockwise ordering of the corner points        *
*                                                                      *
************************************************************************
C
      SUBROUTINE E014 (XI1,XI2,IPAR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NDOFL=16)
      DIMENSION DHELP(NDOFL,NNDER,NNCUBP)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      SAVE /OUTPUT/,/ERRCTL/,/ELEM/,/CHAR/,/CUB/
      SAVE DHELP
C
      SUB='E014'
      IF (ICHECK.GE.998) CALL OTRC('E014  ','10/27/89')
C
      IER=0
      GOTO (10,1,2,3) -IPAR+1
C
C *** Evaluate at point (XI1,XI2)
10    ICUBP0=1
      CALL E014A(XI1,XI2,DHELP,ICUBP0)
      GOTO 4
C
C *** Set number of element
1     IPAR=14
      GOTO 99999
C
C *** Evaluate basis functions on the reference element
C *** for each cubature point
2     DO 11 ICUBP0=1,NCUBP
      CALL E014A(DXI(ICUBP0,1),DXI(ICUBP0,2),DHELP,ICUBP0)
11    CONTINUE
      GOTO 99999
C
C *** Form linear combinations of the values in DHELP
C *** for the actual element
3     ICUBP0=ICUBP
C
4     IF (ICHECK.EQ.0) GOTO 1400
      IF (DETJ.LT.0D0) CALL WERR(-130,'E014  ')
      IF (DETJ.GT.0D0.AND.DETJ.LT.1D-70) CALL WERR(-132,'E014  ')
      IF (IER.NE.0) GOTO 99999
C
1400  B1=BDER(2).OR.BDER(3)
      B2=BDER(4).OR.BDER(5).OR.BDER(6)
      IF (BDER(1)) THEN
       DO 1401 IDFL=1,NDOFL
1401   DBAS(IDFL,1)=DHELP(IDFL,1,ICUBP0)
      ENDIF
      IF (.NOT.(B1.OR.B2)) GOTO 99999
      XJ1=1D0/DETJ
      IF (BDER(2)) THEN
       C1= DJAC(2,2)*XJ1
       C2=-DJAC(2,1)*XJ1
       DO 1402 IDFL=1,NDOFL
1402   DBAS(IDFL,2)=DHELP(IDFL,2,ICUBP0)*C1
     *             +DHELP(IDFL,3,ICUBP0)*C2
      ENDIF
      IF (BDER(3)) THEN
       C1=-DJAC(1,2)*XJ1
       C2= DJAC(1,1)*XJ1
       DO 1403 IDFL=1,NDOFL
1403   DBAS(IDFL,3)=DHELP(IDFL,2,ICUBP0)*C1
     *             +DHELP(IDFL,3,ICUBP0)*C2
      ENDIF
      IF (.NOT.B2) GOTO 99999
      XJ2=XJ1*XJ1
      IF (BDER(4)) THEN
       C1=DJAC(2,2)*DJAC(2,2)*XJ2
       C2=-2D0*DJAC(2,2)*DJAC(2,1)*XJ2
       C3=DJAC(2,1)*DJAC(2,1)*XJ2
       DO 1404 IDFL=1,NDOFL
1404   DBAS(IDFL,4)=DHELP(IDFL,4,ICUBP0)*C1
     *             +DHELP(IDFL,5,ICUBP0)*C2
     *             +DHELP(IDFL,6,ICUBP0)*C3
      ENDIF
      IF (BDER(5)) THEN
       C1=-DJAC(2,2)*DJAC(1,2)*XJ2
       C2=(DJAC(2,2)*DJAC(1,1)+DJAC(2,1)*DJAC(1,2))*XJ2
       C3=-DJAC(2,1)*DJAC(1,1)*XJ2
       DO 1405 IDFL=1,NDOFL
1405   DBAS(IDFL,5)=DHELP(IDFL,4,ICUBP0)*C1
     *             +DHELP(IDFL,5,ICUBP0)*C2
     *             +DHELP(IDFL,6,ICUBP0)*C3
      ENDIF
      IF (BDER(6)) THEN
       C1=DJAC(1,2)*DJAC(1,2)*XJ2
       C2=-2D0*DJAC(1,2)*DJAC(1,1)*XJ2
       C3=DJAC(1,1)*DJAC(1,1)*XJ2
       DO 1406 IDFL=1,NDOFL
1406   DBAS(IDFL,6)=DHELP(IDFL,4,ICUBP0)*C1
     *             +DHELP(IDFL,5,ICUBP0)*C2
     *             +DHELP(IDFL,6,ICUBP0)*C3
      ENDIF
C
99999 END
C
C
C
      SUBROUTINE E014A(X1,X2,DHELP,ICUBP0)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNBAS=21,NNDER=6,NNVE=4,NDOFL=16)
      PARAMETER (Q81=81/256D0,Q3=1D0/3D0,Q243=243D0/256D0,
     *           Q729=729D0/256D0,Q9=1D0/9D0,Q32=2D0/3D0)
      DIMENSION DHELP(NDOFL,NNDER,*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      SAVE /OUTPUT/,/ERRCTL/,/ELEM/
C
      IF (ICHECK.EQ.999) CALL OTRC('E014A ','10/27/89')
C
      IF (BDER(1)) THEN
       DHELP(1,1,ICUBP0) = Q81*(X1+Q3)*(X1-Q3)*(X1-1D0)
     *                        *(X2+Q3)*(X2-Q3)*(X2-1D0)
       DHELP(2,1,ICUBP0) =-Q81*(X1+1D0)*(X1+Q3)*(X1-Q3)
     *                        *(X2+Q3)*(X2-Q3)*(X2-1D0)
       DHELP(3,1,ICUBP0) = Q81*(X1+1D0)*(X1+Q3)*(X1-Q3)
     *                        *(X2+1D0)*(X2+Q3)*(X2-Q3)
       DHELP(4,1,ICUBP0) =-Q81*(X1+Q3)*(X1-Q3)*(X1-1D0)
     *                        *(X2+1D0)*(X2+Q3)*(X2-Q3)
       DHELP(5,1,ICUBP0) =-Q243*(X1+1D0)*(X1-Q3)*(X1-1D0)
     *                         *(X2+Q3)*(X2-Q3)*(X2-1D0)
       DHELP(6,1,ICUBP0) = Q243*(X1+1D0)*(X1+Q3)*(X1-1D0)
     *                         *(X2+Q3)*(X2-Q3)*(X2-1D0)
       DHELP(7,1,ICUBP0) = Q243*(X1+1D0)*(X1+Q3)*(X1-Q3)
     *                         *(X2+1D0)*(X2-Q3)*(X2-1D0)
       DHELP(8,1,ICUBP0) =-Q243*(X1+1D0)*(X1+Q3)*(X1-Q3)
     *                         *(X2+1D0)*(X2+Q3)*(X2-1D0)
       DHELP(9,1,ICUBP0) =-Q243*(X1+1D0)*(X1+Q3)*(X1-1D0)
     *                         *(X2+1D0)*(X2+Q3)*(X2-Q3)
       DHELP(10,1,ICUBP0)= Q243*(X1+1D0)*(X1-Q3)*(X1-1D0)
     *                         *(X2+1D0)*(X2+Q3)*(X2-Q3)
       DHELP(11,1,ICUBP0)= Q243*(X1+Q3)*(X1-Q3)*(X1-1D0)
     *                         *(X2+1D0)*(X2+Q3)*(X2-1D0)
       DHELP(12,1,ICUBP0)=-Q243*(X1+Q3)*(X1-Q3)*(X1-1D0)
     *                         *(X2+1D0)*(X2-Q3)*(X2-1D0)
       DHELP(13,1,ICUBP0)= Q729*(X1+1D0)*(X1-Q3)*(X1-1D0)
     *                         *(X2+1D0)*(X2-Q3)*(X2-1D0)
       DHELP(14,1,ICUBP0)=-Q729*(X1+1D0)*(X1+Q3)*(X1-1D0)
     *                         *(X2+1D0)*(X2-Q3)*(X2-1D0)
       DHELP(15,1,ICUBP0)= Q729*(X1+1D0)*(X1+Q3)*(X1-1D0)
     *                         *(X2+1D0)*(X2+Q3)*(X2-1D0)
       DHELP(16,1,ICUBP0)=-Q729*(X1+1D0)*(X1-Q3)*(X1-1D0)
     *                         *(X2+1D0)*(X2+Q3)*(X2-1D0)
      ENDIF
      IF (BDER(2).OR.BDER(3)) THEN
       DHELP(1,2,ICUBP0) = Q81*((3D0*X1-2D0)*X1-Q9)
     *                        *(X2+Q3)*(X2-Q3)*(X2-1D0)
       DHELP(2,2,ICUBP0) =-Q81*((3D0*X1+2D0)*X1-Q9)
     *                        *(X2+Q3)*(X2-Q3)*(X2-1D0)
       DHELP(3,2,ICUBP0) = Q81*((3D0*X1+2D0)*X1-Q9)
     *                        *(X2+1D0)*(X2+Q3)*(X2-Q3)
       DHELP(4,2,ICUBP0) =-Q81*((3D0*X1-2D0)*X1-Q9)
     *                        *(X2+1D0)*(X2+Q3)*(X2-Q3)
       DHELP(5,2,ICUBP0) =-Q243*((3D0*X1-Q32)*X1-1D0)
     *                         *(X2+Q3)*(X2-Q3)*(X2-1D0)
       DHELP(6,2,ICUBP0) = Q243*((3D0*X1+Q32)*X1-1D0)
     *                         *(X2+Q3)*(X2-Q3)*(X2-1D0)
       DHELP(7,2,ICUBP0) = Q243*((3D0*X1+2D0)*X1-Q9)
     *                         *(X2+1D0)*(X2-Q3)*(X2-1D0)
       DHELP(8,2,ICUBP0) =-Q243*((3D0*X1+2D0)*X1-Q9)
     *                         *(X2+1D0)*(X2+Q3)*(X2-1D0)
       DHELP(9,2,ICUBP0) =-Q243*((3D0*X1+Q32)*X1-1D0)
     *                         *(X2+1D0)*(X2+Q3)*(X2-Q3)
       DHELP(10,2,ICUBP0)= Q243*((3D0*X1-Q32)*X1-1D0)
     *                         *(X2+1D0)*(X2+Q3)*(X2-Q3)
       DHELP(11,2,ICUBP0)= Q243*((3D0*X1-2D0)*X1-Q9)
     *                         *(X2+1D0)*(X2+Q3)*(X2-1D0)
       DHELP(12,2,ICUBP0)=-Q243*((3D0*X1-2D0)*X1-Q9)
     *                         *(X2+1D0)*(X2-Q3)*(X2-1D0)
       DHELP(13,2,ICUBP0)= Q729*((3D0*X1-Q32)*X1-1D0)
     *                         *(X2+1D0)*(X2-Q3)*(X2-1D0)
       DHELP(14,2,ICUBP0)=-Q729*((3D0*X1+Q32)*X1-1D0)
     *                         *(X2+1D0)*(X2-Q3)*(X2-1D0)
       DHELP(15,2,ICUBP0)= Q729*((3D0*X1+Q32)*X1-1D0)
     *                         *(X2+1D0)*(X2+Q3)*(X2-1D0)
       DHELP(16,2,ICUBP0)=-Q729*((3D0*X1-Q32)*X1-1D0)
     *                         *(X2+1D0)*(X2+Q3)*(X2-1D0)
C
       DHELP(1,3,ICUBP0) = Q81*(X1+Q3)*(X1-Q3)*(X1-1D0)
     *                        *((3D0*X2-2D0)*X2-Q9)
       DHELP(2,3,ICUBP0) =-Q81*(X1+1D0)*(X1+Q3)*(X1-Q3)
     *                        *((3D0*X2-2D0)*X2-Q9)
       DHELP(3,3,ICUBP0) = Q81*(X1+1D0)*(X1+Q3)*(X1-Q3)
     *                        *((3D0*X2+2D0)*X2-Q9)
       DHELP(4,3,ICUBP0) =-Q81*(X1+Q3)*(X1-Q3)*(X1-1D0)
     *                        *((3D0*X2+2D0)*X2-Q9)
       DHELP(5,3,ICUBP0) =-Q243*(X1+1D0)*(X1-Q3)*(X1-1D0)
     *                         *((3D0*X2-2D0)*X2-Q9)
       DHELP(6,3,ICUBP0) = Q243*(X1+1D0)*(X1+Q3)*(X1-1D0)
     *                         *((3D0*X2-2D0)*X2-Q9)
       DHELP(7,3,ICUBP0) = Q243*(X1+1D0)*(X1+Q3)*(X1-Q3)
     *                         *((3D0*X2-Q32)*X2-1D0)
       DHELP(8,3,ICUBP0) =-Q243*(X1+1D0)*(X1+Q3)*(X1-Q3)
     *                         *((3D0*X2+Q32)*X2-1D0)
       DHELP(9,3,ICUBP0) =-Q243*(X1+1D0)*(X1+Q3)*(X1-1D0)
     *                         *((3D0*X2+2D0)*X2-Q9)
       DHELP(10,3,ICUBP0)= Q243*(X1+1D0)*(X1-Q3)*(X1-1D0)
     *                         *((3D0*X2+2D0)*X2-Q9)
       DHELP(11,3,ICUBP0)= Q243*(X1+Q3)*(X1-Q3)*(X1-1D0)
     *                         *((3D0*X2+Q32)*X2-1D0)
       DHELP(12,3,ICUBP0)=-Q243*(X1+Q3)*(X1-Q3)*(X1-1D0)
     *                         *((3D0*X2-Q32)*X2-1D0)
       DHELP(13,3,ICUBP0)= Q729*(X1+1D0)*(X1-Q3)*(X1-1D0)
     *                         *((3D0*X2-Q32)*X2-1D0)
       DHELP(14,3,ICUBP0)=-Q729*(X1+1D0)*(X1+Q3)*(X1-1D0)
     *                         *((3D0*X2-Q32)*X2-1D0)
       DHELP(15,3,ICUBP0)= Q729*(X1+1D0)*(X1+Q3)*(X1-1D0)
     *                         *((3D0*X2+Q32)*X2-1D0)
       DHELP(16,3,ICUBP0)=-Q729*(X1+1D0)*(X1-Q3)*(X1-1D0)
     *                         *((3D0*X2+Q32)*X2-1D0)
      ENDIF
C
      IF (BDER(4).OR.BDER(5).OR.BDER(6)) THEN
       DHELP(1,4,ICUBP0) = Q81*(6D0*X1-2D0)
     *                        *(X2+Q3)*(X2-Q3)*(X2-1D0)
       DHELP(2,4,ICUBP0) =-Q81*(6D0*X1+2D0)
     *                        *(X2+Q3)*(X2-Q3)*(X2-1D0)
       DHELP(3,4,ICUBP0) = Q81*(6D0*X1+2D0)
     *                        *(X2+1D0)*(X2+Q3)*(X2-Q3)
       DHELP(4,4,ICUBP0) =-Q81*(6D0*X1-2D0)
     *                        *(X2+1D0)*(X2+Q3)*(X2-Q3)
       DHELP(5,4,ICUBP0) =-Q243*(6D0*X1-Q32)
     *                         *(X2+Q3)*(X2-Q3)*(X2-1D0)
       DHELP(6,4,ICUBP0) = Q243*(6D0*X1+Q32)
     *                         *(X2+Q3)*(X2-Q3)*(X2-1D0)
       DHELP(7,4,ICUBP0) = Q243*(6D0*X1+2D0)
     *                         *(X2+1D0)*(X2-Q3)*(X2-1D0)
       DHELP(8,4,ICUBP0) =-Q243*(6D0*X1+2D0)
     *                         *(X2+1D0)*(X2+Q3)*(X2-1D0)
       DHELP(9,4,ICUBP0) =-Q243*(6D0*X1+Q32)
     *                         *(X2+1D0)*(X2+Q3)*(X2-Q3)
       DHELP(10,4,ICUBP0)= Q243*(6D0*X1-Q32)
     *                         *(X2+1D0)*(X2+Q3)*(X2-Q3)
       DHELP(11,4,ICUBP0)= Q243*(6D0*X1-2D0)
     *                         *(X2+1D0)*(X2+Q3)*(X2-1D0)
       DHELP(12,4,ICUBP0)=-Q243*(6D0*X1-2D0)
     *                         *(X2+1D0)*(X2-Q3)*(X2-1D0)
       DHELP(13,4,ICUBP0)= Q729*(6D0*X1-Q32)
     *                         *(X2+1D0)*(X2-Q3)*(X2-1D0)
       DHELP(14,4,ICUBP0)=-Q729*(6D0*X1+Q32)
     *                         *(X2+1D0)*(X2-Q3)*(X2-1D0)
       DHELP(15,4,ICUBP0)= Q729*(6D0*X1+Q32)
     *                         *(X2+1D0)*(X2+Q3)*(X2-1D0)
       DHELP(16,4,ICUBP0)=-Q729*(6D0*X1-Q32)
     *                         *(X2+1D0)*(X2+Q3)*(X2-1D0)
C
       DHELP(1,5,ICUBP0) = Q81*((3D0*X1-2D0)*X1-Q9)
     *                        *((3D0*X2-2D0)*X2-Q9)
       DHELP(2,5,ICUBP0) =-Q81*((3D0*X1+2D0)*X1-Q9)
     *                        *((3D0*X2-2D0)*X2-Q9)
       DHELP(3,5,ICUBP0) = Q81*((3D0*X1+2D0)*X1-Q9)
     *                        *((3D0*X2+2D0)*X2-Q9)
       DHELP(4,5,ICUBP0) =-Q81*((3D0*X1-2D0)*X1-Q9)
     *                        *((3D0*X2+2D0)*X2-Q9)
       DHELP(5,5,ICUBP0) =-Q243*((3D0*X1-Q32)*X1-1D0)
     *                         *((3D0*X2-2D0)*X2-Q9)
       DHELP(6,5,ICUBP0) = Q243*((3D0*X1+Q32)*X1-1D0)
     *                         *((3D0*X2-2D0)*X2-Q9)
       DHELP(7,5,ICUBP0) = Q243*((3D0*X1+2D0)*X1-Q9)
     *                         *((3D0*X2-Q32)*X2-1D0)
       DHELP(8,5,ICUBP0) =-Q243*((3D0*X1+2D0)*X1-Q9)
     *                         *((3D0*X2+Q32)*X2-1D0)
       DHELP(9,5,ICUBP0) =-Q243*((3D0*X1+Q32)*X1-1D0)
     *                         *((3D0*X2+2D0)*X2-Q9)
       DHELP(10,5,ICUBP0)= Q243*((3D0*X1-Q32)*X1-1D0)
     *                         *((3D0*X2+2D0)*X2-Q9)
       DHELP(11,5,ICUBP0)= Q243*((3D0*X1-2D0)*X1-Q9)
     *                         *((3D0*X2+Q32)*X2-1D0)
       DHELP(12,5,ICUBP0)=-Q243*((3D0*X1-2D0)*X1-Q9)
     *                         *((3D0*X2-Q32)*X2-1D0)
       DHELP(13,5,ICUBP0)= Q729*((3D0*X1-Q32)*X1-1D0)
     *                         *((3D0*X2-Q32)*X2-1D0)
       DHELP(14,5,ICUBP0)=-Q729*((3D0*X1+Q32)*X1-1D0)
     *                         *((3D0*X2-Q32)*X2-1D0)
       DHELP(15,5,ICUBP0)= Q729*((3D0*X1+Q32)*X1-1D0)
     *                         *((3D0*X2+Q32)*X2-1D0)
       DHELP(16,5,ICUBP0)=-Q729*((3D0*X1-Q32)*X1-1D0)
     *                         *((3D0*X2+Q32)*X2-1D0)
C
       DHELP(1,6,ICUBP0) = Q81*(X1+Q3)*(X1-Q3)*(X1-1D0)
     *                        *(6D0*X2-2D0)
       DHELP(2,6,ICUBP0) =-Q81*(X1+1D0)*(X1+Q3)*(X1-Q3)
     *                        *(6D0*X2-2D0)
       DHELP(3,6,ICUBP0) = Q81*(X1+1D0)*(X1+Q3)*(X1-Q3)
     *                        *(6D0*X2+2D0)
       DHELP(4,6,ICUBP0) =-Q81*(X1+Q3)*(X1-Q3)*(X1-1D0)
     *                        *(6D0*X2+2D0)
       DHELP(5,6,ICUBP0) =-Q243*(X1+1D0)*(X1-Q3)*(X1-1D0)
     *                         *(6D0*X2-2D0)
       DHELP(6,6,ICUBP0) = Q243*(X1+1D0)*(X1+Q3)*(X1-1D0)
     *                         *(6D0*X2-2D0)
       DHELP(7,6,ICUBP0) = Q243*(X1+1D0)*(X1+Q3)*(X1-Q3)
     *                         *(6D0*X2-Q32)
       DHELP(8,6,ICUBP0) =-Q243*(X1+1D0)*(X1+Q3)*(X1-Q3)
     *                         *(6D0*X2+Q32)
       DHELP(9,6,ICUBP0) =-Q243*(X1+1D0)*(X1+Q3)*(X1-1D0)
     *                         *(6D0*X2+2D0)
       DHELP(10,6,ICUBP0)= Q243*(X1+1D0)*(X1-Q3)*(X1-1D0)
     *                         *(6D0*X2+2D0)
       DHELP(11,6,ICUBP0)= Q243*(X1+Q3)*(X1-Q3)*(X1-1D0)
     *                         *(6D0*X2+Q32)
       DHELP(12,6,ICUBP0)=-Q243*(X1+Q3)*(X1-Q3)*(X1-1D0)
     *                         *(6D0*X2-Q32)
       DHELP(13,6,ICUBP0)= Q729*(X1+1D0)*(X1-Q3)*(X1-1D0)
     *                         *(6D0*X2-Q32)
       DHELP(14,6,ICUBP0)=-Q729*(X1+1D0)*(X1+Q3)*(X1-1D0)
     *                         *(6D0*X2-Q32)
       DHELP(15,6,ICUBP0)= Q729*(X1+1D0)*(X1+Q3)*(X1-1D0)
     *                         *(6D0*X2+Q32)
       DHELP(16,6,ICUBP0)=-Q729*(X1+1D0)*(X1-Q3)*(X1-1D0)
     *                         *(6D0*X2+Q32)
      ENDIF
C
      END
