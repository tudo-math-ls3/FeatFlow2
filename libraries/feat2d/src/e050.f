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
* E050                                                                 *
*                                                                      *
* Purpose  Calculation of values and derivatives of element  50        *
*          Morley element                                              *
*                                                                      *
* Subroutines/functions called   E050A                                 *
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
* IPAR     I*4    Set to  50  if called with IPAR = -1                 *
* IER      I*4    Error indicator                                      *
*                 -131  Desired derivative not available               *
*                 -132  Triangle has vanishing area                    *
*                 -130  Clockwise ordering of the corner points        *
*                                                                      *
************************************************************************
C
      SUBROUTINE E050(XI1,XI2,XI3,IPAR)
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
      SUB='E050'
      IF (ICHECK.GE.998) CALL OTRC('E050  ','01/02/89')
C
      IER=0
      GOTO (10,1,2,3) -IPAR+1
C
C *** Evaluate at point (XI1,XI2,XI3)
10    ICUBP0=1
      CALL E050A(XI1,XI2,XI3,DHELP,ICUBP0)
      GOTO 4
C
C *** Set number of element
1     IPAR=50
      GOTO 99999
C
C *** Evaluate basis functions on the reference element
C *** for each cubature point
C
2     DO 11 ICUBP0=1,NCUBP
      CALL E050A(DXI(ICUBP0,1),DXI(ICUBP0,2),DXI(ICUBP0,3),
     *            DHELP,ICUBP0)
11    CONTINUE
      GOTO 99999
C
C *** Form linear combination according to actual element coordinates
C
3     ICUBP0=ICUBP
C
4     IF (ICHECK.EQ.0) GOTO 5000
      IF (DETJ.LT.0D0) CALL WERR(-130,'E040  ')
      IF (DETJ.GT.0D0.AND.DETJ.LT.1D-70) CALL WERR(-132,'E040  ')
      IF (IER.NE.0) GOTO 99999
C
5000  XJ1=1D0/DETJ
      AMUL1=1D0
      AMUL2=1D0
      AMUL3=1D0
      IF(KVE(2).GT.KVE(1)) AMUL1=-1D0
      IF(KVE(3).GT.KVE(2)) AMUL2=-1D0
      IF(KVE(1).GT.KVE(3)) AMUL3=-1D0
      A1=1D0/(DJAC(1,1)*DJAC(1,1)+DJAC(2,1)*DJAC(2,1))
      A2=1D0/(DJAC(1,2)*DJAC(1,2)+DJAC(2,2)*DJAC(2,2))
      X32=DJAC(1,2)-DJAC(1,1)
      Y32=DJAC(2,2)-DJAC(2,1)
      A3=1D0/(X32*X32+Y32*Y32)
      A4=DJAC(1,1)*DJAC(1,2)+DJAC(2,1)*DJAC(2,2)
      A5=DJAC(1,1)*X32+DJAC(2,1)*Y32
      A6=DJAC(1,2)*X32+DJAC(2,2)*Y32
      C14= A5*A1
      C16=-A6*A2
      C24=-A4*A1
      C25=-A6*A3
      C35= A5*A3
      C36=-A4*A2
      C44=SQRT(A1)/XJ1
      C55=SQRT(A3)/XJ1
      C66=SQRT(A2)/XJ1
C
C *** Function values
C
      IF (.NOT.BDER(1)) GOTO 201
      P1=DHELP(1,1,ICUBP0)
      P2=DHELP(2,1,ICUBP0)
      P3=DHELP(3,1,ICUBP0)
      P4=DHELP(4,1,ICUBP0)
      P5=DHELP(5,1,ICUBP0)
      P6=DHELP(6,1,ICUBP0)
      DBAS(1,1)=P1+C14*P4+C16*P6
      DBAS(2,1)=P2+C24*P4+C25*P5
      DBAS(3,1)=P3+C35*P5+C36*P6
      DBAS(4,1)=P4*C44*AMUL1
      DBAS(5,1)=P5*C55*AMUL2
      DBAS(6,1)=P6*C66*AMUL3
C
201   IF (.NOT.(BDER(2).OR.BDER(3))) GOTO 203
C
C *** First order derivatives
C
      IF (.NOT.BDER(2)) GOTO 202
      A31= DJAC(2,2)
      A21=-DJAC(2,1)
      P1=-2D0*XI1*(A31+A21)
      P2= 2D0*XI2*A31
      P3= 2D0*XI3*A21
      P4= (2D0*XI3-1D0)*A21
      P5= (1D0-2D0*XI1)*(A31+A21)
      P6= (2D0*XI2-1D0)*A31
      DBAS(1,2)=(P1+C14*P4+C16*P6)*XJ1
      DBAS(2,2)=(P2+C24*P4+C25*P5)*XJ1
      DBAS(3,2)=(P3+C35*P5+C36*P6)*XJ1
      DBAS(4,2)=(P4*C44*AMUL1)*XJ1
      DBAS(5,2)=(P5*C55*AMUL2)*XJ1
      DBAS(6,2)=(P6*C66*AMUL3)*XJ1
202   IF (.NOT.BDER(3)) GOTO 203
      A31=-DJAC(1,2)
      A21= DJAC(1,1)
      P1=-2D0*XI1*(A31+A21)
      P2= 2D0*XI2*A31
      P3= 2D0*XI3*A21
      P4= (2D0*XI3-1D0)*A21
      P5= (1D0-2D0*XI1)*(A31+A21)
      P6= (2D0*XI2-1D0)*A31
      DBAS(1,3)=(P1+C14*P4+C16*P6)*XJ1
      DBAS(2,3)=(P2+C24*P4+C25*P5)*XJ1
      DBAS(3,3)=(P3+C35*P5+C36*P6)*XJ1
      DBAS(4,3)=(P4*C44*AMUL1)*XJ1
      DBAS(5,3)=(P5*C55*AMUL2)*XJ1
      DBAS(6,3)=(P6*C66*AMUL3)*XJ1
C
203   IF(.NOT.(BDER(4).OR.BDER(5).OR.BDER(6))) GOTO 99999
C
C *** Second order derivatives
C
      IF (.NOT.BDER(4)) GOTO 204
      A31= DJAC(2,2)*DJAC(2,2)
      A32=-2D0*DJAC(2,1)*DJAC(2,2)
      A21= DJAC(2,1)*DJAC(2,1)
      P1=2D0*(A31+A32+A21)
      P2=2D0*A31
      P3=2D0*A21
      P4=P3
      P5=P1
      P6=P2
      DBAS(1,4)=(P1+C14*P4+C16*P6)*XJ1*XJ1
      DBAS(2,4)=(P2+C24*P4+C25*P5)*XJ1*XJ1
      DBAS(3,4)=(P3+C35*P5+C36*P6)*XJ1*XJ1
      DBAS(4,4)=(P4*C44*AMUL1)*XJ1*XJ1
      DBAS(5,4)=(P5*C55*AMUL2)*XJ1*XJ1
      DBAS(6,4)=(P6*C66*AMUL3)*XJ1*XJ1
C
204   IF (.NOT.BDER(5)) GOTO 205
      A31=-DJAC(1,2)*DJAC(2,2)
      A32= DJAC(1,2)*DJAC(2,1)+DJAC(1,1)*DJAC(2,2)
      A21=-DJAC(1,1)*DJAC(2,1)
      P1=2D0*(A31+A32+A21)
      P2=2D0*A31
      P3=2D0*A21
      P4=P3
      P5=P1
      P6=P2
      DBAS(1,5)=(P1+C14*P4+C16*P6)*XJ1*XJ1
      DBAS(2,5)=(P2+C24*P4+C25*P5)*XJ1*XJ1
      DBAS(3,5)=(P3+C35*P5+C36*P6)*XJ1*XJ1
      DBAS(4,5)=(P4*C44*AMUL1)*XJ1*XJ1
      DBAS(5,5)=(P5*C55*AMUL2)*XJ1*XJ1
      DBAS(6,5)=(P6*C66*AMUL3)*XJ1*XJ1
205   IF (.NOT.BDER(6)) GOTO 99999
      A31= DJAC(1,2)*DJAC(1,2)
      A32=-2D0*DJAC(1,1)*DJAC(1,2)
      A21= DJAC(1,1)*DJAC(1,1)
      P1=2D0*(A31+A32+A21)
      P2=2D0*A31
      P3=2D0*A21
      P4=P3
      P5=P1
      P6=P2
      DBAS(1,6)=(P1+C14*P4+C16*P6)*XJ1*XJ1
      DBAS(2,6)=(P2+C24*P4+C25*P5)*XJ1*XJ1
      DBAS(3,6)=(P3+C35*P5+C36*P6)*XJ1*XJ1
      DBAS(4,6)=(P4*C44*AMUL1)*XJ1*XJ1
      DBAS(5,6)=(P5*C55*AMUL2)*XJ1*XJ1
      DBAS(6,6)=(P6*C66*AMUL3)*XJ1*XJ1
C
99999 END
C
C
C
      SUBROUTINE E050A(XI1,XI2,XI3,DHELP,ICUBP0)
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
      IF (ICHECK.EQ.999) CALL OTRC('E050A ','01/02/89')
C
      IF (BDER(1)) THEN
       DHELP(1,1,ICUBP0)= XI1*XI1
       DHELP(2,1,ICUBP0)= XI2*XI2
       DHELP(3,1,ICUBP0)= XI3*XI3
       DHELP(4,1,ICUBP0)= XI3*XI3-XI3
       DHELP(5,1,ICUBP0)= XI1*XI1-XI1
       DHELP(6,1,ICUBP0)= XI2*XI2-XI2
      ENDIF
C
      END
