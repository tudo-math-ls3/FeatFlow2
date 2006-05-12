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
* SBV                                                                  *
*                                                                      *
* Purpose  Adjust dimensions of DCORVG, KVERT, KEDGE, KAREA, KADJ      *
*          and KNPR                                                    *
*                                                                      *
* Subroutines/functions called  none                                   *
*                                                                      *
* Version from  12/11/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* IDISP    I*4    =1 Release free space on the arrays                  *
*                    after determination of the new subdivision        *
* For the description of the remaining parameters see SB0              *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
*                                                                      *
************************************************************************
C
      SUBROUTINE SBV(DCORVG,KVERT,KEDGE,KADJ,KNPR,PARX,PARY,PARZ,SEDB,
     *               SADB)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNVE=8,NNAE=6,NNEE=12)
      DIMENSION DCORVG(3,*)
      DIMENSION KVERT(NNVE,*),KADJ(NNAE,*),KEDGE(NNEE,*),KNPR(*)
      DIMENSION KEL(8),KIEL(4),KIELH(4),KIV(2,12),KIAD(4,6)
      DATA KIV /1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
      DATA KIAD/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      EXTERNAL PARX,PARY,PARZ,SEDB,SADB
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='SBV'
      IF (ICHECK.GE.997) CALL OTRC('SBV  ','12/11/89')
C
C
C
      NELOLD=NEL
      NVTOLD=NVT
      NETOLD=NET
      NATOLD=NAT
C
      NVT=NVTOLD+NETOLD
      NAT=0
C
      ICEDGE=0
C
C
      DO 110 IEL=1,NELOLD
C
C Calculate coordinates of new vertices on old edges
C
      DO 120 IED=1,NEE
      IEDGE=KEDGE(IED,IEL)
C
      IF (IEDGE.GT.ICEDGE) THEN
       ICEDGE=ICEDGE+1
      ELSE
       GOTO 120
      ENDIF
C
      IVE1=KIV(1,IED)
      IVE2=KIV(2,IED)
C
      IVT1=KVERT(IVE1,IEL)
      IBCT1=KNPR(IVT1)
      IVT2=KVERT(IVE2,IEL)
      IBCT2=KNPR(IVT2)
C
      IF (IBCT1.EQ.0) THEN
       X1=DCORVG(1,IVT1)
       Y1=DCORVG(2,IVT1)
       Z1=DCORVG(3,IVT1)
      ELSE
       X1=PARX(DCORVG(1,IVT1),DCORVG(2,IVT1),DCORVG(3,IVT1),IBCT1)
       Y1=PARY(DCORVG(1,IVT1),DCORVG(2,IVT1),DCORVG(3,IVT1),IBCT1)
       Z1=PARZ(DCORVG(1,IVT1),DCORVG(2,IVT1),DCORVG(3,IVT1),IBCT1)
      ENDIF
C
      IF (IBCT2.EQ.0) THEN
       X2=DCORVG(1,IVT2)
       Y2=DCORVG(2,IVT2)
       Z2=DCORVG(3,IVT2)
      ELSE
       X2=PARX(DCORVG(1,IVT2),DCORVG(2,IVT2),DCORVG(3,IVT2),IBCT2)
       Y2=PARY(DCORVG(1,IVT2),DCORVG(2,IVT2),DCORVG(3,IVT2),IBCT2)
       Z2=PARZ(DCORVG(1,IVT2),DCORVG(2,IVT2),DCORVG(3,IVT2),IBCT2)
      ENDIF
C
      IF ((IBCT1.EQ.0).OR.(IBCT2.EQ.0).OR.(IBCT1.NE.IBCT2)) THEN
       KNPR(NVTOLD+IEDGE)=0
       DCORVG(1,NVTOLD+IEDGE)=0.5D0*(X1+X2)
       DCORVG(2,NVTOLD+IEDGE)=0.5D0*(Y1+Y2)
       DCORVG(3,NVTOLD+IEDGE)=0.5D0*(Z1+Z2)
      ELSE
       KNPR(NVTOLD+IEDGE)=IBCT1
       CALL SEDB(DCORVG(1,IVT1),DCORVG(2,IVT1),DCORVG(3,IVT1),
     *           DCORVG(1,IVT2),DCORVG(2,IVT2),DCORVG(3,IVT2),
     *           PAR1,PAR2,PAR3,IBCT1)
       DCORVG(1,NVTOLD+IEDGE)=PAR1
       DCORVG(2,NVTOLD+IEDGE)=PAR2
       DCORVG(3,NVTOLD+IEDGE)=PAR3
      ENDIF
C
120   CONTINUE
C
C Actualization of KVERT on old edges
C
      KEL(1)=IEL
      DO 130 I=1,7
130   KEL(I+1)=NEL+I
C
      DO 131 I=1,8
131   KVERT(1,KEL(I))=KVERT(I,IEL)
C
      DO 140 I1=1,4
      I2=MOD(I1,4)+1
      IEDGE=KEDGE(I1,IEL)
      KVERT(2,KEL(I1))=IEDGE+NVTOLD
      KVERT(4,KEL(I2))=IEDGE+NVTOLD
140   CONTINUE
C
      DO 141 I1=1,4
      I2=I1+4
      IEDGE=KEDGE(I1+4,IEL)
      KVERT(5,KEL(I1))=IEDGE+NVTOLD
      KVERT(5,KEL(I2))=IEDGE+NVTOLD
141   CONTINUE
C
      DO 142 I1=1,4
      I2=MOD(I1,4)+1
      IEDGE=KEDGE(I1+8,IEL)
      KVERT(2,KEL(I1+4))=IEDGE+NVTOLD
      KVERT(4,KEL(I2+4))=IEDGE+NVTOLD
142   CONTINUE
C
C Calculation of coordinates and numbers of midpoints of areas and
C of neighbours
C
      DO 150 IAR=1,NAE
      INEIGH=KADJ(IAR,IEL)
      IF ((INEIGH.EQ.0).OR.(INEIGH.GE.IEL)) THEN
C
       NVT=NVT+1
       NAT=NAT+4
C
       IVT1=KVERT(1,KEL(KIAD(1,IAR)))
       IVT2=KVERT(1,KEL(KIAD(2,IAR)))
       IVT3=KVERT(1,KEL(KIAD(3,IAR)))
       IVT4=KVERT(1,KEL(KIAD(4,IAR)))
       IBCT1=KNPR(IVT1)
       IBCT2=KNPR(IVT2)
       IBCT3=KNPR(IVT3)
       IBCT4=KNPR(IVT4)
C
       IF (IBCT1.EQ.0) THEN
        X1=DCORVG(1,IVT1)
        Y1=DCORVG(2,IVT1)
        Z1=DCORVG(3,IVT1)
       ELSE
        X1=PARX(DCORVG(1,IVT1),DCORVG(2,IVT1),DCORVG(3,IVT1),IBCT1)
        Y1=PARY(DCORVG(1,IVT1),DCORVG(2,IVT1),DCORVG(3,IVT1),IBCT1)
        Z1=PARZ(DCORVG(1,IVT1),DCORVG(2,IVT1),DCORVG(3,IVT1),IBCT1)
       ENDIF
C
       IF (IBCT2.EQ.0) THEN
        X2=DCORVG(1,IVT2)
        Y2=DCORVG(2,IVT2)
        Z2=DCORVG(3,IVT2)
       ELSE
        X2=PARX(DCORVG(1,IVT2),DCORVG(2,IVT2),DCORVG(3,IVT2),IBCT2)
        Y2=PARY(DCORVG(1,IVT2),DCORVG(2,IVT2),DCORVG(3,IVT2),IBCT2)
        Z2=PARZ(DCORVG(1,IVT2),DCORVG(2,IVT2),DCORVG(3,IVT2),IBCT2)
       ENDIF
C
       IF (IBCT3.EQ.0) THEN
        X3=DCORVG(1,IVT3)
        Y3=DCORVG(2,IVT3)
        Z3=DCORVG(3,IVT3)
       ELSE
        X3=PARX(DCORVG(1,IVT3),DCORVG(2,IVT3),DCORVG(3,IVT3),IBCT3)
        Y3=PARY(DCORVG(1,IVT3),DCORVG(2,IVT3),DCORVG(3,IVT3),IBCT3)
        Z3=PARZ(DCORVG(1,IVT3),DCORVG(2,IVT3),DCORVG(3,IVT3),IBCT3)
       ENDIF
C
       IF (IBCT4.EQ.0) THEN
        X4=DCORVG(1,IVT4)
        Y4=DCORVG(2,IVT4)
        Z4=DCORVG(3,IVT4)
       ELSE
        X4=PARX(DCORVG(1,IVT4),DCORVG(2,IVT4),DCORVG(3,IVT4),IBCT4)
        Y4=PARY(DCORVG(1,IVT4),DCORVG(2,IVT4),DCORVG(3,IVT4),IBCT4)
        Z4=PARZ(DCORVG(1,IVT4),DCORVG(2,IVT4),DCORVG(3,IVT4),IBCT4)
       ENDIF
C
       IF ((IAR.EQ.1).OR.(IAR.EQ.6)) THEN
        KVERT(3,KEL(KIAD(1,IAR)))=NVT
        KVERT(3,KEL(KIAD(2,IAR)))=NVT
        KVERT(3,KEL(KIAD(3,IAR)))=NVT
        KVERT(3,KEL(KIAD(4,IAR)))=NVT
       ELSE
        KVERT(6,KEL(KIAD(1,IAR)))=NVT
        KVERT(8,KEL(KIAD(2,IAR)))=NVT
        KVERT(8,KEL(KIAD(3,IAR)))=NVT
        KVERT(6,KEL(KIAD(4,IAR)))=NVT
       ENDIF
C
       DCORVG(1,NVT)=0.25D0*(X1+X2+X3+X4)
       DCORVG(2,NVT)=0.25D0*(Y1+Y2+Y3+Y4)
       DCORVG(3,NVT)=0.25D0*(Z1+Z2+Z3+Z4)
C
       IF (IBCT1*IBCT2*IBCT3*IBCT4.EQ.0) THEN
        KNPR(NVT)=0
        DCORVG(1,NVT)=0.25D0*(X1+X2+X3+X4)
        DCORVG(2,NVT)=0.25D0*(Y1+Y2+Y3+Y4)
        DCORVG(3,NVT)=0.25D0*(Z1+Z2+Z3+Z4)
       ELSE
        IF ((IBCT1.EQ.IBCT2).AND.(IBCT3.EQ.IBCT4).AND.
     *      (IBCT1.EQ.IBCT4).AND.(INEIGH.EQ.0)) THEN
         KNPR(NVT)=IBCT1
         CALL SADB(DCORVG(1,IVT1),DCORVG(2,IVT1),DCORVG(3,IVT1),
     *             DCORVG(1,IVT2),DCORVG(2,IVT2),DCORVG(3,IVT2),
     *             DCORVG(1,IVT3),DCORVG(2,IVT3),DCORVG(3,IVT3),
     *             DCORVG(1,IVT4),DCORVG(2,IVT4),DCORVG(3,IVT4),
     *             PAR1,PAR2,PAR3,IBCT1)
         DCORVG(1,NVT)=PAR1
         DCORVG(2,NVT)=PAR2
         DCORVG(3,NVT)=PAR3
        ELSE
         KNPR(NVT)=0
         DCORVG(1,NVT)=0.25D0*(X1+X2+X3+X4)
         DCORVG(2,NVT)=0.25D0*(Y1+Y2+Y3+Y4)
         DCORVG(3,NVT)=0.25D0*(Z1+Z2+Z3+Z4)
        ENDIF
       ENDIF
C
       IF ((IAR.EQ.1).OR.(IAR.EQ.6)) THEN
        KADJ(1,KEL(KIAD(1,IAR)))=INEIGH
        KADJ(1,KEL(KIAD(2,IAR)))=INEIGH
        KADJ(1,KEL(KIAD(3,IAR)))=INEIGH
        KADJ(1,KEL(KIAD(4,IAR)))=INEIGH
       ELSE
        KADJ(2,KEL(KIAD(1,IAR)))=INEIGH
        KADJ(5,KEL(KIAD(2,IAR)))=INEIGH
        KADJ(5,KEL(KIAD(3,IAR)))=INEIGH
        KADJ(2,KEL(KIAD(4,IAR)))=INEIGH
       ENDIF
C
C
      ELSE
C
C
       IF (IEL.EQ.KADJ(1,INEIGH)) THEN
        KIELH(1)=INEIGH
        KIELH(2)=KADJ(3,KIELH(1))
        KIELH(3)=KADJ(3,KIELH(2))
        KIELH(4)=KADJ(3,KIELH(3))
        IARH=1
        GOTO 152
       ENDIF
C
       IF (IEL.EQ.KADJ(1,KADJ(6,INEIGH))) THEN
        KIELH(1)=KADJ(6,INEIGH)
        KIELH(2)=KADJ(3,KIELH(1))
        KIELH(3)=KADJ(3,KIELH(2))
        KIELH(4)=KADJ(3,KIELH(3))
        IARH=6
        GOTO 152
       ENDIF
C
       IELH=INEIGH
       DO 151 I=1,4
       IF (IEL.EQ.KADJ(2,IELH)) THEN
        KIELH(1)=IELH
        KIELH(2)=KADJ(3,KIELH(1))
        KIELH(3)=KADJ(6,KIELH(2))
        KIELH(4)=KADJ(4,KIELH(3))
        IARH=I+1
        GOTO 152
       ENDIF
       IELH=KADJ(3,IELH)
151    CONTINUE
C
152    KIEL(1)=KEL(KIAD(1,IAR))
       KIEL(2)=KEL(KIAD(2,IAR))
       KIEL(3)=KEL(KIAD(3,IAR))
       KIEL(4)=KEL(KIAD(4,IAR))
C
       DO 153 I=1,4
       IF ((IAR.EQ.1).OR.(IAR.EQ.6)) THEN
        IVT1=KVERT(1,KIEL(I))
        IVT2=KVERT(2,KIEL(I))
       ELSE
        IVT1=KVERT(1,KIEL(I))
        IVT2=KVERT(5,KIEL(I))
       ENDIF
C
       DO 154 J=1,4
       IV1=KVERT(1,KIELH(J))
       IV2=KVERT(2,KIELH(J))
       IV3=KVERT(3,KIELH(J))
       IV4=KVERT(4,KIELH(J))
       IV5=KVERT(5,KIELH(J))
       IV6=KVERT(6,KIELH(J))
       IV7=KVERT(7,KIELH(J))
       IV8=KVERT(8,KIELH(J))
       IF (((IV1.EQ.IVT1).OR.(IV2.EQ.IVT1).OR.
     *      (IV3.EQ.IVT1).OR.(IV4.EQ.IVT1).OR.
     *      (IV5.EQ.IVT1).OR.(IV6.EQ.IVT1).OR.
     *      (IV7.EQ.IVT1).OR.(IV8.EQ.IVT1)).AND.
     *     ((IV1.EQ.IVT2).OR.(IV2.EQ.IVT2).OR.
     *      (IV3.EQ.IVT2).OR.(IV4.EQ.IVT2).OR.
     *      (IV5.EQ.IVT2).OR.(IV6.EQ.IVT2).OR.
     *      (IV7.EQ.IVT2).OR.(IV8.EQ.IVT2))) GOTO 155
154     CONTINUE
C
155     IF ((IAR.EQ.1).OR.(IAR.EQ.6)) THEN
         INEI1=1
         IVER1=3
        ELSE
         IF ((I.EQ.1).OR.(I.EQ.4)) THEN
          INEI1=2
          IVER1=6
         ELSE
          INEI1=5
          IVER1=8
         ENDIF
        ENDIF
C
        IF ((IARH.EQ.1).OR.(IARH.EQ.6)) THEN
         INEI2=1
         IVER2=3
        ELSE
         IF ((J.EQ.1).OR.(J.EQ.4)) THEN
          INEI2=2
          IVER2=6
         ELSE
          INEI2=5
          IVER2=8
         ENDIF
        ENDIF
C
        KADJ(INEI1,KIEL(I))=KIELH(J)
        KADJ(INEI2,KIELH(J))=KIEL(I)
C
        KVERT(IVER1,KIEL(I))=KVERT(IVER2,KIELH(J))
C
153     CONTINUE
C
C
      ENDIF
C
150   CONTINUE
C
      DO 160 I1=1,4
      I2=MOD(I1,4)+1
      NAT=NAT+1
      KADJ(3,KEL(I1))=KEL(I2)
160   KADJ(4,KEL(I2))=KEL(I1)
C
      DO 161 I1=1,4
      I2=I1+4
      NAT=NAT+1
      KADJ(6,KEL(I1))=KEL(I2)
161   KADJ(6,KEL(I2))=KEL(I1)
C
      DO 162 I1=1,4
      I2=MOD(I1,4)+1
      NAT=NAT+1
      KADJ(3,KEL(I1+4))=KEL(I2+4)
162   KADJ(4,KEL(I2+4))=KEL(I1+4)
C
      NEL=NEL+7
C
110   CONTINUE
C
C Calculation of coordinates and numbers of centers
C
      DO 200 IEL=1,NELOLD
C
      KEL(1)=IEL
      DO 201 I=1,3
201   KEL(I+1)=KADJ(3,KEL(I))
C
      KEL(5)=KADJ(6,KEL(I))
      DO 202 I=5,7
202   KEL(I+1)=KADJ(3,KEL(I))
C
      NVT=NVT+1
      DO 210 IELH=1,8
210   KVERT(7,KEL(IELH))=NVT
C
      IVT1=KVERT(1,KEL(1))
      IVT2=KVERT(1,KEL(2))
      IVT3=KVERT(1,KEL(3))
      IVT4=KVERT(1,KEL(4))
      IVT5=KVERT(1,KEL(5))
      IVT6=KVERT(1,KEL(6))
      IVT7=KVERT(1,KEL(7))
      IVT8=KVERT(1,KEL(8))
      IBCT1=KNPR(IVT1)
      IBCT2=KNPR(IVT2)
      IBCT3=KNPR(IVT3)
      IBCT4=KNPR(IVT4)
      IBCT5=KNPR(IVT5)
      IBCT6=KNPR(IVT6)
      IBCT7=KNPR(IVT7)
      IBCT8=KNPR(IVT8)
C
      IF (IBCT1.EQ.0) THEN
       X1=DCORVG(1,IVT1)
       Y1=DCORVG(2,IVT1)
       Z1=DCORVG(3,IVT1)
      ELSE
       X1=PARX(DCORVG(1,IVT1),DCORVG(2,IVT1),DCORVG(3,IVT1),IBCT1)
       Y1=PARY(DCORVG(1,IVT1),DCORVG(2,IVT1),DCORVG(3,IVT1),IBCT1)
       Z1=PARZ(DCORVG(1,IVT1),DCORVG(2,IVT1),DCORVG(3,IVT1),IBCT1)
      ENDIF
C
      IF (IBCT2.EQ.0) THEN
       X2=DCORVG(1,IVT2)
       Y2=DCORVG(2,IVT2)
       Z2=DCORVG(3,IVT2)
      ELSE
       X2=PARX(DCORVG(1,IVT2),DCORVG(2,IVT2),DCORVG(3,IVT2),IBCT2)
       Y2=PARY(DCORVG(1,IVT2),DCORVG(2,IVT2),DCORVG(3,IVT2),IBCT2)
       Z2=PARZ(DCORVG(1,IVT2),DCORVG(2,IVT2),DCORVG(3,IVT2),IBCT2)
      ENDIF
C
      IF (IBCT3.EQ.0) THEN
       X3=DCORVG(1,IVT3)
       Y3=DCORVG(2,IVT3)
       Z3=DCORVG(3,IVT3)
      ELSE
       X3=PARX(DCORVG(1,IVT3),DCORVG(2,IVT3),DCORVG(3,IVT3),IBCT3)
       Y3=PARY(DCORVG(1,IVT3),DCORVG(2,IVT3),DCORVG(3,IVT3),IBCT3)
       Z3=PARZ(DCORVG(1,IVT3),DCORVG(2,IVT3),DCORVG(3,IVT3),IBCT3)
      ENDIF
C
      IF (IBCT4.EQ.0) THEN
       X4=DCORVG(1,IVT4)
       Y4=DCORVG(2,IVT4)
       Z4=DCORVG(3,IVT4)
      ELSE
       X4=PARX(DCORVG(1,IVT4),DCORVG(2,IVT4),DCORVG(3,IVT4),IBCT4)
       Y4=PARY(DCORVG(1,IVT4),DCORVG(2,IVT4),DCORVG(3,IVT4),IBCT4)
       Z4=PARZ(DCORVG(1,IVT4),DCORVG(2,IVT4),DCORVG(3,IVT4),IBCT4)
      ENDIF
C
      IF (IBCT5.EQ.0) THEN
       X5=DCORVG(1,IVT5)
       Y5=DCORVG(2,IVT5)
       Z5=DCORVG(3,IVT5)
      ELSE
       X5=PARX(DCORVG(1,IVT5),DCORVG(2,IVT5),DCORVG(3,IVT5),IBCT5)
       Y5=PARY(DCORVG(1,IVT5),DCORVG(2,IVT5),DCORVG(3,IVT5),IBCT5)
       Z5=PARZ(DCORVG(1,IVT5),DCORVG(2,IVT5),DCORVG(3,IVT5),IBCT5)
      ENDIF
C
      IF (IBCT6.EQ.0) THEN
       X6=DCORVG(1,IVT6)
       Y6=DCORVG(2,IVT6)
       Z6=DCORVG(3,IVT6)
      ELSE
       X6=PARX(DCORVG(1,IVT6),DCORVG(2,IVT6),DCORVG(3,IVT6),IBCT6)
       Y6=PARY(DCORVG(1,IVT6),DCORVG(2,IVT6),DCORVG(3,IVT6),IBCT6)
       Z6=PARZ(DCORVG(1,IVT6),DCORVG(2,IVT6),DCORVG(3,IVT6),IBCT6)
      ENDIF
C
      IF (IBCT7.EQ.0) THEN
       X7=DCORVG(1,IVT7)
       Y7=DCORVG(2,IVT7)
       Z7=DCORVG(3,IVT7)
      ELSE
       X7=PARX(DCORVG(1,IVT7),DCORVG(2,IVT7),DCORVG(3,IVT7),IBCT7)
       Y7=PARY(DCORVG(1,IVT7),DCORVG(2,IVT7),DCORVG(3,IVT7),IBCT7)
       Z7=PARZ(DCORVG(1,IVT7),DCORVG(2,IVT7),DCORVG(3,IVT7),IBCT7)
      ENDIF
C
      IF (IBCT8.EQ.0) THEN
       X8=DCORVG(1,IVT8)
       Y8=DCORVG(2,IVT8)
       Z8=DCORVG(3,IVT8)
      ELSE
       X8=PARX(DCORVG(1,IVT8),DCORVG(2,IVT8),DCORVG(3,IVT8),IBCT8)
       Y8=PARY(DCORVG(1,IVT8),DCORVG(2,IVT8),DCORVG(3,IVT8),IBCT8)
       Z8=PARZ(DCORVG(1,IVT8),DCORVG(2,IVT8),DCORVG(3,IVT8),IBCT8)
      ENDIF
C
      KNPR(NVT)=0
      DCORVG(1,NVT)=0.125D0*(X1+X2+X3+X4+X5+X6+X7+X8)
      DCORVG(2,NVT)=0.125D0*(Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8)
      DCORVG(3,NVT)=0.125D0*(Z1+Z2+Z3+Z4+Z5+Z6+Z7+Z8)
C
200   CONTINUE
C
C
C
      END
