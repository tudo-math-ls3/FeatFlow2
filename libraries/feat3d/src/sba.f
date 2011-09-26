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
* SBA                                                                  *
*                                                                      *
* Purpose  Adjust dimensions of KAREA and NAT                          *
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
      SUBROUTINE SBA(KVERT,KNPR,KADJ,DCORVG,KAREA,DCORAG,ISA,IND,
     *               PARX,PARY,PARZ,SADB)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNVE=8,NNAE=6)
      DIMENSION KVERT(NNVE,*),KNPR(*),KADJ(NNAE,*),KAREA(NNAE,*)
      DIMENSION DCORVG(3,*),DCORAG(3,*)
      DIMENSION KIAD(4,6)
      DATA KIAD/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      EXTERNAL PARX,PARY,PARZ,SADB
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='SBA'
      IF (ICHECK.GE.997) CALL OTRC('SBA  ','12/11/89')
C
C
C
      NAT=0
C
      DO 100 IEL=1,NEL
      DO 110 IAR=1,NAE
      INEIGH=KADJ(IAR,IEL)
C
      IF ((INEIGH.EQ.0).OR.(INEIGH.GE.IEL)) THEN
       NAT=NAT+1
       KAREA(IAR,IEL)=NAT
C
       IF (ISA.GE.2) THEN
        IVT1=KVERT(KIAD(1,IAR),IEL)
        IVT2=KVERT(KIAD(2,IAR),IEL)
        IVT3=KVERT(KIAD(3,IAR),IEL)
        IVT4=KVERT(KIAD(4,IAR),IEL)
        INPR1=KNPR(IVT1)
        INPR2=KNPR(IVT2)
        INPR3=KNPR(IVT3)
        INPR4=KNPR(IVT4)
        IF ((INPR1.EQ.INPR2).AND.(INPR1.EQ.INPR3).AND.
     *      (INPR1.EQ.INPR4).AND.(INEIGH.EQ.0)) THEN
         INPR=INPR1
         KNPR(IND+NAT)=INPR
        ELSE
         INPR=0
         KNPR(IND+NAT)=INPR
        ENDIF
       ENDIF
C
       IF (ISA.GE.3) THEN
        IF (INPR1.EQ.0) THEN
         PX1=DCORVG(1,IVT1)
         PY1=DCORVG(2,IVT1)
         PZ1=DCORVG(3,IVT1)
        ELSE
         PX1=PARX(DCORVG(1,IVT1),DCORVG(2,IVT1),DCORVG(3,IVT1),INPR1)
         PY1=PARY(DCORVG(1,IVT1),DCORVG(2,IVT1),DCORVG(3,IVT1),INPR1)
         PZ1=PARZ(DCORVG(1,IVT1),DCORVG(2,IVT1),DCORVG(3,IVT1),INPR1)
        ENDIF
C
        IF (INPR2.EQ.0) THEN
         PX2=DCORVG(1,IVT2)
         PY2=DCORVG(2,IVT2)
         PZ2=DCORVG(3,IVT2)
        ELSE
         PX2=PARX(DCORVG(1,IVT2),DCORVG(2,IVT2),DCORVG(3,IVT2),INPR2)
         PY2=PARY(DCORVG(1,IVT2),DCORVG(2,IVT2),DCORVG(3,IVT2),INPR2)
         PZ2=PARZ(DCORVG(1,IVT2),DCORVG(2,IVT2),DCORVG(3,IVT2),INPR2)
        ENDIF
C
        IF (INPR3.EQ.0) THEN
         PX3=DCORVG(1,IVT3)
         PY3=DCORVG(2,IVT3)
         PZ3=DCORVG(3,IVT3)
        ELSE
         PX3=PARX(DCORVG(1,IVT3),DCORVG(2,IVT3),DCORVG(3,IVT3),INPR3)
         PY3=PARY(DCORVG(1,IVT3),DCORVG(2,IVT3),DCORVG(3,IVT3),INPR3)
         PZ3=PARZ(DCORVG(1,IVT3),DCORVG(2,IVT3),DCORVG(3,IVT3),INPR3)
        ENDIF
C
        IF (INPR4.EQ.0) THEN
         PX4=DCORVG(1,IVT4)
         PY4=DCORVG(2,IVT4)
         PZ4=DCORVG(3,IVT4)
        ELSE
         PX4=PARX(DCORVG(1,IVT4),DCORVG(2,IVT4),DCORVG(3,IVT4),INPR4)
         PY4=PARY(DCORVG(1,IVT4),DCORVG(2,IVT4),DCORVG(3,IVT4),INPR4)
         PZ4=PARZ(DCORVG(1,IVT4),DCORVG(2,IVT4),DCORVG(3,IVT4),INPR4)
        ENDIF
C
        IF (INPR.EQ.0) THEN
         DCORAG(1,NAT)=0.25D0*(PX1+PX2+PX3+PX4)
         DCORAG(2,NAT)=0.25D0*(PY1+PY2+PY3+PY4)
         DCORAG(3,NAT)=0.25D0*(PZ1+PZ2+PZ3+PZ4)
        ELSE
         CALL SADB(DCORVG(1,IVT1),DCORVG(2,IVT1),DCORVG(3,IVT1),
     *             DCORVG(1,IVT2),DCORVG(2,IVT2),DCORVG(3,IVT2),
     *             DCORVG(1,IVT3),DCORVG(2,IVT3),DCORVG(3,IVT3),
     *             DCORVG(1,IVT4),DCORVG(2,IVT4),DCORVG(3,IVT4),
     *             PAR1,PAR2,PAR3,INPR)
         DCORAG(1,NAT)=PARX(PAR1,PAR2,PAR3,INPR)
         DCORAG(2,NAT)=PARY(PAR1,PAR2,PAR3,INPR)
         DCORAG(3,NAT)=PARZ(PAR1,PAR2,PAR3,INPR)
        ENDIF
       ENDIF
C
      ELSE
       DO 120 JAR=1,NAE
       JNEIGH=KADJ(JAR,INEIGH)
       IF (JNEIGH.EQ.IEL) GOTO 121
120    CONTINUE
121    KAREA(IAR,IEL)=KAREA(JAR,INEIGH)
      ENDIF
C
110   CONTINUE
100   CONTINUE
C
C
C
      END
