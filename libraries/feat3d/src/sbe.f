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
* SBE                                                                  *
*                                                                      *
* Purpose  Adjust dimensions of KEDGE and NET                          *
*                                                                      *
* Subroutines/functions called  KGVVEL                                 *
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
      SUBROUTINE SBE(KVERT,KNPR,KVEL,DCORVG,KEDGE,DCOREG,ISE,
     *               PARX,PARY,PARZ,SEDB)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNVE=8,NNEE=12)
      DIMENSION KVERT(NNVE,*),KNPR(*),KVEL(*),KEDGE(NNEE,*)
      DIMENSION DCORVG(3,*),DCOREG(3,*)
      DIMENSION KIV(2,12)
      DATA KIV /1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      EXTERNAL PARX,PARY,PARZ,SEDB
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='SBE'
      IF (ICHECK.GE.997) CALL OTRC('SBE   ','08/29/89')
C
C
C

      NET=0
      DO 100 IEL=1,NEL
      DO 110 IED=1,NEE
C
      IV1=KIV(1,IED)
      IV2=KIV(2,IED)
C
      IVT1=KVERT(IV1,IEL)
      IVT2=KVERT(IV2,IEL)
C
      CALL KGVVEL(KVEL(NVEL*(IVT1-1)+1),KVEL(NVEL*(IVT2-1)+1),
     *            NVEL,IEL,IGVEL)
      IF (IGVEL.GE.IEL) THEN
       
       NET=NET+1
       KEDGE(IED,IEL)=NET
C
       IF (ISE.GE.2) THEN
        INPR1=KNPR(IVT1)
        INPR2=KNPR(IVT2)
        IF (INPR1.EQ.INPR2) THEN
         INPR=INPR1
         KNPR(NVT+NET)=INPR
        ELSE
         INPR=0
         KNPR(NVT+NET)=INPR
        ENDIF
       ENDIF
C
       IF (ISE.GE.3) THEN
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
        IF (INPR.EQ.0) THEN
         DCOREG(1,NET)=0.5D0*(PX1+PX2)
         DCOREG(2,NET)=0.5D0*(PY1+PY2)
         DCOREG(3,NET)=0.5D0*(PZ1+PZ2)
        ELSE
         CALL SEDB(DCORVG(1,IVT1),DCORVG(2,IVT1),DCORVG(3,IVT1),
     *             DCORVG(1,IVT2),DCORVG(2,IVT2),DCORVG(3,IVT2),
     *             PAR1,PAR2,PAR3,INPR)
         DCOREG(1,NET)=PARX(PAR1,PAR2,PAR3,INPR)
         DCOREG(2,NET)=PARY(PAR1,PAR2,PAR3,INPR)
         DCOREG(3,NET)=PARZ(PAR1,PAR2,PAR3,INPR)
        ENDIF        
       ENDIF
C
      ELSE
       DO 120 IEDH=1,NEE
       IVH1=KIV(1,IEDH)
       IVH2=KIV(2,IEDH)
       IVTH1=KVERT(IVH1,IGVEL)
       IVTH2=KVERT(IVH2,IGVEL)
       IF (((IVTH1.EQ.IVT1).AND.(IVTH2.EQ.IVT2)).OR.
     *     ((IVTH1.EQ.IVT2).AND.(IVTH2.EQ.IVT1))) GOTO 121
120    CONTINUE
121    KEDGE(IED,IEL)=KEDGE(IEDH,IGVEL)
      ENDIF
C
110   CONTINUE
100   CONTINUE
C
C
C
      END
C
C
C
      SUBROUTINE KGVVEL(KVL1,KVL2,NVEL,IEL,IGVEL)
      DIMENSION KVL1(*),KVL2(*)
C
C**** IGVEL = KLEINSTES GEMEINSAMES ELEMENT VON KVL1 UND KVL2
      IGVEL=IEL
      DO 10 I=1,NVEL
      IVL1=KVL1(I)
      IF (IVL1.EQ.0) GOTO 99999
C
      DO 20 J=1,NVEL
      IVL2=KVL2(J)
      IF (IVL2.EQ.0) GOTO 10
      IF ((IVL1.EQ.IVL2).AND.(IVL1.LT.IGVEL)) THEN
       IGVEL=IVL1
       GOTO 10
      ENDIF
20    CONTINUE
C
10    CONTINUE
C
C
99999 END
