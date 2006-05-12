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
* SBCA                                                                 *
*                                                                      *
* Purpose  Determine KADJ from coarse grid                             *
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
      SUBROUTINE XSBCA(IDISP)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNVE=8,NNAE=6,NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/,/TRIAA/
C
      SUB='XSBCA'
      IF (ICHECK.GE.997) CALL OTRC('XSBCA ','12/11/89')
C
c      CALL ZLEN(LADJ,ILEN)
c      IF ((ILEN.LT.NNAE*NEL).AND.(ILEN.GT.0))
c     *    CALL ZDISP(0,LADJ,'KADJ  ')
      IF (LADJ.NE.0) THEN
       CALL ZDISP(0,LADJ,'KADJ  ')
       IF (IER.NE.0) GOTO 99999
       LADJ=0
      ENDIF
C
      CALL ZNEW(NNAE*NEL,3,LADJ,'KADJ  ')
      IF (IER.NE.0) GOTO 99999
C
      DO 1 I=1,NNAE*NEL
1     KWORK(L(LADJ)+I-1)=-1
C
      CALL SBCA(KWORK(L(LVERT)),KWORK(L(LADJ)))
      IF (IER.NE.0) GOTO 99999
C
      IF (IDISP.EQ.1) CALL ZDISP(NNAE*NEL,LADJ,'KADJ  ')
C
C
99999 END
C
C
C
      SUBROUTINE SBCA(KVERT,KADJ)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNVE=8,NNAE=6)
      DIMENSION KVERT(NNVE,*),KADJ(NNAE,*)
      DIMENSION KIAD(4,6)
      DATA KIAD/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='SBCA'
      IF (ICHECK.GE.997) CALL OTRC('SBCA  ','12/11/89')
C
C
C

      NAT=0
C
      DO 100 IEL=1,NEL
      DO 110 IAR=1,NAE
C
      IF (KADJ(IAR,IEL).GE.0) GOTO 110
C
      NAT=NAT+1
      IVE1=KIAD(1,IAR)
      IVE2=KIAD(2,IAR)
      IVE3=KIAD(3,IAR)
      IVE4=KIAD(4,IAR)
C
      IVT1=KVERT(IVE1,IEL)
      IVT2=KVERT(IVE2,IEL)
      IVT3=KVERT(IVE3,IEL)
      IVT4=KVERT(IVE4,IEL)
C
      DO 120 JEL=1,NEL
      IF (JEL.EQ.IEL) GOTO 120
      DO 130 JAR=1,NAE
C
      JVE1=KIAD(1,JAR)
      JVE2=KIAD(2,JAR)
      JVE3=KIAD(3,JAR)
      JVE4=KIAD(4,JAR)
C
      JVT1=KVERT(JVE1,JEL)
      JVT2=KVERT(JVE2,JEL)
      JVT3=KVERT(JVE3,JEL)
      JVT4=KVERT(JVE4,JEL)
C
      IF (((JVT1.EQ.IVT1).OR.(JVT2.EQ.IVT1).OR.
     *     (JVT3.EQ.IVT1).OR.(JVT4.EQ.IVT1)).AND.
     *    ((JVT1.EQ.IVT2).OR.(JVT2.EQ.IVT2).OR.
     *     (JVT3.EQ.IVT2).OR.(JVT4.EQ.IVT2)).AND.
     *    ((JVT1.EQ.IVT3).OR.(JVT2.EQ.IVT3).OR.
     *     (JVT3.EQ.IVT3).OR.(JVT4.EQ.IVT3)).AND.
     *    ((JVT1.EQ.IVT4).OR.(JVT2.EQ.IVT4).OR.
     *     (JVT3.EQ.IVT4).OR.(JVT4.EQ.IVT4))) THEN
       KADJ(IAR,IEL)=JEL
       KADJ(JAR,JEL)=IEL
       GOTO 110
      ENDIF
C
130   CONTINUE
120   CONTINUE
C
      KADJ(IAR,IEL)=0
C
110   CONTINUE
100   CONTINUE
C
C
C
      END
