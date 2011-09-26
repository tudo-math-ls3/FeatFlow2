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
* SBCB                                                                 *
*                                                                      *
* Purpose  Determine KNPR from movie grid                              *
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
      SUBROUTINE XSBCB
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNARR=299)
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
      SUB='XSBCB'
      IF (ICHECK.GE.997) CALL OTRC('XSBCA ','12/11/89')
C
      CALL ZLEN(LNPR,ILEN)
      IF ((ILEN.LT.NVT).AND.(ILEN.GT.0)) CALL ZDISP(0,LNPR,'KNPR  ')
      CALL ZNEW(NVT,3,LNPR,'KNPR  ')
      IF (IER.NE.0) GOTO 99999
C
      CALL SBCB(KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LVEL)),
     *          KWORK(L(LNPR)))
      IF (IER.NE.0) GOTO 99999
C
C
C
99999 END
C
C
C
      SUBROUTINE SBCB(KVERT,KADJ,KVEL,KNPR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNVE=8,NNAE=6)
      DIMENSION KVERT(NNVE,*),KADJ(NNAE,*),KVEL(*),KNPR(*)
      DIMENSION KIAD(4,6),KHIA(3)
      DATA KIAD/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='SBCB'
      IF (ICHECK.GE.997) CALL OTRC('SBCB  ','12/11/89')
C
C
C
      IF (NBCT.EQ.1) THEN
C
       DO 100 IVT=1,NVT
       IEL=KVEL(NVEL*(IVT-1)+1)
C
       IARC=0
       DO 110 IAR=1,NAE
       IV1=KVERT(KIAD(1,IAR),IEL)
       IV2=KVERT(KIAD(2,IAR),IEL)
       IV3=KVERT(KIAD(3,IAR),IEL)
       IV4=KVERT(KIAD(4,IAR),IEL)
C
       IF ((IV1.EQ.IVT).OR.(IV2.EQ.IVT).OR.(IV3.EQ.IVT).OR.
     *     (IV4.EQ.IVT)) THEN
        IARC=IARC+1
        KHIA(IARC)=IAR
        IF (IARC.EQ.3) GOTO 111
       ENDIF
C
110    CONTINUE
C
111    IF ((KADJ(KHIA(1),IEL).EQ.0).OR.(KADJ(KHIA(2),IEL).EQ.0).OR.
     *     (KADJ(KHIA(3),IEL).EQ.0)) KNPR(IVT)=NBCT
C
100    CONTINUE
C
       GOTO 99999
      ENDIF
C
C
C
99999 END
