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
* SBEEL                                                                *
*                                                                      *
* Purpose  Determine KEEL and NEEL                                     *
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
      SUBROUTINE SBEEL(KEDGE,KEEL,IPAR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNEE=12)
      DIMENSION KEDGE(NNEE,*),KEEL(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='SBEEL'
      IF (ICHECK.GE.997) CALL OTRC('SBEEL ','12/11/89')
C
C
C
      IF (IPAR.EQ.0) THEN
       NEEL=0
       DO 100 IET=1,NET
100    KEEL(IET)=0
C
       DO 110 IEL=1,NEL
       DO 120 IEE=1,NEE
       IET=KEDGE(IEE,IEL)
       KEEL(IET)=KEEL(IET)+1
       NEEL=MAX0(NEEL,KEEL(IET))
120    CONTINUE
110    CONTINUE
C
       GOTO 99999
      ENDIF
C
C
      IF (IPAR.EQ.1) THEN
       DO 200 IEEL=1,NEEL*NET
200    KEEL(IEEL)=0
C
       DO 210 IEL=1,NEL
       DO 220 IEE=1,NEE
       IET=KEDGE(IEE,IEL)
       DO 230 IEEL=1,NEEL
230    IF (KEEL(NEEL*(IET-1)+IEEL).EQ.0) GOTO 231
231    KEEL(NEEL*(IET-1)+IEEL)=IEL
220    CONTINUE
210    CONTINUE
C
C sortieren von KEEL
C
       DO 250 IET=1,NET
       DO 251 IEEL=1,NEEL
       MIN=KEEL(NEEL*(IET-1)+IEEL)
       IF (MIN.EQ.0) GOTO 250
       DO 252 JEEL=IEEL+1,NEEL
       IF (KEEL(NEEL*(IET-1)+JEEL).EQ.0) GOTO 250
       IF (KEEL(NEEL*(IET-1)+JEEL).LT.MIN) THEN
        KEEL(NEEL*(IET-1)+IEEL)=KEEL(NEEL*(IET-1)+JEEL)
        KEEL(NEEL*(IET-1)+JEEL)=MIN
        MIN=KEEL(NEEL*(IET-1)+IEEL)
       ENDIF
252    CONTINUE
251    CONTINUE
250    CONTINUE
C
       GOTO 99999
      ENDIF
C
C
C
99999 END
