************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT3D  (Release 1.1)               *
*                                                                      *
* Authors: J. Harig, S. Turek                                          *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                     *
************************************************************************
*                                                                      *
* SBVEL                                                                *
*                                                                      *
* Purpose  Determine KVEL and NVEL                                     *
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
      SUBROUTINE SBVEL(KVERT,KVEL,IPAR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNVE=8)
      DIMENSION KVERT(NNVE,*),KVEL(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='SBVEL'
      IF (ICHECK.GE.997) CALL OTRC('SBVEL ','12/11/89')
C
C
C
      IF (IPAR.EQ.0) THEN
       NVEL=0
       DO 100 IVT=1,NVT
100    KVEL(IVT)=0
C
       DO 110 IEL=1,NEL
       DO 120 IVE=1,NVE
       IVT=KVERT(IVE,IEL)
       KVEL(IVT)=KVEL(IVT)+1
       NVEL=MAX0(NVEL,KVEL(IVT))
120    CONTINUE
110    CONTINUE
C
       GOTO 99999
      ENDIF
C
C
      IF (IPAR.EQ.1) THEN
       DO 200 IVEL=1,NVEL*NVT
200    KVEL(IVEL)=0
C
       DO 210 IEL=1,NEL
       DO 220 IVE=1,NVE
       IVT=KVERT(IVE,IEL)
       DO 230 IVEL=1,NVEL
230    IF (KVEL(NVEL*(IVT-1)+IVEL).EQ.0) GOTO 231
231    KVEL(NVEL*(IVT-1)+IVEL)=IEL
220    CONTINUE
210    CONTINUE
C
C sortieren von KVEL
C
       DO 250 IVT=1,NVT
       DO 251 IVEL=1,NVEL
       MIN=KVEL(NVEL*(IVT-1)+IVEL)
       IF (MIN.EQ.0) GOTO 250
       DO 252 JVEL=IVEL+1,NVEL
       IF (KVEL(NVEL*(IVT-1)+JVEL).EQ.0) GOTO 250
       IF (KVEL(NVEL*(IVT-1)+JVEL).LT.MIN) THEN
        KVEL(NVEL*(IVT-1)+IVEL)=KVEL(NVEL*(IVT-1)+JVEL)
        KVEL(NVEL*(IVT-1)+JVEL)=MIN
        MIN=KVEL(NVEL*(IVT-1)+IVEL)
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
