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
* XORSC                                                                *
*                                                                      *
* Purpose  Allocate arrays and call of subroutine ORSC to get          *
*          coarse grid from unit MFILE                                 *
*                                                                      *
* Subroutines/functions called  OF0, ORSC, SBD01, ZLEN, ZNEW, ZDISP    *
*                                                                      *
* Version from  08/15/90                                               *
*                                                                      *
* INPUT   TYPE                                                         *
* -----   ----                                                         *
* MFILE   I*4    Unit number of external file                          *
* CFILE   C*(*)  Filename                                              *
*                                                                      *
* OUTPUT  TYPE                                                         *
* ------  ----                                                         *
* MFILE   I*4    Unit number used                                      *
* CFILE   C*(*)  Filename used (length not larger than for input)      *
* IER     I*4    Error indicator                                       *
*                -110  MFILE  exceeds maximum range 16,...,80          *
*                -110  Read error on unit  MFILE                       *
*                -112  File  CFILE  could not be opened                *
*                                                                      *
************************************************************************
*                                                                      *
* ORSC                                                                 *
*                                                                      *
* Purpose  Read coarse grid data from unit MFILE                       *
*                                                                      *
* Subroutines/functions called  none                                   *
*                                                                      *
* Version from  08/15/90                                               *
*                                                                      *
* INPUT   TYPE                                                         *
* -----   ----                                                         *
* DCORVG  R*8    Arrays on external file to be read in                 *
* KVERT   I*4                                                          *
* KADJ    I*4                                                          *
* KNPR    I*4                                                          *
* KMM     I*4                                                          *
* MFILE   I*4    Unit number of external file                          *
*                                                                      *
* OUTPUT  TYPE                                                         *
* ------  ----                                                         *
* IER     I*4    Error indicator                                       *
*                -110  Read error on unit  MFILE                       *
* The remaining output quantities are contained in COMMON block /TRIAD/*
*                                                                      *
************************************************************************
C
      SUBROUTINE XORSC(MFILE,CFILE)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CFILE*(*)
C
      PARAMETER (NNVE=4,NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/,/TRIAA/
C
      SUB='XORSC'
      IF (ICHECK.GE.997) CALL OTRC('XORSC ','08/15/90')
      IER=0
C
C *** Open I/O-file
C
      CALL OF0(MFILE,CFILE,1)
      IF (IER.NE.0) GOTO 99999
      REWIND (MFILE)
C
      READ (MFILE,*,ERR=99997,END=99997)
      READ (MFILE,*,ERR=99997,END=99997)
      READ (MFILE,*,ERR=99997,END=99997) NEL,NVT,NMT,NVE,NBCT
      READ (MFILE,*,ERR=99997,END=99997)
C
      CALL ZLEN(LCORVG,JLEN)
      IF (JLEN.LT.0) GOTO 99997
      IF (JLEN.LT.2*NVT) THEN
       IF (JLEN.GT.0) CALL ZDISP(0,LCORVG,'DCORVG')
       CALL ZNEW(2*NVT,-1,LCORVG,'DCORVG')
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
      CALL ZLEN(LVERT,JLEN)
      IF (JLEN.LT.0) GOTO 99997
      IF (JLEN.LT.NNVE*NEL) THEN
       IF (JLEN.GT.0) CALL ZDISP(0,LVERT,'KVERT ')
       CALL ZNEW(NNVE*NEL,3,LVERT,'KVERT ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
      CALL ZLEN(LNPR,JLEN)
      IF (JLEN.LT.0) GOTO 99997
      IF (JLEN.LT.NVT+NMT) THEN
       IF (JLEN.GT.0) CALL ZDISP(0,LNPR,'KNPR  ')
       CALL ZNEW(NVT+NMT,-3,LNPR,'KNPR  ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
      CALL ZLEN(LMM,JLEN)
      IF (JLEN.LT.0) GOTO 99997
      IF (JLEN.LT.2*NBCT) THEN
       IF (JLEN.GT.0) CALL ZDISP(0,LMM,'KMM   ')
       CALL ZNEW(2*NBCT,-3,LMM,'KMM   ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
      CALL ORSC(DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *          KWORK(L(LNPR)),KWORK(L(LMM)),MFILE)
      IF (IER.NE.0) GOTO 99999
C
      IF (LADJ.GT.0) THEN
       CALL ZDISP(0,LADJ,'KADJ  ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
      CALL XS2A
C
      CALL ZLEN(LVBD,JLEN)
      IF (JLEN.LT.0) GOTO 99997
      IF (JLEN.LT.NVBD) THEN
       IF (JLEN.GT.0) CALL ZDISP(0,LVBD,'KVBD  ')
       CALL ZNEW(NVBD,-3,LVBD,'DVBD  ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
      CALL ZLEN(LVBDP,JLEN)
      IF (JLEN.LT.0) GOTO 99997
      IF (JLEN.LT.NVBD) THEN
       IF (JLEN.GT.0) CALL ZDISP(0,LVBDP,'DVBDP ')
       CALL ZNEW(NVBD,-1,LVBDP,'DVBDP ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
      IF (NMT.GT.0) THEN
       CALL ZLEN(LMBDP,JLEN)
       IF (JLEN.LT.0) GOTO 99997
       IF (JLEN.LT.NVBD) THEN
        IF (JLEN.GT.0) CALL ZDISP(0,LMBDP,'DMBDP ')
        CALL ZNEW(NVBD,-1,LMBDP,'DMBDP ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
      ENDIF
C
      CALL SBD01(DWORK(L(LCORVG)),KWORK(L(LNPR)),
     *           KWORK(L(LVBD)),DWORK(L(LVBDP)),DWORK(L(LMBDP)))
C
      GOTO 99999
C      
99997 WRITE (CPARAM,'(I15)') MFILE
      CALL WERR(-110,'XORSC ')
      GOTO 99999
C
99999 END
C
C
C
      SUBROUTINE ORSC(DCORVG,KVERT,KNPR,KMM,MFILE)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNVE=4)
      DIMENSION DCORVG(2,*)
      DIMENSION KVERT(NNVE,*),KNPR(*),KMM(2,*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='ORSC'
      IF (ICHECK.GE.997) CALL OTRC('ORSC  ','08/15/90')
C
      READ (MFILE,*,ERR=99997,END=99997)
     *     ((DCORVG(IDIM,IVT),IDIM=1,2),IVT=1,NVT)
      READ (MFILE,*,ERR=99997,END=99997)
      READ (MFILE,*,ERR=99997,END=99997)
     *     ((KVERT(IVE,IEL),IVE=1,NVE),IEL=1,NEL)
      READ (MFILE,*,ERR=99997,END=99997)
      READ (MFILE,*,ERR=99997,END=99997)
     *     (KNPR(IVMT),IVMT=1,NVT+NMT)
      READ (MFILE,*,ERR=99997,END=99997)
      READ (MFILE,*,ERR=99997,END=99997)
     *     ((KMM(I,IBCT),I=1,2),IBCT=1,NBCT)
C
      NVBD=0
      DO 10 IVT=1,NVT
10    IF (KNPR(IVT).GT.0) NVBD=NVBD+1
C
      GOTO 99999
C
99997 WRITE (CPARAM,'(I15)') MFILE
      CALL WERR(-110,'ORSC  ')
C
99999 END
