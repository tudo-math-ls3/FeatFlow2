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
* XORSC                                                                *
*                                                                      *
* Purpose  Allocate arrays and call of subroutine ORSC to get          *
*          coarse grid from unit MFILE                                 *
*                                                                      *
* Subroutines/functions called  OF0, ORSC, ZLEN, ZNEW, ZDISP           *
*                                                                      *
* Version from  12/11/89                                               *
*                                                                      *
* INPUT   TYPE                                                         *
* -----   ----                                                         *
* MFILE   I*4    Unit number of external file                          *
* CFILE   C*15   Filename                                              *
*                                                                      *
* OUTPUT  TYPE                                                         *
* ------  ----                                                         *
* MFILE   I*4    Unit number used                                      *
* CFILE   C*15   Filename used                                         *
* IER     I*4    Error indicator                                       *
*                -110  MFILE  exceeds maximum range 16,...,80          *
*                -110  Read error on unit  MFILE                       *
*                -112  File  CFILE  could not be opened                *
*                                                                      *
************************************************************************
C
      SUBROUTINE XORSC(MFILE,CFILE)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CFILE*15,BLANK*15
C
      PARAMETER (NNVE=8,NNEE=12,NNAE=6,NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/,/TRIAA/
      DATA BLANK/'              '/
C
      SUB='XORSC'
      IF (ICHECK.GE.997) CALL OTRC('XORSC ','12/11/89')
      IER=0
C
C *** Open I/O-file
C
      CALL OF0(MFILE,CFILE,1)
      IF (IER.NE.0) GOTO 99999
C
      READ (MFILE,*,ERR=99997,END=99997)
      READ (MFILE,*,ERR=99997,END=99997)
      READ (MFILE,*,ERR=99997,END=99997) NEL,NVT,NBCT,NVE,NEE,NAE
      READ (MFILE,*,ERR=99997,END=99997)
C
      CALL ZLEN(LCORVG,JLEN)
      IF (JLEN.LT.0) GOTO 99997
      IF (JLEN.LT.3*NVT) THEN
       IF (JLEN.GT.0) CALL ZDISP(0,LCORVG,'DCORVG')
       CALL ZNEW(3*NVT,-1,LCORVG,'DCORVG')
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
      CALL ZLEN(LVERT,JLEN)
      IF (JLEN.LT.0) GOTO 99997
      IF (JLEN.LT.NNVE*NEL) THEN
       IF (JLEN.GT.0) CALL ZDISP(0,LVERT,'KVERT ')
       CALL ZNEW(NNVE*NEL,-3,LVERT,'KVERT ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
      CALL ZLEN(LNPR,JLEN)
      IF (JLEN.LT.0) GOTO 99997
      IF (JLEN.LT.NVT) THEN
       IF (JLEN.GT.0) CALL ZDISP(0,LNPR,'KNPR  ')
       CALL ZNEW(NVT,-3,LNPR,'KNPR  ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
      CALL ORSC(DWORK(L(LCORVG)),KWORK(L(LVERT)),KWORK(L(LNPR)),
     *           MFILE)
      IF (IER.EQ.0) GOTO 99999
C
C
99997 WRITE (CPARAM,'(I15)') MFILE
      CALL WERR(-110,'XORSC ')
      GOTO 99999
C
C
99999 END
C
C
C
      SUBROUTINE ORSC(DCORVG,KVERT,KNPR,MFILE)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNVE=8)
      DIMENSION DCORVG(3,*)
      DIMENSION KVERT(NNVE,*),KNPR(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='ORSC'
      IF (ICHECK.GE.997) CALL OTRC('ORSC  ','12/11/89')
C
      READ (MFILE,*,ERR=99997,END=99997)
     *     ((DCORVG(IDIM,IVT),IDIM=1,3),IVT=1,NVT)
      READ (MFILE,*,ERR=99997,END=99997)
      READ (MFILE,*,ERR=99997,END=99997)
     *     ((KVERT(IVE,IEL),IVE=1,NVE),IEL=1,NEL)
      READ (MFILE,*,ERR=99997,END=99997)
      READ (MFILE,*,ERR=99997,END=99997)
     *     (KNPR(IVT),IVT=1,NVT)
      GOTO 99999
C
99997 WRITE (CPARAM,'(I15)') MFILE
      CALL WERR(-110,'ORSC  ')
C
99999 END
