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
* XOWS                                                                 *
*                                                                      *
* Purpose  Write subdivision onto unit MFILE                           *
*                                                                      *
* Subroutines/functions called  XOWA                                   *
*                                                                      *
* Version from  12/11/89                                               *
*                                                                      *
* INPUT   TYPE                                                         *
* -----   ----                                                         *
* MFILE   I*4    Unit number of external file                          *
* CFILE   C*(*)  Filename                                              *
* IFMT    I*4    0  Format free output                                 *
*                1  Use format FMT(ITYPE) depending on datatype        *
*                                                                      *
* OUTPUT  TYPE                                                         *
* ------  ----                                                         *
* MFILE   I*4    Unit number used                                      *
* CFILE   C*(*)  Filename used (length not larger than for input)      *
* IER     I*4    Error indicator                                       *
*                -111  MFILE  exceeds maximum range 16,...,80          *
*                -111  Write error on unit  MFILE                      *
*                -112  File  CFILE  could not be opened                *
*                                                                      *
************************************************************************
C
      SUBROUTINE XOWS(MFILE,CFILE,IFMT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CFILE*(*)
C
      PARAMETER (NNARR=299)
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
      SUB='XOWS'
      IF (ICHECK.GE.997) CALL OTRC('XOWS  ','12/11/89')
      IER=0
C
C *** Open I/O-file
C
      CALL OF0(MFILE,CFILE,IFMT)
      IF (IER.NE.0) GOTO 99999
C
      IF (IFMT.EQ.1) THEN
       WRITE (MFILE,'(5I15)',ERR=99997)
     *                         LCORMG,LMID,LADJ,LVEL,LMEL
       WRITE (MFILE,'(5I15)',ERR=99997)
     *                         LVBD,LEBD,LBCT,LVBDP,LMBDP
       WRITE (MFILE,'(7I15)',ERR=99997)
     *                         NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      ELSE
       WRITE (MFILE,ERR=99997) LCORMG,LMID,LADJ,LVEL,LMEL
       WRITE (MFILE,ERR=99997) LVBD,LEBD,LBCT,LVBDP,LMBDP
       WRITE (MFILE,ERR=99997) NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      ENDIF
C
      CALL XOWA(LCORVG,'DCORVG',MFILE,CFILE,IFMT)
      IF (IER.NE.0) GOTO 99997
      IF (LCORMG.GT.0) THEN
       CALL XOWA(LCORMG,'DCORMG',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
      CALL XOWA(LVERT,'KVERT ',MFILE,CFILE,IFMT)
      IF (IER.NE.0) GOTO 99997
      IF (LMID.GT.0) THEN
       CALL XOWA(LMID,'KMID  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
      IF (LADJ.GT.0) THEN
       CALL XOWA(LADJ,'KADJ  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
      IF (LVEL.GT.0) THEN
       CALL XOWA(LVEL,'KVEL  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
      IF (LMEL.GT.0) THEN
       CALL XOWA(LMEL,'KMEL  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
      IF (LVBD.GT.0) THEN
       CALL XOWA(LVBD,'KVBD  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
      IF (LEBD.GT.0) THEN
       CALL XOWA(LEBD,'KEBD  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
      IF (LBCT.GT.0) THEN
       CALL XOWA(LBCT,'KBCT  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
      IF (LVBDP.GT.0) THEN
       CALL XOWA(LVBDP,'DVBDP ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
      IF (LMBDP.GT.0) THEN
       CALL XOWA(LMBDP,'DMBDP ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
      CALL XOWA(LNPR,'KNPR  ',MFILE,CFILE,IFMT)
      IF (IER.NE.0) GOTO 99997
      CALL XOWA(LMM,'KMM   ',MFILE,CFILE,IFMT)
      IF (IER.EQ.0) GOTO 99999
C
99997 WRITE (CPARAM,'(I15)') MFILE
      CALL WERR(-111,'XOWS  ')
C
99999 END
