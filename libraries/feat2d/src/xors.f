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
* XORS                                                                 *
*                                                                      *
* Purpose  Read subdivision from unit MFILE                            *
*          Data are assumed to be written by XOWS                      *
*                                                                      *
* Subroutines/functions called  OF0, XORA                              *
*                                                                      *
* Version from  12/11/89                                               *
*                                                                      *
* INPUT   TYPE                                                         *
* -----   ----                                                         *
* MFILE   I*4    Unit number of external file                          *
* CFILE   C*(*)  Filename                                              *
* IFMT    I*4    0  Format free input                                  *
*                1  Formatted input                                    *
*                                                                      *
* OUTPUT  TYPE                                                         *
* ------  ----                                                         *
* MFILE   I*4    Unit number used                                      *
* CFILE   C*(*)  Filename used (length not larger than for input)      *
* IER     I*4    Error indicator                                       *
*                -110  MFILE  exceeds maximum range 16,...,80          *
*                -110  Read error on unit MFILE                        *
*                -112  File CFILE could not be opened                  *
* The remaining output quantities are contained in                     *
* COMMON blocks /TRIAA/ and /TRIAD/                                    *
*                                                                      *
************************************************************************
C
      SUBROUTINE XORS(MFILE,CFILE,IFMT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CFILE*(*)
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
      SUB='XORS'
      IF (ICHECK.GE.997) CALL OTRC('XORS  ','12/11/89')
      IER=0
C
C *** Open I/O-file
C
      CALL OF0(MFILE,CFILE,IFMT)
      IF (IER.NE.0) GOTO 99999
C
      IF (IFMT.EQ.1) THEN
       READ (MFILE,'(5I15)',ERR=99997,END=99997)
     *                                  LCM0,LMID0,LADJ0,LVEL0,LMEL0
       READ (MFILE,'(5I15)',ERR=99997,END=99997)
     *                         LVBD0,LEBD0,LBCT0,LVBDP0,LMBDP0
       READ (MFILE,'(7I15)',ERR=99997,END=99997)
     *                                  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      ELSE
       READ (MFILE,ERR=99997,END=99997) LCM0,LMID0,LADJ0,LVEL0,LMEL0
       READ (MFILE,ERR=99997,END=99997) LVBD0,LEBD0,LBCT0,LVBDP0,LMBDP0
       READ (MFILE,ERR=99997,END=99997) NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      ENDIF
C
      CALL XORA(LCORVG,'DCORVG',MFILE,CFILE,IFMT)
      IF (IER.NE.0) GOTO 99997
      IF (LCM0.GT.0) THEN
       CALL XORA(LCORMG,'DCORMG',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
      CALL XORA(LVERT,'KVERT ',MFILE,CFILE,IFMT)
      IF (IER.NE.0) GOTO 99997
      IF (LMID0.GT.0) THEN
       CALL XORA(LMID,'KMID  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
      IF (LADJ0.GT.0) THEN
       CALL XORA(LADJ,'KADJ  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
      IF (LVEL0.GT.0) THEN
       CALL XORA(LVEL,'KVEL  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
      IF (LMEL0.GT.0) THEN
       CALL XORA(LMEL,'KMEL  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
      IF (LVBD0.GT.0) THEN
       CALL XORA(LVBD,'KVBD  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
      IF (LEBD0.GT.0) THEN
       CALL XORA(LEBD,'KEBD  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
      IF (LBCT0.GT.0) THEN
       CALL XORA(LBCT,'KBCT  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
      IF (LVBDP0.GT.0) THEN
       CALL XORA(LVBDP,'DVBDP ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
      IF (LMBDP0.GT.0) THEN
       CALL XORA(LMBDP,'DMBDP ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
      CALL XORA(LNPR,'KNPR  ',MFILE,CFILE,IFMT)
      IF (IER.NE.0) GOTO 99997
      CALL XORA(LMM,'KMM   ',MFILE,CFILE,IFMT)
      IF (IER.EQ.0) GOTO 99999
C
99997 WRITE (CPARAM,'(I15)') MFILE
      CALL WERR(-110,'XORS  ')
C
99999 END
