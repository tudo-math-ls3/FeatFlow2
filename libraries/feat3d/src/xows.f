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
* XOWS                                                                 *
*                                                                      *
* Purpose  Write subdivision onto unit MFILE                           *
*                                                                      *
* Subroutines/functions called  OF0, XOWA                              *
*                                                                      *
* Version from  12/11/89                                               *
*                                                                      *
* INPUT   TYPE                                                         *
* -----   ----                                                         *
* MFILE   I*4    Unit number of external file                          *
* CFILE   C*15   Filename                                              *
* IFMT    I*4    0  Format free output                                 *
*                1  Use format FMT(ITYPE) depending on datatype        *
*                                                                      *
* OUTPUT  TYPE                                                         *
* ------  ----                                                         *
* MFILE   I*4    Unit number used                                      *
* CFILE   C*15   Filename used                                         *
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
      CHARACTER CFILE*15
      PARAMETER (NNARR=299)
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
       WRITE (MFILE,'(3I8)',ERR=99997)
     *                         LCORVG,LCORMG,LCORAG
       WRITE (MFILE,'(5I8)',ERR=99997)
     *                         LVERT,LEDGE,LAREA,LADJ,LNPR
       WRITE (MFILE,'(9I8)',ERR=99997)
     *                         LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,
     *                         LEVE,LAVE
       WRITE (MFILE,'(4I8)',ERR=99997)
     *                         LBCT,LVBD,LEBD,LABD
       WRITE (MFILE,'(4I8)',ERR=99997)
     *                         NEL,NVT,NET,NAT
       WRITE (MFILE,'(3I8)',ERR=99997)
     *                         NVE,NEE,NAE
       WRITE (MFILE,'(5I8)',ERR=99997)
     *                         NVEL,NEEL,NVED,NVAR,NEAR
       WRITE (MFILE,'(4I8)',ERR=99997)
     *                         NBCT,NVBD,NEBD,NABD
      ELSE
       WRITE (MFILE,ERR=99997) LCORVG,LCORMG,LCORAG
       WRITE (MFILE,ERR=99997) LVERT,LEDGE,LAREA,LADJ,LNPR
       WRITE (MFILE,ERR=99997) LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,
     *                         LEVE,LAVE
       WRITE (MFILE,ERR=99997) LBCT,LVBD,LEBD,LABD
       WRITE (MFILE,ERR=99997) NEL,NVT,NET,NAT
       WRITE (MFILE,ERR=99997) NVE,NEE,NAE
       WRITE (MFILE,ERR=99997) NVEL,NEEL,NVED,NVAR,NEAR
       WRITE (MFILE,ERR=99997) NBCT,NVBD,NEBD,NABD
      ENDIF
C
      IF (LCORVG.GT.0) THEN
       CALL XOWA(LCORVG,'DCORVG',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LCORMG.GT.0) THEN
       CALL XOWA(LCORMG,'DCOREG',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LCORAG.GT.0) THEN
       CALL XOWA(LCORAG,'DCORAG',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LVERT.GT.0) THEN
       CALL XOWA(LVERT,'KVERT ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LEDGE.GT.0) THEN
       CALL XOWA(LEDGE,'KEDGE ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LAREA.GT.0) THEN
       CALL XOWA(LAREA,'KAREA ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LADJ.GT.0) THEN
       CALL XOWA(LADJ,'KADJ  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LNPR.GT.0) THEN
       CALL XOWA(LNPR,'KNPR  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LVEL.GT.0) THEN
       CALL XOWA(LVEL,'KVEL  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LEEL.GT.0) THEN
       CALL XOWA(LEEL,'KEEL  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LAEL.GT.0) THEN
       CALL XOWA(LAEL,'KAEL  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LVED.GT.0) THEN
       CALL XOWA(LVED,'KVED  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LAED.GT.0) THEN
       CALL XOWA(LAED,'KAED  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LVAR.GT.0) THEN
       CALL XOWA(LVAR,'KVAR  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LEAR.GT.0) THEN
       CALL XOWA(LEAR,'KEAR  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LEVE.GT.0) THEN
       CALL XOWA(LEVE,'KEVE  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LAVE.GT.0) THEN
       CALL XOWA(LAVE,'KAVE  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LBCT.GT.0) THEN
       CALL XOWA(LBCT,'KBCT  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LVBD.GT.0) THEN
       CALL XOWA(LVBD,'KVBD  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LEBD.GT.0) THEN
       CALL XOWA(LEBD,'KEBD  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LABD.GT.0) THEN
       CALL XOWA(LABD,'KABD  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
      GOTO 99999
C
99997 WRITE (CPARAM,'(I15)') MFILE
      CALL WERR(-111,'XOWS  ')
C
99999 END
