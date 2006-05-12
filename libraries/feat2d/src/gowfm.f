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
* XGOWFM, GOWFM1, GOWFM2                                               *
*                                                                      *
* Purpose  Write array as MOVIE.BYU function file onto unit MFILE      *   
*                                                                      *
* Subroutines/functions called  OF0                                    *
*                                                                      *
* Version from  08/16/90                                               *
*                                                                      *
* INPUT   TYPE                                                         *
* -----   ----                                                         *
* LNR     I*4    Number of array                                       *
* MFILE   I*4    Unit number of external file                          *
* CFILE   C*(*)  Filename                                              *
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
      SUBROUTINE XGOWFM(LNR,MFILE,CFILE)
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
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAA/
C
      SUB='XGOWFM'
      IF (ICHECK.GE.997) CALL OTRC('XGOWFM','08/16/90')
      IER=0
C
C *** Open I/O-file
C
      CALL OF0(MFILE,CFILE,1)
      IF (IER.NE.0) GOTO 99999
C
      CALL ZTYPE(LNR,ITYPE)
      CALL ZLEN(LNR,ILEN)
      IF (ITYPE.EQ.1) THEN
       CALL GOWFM1(DWORK(L(LNR)),ILEN,MFILE)
      ELSE IF (ITYPE.EQ.2) THEN
       CALL GOWFM2(VWORK(L(LNR)),ILEN,MFILE)
      ENDIF
C
99999 END
C
C
C
      SUBROUTINE GOWFM1(DX,NX,MFILE)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      DIMENSION DX(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /ERRCTL/,/CHAR/
C
      SUB='GOWFM1'
      IF (ICHECK.GE.997) CALL OTRC('GOWFM1','08/16/90')
      IER=0
C
      WRITE (MFILE,'(6E12.5)') (REAL(DX(IX)),IX=1,NX)
C
99999 END
C
C
C
      SUBROUTINE GOWFM2(VX,NX,MFILE)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      DIMENSION VX(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /ERRCTL/,/CHAR/
C
      SUB='GOWFM2'
      IF (ICHECK.GE.997) CALL OTRC('GOWFM2','08/16/90')
      IER=0
C
      WRITE (MFILE,'(6E12.5)') (VX(IX),IX=1,NX)
C
99999 END
