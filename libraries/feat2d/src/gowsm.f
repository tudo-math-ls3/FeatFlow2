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
* XGOWSM, GOWSM                                                        *
*                                                                      *
* Purpose  Write subdivision in MOVIE.BYU format onto unit MFILE       *   
*                                                                      *
* Subroutines/functions called  OF0                                    *
*                                                                      *
* Version from  08/16/90                                               *
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
*                -111  MFILE  exceeds maximum range 16,...,80          *
*                -111  Write error on unit  MFILE                      *
*                -112  File  CFILE  could not be opened                *
*                                                                      *
************************************************************************
C
      SUBROUTINE XGOWSM(MFILE,CFILE)
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
      SUB='XGOWSM'
      IF (ICHECK.GE.997) CALL OTRC('XGOWSM','08/16/90')
      IER=0
C
C *** Open I/O-file
C
      CALL OF0(MFILE,CFILE,1)
      IF (IER.NE.0) GOTO 99999
C
      CALL GOWSM(DWORK(L(LCORVG)),KWORK(L(LVERT)),MFILE)
C
99999 END
C
C
C
      SUBROUTINE GOWSM(DCORVG,KVERT,MFILE)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNVE=4)
      DIMENSION DCORVG(2,*),KVERT(NNVE,*),VAUX(6),KAUX(10)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      SAVE VAUX,KAUX   
      SAVE /ERRCTL/,/CHAR/,/TRIAD/
      DATA VAUX/6*0./
C
      SUB='GOWSM'
      IF (ICHECK.GE.997) CALL OTRC('GOWSM ','08/16/90')
      IER=0
C
      NP=1
      NJ=NVT
      NPT=NEL
      NEDGE=NEL*NVE
C
      WRITE (MFILE,10001,ERR=99997) NP,NJ,NPT,NEDGE
      WRITE (MFILE,10001,ERR=99997) 1,NEL
C
      DO 1 IEL=1,NEL
      KVERT(NVE,IEL)=-KVERT(NVE,IEL)
1     CONTINUE
C
      DO 2 IVT=2,NVT,2
      VAUX(1)=REAL(DCORVG(1,IVT-1))
      VAUX(2)=REAL(DCORVG(2,IVT-1))
      VAUX(4)=REAL(DCORVG(1,IVT))
      VAUX(5)=REAL(DCORVG(2,IVT))
2     WRITE (MFILE,10002,ERR=99997) (VAUX(I),I=1,6)
      IF (MOD(NVT,2).NE.0) THEN
       VAUX(1)=REAL(DCORVG(1,NVT))
       VAUX(2)=REAL(DCORVG(2,NVT))
       WRITE (MFILE,10002,ERR=99997) (VAUX(I),I=1,3)
      ENDIF
C
      IP=1
      DO 3 IEL=1,NEL
      DO 3 IVE=1,NVE
      KAUX(IP)=KVERT(IVE,IEL)
      IF (IP.EQ.10) THEN
       IP=1
       WRITE (MFILE,10001,ERR=99997) (KAUX(I),I=1,10)
      ELSE
       IP=IP+1
      ENDIF
3     CONTINUE
      IF (IP.GT.1) WRITE (MFILE,10001,ERR=99997) (KAUX(I),I=1,IP-1)
C
      DO 4 IEL=1,NEL
      KVERT(NVE,IEL)=-KVERT(NVE,IEL)
4     CONTINUE
C
      GOTO 99999
C
99997 WRITE (CPARAM,'(I15)') MFILE
      CALL WERR(-111,'XOWSM ')
C
10001 FORMAT(10I8)
10002 FORMAT(6E12.5)
C
99999 END
