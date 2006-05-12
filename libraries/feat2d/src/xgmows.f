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
* XMOWS                                                                *
*                                                                      *
* Purpose  Write multiple triangulations (multigrid version)           *
*          in MOVIE.BYU format                                         *
*                                                                      *
* Subroutines/functions called  XGOWSM                                 *
*                                                                      *
* Version from  08/25/90                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NLEV     I*4    Vector of numbers of the arrays                      *
* MFILE    I*4    Unit number used for output                          *
* CCFILE   C*(*)  Array of filenames used for output                   *
* Meshes on COMMON /MGTRD/ and /MGTRA/                                 *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* IER      I*4    Error indicator                                      *
*                 Set by OF0 or OWA                                    *
*                                                                      *
************************************************************************
C
      SUBROUTINE XGMOWS(MFILE,CCFILE)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER*(*) CCFILE(*)
C
      PARAMETER (NNARR=299,NNAB=21,NNDER=6,NNLEV=9)
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /TRIAA/,/TRIAD/,/MGTRD/,/MGTRA/,/MGPAR/,/ERRCTL/,/CHAR/
C
      SUB='XGMOWS'
      IF (ICHECK.GE.997) CALL OTRC('XGMOWS','08/25/90')
      IER=0
C
      DO 10 ILEV=1,NLEV
C
      NEL =KNEL(ILEV)
      NVT =KNVT(ILEV)
      NMT =KNMT(ILEV)
      NVEL=KNVEL(ILEV)
      NVBD=KNVBD(ILEV)
C
      LCORVG=KLCVG(ILEV)
      LCORMG=KLCMG(ILEV)
      LVERT =KLVERT(ILEV)
      LMID  =KLMID(ILEV)
      LADJ  =KLADJ(ILEV)
      LVEL  =KLVEL(ILEV)
      LMEL  =KLMEL(ILEV)
      LNPR  =KLNPR(ILEV)
      LMM   =KLMM(ILEV)
      LVBD  =KLVBD(ILEV)
      LEBD  =KLEBD(ILEV)
      LBCT  =KLBCT(ILEV)
      LVBDP =KLVBDP(ILEV)
      LMBDP =KLMBDP(ILEV)
C
      CALL XGOWSM(MFILE,CCFILE(ILEV))
      IF (IER.NE.0) GOTO 99999
      CLOSE(MFILE)
C     
10    CONTINUE   
C     
99999 END
      
      
