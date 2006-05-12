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
*                                                                      *
* Subroutines/functions called  XOWS                                   *
*                                                                      *
* Version from  08/25/90                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NLEV     I*4    Vector of numbers of the arrays                      *
*                 New arrays are allocated on DWORK for LA(IBLOC)=0    *
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
      SUBROUTINE XMOWS(MFILE,CCFILE,IFMT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER*(*) CCFILE(*)
C
      PARAMETER (NNLEV=9)
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
      SAVE /TRIAA/,/TRIAD/,/MGTRD/,/MGTRA/,/MGPAR/,/ERRCTL/,/CHAR/
C
      SUB='XMOWS'
      IF (ICHECK.GE.997) CALL OTRC('XMOWS ','08/25/90')
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
      CALL XOWS(MFILE,CCFILE(ILEV),IFMT)
      IF (IER.NE.0) GOTO 99999
      CLOSE(MFILE)
C     
10    CONTINUE   
C     
99999 END
      
      
