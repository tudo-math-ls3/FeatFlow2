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
* XMAPn                                                                *
*                                                                      *
* Purpose  Calculate pointer vectors                                   *
*          (multigrid version)                                         *
*          Successive call of XAPn                                     *
*                                                                      *
* Subroutines/functions called  XAPn                                   *
*                                                                      *
* Version from  02/18/91                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NLEV     I*4    Number of levels                                     *
* ELE      SUBR                                                        *
* ISYMM    I*4                                                         *
* Meshes on COMMON /MGTRD/ and /MGTRA/                                 *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* KNA      I*4    Number of entries for each level                     *
* KNEQ     I*4    Number of unknowns for each level                    *
* IER      I*4    Error indicator                                      *
*                 Set by ZNEW                                          *
*                                                                      *
************************************************************************
C
      SUBROUTINE XMAP3(KLDIA,KLDIAS,KNDIA,KNA,KNEQ,ELE,ISYMM)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNARR=299,NNLEV=9)
      DIMENSION KLDIA(*),KLDIAS(*),KNDIA(*),KNA(*),KNEQ(*)
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
      SUB='XMAP3 '
      IF (ICHECK.GE.997) CALL OTRC('XMAP3 ','02/18/91')
      IER=0
C
      DO 10 ILEV=NLMIN,NLMAX
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
      CALL XAP3(KLDIA(ILEV),KLDIAS(ILEV),KNDIA(ILEV),KNA(ILEV),
     *          KNEQ(ILEV),ELE,ISYMM)
      IF (IER.NE.0) GOTO 99999
C     
10    CONTINUE
C     
99999 END
C
C
C
C
      SUBROUTINE XMAP7(KLCOL,KLLD,KNA,KNEQ,ELE,ISYMM)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNARR=299,NNLEV=9)
      DIMENSION KLCOL(*),KLLD(*),KNA(*),KNEQ(*)
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
      SUB='XMAP7 '
      IF (ICHECK.GE.997) CALL OTRC('XMAP7 ','08/25/90')
      IER=0
C
      DO 10 ILEV=NLMIN,NLMAX
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
      CALL XAP7(KLCOL(ILEV),KLLD(ILEV),KNA(ILEV),KNEQ(ILEV),ELE,ISYMM)
      IF (IER.NE.0) GOTO 99999
C     
10    CONTINUE
C     
99999 END
      
      
