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
* XMVA0                                                                *
*                                                                      *
* Purpose  Calculate matrices corresponding to linear forms            *
*          (multigrid version)                                         *
*          Successive call of XVA0                                     *
*                                                                      *
* Subroutines/functions called  XVA0                                   *
*                                                                      *
* Version from  11/22/90                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NLEV     I*4    Number of levels                                     *                                                 * NLMAX    I*4    for which vectors are needed                         *                                                 * KLB      I*4    NBLOC*NLEV numbers of vectors                        *
* NBLOC    I*4                                                         *
* ICLEAR   I*4                                                         *
* ELE      SUBR                                                        *
* COEFF    SUBR                                                        *
* BCON     L*4                                                         *
* KB       I*4    See parameters of XVA0 and VA0                       *
* KBN      I*4                                                         *
* ICUB     I*4                                                         *
* CARR     C*6    Names of vectors for all levels                      *
* BSNGL    L*4    NBLOC*NLEV logical values                            *
*                 .TRUE. meand that corresponding array is converted   *
*                 to single precision after completion                 *
* Meshes on COMMON /MGTRD/ and /MGTRA/                                 *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* IER      I*4    Error indicator                                      *
*                                                                      *
************************************************************************
C
      SUBROUTINE XMVA0(KLB,KNEQ,NBLOC,ICLEAR,ELE,
     *                 COEFF,BCON,KB,KBN,ICUB,CARR,BSNGL)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CARR*6
C
      PARAMETER (NNARR=299,NNAB=21,NNLEV=9)
      DIMENSION KLB(NBLOC,*),KNEQ(*),KB(2,NNAB,*),KBN(*),BCON(*)
      DIMENSION CARR(NBLOC,*),BSNGL(NBLOC,*)
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /ERRCTL/ IER,ICHECK
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
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /TRIAA/,/TRIAD/,/MGTRD/,/MGTRA/,/MGPAR/,/ERRCTL/,/CHAR/
C
      SUB='XMVA0'
      IF (ICHECK.GE.997) CALL OTRC('XMVA0 ','11/22/90')
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
      CALL XVA0(KLB(1,ILEV),KNEQ(ILEV),NBLOC,ICLEAR,ELE,
     *          COEFF,BCON,KB,KBN,ICUB,CARR(1,ILEV),BSNGL(1,ILEV))
      IF (IER.NE.0) GOTO 99999
C     
10    CONTINUE
C     
99999 END
      
      
