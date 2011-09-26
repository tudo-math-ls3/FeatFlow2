************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.3)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek                     *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
* based on the 2-D routine, modified by P.Schreiber                 *
* XMAB07                                                               *
*                                                                      *
* Purpose  Calculate matrices corresponding to bilinear forms          *
*          (multigrid version)                                         *
*          Successive call of XAB07                                    *
*                                                                      *
* Subroutines/functions called  XAB07                                  *
*                                                                      *
* Version from  12/04/90                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NLEV     I*4    Number of levels                                     *
*                 for which matrices are generated                     *
* KLA      I*4    NBLOC*NLEV numbers of arrays                         *
* KLCOL    I*4    NLEV numbers of pointer vectors                      * 
* KLLD     I*4                                                         *
* NBLOC    I*4                                                         *
* ICLEAR   I*4                                                         *
* ELE      SUBR                                                        *
* COEFF    SUBR                                                        *
* BCON     L*4                                                         *
* KAB      I*4    See parameters of XAB07 and AB07                     *
* KABN     I*4                                                         *
* ICUB     I*4                                                         *
* ISYMM    I*4                                                         *
* CARR     C*6    Names of matrices for all levels                     *
* BSNGL    L*4    NBLOC*NLEV logical values                            *
*                 .TRUE. meand that corresponding array is converted   *
*                 to single precision after completion                 *
* Meshes on COMMON /MGTRD/ and /MGTRA/                                 *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* KLA      I*4    NBLOC*NLEV numbers of arrays on DA                   *
*                                                                      *
************************************************************************
C
      SUBROUTINE XMAB07(KLA,KLCOL,KLLD,KNA,KNEQ,NBLCA,ICLR,ELE,
     *                  COEFF,BCON,KAB,KABN,ICUB,ISYMM,ILINT,
     *                  BSNGL,CARR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CARR*6
C
      PARAMETER (NNARR=299,NNAB=21,NNLEV=9)
      DIMENSION KLA(NBLCA,*),KLCOL(*),KLLD(*)
      DIMENSION KNA(*),KNEQ(*),KAB(2,NNAB,*),KABN(*),BCON(*)
      DIMENSION CARR(NBLCA,*),BSNGL(NBLCA,*)
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /TRIAA/,/TRIAD/,/MGTRD/,/MGTRA/,/MGPAR/,/ERRCTL/,/CHAR/
C
      SUB='XMAB07'
      IF (ICHECK.GE.997) CALL OTRC('XMAB07','12/04/90')
      IER=0
C
C
      DO 10 ILEV=NLMIN,NLMAX
C
      NEL =KNEL(ILEV)    
      NVT =KNVT(ILEV)
      NET =KNET(ILEV)
      NAT =KNAT(ILEV)
      NVE =KNVE(ILEV)
      NEE =KNEE(ILEV)
      NAE =KNAE(ILEV)
      NVEL =KNVEL(ILEV)
      NEEL =KNEEL(ILEV)
      NVED =KNVED(ILEV)
      NVAR =KNVAR(ILEV)
      NEAR =KNEAR(ILEV)
      NBCT =KNBCT(ILEV)
      NVBD =KNVBD(ILEV)
      NEBD =KNEBD(ILEV)
      NABD =KNABD(ILEV)
C
      LCORVG =KLCVG(ILEV) 
      LCORMG =KLCMG(ILEV) 
      LCORAG =KLCAG(ILEV) 
      LVERT  =KLVERT(ILEV)
      LEDGE  =KLEDGE(ILEV) 
      LAREA  =KLAREA(ILEV) 
      LADJ   =KLADJ(ILEV) 
      LVEL   =KLVEL(ILEV) 
      LEEL   =KLEEL(ILEV)
      LAEL   =KLAEL(ILEV)   
      LVED   =KLVED(ILEV) 
      LAED   =KLAED(ILEV) 
      LVAR   =KLVAR(ILEV) 
      LEAR   =KLEAR(ILEV) 
      LEVE   =KLEVE(ILEV) 
      LAVE   =KLAVE(ILEV) 
      LNPR   =KLNPR(ILEV) 
      LBCT   =KLBCT(ILEV) 
      LVBD   =KLVBD(ILEV) 
      LEBD   =KLEBD(ILEV) 
      LABD   =KLABD(ILEV) 
C
      CALL XAB07(KLA(1,ILEV),KLCOL(ILEV),KLLD(ILEV),KNA(ILEV),
     *           KNEQ(ILEV),NBLCA,ICLR,ELE,COEFF,BCON,
     *           KAB,KABN,ICUB,ISYMM,ILINT,
     *           BSNGL,CARR(1,ILEV))
      IF (IER.NE.0) GOTO 99999
C
10    CONTINUE
C     
99999 END
