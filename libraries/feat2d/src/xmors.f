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
* XMORS                                                                *
*                                                                      *
* Purpose  Write multiple triangulations (multigrid version)           *
*                                                                      *
* Subroutines/functions called  XORS                                   *
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
      SUBROUTINE XMORS(MFILE,CCFILE,IFMT)
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
      SUB='XMORS'
      IF (ICHECK.GE.997) CALL OTRC('XMORS ','08/25/90')
      IER=0
C
      CALL XMSCL
      IF (IER.NE.0) GOTO 99999
      CALL XSCL
      IF (IER.NE.0) GOTO 99999
C      
      DO 10 ILEV=1,NLEV
C
      NEL =0
      NVT =0
      NMT =0
      NVEL=0
      NVBD=0
C
      LCORVG=0
      LCORMG=0
      LVERT =0
      LMID  =0
      LADJ  =0
      LVEL  =0
      LMEL  =0
      LNPR  =0
      LMM   =0
      LVBD  =0
      LEBD  =0
      LBCT  =0
      LVBDP =0
      LMBDP =0
C
      CALL XORS(MFILE,CCFILE(ILEV),IFMT)
      IF (IER.NE.0) GOTO 99999
      CLOSE(MFILE)
C     
      KNEL(ILEV) =NEL
      KNVT(ILEV) =NVT
      KNMT(ILEV) =NMT
      KNVEL(ILEV)=NVEL
      KNVBD(ILEV)=NVBD
C
      KLCVG(ILEV) =LCORVG
      KLCMG(ILEV) =LCORMG
      KLVERT(ILEV)=LVERT
      KLMID(ILEV) =LMID
      KLADJ(ILEV) =LADJ
      KLVEL(ILEV) =LVEL
      KLMEL(ILEV) =LMEL
      KLNPR(ILEV) =LNPR
      KLMM(ILEV)  =LMM  
      KLVBD(ILEV) =LVBD
      KLEBD(ILEV) =LEBD 
      KLBCT(ILEV) =LBCT
      KLVBDP(ILEV)=LVBDP 
      KLMBDP(ILEV)=LMBDP
C
10    CONTINUE   
C     
99999 END
      
      
