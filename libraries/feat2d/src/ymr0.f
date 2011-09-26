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
* YMR0nn                                                               *
*                                                                      *
* Purpose  CALL MR0nn                                                  *
*          Multigrid restriction routines (Version 0)                  *
*                                                                      *
* Subroutines/functions called  MR0nn                                  *
*                                                                      *
* Version from  08/25/90                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DX2      R*8    Fine grid vector                                     *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX2      R*8    Coarse grid vector                                   *
*                                                                      *
************************************************************************
C
      SUBROUTINE YMR001(DX2,DX1)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNARR=299,NNLEV=9)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION DX1(*),DX2(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /ERRCTL/,/CHAR/,/MGTRD/,/MGTRA/,/MGPAR/
C
      SUB='YMR001'
      IF (ICHECK.GE.998) CALL OTRC('YMR001','08/25/90')
      IER=0
C
C
      LLV2=L(KLVERT(ILEV+1))
      LLA2=L(KLADJ(ILEV+1))
      NVT1=KNVT(ILEV)
      NEL2=KNEL(ILEV+1)
      NEL1=KNEL(ILEV)
C
      CALL MR001(DX2,DX1,KWORK(LLV2),KWORK(LLA2),NVT1,NEL2,NEL1)
C
99999 END
C
C
C
C
C
      SUBROUTINE YMR011(DX2,DX1)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNARR=299,NNLEV=9)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION DX1(*),DX2(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /ERRCTL/,/CHAR/,/MGTRA/,/MGPAR/
C
      SUB='YMR011'
      IF (ICHECK.GE.998) CALL OTRC('YMR011','08/25/90')
      IER=0
C
      LLV2=L(KLVERT(ILEV+1))
      LLA2=L(KLADJ(ILEV+1))
      NVT1=KNVT(ILEV)
      NEL2=KNEL(ILEV+1)
C
      CALL MR011(DX2,DX1,KWORK(LLV2),KWORK(LLA2),NVT1,NEL2)
C
99999 END
