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
* YMP1nn                                                               *
*                                                                      *
* Purpose  CALL MP1nn                                                  *
*          Multigrid prolongation routines (Version 1)                 *
*                                                                      *
* Subroutines/functions called  MP1nn                                  *
*                                                                      *
* Version from  04/12/91                                               *
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
      SUBROUTINE YMP111(DX1,DX2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNARR=299,NNLEV=9)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION DX1(*),DX2(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /MACROA/ LMACVG,LMACMG,LMAVT,LMAMID,LMAADJ,LMAVEL,LMAMEL,
     *                LMANPR,LMAMM,LMAVBD,LMAEBD,LMABCT,LMAVBP,LMAMBP,
     *                LMAVE
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
      SAVE /ERRCTL/,/CHAR/,/MACROA/,/MGTRD/,/MGTRA/,/MGPAR/
C
      SUB='YMP111'
      IF (ICHECK.GE.998) CALL OTRC('YMP111','04/12/91')
      IER=0
C
      LLV1  =L(KLVERT(ILEV-1))
      LLV2  =L(KLVERT(ILEV))
      LLMAVT=L(LMAVT)
      LLMADJ=L(LMAVT)
      LLMAVE=L(LMAVE)
      NEL1=KNEL(ILEV-1)
C
      CALL MP111(DX1,DX2,NEL1,KWORK(LLV1),KWORK(LLV2),
     *           KWORK(LLMAVT),KWORK(LLMADJ),KWORK(LLMAVE))
C
99999 END
