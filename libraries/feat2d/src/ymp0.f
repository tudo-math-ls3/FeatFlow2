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
* YMP0nn                                                               *
*                                                                      *
* Purpose  CALL MP0nn                                                  *
*          Multigrid prolongation routines (Version 0)                 *
*                                                                      *
* Subroutines/functions called  MP0nn                                  *
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
      SUBROUTINE YMP001(DX1,DX2)
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
      SUB='YMP001'
      IF (ICHECK.GE.998) CALL OTRC('YMP001','08/25/90')
      IER=0
C
C
      LLV1=L(KLVERT(ILEV-1))
      LLV2=L(KLVERT(ILEV))
      LLA1=L(KLADJ(ILEV-1))
      NVT1=KNVT(ILEV-1)
      NEL1=KNEL(ILEV-1)
C
      CALL MP001(DX1,DX2,KWORK(LLV1),KWORK(LLV2),KWORK(LLA1),NVT1,NEL1)
C
99999 END
C
C
C
C
C
      SUBROUTINE YMP011(DX1,DX2)
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
      SUB='YMP011'
      IF (ICHECK.GE.998) CALL OTRC('YMP011','08/25/90')
      IER=0
C
C
      LLV1=L(KLVERT(ILEV-1))
      LLV2=L(KLVERT(ILEV))
      LLA1=L(KLADJ(ILEV-1))
      LLA2=L(KLADJ(ILEV))
      NVT1=KNVT(ILEV-1)
      NEL1=KNEL(ILEV-1)
C
      CALL MP011(DX1,DX2,KWORK(LLV1),KWORK(LLV2),KWORK(LLA1),
     *           KWORK(LLA2),NVT1,NEL1)
C
99999 END
