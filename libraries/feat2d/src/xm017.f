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
* XM017                                                                *
*                                                                      *
* Purpose  Save information on system matrices for all levels          *
*          Call XM010                                                  *
*                                                                      *
* Subroutines/functions called  XM010                                  *
*                                                                      *
* Version from  08/25/90                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* KLA      I*4    Numbers of arrys characterisong the matrix A         *
* KLCOL    I*4    for every level between 1,...,NLMAX                  *
* KLLD     I*4                                                         *
* The remaining parameters are those of XM010                          *
*                                                                      *
*                                                                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
*                                                                      *
************************************************************************
C
      SUBROUTINE XM017(KLA,KLCOL,KLLD,LX,LB,KNEQ,NIT,ITE,
     *                 EPS,YAX,YPROL,YREST,YPRSM,YPOSM,YEX,YBC,IDISP)
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B) 
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNLEV=9,NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION KLA(*),KLCOL(*),KLLD(*)
      DIMENSION KNEQ(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL YAX,YPROL,YREST,YPRSM,YPOSM,YEX,YBC
      SAVE /ERRCTL/,/CHAR/,/MGPAR/,/XYPAR/
C
      SUB='XM017 '
      IF (ICHECK.GE.997) CALL OTRC('XM017 ','08/25/90')
      IER=0
C
      IP=100     
      DO 1 ILEV=1,NLMAX
      KXYPAR(IP+1)=L(KLA(ILEV))
      KXYPAR(IP+2)=L(KLCOL(ILEV))
      KXYPAR(IP+3)=L(KLLD(ILEV))
      IP=IP+3
1     CONTINUE
C
      CALL XM010(LX,LB,KNEQ,NIT,ITE,EPS,
     *           YAX,YPROL,YREST,YPRSM,YPOSM,YEX,YBC,IDISP,0)
C
      END
