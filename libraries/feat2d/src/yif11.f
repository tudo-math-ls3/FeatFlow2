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
* YIF11n                                                               *
*                                                                      *
* Purpose  Preconditioning by ILU decomposition                        *
*          Matrix stored in technique  n  (see Reference Manual)       *
*          Double/double precision version                             *
*                                                                      *
* Subroutines/functions called   IF11n                                 *
*                                                                      *
* Version from  12/02/89                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DX       R*8    Vector                                               *
* NEQ      I*4    Number of equations (length of DX)                   *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8    Resulting vector                                     *
*                                                                      *
************************************************************************
C
      SUBROUTINE YIF117(DX,NEQ)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION DX(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
      COMMON /ERRCTL/ IER,ICHECK
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /XYPAR/,/ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('YIF117','12/02/89')
C     
      CALL IF117(DWORK(KXYPAR(4)),KWORK(KXYPAR(2)),
     *           KWORK(KXYPAR(3)),DX,NEQ)
      END
