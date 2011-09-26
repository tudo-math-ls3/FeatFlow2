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
* YID12n                                                               *
*                                                                      *
* Purpose  SSOR Preconditioning                                        *
*          Matrix stored in technique  n  (see Reference Manual)       *
*          Single/Single precision version                             *
*                                                                      *
* Subroutines/functions called   ID12n                                 *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* VX       R*4    Vector                                               *
* NEQ      I*4    Number of equations (length of VX)                   *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* VX       R*4    Resulting vector                                     *
*                                                                      *
************************************************************************
C
      SUBROUTINE YID127(VX,NEQ)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION VX(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /XYPAR/,/ERRCTL/
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      IF (ICHECK.EQ.999) CALL OTRC('YID127','01/02/89')
C
      CALL ID127(VWORK(KXYPAR(1)),KWORK(KXYPAR(2)),
     *           KWORK(KXYPAR(3)),VX,NEQ,DXYPAR(1))
      END
C
C
C
      SUBROUTINE YID128(VX,NEQ)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION VX(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /XYPAR/,/ERRCTL/
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      IF (ICHECK.EQ.999) CALL OTRC('YID128','01/02/89')
C
      CALL ID128(VWORK(KXYPAR(1)),KWORK(KXYPAR(2)),
     *           KWORK(KXYPAR(3)),VX,NEQ,DXYPAR(1))
      END
C
C
C
      SUBROUTINE YID12A(VX,NEQ)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION VX(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /XYPAR/,/ERRCTL/
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      IF (ICHECK.EQ.999) CALL OTRC('YID12A','01/02/89')
C
      CALL ID12A(VWORK(KXYPAR(1)),KWORK(KXYPAR(2)),
     *           KWORK(KXYPAR(3)),KWORK(KXYPAR(4)),VX,NEQ,DXYPAR(1))
      END
