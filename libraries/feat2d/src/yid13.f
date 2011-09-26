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
* YID13n                                                               *
*                                                                      *
* Purpose  SSOR Preconditioning                                        *
*          Matrix stored in technique  n  (see Reference Manual)       *
*          Single/double precision version                             *
*                                                                      *
* Subroutines/functions called   ID13n                                 *
*                                                                      *
* Version from  01/02/89                                               *
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
      SUBROUTINE YID137(DX,NEQ)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION DX(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /XYPAR/,/ERRCTL/
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      IF (ICHECK.EQ.999) CALL OTRC('YID137','01/02/89')
C
      CALL ID137(VWORK(KXYPAR(1)),KWORK(KXYPAR(2)),
     *           KWORK(KXYPAR(3)),DX,NEQ,DXYPAR(1))
      END
C
C
C
      SUBROUTINE YID138(DX,NEQ)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION DX(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /XYPAR/,/ERRCTL/
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      IF (ICHECK.EQ.999) CALL OTRC('YID138','01/02/89')
C
      CALL ID138(VWORK(KXYPAR(1)),KWORK(KXYPAR(2)),
     *           KWORK(KXYPAR(3)),DX,NEQ,DXYPAR(1))
      END
C
C
C
      SUBROUTINE YID13A(DX,NEQ)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION DX(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /XYPAR/,/ERRCTL/
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      IF (ICHECK.EQ.999) CALL OTRC('YID13A','01/02/89')
C
      CALL ID13A(VWORK(KXYPAR(1)),KWORK(KXYPAR(2)),
     *           KWORK(KXYPAR(3)),KWORK(KXYPAR(4)),DX,NEQ,DXYPAR(1))
      END
