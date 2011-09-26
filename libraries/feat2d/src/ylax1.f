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
* YLAX1n   (corresponding Call subroutines used, e.g., in XIE010)      *
*                                                                      *
* Purpose  Matrix vector product                                       *
*          Matrix stored in technique  n  (see Reference Manual)       *
*          Double/double precision version                             *
*                                                                      *
* Subroutines/functions called   LAX1n                                 *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DX       R*8    Vector                                               *
* NEQ      I*4    Number of equations (length of DX, DAX)              *
* A1,A2    R*8    DAX := A1*DA*DX + A2*DAX                             *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DAX      R*8    Resulting vector                                     *
*                                                                      *
************************************************************************
C
      SUBROUTINE YLAX13(DX,DAX,NEQ,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION DX(*),DAX(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
      COMMON /ERRCTL/ IER,ICHECK
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /XYPAR/,/ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('YLAX13','11/29/90')
C
      CALL LAX13(DWORK(KXYPAR(1)),KWORK(KXYPAR(2)),
     *           KWORK(KXYPAR(3)),KXYPAR(4),NEQ,DX,DAX,A1,A2)
      END
C
C
C
      SUBROUTINE YLAX14(DX,DAX,NEQ,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION DX(*),DAX(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
      COMMON /ERRCTL/ IER,ICHECK
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /XYPAR/,/ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('YLAX14','01/07/91')
C
      CALL LAX14(DWORK(KXYPAR(1)),KWORK(KXYPAR(2)),
     *           KWORK(KXYPAR(3)),KXYPAR(4),NEQ,DX,DAX,A1,A2)
      END
C
C
C
      SUBROUTINE YLAX17(DX,DAX,NEQ,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION DX(*),DAX(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
      COMMON /ERRCTL/ IER,ICHECK
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /XYPAR/,/ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('YLAX17','01/02/89')
C
      CALL LAX17(DWORK(KXYPAR(1)),KWORK(KXYPAR(2)),
     *           KWORK(KXYPAR(3)),NEQ,DX,DAX,A1,A2)
C
      END
C
C
C
      SUBROUTINE YLAX18(DX,DAX,NEQ,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION DX(*),DAX(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /XYPAR/,/ERRCTL/
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      IF (ICHECK.EQ.999) CALL OTRC('YLAX18','01/02/89')
C
      CALL LAX18(DWORK(KXYPAR(1)),KWORK(KXYPAR(2)),
     *           KWORK(KXYPAR(3)),NEQ,DX,DAX,A1,A2)
C
      END
C
C
C
      SUBROUTINE YLAX19(DX,DAX,NEQ,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION DX(*),DAX(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
      COMMON /ERRCTL/ IER,ICHECK
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /XYPAR/,/ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('YLAX19','01/02/89')
C
      CALL LAX19(DWORK(KXYPAR(1)),KWORK(KXYPAR(2)),
     *           KWORK(KXYPAR(3)),NEQ,DX,DAX,A1,A2)
C
      END
C
C
C
      SUBROUTINE YLAX1A(DX,DAX,NEQ,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION DX(*),DAX(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
      COMMON /ERRCTL/ IER,ICHECK
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /XYPAR/,/ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('YLAX1A','01/02/89')
C
      CALL LAX1A(DWORK(KXYPAR(1)),KWORK(KXYPAR(2)),
     *           KWORK(KXYPAR(3)),KWORK(KXYPAR(4)),
     *           NEQ,DX,DAX,A1,A2)
C
      END
