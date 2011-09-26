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
* YLTX2n   (corresponding Call subroutines used, e.g., in XIE020)      *
*                                                                      *
* Purpose  Matrix vector product with the transposed matrix            *
*          Matrix stored in technique  n  (see Reference Manual)       *
*          Single/Single precision version                             *
*                                                                      *
* Subroutines/functions called   LTX2n                                 *
*                                                                      *
* Version from  11/15/89                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* VX       R*4    Vector                                               *
* NEQ      I*4    Number of equations (length of VX, VTX)              *
*                 For storage technique 9: NEQ1 (length of VX) and     *
*                 NEQ2 (length of VTX) are needed                      *
* A1,A2    R*8    VTX := A1*VAT*VX + A2*VTX                            *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* VTX      R*4    Resulting vector                                     *
*                                                                      *
************************************************************************
C
      SUBROUTINE YLTX23(VX,VTX,NEQ,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION VX(*),VTX(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
      COMMON /ERRCTL/ IER,ICHECK
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /XYPAR/,/ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('YLTX23','11/29/90')
C
      CALL LTX23(VWORK(KXYPAR(1)),KWORK(KXYPAR(2)),
     *           KWORK(KXYPAR(3)),KXYPAR(4),NEQ,VX,VTX,A1,A2)
      END
C
C
C
      SUBROUTINE YLTX27(VX,VTX,NEQ,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION VX(*),VTX(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
      COMMON /ERRCTL/ IER,ICHECK
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /XYPAR/,/ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('YLTX27','01/02/89')
C
      CALL LTX27(VWORK(KXYPAR(1)),KWORK(KXYPAR(2)),
     *           KWORK(KXYPAR(3)),NEQ,VX,VTX,A1,A2)
      END
C
C
C
      SUBROUTINE YLTX29(VX,VTX,NEQ1,NEQ2,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION VX(*),VTX(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
      COMMON /ERRCTL/ IER,ICHECK
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /XYPAR/,/ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('YLTX29','11/15/89')
C
      CALL LTX29(VWORK(KXYPAR(1)),KWORK(KXYPAR(2)),
     *           KWORK(KXYPAR(3)),NEQ1,NEQ2,VX,VTX,A1,A2)
      END
C
C
C
      SUBROUTINE YLTX2A(VX,VTX,NEQ,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION VX(*),VTX(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
      COMMON /ERRCTL/ IER,ICHECK
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /XYPAR/,/ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('YLTX2A','01/02/89')
C
      CALL LTX2A(VWORK(KXYPAR(1)),KWORK(KXYPAR(2)),
     *           KWORK(KXYPAR(3)),KWORK(KXYPAR(4)),
     *           NEQ,VX,VTX,A1,A2)
      END
