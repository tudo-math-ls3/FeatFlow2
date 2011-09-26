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
* YLTX3n   (corresponding Call subroutines used, e.g., in XIE010)      *
*                                                                      *
* Purpose  Matrix vector product with the transposed matrix            *
*          Matrix stored in technique  n  (see Reference Manual)       *
*          Single/double precision version                             *
*                                                                      *
* Subroutines/functions called   LTX3n                                 *
*                                                                      *
* Version from  11/15/89                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DX       R*8    Vector                                               *
* NEQ      I*4    Number of equations (length of DX, DTX)              *
*                 For storage technique 9: NEQ1 (length of DX) and     *
*                 NEQ2 (length of DTX) are needed                      *
* A1,A2    R*8    DTX := A1*VAT*DX + A2*DTX                            *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DTX      R*8    Resulting vector                                     *
*                                                                      *
************************************************************************
C
      SUBROUTINE YLTX33(DX,DTX,NEQ,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION DX(*),DTX(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)      
      COMMON /ERRCTL/ IER,ICHECK
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /XYPAR/,/ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('YLTX33','11/29/90')
C
      CALL LTX33(VWORK(KXYPAR(1)),KWORK(KXYPAR(2)),
     *           KWORK(KXYPAR(3)),KXYPAR(4),NEQ,DX,DTX,A1,A2)
      END
C
C
C
      SUBROUTINE YLTX37(DX,DTX,NEQ,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION DX(*),DTX(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)      
      COMMON /ERRCTL/ IER,ICHECK
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /XYPAR/,/ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('YLTX37','01/02/89')
C
      CALL LTX37(VWORK(KXYPAR(1)),KWORK(KXYPAR(2)),
     *           KWORK(KXYPAR(3)),NEQ,DX,DTX,A1,A2)
      END
C
C
C
      SUBROUTINE YLTX39(DX,DTX,NEQ1,NEQ2,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION DX(*),DTX(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
      COMMON /ERRCTL/ IER,ICHECK
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /XYPAR/,/ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('YLTX39','11/15/89')
C
      CALL LTX39(VWORK(KXYPAR(1)),KWORK(KXYPAR(2)),
     *           KWORK(KXYPAR(3)),NEQ1,NEQ2,DX,DTX,A1,A2)
      END
C
C
C
      SUBROUTINE YLTX3A(DX,DTX,NEQ,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION DX(*),DTX(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
      COMMON /ERRCTL/ IER,ICHECK
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /XYPAR/,/ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('YLTX3A','01/02/89')
C
      CALL LTX3A(VWORK(KXYPAR(1)),KWORK(KXYPAR(2)),
     *           KWORK(KXYPAR(3)),KWORK(KXYPAR(4)),
     *           NEQ,DX,DTX,A1,A2)
      END
