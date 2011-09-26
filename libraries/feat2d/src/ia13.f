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
* IA13n                                                                *
*                                                                      *
* Purpose  Preconditioning by scaling                                  *
*          Matrix stored in technique  n  (see Reference Manual)       *
*          Single/double precision version                             *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* VA       R*4    Matrix stored in technique  n                        *
* KLD      I*4    Pointer vectors corresponding to the                 *
* KOP      I*4    storage technique                                    *
* DX       R*8    Vector                                               *
* NEQ      I*4    Number of equations (length of DX)                   *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8    Resulting vector                                     *
*                                                                      *
************************************************************************
C
      SUBROUTINE IA133(VA,DX,NEQ)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),DX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('IA133 ','11/19/90')
C
      DO 1 IEQ=1,NEQ
1     DX(IEQ)=DX(IEQ)/DBLE(VA(IEQ))
C
      END
C
C
C
      SUBROUTINE IA137(VA,KLD,DX,NEQ)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KLD(*),DX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('IA137 ','01/02/89')
C
      DO 1 IEQ=1,NEQ
1     DX(IEQ)=DX(IEQ)/DBLE(VA(KLD(IEQ)))
C
      END
C
C
C
      SUBROUTINE IA13A(VA,KLD,KOP,DX,NEQ)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KLD(*),DX(*),KOP(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('IA13A ','01/02/89')
C
      DO 1 IEQ=1,NEQ
1     DX(IEQ)=DX(IEQ)/DBLE(VA(KLD(KOP(IEQ))))
C
      END
