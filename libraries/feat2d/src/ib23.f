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
* IB23n                                                                *
*                                                                      *
* Purpose  Smoothing using Gauss-Seidel-method                         *
*          Single/double precision version                             *
*                                                                      *
* Subroutines/functions called  none                                   *
*                                                                      *
* Version from  01/09/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* VA       R*4    Matrix                                               *
* KCOL     I*4    Pointer vectors for the matrix VA corresponding to   *
* KLD      I*4    the storage technique n                              *
* KOP      I*4                                                         *
* DX       R*8    Starting vector                                      *
* DB       R*8    Right hand side                                      *
* NEQ      I*4    Number of equations                                  *
* NIT      I*4    Number of smoothing steps                            *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8    Smoothed vector                                      *
*                                                                      *
************************************************************************
C
      SUBROUTINE IB237(VA,KCOL,KLD,DX,DB,NEQ,NIT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),DX(*),DB(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('IB237 ','01/09/89')
C
      DO 1 ITE=1,NIT
      DO 2 IEQ=1,NEQ
      AUX=DB(IEQ)
      DO 3 ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
3     AUX=AUX-DBLE(VA(ICOL))*DX(KCOL(ICOL))
      AUX=AUX/DBLE(VA(KLD(IEQ)))
2     DX(IEQ)=AUX
1     CONTINUE
C
      END
C
C
C
      SUBROUTINE IB23A(VA,KCOL,KLD,KOP,DX,DB,NEQ,NIT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),KOP(*),DX(*),DB(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('IB23A ','01/09/89')
C
      DO 1 ITE=1,NIT
      DO 2 IEQ=1,NEQ
      AUX=DB(IEQ)
      JOP=KOP(IEQ)
      DO 3 ICOL=KLD(JOP)+1,KLD(JOP+1)-1
3     AUX=AUX-DBLE(VA(ICOL))*DX(KCOL(ICOL)+IEQ)
      AUX=AUX/DBLE(VA(KLD(JOP)))
2     DX(IEQ)=AUX
1     CONTINUE
C
      END
