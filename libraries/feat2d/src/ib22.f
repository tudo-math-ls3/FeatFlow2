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
* IB22n                                                                *
*                                                                      *
* Purpose  Smoothing using Gauss-Seidel-method                         *
*          Single/single precision version                             *
*                                                                      *
* Subroutines/functions called  none                                   *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* VA       R*4    Matrix                                               *
* KCOL     I*4    Pointer vectors for the matrix VA corresponding to   *
* KLD      I*4    the storage technique n                              *
* KOP      I*4                                                         *
* VX       R*4    Starting vector                                      *
* VB       R*4    Right hand side                                      *
* NEQ      I*4    Number of equations                                  *
* NIT      I*4    Number of smoothing steps                            *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* VX       R*4    Smoothed vector                                      *
*                                                                      *
************************************************************************
C
      SUBROUTINE IB227(VA,KCOL,KLD,VX,VB,NEQ,NIT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),VX(*),VB(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('IB227 ','01/02/89')
C
      DO 1 ITE=1,NIT
      DO 2 IEQ=1,NEQ
      AUX=VB(IEQ)
      DO 3 ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
3     AUX=AUX-VA(ICOL)*VX(KCOL(ICOL))
      AUX=AUX/VA(KLD(IEQ))
2     VX(IEQ)=AUX
1     CONTINUE
C
      END
C
C
C
      SUBROUTINE IB22A(VA,KCOL,KLD,KOP,VX,VB,NEQ,NIT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),KOP(*),VX(*),VB(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('IB22A ','01/02/89')
C
      DO 1 ITE=1,NIT
      DO 2 IEQ=1,NEQ
      AUX=VB(IEQ)
      JOP=KOP(IEQ)
      DO 3 ICOL=KLD(JOP)+1,KLD(JOP+1)-1
3     AUX=AUX-VA(ICOL)*VX(KCOL(ICOL)+IEQ)
      AUX=AUX/VA(KLD(JOP))
2     VX(IEQ)=AUX
1     CONTINUE
C
      END
