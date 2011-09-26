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
* IC23n                                                                *
*                                                                      *
* Purpose  Smoothing using SOR-method                                  *
*          Single/double precision version                             *
*                                                                      *
* Subroutines/functions called  none                                   *
*                                                                      *
* Version from  10/27/89                                               *
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
* NIT      I*4    Maximum number of smoothing steps                    *
* OMEGA    R*8    Relaxation parameter                                 *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8    Smoothed vector                                      *
*                                                                      *
************************************************************************
C
      SUBROUTINE IC237(VA,KCOL,KLD,DX,DB,NEQ,NIT,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),DX(*),DB(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('IC237 ','10/27/89')
C
      DO 1 ITE=1,NIT
      DO 2 IEQ=1,NEQ
      AUX=DB(IEQ)
      DO 3 ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
3     AUX=AUX-DBLE(VA(ICOL))*DX(KCOL(ICOL))
      AUX=OMEGA*(AUX/DBLE(VA(KLD(IEQ)))-DX(IEQ))+DX(IEQ)
2     DX(IEQ)=AUX
1     CONTINUE
C
      END
C
C
C
      SUBROUTINE IC23A(VA,KCOL,KLD,KOP,DX,DB,NEQ,NIT,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),KOP(*),DX(*),DB(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.997) CALL OTRC('IC23A ','10/27/89')
C
      DO 1 ITE=1,NIT
      DO 2 IEQ=1,NEQ
      AUX=DB(IEQ)
      JOP=KOP(IEQ)
      DO 3 ICOL=KLD(JOP)+1,KLD(JOP+1)-1
3     AUX=AUX-DBLE(VA(ICOL))*DX(KCOL(ICOL)+IEQ)
      AUX=OMEGA*(AUX/DBLE(VA(KLD(JOP)))-DX(IEQ))+DX(IEQ)
2     DX(IEQ)=AUX
1     CONTINUE
C
      END
