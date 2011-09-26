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
* IA22n                                                                *
*                                                                      *
* Purpose  Smoothing using (damped) Jacobi-method                      *
*          Single/single precision version                             *
*                                                                      *
* Subroutines/functions called  LCP2, LLC2, LVM2                       *
*                                                                      *
* Version from  10/26/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* VA       R*4    Matrix                                               *
* KCOL     I*4    Pointer vectors for the matrix VA corresponding to   *
* KLD      I*4    the storage technique n                              *
* KOP      I*4                                                         *
* VX       R*4    Starting vector                                      *
* VB       R*4    Right hand side                                      *
* VD       R*4    Help vector                                          *
* NEQ      I*4    Number of equations                                  *
* NIT      I*4    Number of smoothing steps                            *
* OMEGA    R*8    Damping factor                                       *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* VX       R*4    Smoothed vector                                      *
*                                                                      *
************************************************************************
C
      SUBROUTINE IA223(VA,KDIA,KDIAS,NDIA,VX,VB,VD,NEQ,NIT,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KDIA(*),KDIAS(*),VX(*),VB(*),VD(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('IA223 ','11/19/90')
C
      DO 1 IEQ=1,NEQ
1     VA(IEQ)=1D0/VA(IEQ)
C
      DO 2 ITE=1,NIT
      CALL LCP2(VB,VD,NEQ)
      DO 2 IDIA=2,NDIA
      J1=KDIA(IDIA)
      IF (J1.GT.0) THEN
       I1=1
       NEQ1=NEQ-J1
      ELSE
       I1=1-J1
       NEQ1=NEQ+J1
      ENDIF
      J0=KDIAS(IDIA)-I1
      CALL LVM2(VX(I1+J1),VA(I1+J0),VD(I1),NEQ1,-1D0,1D0)
3     CONTINUE
      CALL LVM2(VD,VA,VD,NEQ,1D0,0D0)
      CALL LLC2(VD,VX,NEQ,1D0-OMEGA,OMEGA)
2     CONTINUE
C
      DO 10 IEQ=1,NEQ
10    VA(IEQ)=1D0/VA(IEQ)
C
      END
C
C
C
      SUBROUTINE IA227(VA,KCOL,KLD,VX,VB,VD,NEQ,NIT,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),VX(*),VB(*),VD(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('IA227 ','10/26/89')
C
      DO 1 ITE=1,NIT
      CALL LCP2(VB,VD,NEQ)
      DO 2 IEQ=1,NEQ
      DO 3 ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
3     VD(IEQ)=VD(IEQ)-VA(ICOL)*VX(KCOL(ICOL))
2     CONTINUE
      DO 4 IEQ=1,NEQ
4     VX(IEQ)=(1D0-OMEGA)*VX(IEQ)+OMEGA*VD(IEQ)/VA(KLD(IEQ))
1     CONTINUE
C
      END
C
C
C
      SUBROUTINE IA22A(VA,KCOL,KLD,KOP,VX,VB,VD,NEQ,NIT,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),KOP(*),VX(*),VB(*),VD(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('IA22A ','10/26/89')
C
      DO 1 ITE=1,NIT
      CALL LCP2(VB,VD,NEQ)
      DO 2 IEQ=1,NEQ
      JOP=KOP(IEQ)
      DO 3 ICOL=KLD(JOP)+1,KLD(JOP+1)-1
3     VD(IEQ)=VD(IEQ)-VA(ICOL)*VX(KCOL(ICOL)+IEQ)
2     CONTINUE
      DO 4 IEQ=1,NEQ
4     VX(IEQ)=(1D0-OMEGA)*VX(IEQ)+OMEGA*VD(IEQ)/VA(KLD(KOP(IEQ)))
1     CONTINUE
C
      END
