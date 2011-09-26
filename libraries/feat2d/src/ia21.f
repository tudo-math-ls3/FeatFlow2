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
* IA21n                                                                *
*                                                                      *
* Purpose  Smoothing using (damped) Jacobi-method                      *
*          Double/double precision version                             *
*                                                                      *
* Subroutines/functions called  LCP1, LLC1, LVM1                       *
*                                                                      *
* Version from  10/26/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DA       R*8    Matrix                                               *
* KCOL     I*4    Pointer vectors for the matrix DA corresponding to   *
* KLD      I*4    the storage technique n                              *
* KOP      I*4                                                         *
* DX       R*8    Starting vector                                      *
* DB       R*8    Right hand side                                      *
* DD       R*8    Help vector                                          *
* NEQ      I*4    Number of equations                                  *
* NIT      I*4    Number of smoothing steps                            *
* OMEGA    R*8    Damping factor                                       *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8    Smoothed vector                                      *
*                                                                      *
************************************************************************
C
      SUBROUTINE IA213(DA,KDIA,KDIAS,NDIA,DX,DB,DD,NEQ,NIT,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),KDIA(*),KDIAS(*),DX(*),DB(*),DD(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('IA213 ','11/19/90')
C
      DO 1 IEQ=1,NEQ
1     DA(IEQ)=1D0/DA(IEQ)
C
      DO 2 ITE=1,NIT
      CALL LCP1(DB,DD,NEQ)
      DO 3 IDIA=2,NDIA
      J1=KDIA(IDIA)
      IF (J1.GT.0) THEN
       I1=1
       NEQ1=NEQ-J1
      ELSE
       I1=1-J1
       NEQ1=NEQ+J1
      ENDIF
      J0=KDIAS(IDIA)-I1
      CALL LVM1(DX(I1+J1),DA(I1+J0),DD(I1),NEQ1,-1D0,1D0)
3     CONTINUE
      CALL LVM1(DD,DA,DD,NEQ,1D0,0D0)
      CALL LLC1(DD,DX,NEQ,1D0-OMEGA,OMEGA)
2     CONTINUE
C
      DO 10 IEQ=1,NEQ
10    DA(IEQ)=1D0/DA(IEQ)
C
      END
C
C
C
      SUBROUTINE IA217(DA,KCOL,KLD,DX,DB,DD,NEQ,NIT,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),KCOL(*),KLD(*),DX(*),DB(*),DD(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('IA217 ','10/26/89')
C
      DO 1 ITE=1,NIT
      CALL LCP1(DB,DD,NEQ)
      DO 2 IEQ=1,NEQ
      DO 3 ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
3     DD(IEQ)=DD(IEQ)-DA(ICOL)*DX(KCOL(ICOL))
2     CONTINUE
      DO 4 IEQ=1,NEQ
4     DX(IEQ)=(1D0-OMEGA)*DX(IEQ)+OMEGA*DD(IEQ)/DA(KLD(IEQ))
1     CONTINUE
C
      END
C
C
C
      SUBROUTINE IA21A(DA,KCOL,KLD,KOP,DX,DB,DD,NEQ,NIT,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),KCOL(*),KLD(*),KOP(*),DX(*),DB(*),DD(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('IA21A ','10/26/89')
C
      DO 1 ITE=1,NIT
      CALL LCP1(DB,DD,NEQ)
      DO 2 IEQ=1,NEQ
      JOP=KOP(IEQ)
      DO 3 ICOL=KLD(JOP)+1,KLD(JOP+1)-1
3     DD(IEQ)=DD(IEQ)-DA(ICOL)*DX(KCOL(ICOL)+IEQ)
2     CONTINUE
      DO 4 IEQ=1,NEQ
4     DX(IEQ)=(1D0-OMEGA)*DX(IEQ)+OMEGA*DD(IEQ)/DA(KLD(KOP(IEQ)))
1     CONTINUE
C
      END
