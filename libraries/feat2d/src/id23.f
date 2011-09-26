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
* ID23n                                                                *
*                                                                      *
* Purpose  Smoothing using SSOR Iteration                              *
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
* KCOL     I*4    storage technique                                    *
* KOP      I*4                                                         *
* DX       R*8    Vector                                               *
* DB       R*8    Right hand side                                      *
* NEQ      I*4    Number of equations (length of DX)                   *
* NIT      I*4    Number of smoothing steps                            *
* OMEGA    R*8    Relaxation parameter                                 *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8    Smoothed vector                                      *
*                                                                      *
************************************************************************
C
      SUBROUTINE ID237(VA,KCOL,KLD,DX,DB,NEQ,NIT,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),DX(*),DB(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('ID237 ','01/02/89')
C
      DO 5 ITE=1,NIT
      DO 1 IEQ=1,NEQ
      AUX=0D0
      ILD=KLD(IEQ)
      DO 4 ICOL=ILD+1,KLD(IEQ+1)-1
      IF (KCOL(ICOL).GE.IEQ) GOTO 1
4     AUX=AUX+VA(ICOL)*DX(KCOL(ICOL))
1     DX(IEQ)=(DX(IEQ)-AUX*OMEGA)/VA(ILD)
C
      DO 10 IEQ=NEQ-1,1,-1
      AUX=0D0
      ILD=KLD(IEQ)
      DO 12 ICOL=ILD+1,KLD(IEQ+1)-1
      IF (KCOL(ICOL).LE.IEQ) GOTO 12
      AUX=AUX+VA(ICOL)*DX(KCOL(ICOL))
12    CONTINUE
10    DX(IEQ)=DX(IEQ)-AUX*OMEGA/VA(ILD)
5     CONTINUE
C
      END
C
C
C
      SUBROUTINE ID238(VA,KCOL,KLD,DX,DB,NEQ,NIT,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),DX(*),DB(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('ID238 ','01/02/89')
C
      DO 5 ITE=1,NIT
      DO 1 IEQ=1,NEQ
      ILD=KLD(IEQ)
      DX(IEQ)=DX(IEQ)*OMEGA/VA(ILD)
      DO 4 ICOL=ILD+1,KLD(IEQ+1)-1
4     DX(KCOL(ICOL))=DX(KCOL(ICOL))-VA(ICOL)*DX(IEQ)
1     CONTINUE
C
      DO 10 IEQ=NEQ-1,1,-1
      AUX=0D0
      ILD=KLD(IEQ)
      DO 12 ICOL=ILD+1,KLD(IEQ+1)-1
12    AUX=AUX+VA(ICOL)*DX(KCOL(ICOL))
10    DX(IEQ)=DX(IEQ)-AUX*OMEGA/VA(ILD)
5     CONTINUE
C
      END
C
C
C
      SUBROUTINE ID23A(VA,KCOL,KLD,KOP,DX,DB,NEQ,NIT,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),DX(*),DB(*),KOP(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('ID23A ','01/02/89')
C
      DO 5 ITE=1,NIT
      DO 1 IEQ=1,NEQ
      AUX=0D0
      IOP=KOP(IEQ)
      ILD=KLD(IOP)
      DO 4 ICOL=ILD+1,KLD(IOP+1)-1
      IF (KCOL(ICOL).GT.0) GOTO 1
4     AUX=AUX+VA(ICOL)*DX(KCOL(ICOL)+IEQ)
1     DX(IEQ)=(DX(IEQ)-AUX*OMEGA)/VA(ILD)
C
      DO 10 IEQ=NEQ-1,1,-1
      AUX=0D0
      IOP=KOP(IEQ)
      ILD=KLD(IOP)
      DO 12 ICOL=KLD(IOP+1)-1,ILD+1,-1
      IF (KCOL(ICOL).LE.0) GOTO 10
      AUX=AUX+VA(ICOL)*DX(KCOL(ICOL)+IEQ)
12    CONTINUE
10    DX(IEQ)=DX(IEQ)-AUX*OMEGA/VA(ILD)
5     CONTINUE
C
      END
