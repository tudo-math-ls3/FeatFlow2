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
* ID12n                                                                *
*                                                                      *
* Purpose  SSOR Preconditioning                                        *
*          Matrix stored in technique  n  (see Reference Manual)       *
*          Single/Single precision version                             *
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
* VX       R*4    Vector                                               *
* NEQ      I*4    Number of equations (length of VX)                   *
* OMEGA    R*8    Relaxation parameter                                 *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* VX       R*4    Resulting vector                                     *
*                                                                      *
************************************************************************
C
      SUBROUTINE ID127(VA,KCOL,KLD,VX,NEQ,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),VX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('ID127 ','01/02/89')
C
      DO 1 IEQ=1,NEQ
      AUX=0D0
      ILD=KLD(IEQ)
      DO 4 ICOL=ILD+1,KLD(IEQ+1)-1
      IF (KCOL(ICOL).GE.IEQ) GOTO 1
4     AUX=AUX+VA(ICOL)*VX(KCOL(ICOL))
1     VX(IEQ)=(VX(IEQ)-AUX*OMEGA)/VA(ILD)
C
      DO 10 IEQ=NEQ-1,1,-1
      AUX=0D0
      ILD=KLD(IEQ)
      DO 12 ICOL=ILD+1,KLD(IEQ+1)-1
      IF (KCOL(ICOL).LE.IEQ) GOTO 12
      AUX=AUX+VA(ICOL)*VX(KCOL(ICOL))
12    CONTINUE
10    VX(IEQ)=VX(IEQ)-AUX*OMEGA/VA(ILD)
C
      END
C
C
C
      SUBROUTINE ID128(VA,KCOL,KLD,VX,NEQ,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),VX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('ID128 ','01/02/89')
C
      DO 1 IEQ=1,NEQ
      ILD=KLD(IEQ)
      VX(IEQ)=VX(IEQ)*OMEGA/VA(ILD)
      DO 4 ICOL=ILD+1,KLD(IEQ+1)-1
4     VX(KCOL(ICOL))=VX(KCOL(ICOL))-VA(ICOL)*VX(IEQ)
1     CONTINUE
C
      DO 10 IEQ=NEQ-1,1,-1
      AUX=0D0
      ILD=KLD(IEQ)
      DO 12 ICOL=ILD+1,KLD(IEQ+1)-1
12    AUX=AUX+VA(ICOL)*VX(KCOL(ICOL))
10    VX(IEQ)=VX(IEQ)-AUX*OMEGA/VA(ILD)
C
      END
C
C
C
      SUBROUTINE ID12A(VA,KCOL,KLD,KOP,VX,NEQ,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),VX(*),KOP(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('ID12A ','01/02/89')
C
      DO 1 IEQ=1,NEQ
      AUX=0D0
      IOP=KOP(IEQ)
      ILD=KLD(IOP)
      DO 4 ICOL=ILD+1,KLD(IOP+1)-1
      IF (KCOL(ICOL).GT.0) GOTO 1
4     AUX=AUX+VA(ICOL)*VX(KCOL(ICOL)+IEQ)
1     VX(IEQ)=(VX(IEQ)-AUX*OMEGA)/VA(ILD)
C
      DO 10 IEQ=NEQ-1,1,-1
      AUX=0D0
      IOP=KOP(IEQ)
      ILD=KLD(IOP)
      DO 12 ICOL=KLD(IOP+1)-1,ILD+1,-1
      IF (KCOL(ICOL).LE.0) GOTO 10
      AUX=AUX+VA(ICOL)*VX(KCOL(ICOL)+IEQ)
12    CONTINUE
10    VX(IEQ)=VX(IEQ)-AUX*OMEGA/VA(ILD)
C
      END
