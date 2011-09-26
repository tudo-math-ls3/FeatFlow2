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
* IF13n                                                                *
*                                                                      *
* Purpose  Preconditioning by ILU decomposition                        *
*          Matrix stored in technique  n  (see Reference Manual)       *
*          Single/double precision version                             *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  12/02/89                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* VA       R*4    Matrix stored in technique  n                        *
* KCOL     I*4    Pointer vectors corresponding to the                 *
* KLD      I*4    storage technique                                    *
* KOP      I*4                                                         *
* DX       R*8    Vector                                               *
* NEQ      I*4    Number of equations (length of DX)                   *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8    Resulting vector                                     *
*                                                                      *
************************************************************************
C
      SUBROUTINE IF137(VA,KCOL,KLD,DX,NEQ)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),DX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('IF137 ','12/02/89')
C
      DO 501 IEQ=2,NEQ
      AUX=0D0
      DO 502 JCOL=KLD(IEQ)+1,KLD(IEQ+1)-1
      ICOL=KCOL(JCOL)
      IF (ICOL.GE.IEQ) GOTO 501
502   AUX=AUX+DBLE(VA(JCOL))*DX(ICOL)
501   DX(IEQ)=DX(IEQ)-AUX
C
      DX(NEQ)=DX(NEQ)/DBLE(VA(KLD(NEQ)))
      DO 503 IEQ=NEQ-1,1,-1
      ILD=KLD(IEQ)
      AUX=0D0
      DO 504 JCOL=KLD(IEQ+1)-1,ILD+1,-1
      ICOL=KCOL(JCOL)
      IF (ICOL.LE.IEQ) GOTO 503
      AUX=AUX+DBLE(VA(JCOL))*DX(ICOL)
504   CONTINUE
503   DX(IEQ)=(DX(IEQ)-AUX)/DBLE(VA(ILD))
C
      END
