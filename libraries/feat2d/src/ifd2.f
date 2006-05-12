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
* IFD2n                                                                *
*                                                                      *
* Purpose  Calculate ILU decomposition of a matrix VA                  *
*          Matrix stored in technique  n  (see Reference Manual)       *
*          Single precision version                                    *
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
* NEQ      I*4    Number of equations                                  *
* ILU      I*4                                                         *
* ALPHA    R*8                                                         *
* TOL      R*8                                                         *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* VA       R*8    Resulting matrix                                     *
*                                                                      *
************************************************************************
C
      SUBROUTINE IFD27(VA,KCOL,KLD,NEQ,ILU,ALPHA,TOL)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VA(*),KCOL(*),KLD(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /ERRCTL/,/CHAR/
C
      SUB='IFD27'
      IF (ICHECK.GE.997) CALL OTRC('IFD27 ','12/02/89')
C
      IF (ILU.EQ.1) THEN
C
C *** Constant shift for ILU
       ALPHA1=1D0/(1D0+ALPHA)
       DO 11 IEQ=1,NEQ
       DO 11 ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
11     VA(ICOL)=VA(ICOL)*ALPHA1
C
      ELSE
C
C *** Constant shift for MILU
       ALPHA1=1D0+ALPHA
       DO 21 IEQ=1,NEQ
       ILD=KLD(IEQ)
21     VA(ILD)=VA(ILD)*ALPHA1
C
      ENDIF
C
      BMILU=ILU.EQ.2
      ILD=1
C
C *** Incomplete GAUSS elimination
      DO 110 IEQ=2,NEQ
      A=VA(ILD)
C
      IF (ABS(A).LT.TOL) THEN
       VA(ILD)=SIGN(TOL,A)
       IF (IER.EQ.0) IER=IEQ-1
      ENDIF
C
      ILD=KLD(IEQ)
      IF (KCOL(ILD+1).GE.IEQ) GOTO 110
C
      DO 120 JCOL=ILD+1,KLD(IEQ+1)-1
      ICOL=KCOL(JCOL)
      IF (ICOL.GE.IEQ) GOTO 110
      JLD=KLD(ICOL)
      AIJ=VA(JCOL)/VA(JLD)
C *** AIJ is the entry of the elimination matrix at (IEQ,ICOL)
      VA(JCOL)=AIJ
C
      IP=KLD(IEQ+1)-1
      IF (BMILU) THEN
C
C *** Modified incomplete LU decomposition
C
C *** Loop over subdiagonal entries in line ICOL
      DO 130 JC=KLD(ICOL+1)-1,JLD,-1
      JCOL0=KCOL(JC)
C *** First check for diagonal entry - diagonals are stored sep4rately
      IF (JCOL0-IEQ) 131,135,132
131   IF (JCOL0.LE.ICOL) GOTO 120
C *** Look for an entry at position (ICOL,JCOL0)
132   IF (KCOL(IP)-JCOL0) 135,133,134
C *** Off diagonal entry (ICOL,JCOL0) found
133   VA(IP)=VA(IP)-AIJ*VA(JC)
      IP=IP-1
      GOTO 130
C *** Entry not yet found
134   IP=IP-1
C *** Continue
      IF (IP.GT.ILD) GOTO 132
C *** Insert at diagonal position
135   VA(ILD)=VA(ILD)-AIJ*VA(JC)
130   CONTINUE
C
      ELSE
C
C *** Shifted incomplete LU decomposition
C
C *** Loop over subdiagonal entries in line ICOL
      DO 140 JC=KLD(ICOL+1)-1,JLD,-1
      JCOL0=KCOL(JC)
C *** First check for diagonal entry - diagonals are stored separately
      IF (JCOL0-IEQ) 141,145,142
141   IF (JCOL0.LE.ICOL) GOTO 120
C *** Look for an entry at position (ICOL,JCOL0)
142   IF (KCOL(IP)-JCOL0) 140,143,144
C *** Off diagonal entry (ICOL,JCOL0) found
143   VA(IP)=VA(IP)-AIJ*VA(JC)
      IP=IP-1
      GOTO 140
C *** Entry not yet found - continue
144   IP=IP-1
      GOTO 142
C *** Insert at diagonal position
145   VA(ILD)=VA(ILD)-AIJ*VA(JC)
140   CONTINUE
C
      ENDIF
C
120   CONTINUE
110   CONTINUE
C
      END
