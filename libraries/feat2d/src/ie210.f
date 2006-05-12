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
* IE210                                                                *
*                                                                      *
* Purpose  Smoothing using a preconditioned conjugate gradient method  *
*          Double precision version                                    *
*                                                                      *
* Subroutines/functions called  LSP1 , LLC1 , LL21                     *
*                                                                      *
* Version from  10/30/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DX       R*8    Starting vector                                      *
* DB       R*8    Right hand side                                      *
* NEQ      I*4    Number of equations                                  *
* NIT      I*4    Number of smoothing steps                            *
* DAX0     SUBR   EXTERNAL Subroutine DAX0(DX,DAX,NEQ,A1,A2)           *
*                 Results  DAX := A1 * A * DX + A2 * DAX               *
* DCG0C    SUBR   EXTERNAL Subroutine DCG0C(DG,NEQ)                    *
*                 Results  DG := C**(-1) * DG  for the precondioning   *
*                 matrix C                                             *
* BNOCON   L*4    .TRUE.  No Preconditioning                           *
* DR,DD    R*8    Workspace vectors of length NEQ                      *
* DD1,DG   R*8    For BNOCON , DG must be replaced by DR               *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8    Smoothed vector                                      *
*                                                                      *
************************************************************************
C
      SUBROUTINE IE210(DX,DB,NEQ,NIT,DAX0,DCG0C,BNOCON,DR,DD,DD1,DG)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION DX(*),DB(*),DR(*),DG(*),DD(*),DD1(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /ERRCTL/,/CHAR/
C
      SUB='IE210'
      IF (ICHECK.GE.998) CALL OTRC('IE210 ','10/30/89')
C
C *** Initialization
      CALL DAX0(DX,DR,NEQ,1D0,0D0)
      CALL LLC1(DB,DR,NEQ,-1D0,1D0)
      CALL LL21(DR,NEQ,RES)
C
      IF (BNOCON) THEN
       SIGMA0=RES*RES
      ELSE
       CALL LCP1(DR,DG,NEQ)
       CALL DCG0C(DG,NEQ)
       CALL LSP1(DR,DG,NEQ,SIGMA0)
      ENDIF
C
      CALL LLC1(DG,DD,NEQ,-1D0,0D0)
C
C *** Iterative correction
      DO 1 ITE=1,NIT
C
      CALL DAX0(DD,DD1,NEQ,1D0,0D0)
      CALL LSP1(DD,DD1,NEQ,ALPHA)
      ALPHA=SIGMA0/ALPHA
      CALL LLC1(DD,DX,NEQ,ALPHA,1D0)
      CALL LLC1(DD1,DR,NEQ,ALPHA,1D0)
C
      CALL LL21(DR,NEQ,FR)
C
      IF (BNOCON) THEN
       SIGMA1=FR*FR
      ELSE
       CALL LCP1(DR,DG,NEQ)
       CALL DCG0C(DG,NEQ)
       CALL LSP1(DR,DG,NEQ,SIGMA1)
      ENDIF
C
      GAMMA=SIGMA1/SIGMA0
      SIGMA0=SIGMA1
      CALL LLC1(DG,DD,NEQ,-1D0,GAMMA)
C
1     CONTINUE
C
      END
