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
* IE23n                                                                *
*                                                                      *
* Purpose  Smoothing using a preconditioned conjugate gradient method  *
*          Single/double precision version                             *
*                                                                      *
* Subroutines/functions called  LSP1 , LLC1 , LL21 , LAX3n , ID13n     *
*                                                                      *
* Version from  10/30/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* VA       R*4    Matrix in storage technique n                        *
* KCOL     I*4                                                         *
* KLD      I*4                                                         *
* KOP      I*4                                                         *
* DX       R*8    Starting vector                                      *
* DB       R*8    Right hand side                                      *
* NEQ      I*4    Number of smoothing steps                            *
* OMEGA    R*8    0 <  OMEGA     No Preconditioning                    *
*                 0 >= OMEGA  -  SSOR-Preconditioning                  *
* DR,DD    R*8    Workspace vectors of length NEQ                      *
* DD1,DG   R*8    For OMEGA < 0 , DG must be replaced by DR            *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8    Smoothed vector                                      *
*                                                                      *
************************************************************************
C
      SUBROUTINE IE237(VA,KCOL,KLD,DX,DB,NEQ,NIT,OMEGA,DR,DD,DD1,DG)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VA(*),KCOL(*),KLD(*)
      DIMENSION DX(*),DB(*),DR(*),DG(*),DD(*),DD1(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /ERRCTL/,/CHAR/
C
      SUB='IE237'
      IF (ICHECK.GE.998) CALL OTRC('IE237 ','10/30/89')
C
      BNOCON=OMEGA.LT.0D0
C
C *** Initialization
      CALL LAX37(VA,KCOL,KLD,NEQ,DX,DR,1D0,0D0)
      CALL LLC1(DB,DR,NEQ,-1D0,1D0)
      CALL LL21(DR,NEQ,RES)
C
      IF (BNOCON) THEN
       SIGMA0=RES*RES
      ELSE
       CALL LCP1(DR,DG,NEQ)
       CALL ID137(VA,KCOL,KLD,DG,NEQ,OMEGA)
       CALL LSP1(DR,DG,NEQ,SIGMA0)
      ENDIF
C
      CALL LLC1(DG,DD,NEQ,-1D0,0D0)
C
C *** Iterative correction
      DO 1 ITE=1,NIT
C
      CALL LAX37(VA,KCOL,KLD,NEQ,DD,DD1,1D0,0D0)
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
       CALL ID137(VA,KCOL,KLD,DG,NEQ,OMEGA)
       CALL LSP1(DR,DG,NEQ,SIGMA1)
      ENDIF
C
      GAMMA=SIGMA1/SIGMA0
      SIGMA0=SIGMA1
      CALL LLC1(DG,DD,NEQ,-1D0,GAMMA)
1     CONTINUE
C
      END
C
C
C
      SUBROUTINE IE238(VA,KCOL,KLD,DX,DB,NEQ,NIT,OMEGA,DR,DD,DD1,DG)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VA(*),KCOL(*),KLD(*)
      DIMENSION DX(*),DB(*),DR(*),DG(*),DD(*),DD1(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /ERRCTL/,/CHAR/
C
      SUB='IE238'
      IF (ICHECK.GE.998) CALL OTRC('IE238 ','10/30/89')
C
      BNOCON=OMEGA.LT.0D0
C
C *** Initialization
      CALL LAX38(VA,KCOL,KLD,NEQ,DX,DR,1D0,0D0)
      CALL LLC1(DB,DR,NEQ,-1D0,1D0)
      CALL LL21(DR,NEQ,RES)
C
      IF (BNOCON) THEN
       SIGMA0=RES*RES
      ELSE
       CALL LCP1(DR,DG,NEQ)
       CALL ID138(VA,KCOL,KLD,DG,NEQ,OMEGA)
       CALL LSP1(DR,DG,NEQ,SIGMA0)
      ENDIF
C
      CALL LLC1(DG,DD,NEQ,-1D0,0D0)
C
C *** Iterative correction
      DO 1 ITE=1,NIT
C
      CALL LAX38(VA,KCOL,KLD,NEQ,DD,DD1,1D0,0D0)
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
       CALL ID138(VA,KCOL,KLD,DG,NEQ,OMEGA)
       CALL LSP1(DR,DG,NEQ,SIGMA1)
      ENDIF
C
      GAMMA=SIGMA1/SIGMA0
      SIGMA0=SIGMA1
      CALL LLC1(DG,DD,NEQ,-1D0,GAMMA)
1     CONTINUE
C
      END
C
C
C
      SUBROUTINE IE23A(VA,KCOL,KLD,KOP,DX,DB,NEQ,NIT,OMEGA,DR,DD,DD1,DG)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VA(*),KCOL(*),KLD(*),KOP(*)
      DIMENSION DX(*),DB(*),DR(*),DG(*),DD(*),DD1(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /ERRCTL/,/CHAR/
C
      SUB='IE23A'
      IF (ICHECK.GE.998) CALL OTRC('IE23A ','10/30/89')
C
      BNOCON=OMEGA.LT.0D0
C
C *** Initialization
      CALL LAX3A(VA,KCOL,KLD,KOP,NEQ,DX,DR,1D0,0D0)
      CALL LLC1(DB,DR,NEQ,-1D0,1D0)
      CALL LL21(DR,NEQ,RES)
C
      IF (BNOCON) THEN
       SIGMA0=RES*RES
      ELSE
       CALL LCP1(DR,DG,NEQ)
       CALL ID13A(VA,KCOL,KLD,KOP,DG,NEQ,OMEGA)
       CALL LSP1(DR,DG,NEQ,SIGMA0)
      ENDIF
C
      CALL LLC1(DG,DD,NEQ,-1D0,0D0)
C
C *** Iterative correction
      DO 1 ITE=1,NIT
C
      CALL LAX3A(VA,KCOL,KLD,KOP,NEQ,DD,DD1,1D0,0D0)
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
       CALL ID13A(VA,KCOL,KLD,KOP,DG,NEQ,OMEGA)
       CALL LSP1(DR,DG,NEQ,SIGMA1)
      ENDIF
C
      GAMMA=SIGMA1/SIGMA0
      SIGMA0=SIGMA1
      CALL LLC1(DG,DD,NEQ,-1D0,GAMMA)
1     CONTINUE
C
      END
