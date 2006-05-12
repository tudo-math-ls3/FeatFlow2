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
* IE010                                                                *
*                                                                      *
* Purpose  Solution of a linear system  A*X = B  using                 *
*          a preconditioned conjugate gradient method                  *
*          Double precision version                                    *
*                                                                      *
* Subroutines/functions called  LSP1 , LLC1 , LL21 , LLI1, LCL1        *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DX       R*8    Starting vector                                      *
* DB       R*8    Right hand side                                      *
* NEQ      I*4    Number of equations                                  *
* ITE      I*4    Minimum number of iterations (used for IREL=1 only)  *
* NIT      I*4    Maximum number of iterations                         *
* EPS      R*8    Desired precision                                    *
*                 IREL=0: Stop if !!RES!! < EPS                        *
*                 IREL=1: Stop if !!RES!!/!!RES0!! < EPS               *
*                         and a minimum of ITE iterations is performed *
* DAX0     SUBR   EXTERNAL Subroutine DAX0(DX,DAX,NEQ,A1,A2)           *
*                 Results  DAX := A1 * A * DX + A2 * DAX               *
* DCG0C    SUBR   EXTERNAL Subroutine DCG0C(DG,NEQ)                    *
*                 Results  DG := C**(-1) * DG  for the precondioning   *
*                 matrix C                                             *
* BNOCON   L*4    .TRUE.   No Preconditioning                          *
* DR,DD    R*8    Workspace vectors of length NEQ                      *
* DD1,DG   R*8    For BNOCON , DG must be replaced by DR               *
* IREL     I*4    See above                                            *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8    Solution vector                                      *
* ITE      I*4    Number of iterations                                 *
* IER      I*4    Error indicator                                      *
*                 +1  Precision EPS not achieved after NIT iterations  *
* RHO      R*8    Convergence rate                                     *
*                                                                      *
* Changes:                                                             *
* --------                                                             *
* 15.10.2004 - Parameter RHO added                                     *
*                                                                      *
************************************************************************
C
      SUBROUTINE IE010(DX,DB,NEQ,NIT,ITE,EPS,DAX0,DCG0C,BNOCON,
     *                 DR,DD,DD1,DG,IREL,RHO)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION DX(*),DB(*),DR(*),DG(*),DD(*),DD1(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
      DATA FR/0D0/
C
      SUB='IE010'
      IF (ICHECK.GE.997) CALL OTRC('IE010 ','01/02/89')
C
      BMSG2=M.GE.2.OR.MT.GE.2
      NIT0=MAX(ITE,0)
      BREL=IREL.EQ.1
C
      IF (ICHECK.GT.0) THEN
       CALL LLI1(DB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL1(DX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IE010 ')
        GOTO 99999
       ENDIF
      ENDIF
C
C *** Initialization
      CALL DAX0(DX,DR,NEQ,1D0,0D0)
      CALL LLC1(DB,DR,NEQ,-1D0,1D0)
      CALL LL21(DR,NEQ,RES)
      IF (RES.LE.EPS.AND..NOT.BREL) THEN
       ITE=0
       FR=RES
       GOTO 200
      ENDIF
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
      DO 100 ITE=1,NIT
C
      CALL DAX0(DD,DD1,NEQ,1D0,0D0)
      CALL LSP1(DD,DD1,NEQ,ALPHA)
      ALPHA=SIGMA0/ALPHA
      CALL LLC1(DD,DX,NEQ,ALPHA,1D0)
      CALL LLC1(DD1,DR,NEQ,ALPHA,1D0)
C
      CALL LL21(DR,NEQ,FR)
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,FR
       CALL OMSG(73,'IE010 ')
      ENDIF
      IF (BREL) THEN
       IF (FR.LE.RES*EPS.AND.ITE.GE.NIT0) GOTO 200
      ELSE
       IF (FR.LE.EPS) GOTO 200
      ENDIF
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
100   CONTINUE
C
      WRITE (CPARAM,'(I15,2D25.16)') NIT,FR,RES
      CALL OMSG(71,'IE010 ')
      CALL OMSG(72,'IE010 ')
C
      IF (RES.GE.1D-70) THEN
       CAPPA=(FR/RES)**(1D0/NIT)
      ELSE
       CAPPA=0D0
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE010 ')
C
      IER=1
      GOTO 99999
C
200   IER=0
      IF (RES.GE.1D-70) RES=FR/RES
      WRITE (CPARAM,'(I15,2D25.16)') ITE,FR,RES
      CALL OMSG(72,'IE010 ')
C
      IF (ITE.EQ.0) THEN
       CAPPA=0D0
      ELSE
       CAPPA=RES**(1D0/ITE)
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE010 ')
C

C Save the convergence rate CAPPA to RHO and exit:
99999 RHO=CAPPA
      END
