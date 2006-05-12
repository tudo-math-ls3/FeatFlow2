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
* IE33n                                                                *
*                                                                      *
* Purpose  Solution of a linear system  A*X = B  using                 *
*          a preconditioned conjugate gradient method                  *
*          Single/double precision version                             *
*                                                                      *
* Subroutines/functions called  LSP1 , LLC1 , LL21 , LLI1 , LCL1 ,     *
*                               LAX3n , ID13n                          *
*                                                                      *
* Version from  11/12/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* VA       R*4    Matrix in storage technique n                        *
* KCOL     I*4                                                         *
* KLD      I*4                                                         *
* KOP      I*4                                                         *
* DX       R*8    Starting vector                                      *
* DB       R*8    Right hand side                                      *
* NEQ      I*4    Number of equations                                  *
* ITE      I*4    Minimum number of iterations                         *
* NIT      I*4    Maximum number of iterations                         *
* EPS      R*8    Desired precision, stop if !!RES!!/!!RES0!! < EPS    *
* OMEGA    R*8    0 <  OMEGA     No Preconditioning                    *
*                 0 >= OMEGA  -  SSOR Preconditioning                  *
* DR,DD    R*8    Workspace vectors of length NEQ                      *
* DD1,DG   R*8    For OMEGA < 0 , DG must be replaced by DR            *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8    Solution vector                                      *
* ITE      I*4    Number of iterations                                 *
* IER      I*4    Error indicator                                      *
*                 +1  Precision EPS not achieved after NIT iterations  *
*                                                                      *
************************************************************************
C
      SUBROUTINE IE333(VA,KDIA,KDIAS,NDIA,DX,DB,NEQ,NIT,ITE,EPS,OMEGA,
     *                 DR,DD,DD1,DG)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VA(*),KDIA(*),KDIAS(*)
      DIMENSION DX(*),DB(*),DR(*),DG(*),DD(*),DD1(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
      DATA FR/0D0/
C
      SUB='IE333'
      IF (ICHECK.GE.997) CALL OTRC('IE333 ','02/01/91')
C
      BMSG2=M.GE.2.OR.MT.GE.2
      NIT0=MAX(ITE,0)
C
      BNOCON=OMEGA.LT.0D0
C
      IF (ICHECK.GT.0) THEN
       CALL LLI1(DB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL1(DX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IE333 ')
        GOTO 99999
       ENDIF
      ENDIF
C
C *** Initialization
      CALL LAX33(VA,KDIA,KDIAS,NDIA,NEQ,DX,DR,1D0,0D0)
      CALL LLC1(DB,DR,NEQ,-1D0,1D0)
      CALL LL21(DR,NEQ,RES)
      IF (RES.LE.EPS) THEN
       ITE=0
       FR=RES
       GOTO 200
      ENDIF
C
      IF (BNOCON) THEN
       SIGMA0=RES*RES
      ELSE
       CALL LCP1(DR,DG,NEQ)
       CALL IA133(VA,DG,NEQ)
       CALL LSP1(DR,DG,NEQ,SIGMA0)
      ENDIF
C
      CALL LLC1(DG,DD,NEQ,-1D0,0D0)
C
C *** Iterative correction
      DO 100 ITE=1,NIT
C
      CALL LAX33(VA,KDIA,KDIAS,NDIA,NEQ,DD,DD1,1D0,0D0)
      CALL LSP1(DD,DD1,NEQ,ALPHA)
      ALPHA=SIGMA0/ALPHA
      CALL LLC1(DD,DX,NEQ,ALPHA,1D0)
      CALL LLC1(DD1,DR,NEQ,ALPHA,1D0)
C
      CALL LL21(DR,NEQ,FR)
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,FR
       CALL OMSG(73,'IE333 ')
      ENDIF
      IF (FR.LE.RES*EPS.AND.ITE.GE.NIT0) GOTO 200
C
      IF (BNOCON) THEN
       SIGMA1=FR*FR
      ELSE
       CALL LCP1(DR,DG,NEQ)
       CALL IA133(VA,DG,NEQ)
       CALL LSP1(DR,DG,NEQ,SIGMA1)
      ENDIF
C
      GAMMA=SIGMA1/SIGMA0
      SIGMA0=SIGMA1
      CALL LLC1(DG,DD,NEQ,-1D0,GAMMA)
100   CONTINUE
C
      WRITE (CPARAM,'(I15,2D25.16)') NIT,FR,RES
      CALL OMSG(71,'IE333 ')
      CALL OMSG(72,'IE333 ')
C
      IF (RES.GE.1D-70) THEN
       CAPPA=(FR/RES)**(1D0/NIT)
      ELSE
       CAPPA=0D0
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE333 ')
C
      IER=1
      GOTO 99999
C
200   IER=0
      IF (RES.GE.1D-70) RES=FR/RES
      WRITE (CPARAM,'(I15,2D25.16)') ITE,FR,RES
      CALL OMSG(72,'IE333 ')
C
      IF (ITE.EQ.0) THEN
       CAPPA=0D0
      ELSE
       CAPPA=RES**(1D0/ITE)
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE333 ')
C
99999 END
C
C
C
      SUBROUTINE IE334(VA,KDIA,KDIAS,NDIA,DX,DB,NEQ,NIT,ITE,EPS,OMEGA,
     *                 DR,DD,DD1,DG)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VA(*),KDIA(*),KDIAS(*)
      DIMENSION DX(*),DB(*),DR(*),DG(*),DD(*),DD1(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
      DATA FR/0D0/
C
      SUB='IE334'
      IF (ICHECK.GE.997) CALL OTRC('IE334 ','02/01/91')
C
      BMSG2=M.GE.2.OR.MT.GE.2
      NIT0=MAX(ITE,0)
C
      BNOCON=OMEGA.LT.0D0
C
      IF (ICHECK.GT.0) THEN
       CALL LLI1(DB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL1(DX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IE334 ')
        GOTO 99999
       ENDIF
      ENDIF
C
C *** Initialization
      CALL LAX34(VA,KDIA,KDIAS,NDIA,NEQ,DX,DR,1D0,0D0)
      CALL LLC1(DB,DR,NEQ,-1D0,1D0)
      CALL LL21(DR,NEQ,RES)
      IF (RES.LE.EPS) THEN
       ITE=0
       FR=RES
       GOTO 200
      ENDIF
C
      IF (BNOCON) THEN
       SIGMA0=RES*RES
      ELSE
       CALL LCP1(DR,DG,NEQ)
       CALL IA133(VA,DG,NEQ)
       CALL LSP1(DR,DG,NEQ,SIGMA0)
      ENDIF
C
      CALL LLC1(DG,DD,NEQ,-1D0,0D0)
C
C *** Iterative correction
      DO 100 ITE=1,NIT
C
      CALL LAX34(VA,KDIA,KDIAS,NDIA,NEQ,DD,DD1,1D0,0D0)
      CALL LSP1(DD,DD1,NEQ,ALPHA)
      ALPHA=SIGMA0/ALPHA
      CALL LLC1(DD,DX,NEQ,ALPHA,1D0)
      CALL LLC1(DD1,DR,NEQ,ALPHA,1D0)
C
      CALL LL21(DR,NEQ,FR)
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,FR
       CALL OMSG(73,'IE334 ')
      ENDIF
      IF (FR.LE.RES*EPS.AND.ITE.GE.NIT0) GOTO 200
C
      IF (BNOCON) THEN
       SIGMA1=FR*FR
      ELSE
       CALL LCP1(DR,DG,NEQ)
       CALL IA133(VA,DG,NEQ)
       CALL LSP1(DR,DG,NEQ,SIGMA1)
      ENDIF
C
      GAMMA=SIGMA1/SIGMA0
      SIGMA0=SIGMA1
      CALL LLC1(DG,DD,NEQ,-1D0,GAMMA)
100   CONTINUE
C
      WRITE (CPARAM,'(I15,2D25.16)') NIT,FR,RES
      CALL OMSG(71,'IE334 ')
      CALL OMSG(72,'IE334 ')
C
      IF (RES.GE.1D-70) THEN
       CAPPA=(FR/RES)**(1D0/NIT)
      ELSE
       CAPPA=0D0
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE334 ')
C
      IER=1
      GOTO 99999
C
200   IER=0
      IF (RES.GE.1D-70) RES=FR/RES
      WRITE (CPARAM,'(I15,2D25.16)') ITE,FR,RES
      CALL OMSG(72,'IE334 ')
C
      IF (ITE.EQ.0) THEN
       CAPPA=0D0
      ELSE
       CAPPA=RES**(1D0/ITE)
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE334 ')
C
99999 END
C
C
C
      SUBROUTINE IE337(VA,KCOL,KLD,DX,DB,NEQ,NIT,ITE,EPS,OMEGA,
     *                 DR,DD,DD1,DG)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VA(*),KCOL(*),KLD(*)
      DIMENSION DX(*),DB(*),DR(*),DG(*),DD(*),DD1(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
      DATA FR/0D0/
C
      SUB='IE337'
      IF (ICHECK.GE.997) CALL OTRC('IE337 ','11/12/89')
C
      BMSG2=M.GE.2.OR.MT.GE.2
      NIT0=MAX(ITE,0)
C
      BNOCON=OMEGA.LT.0D0
C
      IF (ICHECK.GT.0) THEN
       CALL LLI1(DB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL1(DX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IE337 ')
        GOTO 99999
       ENDIF
      ENDIF
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
      DO 100 ITE=1,NIT
C
      CALL LAX37(VA,KCOL,KLD,NEQ,DD,DD1,1D0,0D0)
      CALL LSP1(DD,DD1,NEQ,ALPHA)
      ALPHA=SIGMA0/ALPHA
      CALL LLC1(DD,DX,NEQ,ALPHA,1D0)
      CALL LLC1(DD1,DR,NEQ,ALPHA,1D0)
C
      CALL LL21(DR,NEQ,FR)
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,FR
       CALL OMSG(73,'IE337 ')
      ENDIF
      IF (FR.LE.RES*EPS.AND.ITE.GE.NIT0) GOTO 200
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
100   CONTINUE
C
      WRITE (CPARAM,'(I15,2D25.16)') NIT,FR,RES
      CALL OMSG(71,'IE337 ')
      CALL OMSG(72,'IE337 ')
C
      IF (RES.GE.1D-70) THEN
       CAPPA=(FR/RES)**(1D0/NIT)
      ELSE
       CAPPA=0D0
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE337 ')
C
      IER=1
      GOTO 99999
C
200   IER=0
      IF (RES.GE.1D-70) RES=FR/RES
      WRITE (CPARAM,'(I15,2D25.16)') ITE,FR,RES
      CALL OMSG(72,'IE337 ')
C
      IF (ITE.EQ.0) THEN
       CAPPA=0D0
      ELSE
       CAPPA=RES**(1D0/ITE)
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE337 ')
C
99999 END
C
C
C
      SUBROUTINE IE338(VA,KCOL,KLD,DX,DB,NEQ,NIT,ITE,EPS,OMEGA,
     *                 DR,DD,DD1,DG)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VA(*),KCOL(*),KLD(*)
      DIMENSION DX(*),DB(*),DR(*),DG(*),DD(*),DD1(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
      DATA FR/0D0/
C
      SUB='IE338'
      IF (ICHECK.GE.997) CALL OTRC('IE338 ','11/12/89')
C
      BMSG2=M.GE.2.OR.MT.GE.2
      NIT0=MAX(ITE,0)
C
      BNOCON=OMEGA.LT.0D0
C
      IF (ICHECK.GT.0) THEN
       CALL LLI1(DB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL1(DX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IE338 ')
        GOTO 99999
       ENDIF
      ENDIF
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
      DO 100 ITE=1,NIT
C
      CALL LAX38(VA,KCOL,KLD,NEQ,DD,DD1,1D0,0D0)
      CALL LSP1(DD,DD1,NEQ,ALPHA)
      ALPHA=SIGMA0/ALPHA
      CALL LLC1(DD,DX,NEQ,ALPHA,1D0)
      CALL LLC1(DD1,DR,NEQ,ALPHA,1D0)
C
      CALL LL21(DR,NEQ,FR)
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,FR
       CALL OMSG(73,'IE338 ')
      ENDIF
      IF (FR.LE.RES*EPS.AND.ITE.GE.NIT0) GOTO 200
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
100   CONTINUE
C
      WRITE (CPARAM,'(I15,2D25.16)') NIT,FR,RES
      CALL OMSG(71,'IE338 ')
      CALL OMSG(72,'IE338 ')
C
      IF (RES.GE.1D-70) THEN
       CAPPA=(FR/RES)**(1D0/NIT)
      ELSE
       CAPPA=0D0
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE338 ')
C
      IER=1
      GOTO 99999
C
200   IER=0
      IF (RES.GE.1D-70) RES=FR/RES
      WRITE (CPARAM,'(I15,2D25.16)') ITE,FR,RES
      CALL OMSG(72,'IE338 ')
C
      IF (ITE.EQ.0) THEN
       CAPPA=0D0
      ELSE
       CAPPA=RES**(1D0/ITE)
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE338 ')
C
99999 END
C
C
C
      SUBROUTINE IE33A(VA,KCOL,KLD,KOP,DX,DB,NEQ,NIT,ITE,EPS,OMEGA,
     *                 DR,DD,DD1,DG)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VA(*),KCOL(*),KLD(*),KOP(*)
      DIMENSION DX(*),DB(*),DR(*),DG(*),DD(*),DD1(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
      DATA FR/0D0/
C
      SUB='IE33A'
      IF (ICHECK.GE.997) CALL OTRC('IE33A ','11/12/89')
C
      BMSG2=M.GE.2.OR.MT.GE.2
      NIT0=MAX(ITE,0)
C
      BNOCON=OMEGA.LT.0D0
C
      IF (ICHECK.GT.0) THEN
       CALL LLI1(DB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL1(DX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IE33A ')
        GOTO 99999
       ENDIF
      ENDIF
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
      DO 100 ITE=1,NIT
C
      CALL LAX3A(VA,KCOL,KLD,KOP,NEQ,DD,DD1,1D0,0D0)
      CALL LSP1(DD,DD1,NEQ,ALPHA)
      ALPHA=SIGMA0/ALPHA
      CALL LLC1(DD,DX,NEQ,ALPHA,1D0)
      CALL LLC1(DD1,DR,NEQ,ALPHA,1D0)
C
      CALL LL21(DR,NEQ,FR)
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,FR
       CALL OMSG(73,'IE33A ')
      ENDIF
      IF (FR.LE.RES*EPS.AND.ITE.GE.NIT0) GOTO 200
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
100   CONTINUE
C
      WRITE (CPARAM,'(I15,2D25.16)') NIT,FR,RES
      CALL OMSG(71,'IE33A ')
      CALL OMSG(72,'IE33A ')
C
      IF (RES.GE.1D-70) THEN
       CAPPA=(FR/RES)**(1D0/NIT)
      ELSE
       CAPPA=0D0
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE33A ')
C
      IER=1
      GOTO 99999
C
200   IER=0
      IF (RES.GE.1D-70) RES=FR/RES
      WRITE (CPARAM,'(I15,2D25.16)') ITE,FR,RES
      CALL OMSG(72,'IE33A ')
C
      IF (ITE.EQ.0) THEN
       CAPPA=0D0
      ELSE
       CAPPA=RES**(1D0/ITE)
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE33A ')
C
99999 END
