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
* IE31n                                                                *
*                                                                      *
* Purpose  Solution of a linear system  A*X = B  using                 *
*          a preconditioned conjugate gradient method                  *
*          Double/double precision version                             *
*                                                                      *
* Subroutines/functions called  LSP1 , LLC1 , LL21 , LLI1 , LCL1 ,     *
*                               LAX1n , ID11n                          *
*                                                                      *
* Version from  11/12/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DA       R*8    Matrix in storage technique n                        *
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
      SUBROUTINE IE313(DA,KDIA,KDIAS,NDIA,DX,DB,NEQ,NIT,ITE,EPS,OMEGA,
     *                 DR,DD,DD1,DG)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION DA(*),KDIA(*),KDIAS(*)
      DIMENSION DX(*),DB(*),DR(*),DG(*),DD(*),DD1(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
      DATA FR/0D0/
C
      SUB='IE313'
      IF (ICHECK.GE.997) CALL OTRC('IE313 ','02/01/91')
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
        IF (BMSG2) CALL OMSG(70,'IE313 ')
        GOTO 99999
       ENDIF
      ENDIF
C
C *** Initialization
      CALL LAX13(DA,KDIA,KDIAS,NDIA,NEQ,DX,DR,1D0,0D0)
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
       CALL IA113(DA,DG,NEQ)
       CALL LSP1(DR,DG,NEQ,SIGMA0)
      ENDIF
C
      CALL LLC1(DG,DD,NEQ,-1D0,0D0)
C
C *** Iterative correction
      DO 100 ITE=1,NIT
C
      CALL LAX13(DA,KDIA,KDIAS,NDIA,NEQ,DD,DD1,1D0,0D0)
      CALL LSP1(DD,DD1,NEQ,ALPHA)
      ALPHA=SIGMA0/ALPHA
      CALL LLC1(DD,DX,NEQ,ALPHA,1D0)
      CALL LLC1(DD1,DR,NEQ,ALPHA,1D0)
C
      CALL LL21(DR,NEQ,FR)
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,FR
       CALL OMSG(73,'IE313 ')
      ENDIF
      IF (FR.LE.RES*EPS.AND.ITE.GE.NIT0) GOTO 200
C
      IF (BNOCON) THEN
       SIGMA1=FR*FR
      ELSE
       CALL LCP1(DR,DG,NEQ)
       CALL IA113(DA,DG,NEQ)
       CALL LSP1(DR,DG,NEQ,SIGMA1)
      ENDIF
C
      GAMMA=SIGMA1/SIGMA0
      SIGMA0=SIGMA1
      CALL LLC1(DG,DD,NEQ,-1D0,GAMMA)
100   CONTINUE
C
      WRITE (CPARAM,'(I15,2D25.16)') NIT,FR,RES
      CALL OMSG(71,'IE313 ')
      CALL OMSG(72,'IE313 ')
C
      IF (RES.GE.1D-70) THEN
       CAPPA=(FR/RES)**(1D0/NIT)
      ELSE
       CAPPA=0D0
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE313 ')
C
      IER=1
      GOTO 99999
C
200   IER=0
      IF (RES.GE.1D-70) RES=FR/RES
      WRITE (CPARAM,'(I15,2D25.16)') ITE,FR,RES
      CALL OMSG(72,'IE313 ')
C
      IF (ITE.EQ.0) THEN
       CAPPA=0D0
      ELSE
       CAPPA=RES**(1D0/ITE)
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE313 ')
C
99999 END
C
C
C
      SUBROUTINE IE314(DA,KDIA,KDIAS,NDIA,DX,DB,NEQ,NIT,ITE,EPS,OMEGA,
     *                 DR,DD,DD1,DG)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION DA(*),KDIA(*),KDIAS(*)
      DIMENSION DX(*),DB(*),DR(*),DG(*),DD(*),DD1(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
      DATA FR/0D0/
C
      SUB='IE314'
      IF (ICHECK.GE.997) CALL OTRC('IE314 ','02/01/91')
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
        IF (BMSG2) CALL OMSG(70,'IE314 ')
        GOTO 99999
       ENDIF
      ENDIF
C
C *** Initialization
      CALL LAX14(DA,KDIA,KDIAS,NDIA,NEQ,DX,DR,1D0,0D0)
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
       CALL IA113(DA,DG,NEQ)
       CALL LSP1(DR,DG,NEQ,SIGMA0)
      ENDIF
C
      CALL LLC1(DG,DD,NEQ,-1D0,0D0)
C
C *** Iterative correction
      DO 100 ITE=1,NIT
C
      CALL LAX14(DA,KDIA,KDIAS,NDIA,NEQ,DD,DD1,1D0,0D0)
      CALL LSP1(DD,DD1,NEQ,ALPHA)
      ALPHA=SIGMA0/ALPHA
      CALL LLC1(DD,DX,NEQ,ALPHA,1D0)
      CALL LLC1(DD1,DR,NEQ,ALPHA,1D0)
C
      CALL LL21(DR,NEQ,FR)
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,FR
       CALL OMSG(73,'IE314 ')
      ENDIF
      IF (FR.LE.RES*EPS.AND.ITE.GE.NIT0) GOTO 200
C
      IF (BNOCON) THEN
       SIGMA1=FR*FR
      ELSE
       CALL LCP1(DR,DG,NEQ)
       CALL IA113(DA,DG,NEQ)
       CALL LSP1(DR,DG,NEQ,SIGMA1)
      ENDIF
C
      GAMMA=SIGMA1/SIGMA0
      SIGMA0=SIGMA1
      CALL LLC1(DG,DD,NEQ,-1D0,GAMMA)
100   CONTINUE
C
      WRITE (CPARAM,'(I15,2D25.16)') NIT,FR,RES
      CALL OMSG(71,'IE314 ')
      CALL OMSG(72,'IE314 ')
C
      IF (RES.GE.1D-70) THEN
       CAPPA=(FR/RES)**(1D0/NIT)
      ELSE
       CAPPA=0D0
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE314 ')
C
      IER=1
      GOTO 99999
C
200   IER=0
      IF (RES.GE.1D-70) RES=FR/RES
      WRITE (CPARAM,'(I15,2D25.16)') ITE,FR,RES
      CALL OMSG(72,'IE314 ')
C
      IF (ITE.EQ.0) THEN
       CAPPA=0D0
      ELSE
       CAPPA=RES**(1D0/ITE)
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE314 ')
C
99999 END
C
C
C
      SUBROUTINE IE317(DA,KCOL,KLD,DX,DB,NEQ,NIT,ITE,EPS,OMEGA,
     *                 DR,DD,DD1,DG)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION DA(*),KCOL(*),KLD(*)
      DIMENSION DX(*),DB(*),DR(*),DG(*),DD(*),DD1(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
      DATA FR/0D0/
C
      SUB='IE317'
      IF (ICHECK.GE.997) CALL OTRC('IE317 ','11/12/89')
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
        IF (BMSG2) CALL OMSG(70,'IE317 ')
        GOTO 99999
       ENDIF
      ENDIF
C
C *** Initialization
      CALL LAX17(DA,KCOL,KLD,NEQ,DX,DR,1D0,0D0)
      CALL LLC1(DB,DR,NEQ,-1D0,1D0)
      CALL LL21(DR,NEQ,RES)
C
      IF (BNOCON) THEN
       SIGMA0=RES*RES
      ELSE
       CALL LCP1(DR,DG,NEQ)
       CALL ID117(DA,KCOL,KLD,DG,NEQ,OMEGA)
       CALL LSP1(DR,DG,NEQ,SIGMA0)
      ENDIF
C
      CALL LLC1(DG,DD,NEQ,-1D0,0D0)
C
C *** Iterative correction
      DO 100 ITE=1,NIT
C
      CALL LAX17(DA,KCOL,KLD,NEQ,DD,DD1,1D0,0D0)
      CALL LSP1(DD,DD1,NEQ,ALPHA)
      ALPHA=SIGMA0/ALPHA
      CALL LLC1(DD,DX,NEQ,ALPHA,1D0)
      CALL LLC1(DD1,DR,NEQ,ALPHA,1D0)
C
      CALL LL21(DR,NEQ,FR)
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,FR
       CALL OMSG(73,'IE317 ')
      ENDIF
      IF (FR.LE.RES*EPS.AND.ITE.GE.NIT0) GOTO 200
C
      IF (BNOCON) THEN
       SIGMA1=FR*FR
      ELSE
       CALL LCP1(DR,DG,NEQ)
       CALL ID117(DA,KCOL,KLD,DG,NEQ,OMEGA)
       CALL LSP1(DR,DG,NEQ,SIGMA1)
      ENDIF
C
      GAMMA=SIGMA1/SIGMA0
      SIGMA0=SIGMA1
      CALL LLC1(DG,DD,NEQ,-1D0,GAMMA)
100   CONTINUE
C
      WRITE (CPARAM,'(I15,2D25.16)') NIT,FR,RES
      CALL OMSG(71,'IE317 ')
      CALL OMSG(72,'IE317 ')
C
      IF (RES.GE.1D-70) THEN
       CAPPA=(FR/RES)**(1D0/NIT)
      ELSE
       CAPPA=0D0
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE317 ')
C
      IER=1
      GOTO 99999
C
200   IER=0
      IF (RES.GE.1D-70) RES=FR/RES
      WRITE (CPARAM,'(I15,2D25.16)') ITE,FR,RES
      CALL OMSG(72,'IE317 ')
C
      IF (ITE.EQ.0) THEN
       CAPPA=0D0
      ELSE
       CAPPA=RES**(1D0/ITE)
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE317 ')
C
99999 END
C
C
C
      SUBROUTINE IE318(DA,KCOL,KLD,DX,DB,NEQ,NIT,ITE,EPS,OMEGA,
     *                 DR,DD,DD1,DG)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION DA(*),KCOL(*),KLD(*)
      DIMENSION DX(*),DB(*),DR(*),DG(*),DD(*),DD1(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
      DATA FR/0D0/
C
      SUB='IE318'
      IF (ICHECK.GE.997) CALL OTRC('IE318 ','11/12/89')
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
        IF (BMSG2) CALL OMSG(70,'IE318 ')
        GOTO 99999
       ENDIF
      ENDIF
C
C *** Initialization
      CALL LAX18(DA,KCOL,KLD,NEQ,DX,DR,1D0,0D0)
      CALL LLC1(DB,DR,NEQ,-1D0,1D0)
      CALL LL21(DR,NEQ,RES)
C
      IF (BNOCON) THEN
       SIGMA0=RES*RES
      ELSE
       CALL LCP1(DR,DG,NEQ)
       CALL ID118(DA,KCOL,KLD,DG,NEQ,OMEGA)
       CALL LSP1(DR,DG,NEQ,SIGMA0)
      ENDIF
C
      CALL LLC1(DG,DD,NEQ,-1D0,0D0)
C
C *** Iterative correction
      DO 100 ITE=1,NIT
C
      CALL LAX18(DA,KCOL,KLD,NEQ,DD,DD1,1D0,0D0)
      CALL LSP1(DD,DD1,NEQ,ALPHA)
      ALPHA=SIGMA0/ALPHA
      CALL LLC1(DD,DX,NEQ,ALPHA,1D0)
      CALL LLC1(DD1,DR,NEQ,ALPHA,1D0)
C
      CALL LL21(DR,NEQ,FR)
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,FR
       CALL OMSG(73,'IE318 ')
      ENDIF
      IF (FR.LE.RES*EPS.AND.ITE.GE.NIT0) GOTO 200
C
      IF (BNOCON) THEN
       SIGMA1=FR*FR
      ELSE
       CALL LCP1(DR,DG,NEQ)
       CALL ID118(DA,KCOL,KLD,DG,NEQ,OMEGA)
       CALL LSP1(DR,DG,NEQ,SIGMA1)
      ENDIF
C
      GAMMA=SIGMA1/SIGMA0
      SIGMA0=SIGMA1
      CALL LLC1(DG,DD,NEQ,-1D0,GAMMA)
100   CONTINUE
C
      WRITE (CPARAM,'(I15,2D25.16)') NIT,FR,RES
      CALL OMSG(71,'IE318 ')
      CALL OMSG(72,'IE318 ')
C
      IF (RES.GE.1D-70) THEN
       CAPPA=(FR/RES)**(1D0/NIT)
      ELSE
       CAPPA=0D0
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE318 ')
C
      IER=1
      GOTO 99999
C
200   IER=0
      IF (RES.GE.1D-70) RES=FR/RES
      WRITE (CPARAM,'(I15,2D25.16)') ITE,FR,RES
      CALL OMSG(72,'IE318 ')
C
      IF (ITE.EQ.0) THEN
       CAPPA=0D0
      ELSE
       CAPPA=RES**(1D0/ITE)
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE318 ')
C
99999 END
C
C
C
      SUBROUTINE IE31A(DA,KCOL,KLD,KOP,DX,DB,NEQ,NIT,ITE,EPS,OMEGA,
     *                 DR,DD,DD1,DG)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION DA(*),KCOL(*),KLD(*),KOP(*)
      DIMENSION DX(*),DB(*),DR(*),DG(*),DD(*),DD1(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
      DATA FR/0D0/
C
      SUB='IE31A'
      IF (ICHECK.GE.997) CALL OTRC('IE31A ','11/12/89')
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
        IF (BMSG2) CALL OMSG(70,'IE31A ')
        GOTO 99999
       ENDIF
      ENDIF
C
C *** Initialization
      CALL LAX1A(DA,KCOL,KLD,KOP,NEQ,DX,DR,1D0,0D0)
      CALL LLC1(DB,DR,NEQ,-1D0,1D0)
      CALL LL21(DR,NEQ,RES)
C
      IF (BNOCON) THEN
       SIGMA0=RES*RES
      ELSE
       CALL LCP1(DR,DG,NEQ)
       CALL ID11A(DA,KCOL,KLD,KOP,DG,NEQ,OMEGA)
       CALL LSP1(DR,DG,NEQ,SIGMA0)
      ENDIF
C
      CALL LLC1(DG,DD,NEQ,-1D0,0D0)
C
C *** Iterative correction
      DO 100 ITE=1,NIT
C
      CALL LAX1A(DA,KCOL,KLD,KOP,NEQ,DD,DD1,1D0,0D0)
      CALL LSP1(DD,DD1,NEQ,ALPHA)
      ALPHA=SIGMA0/ALPHA
      CALL LLC1(DD,DX,NEQ,ALPHA,1D0)
      CALL LLC1(DD1,DR,NEQ,ALPHA,1D0)
C
      CALL LL21(DR,NEQ,FR)
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,FR
       CALL OMSG(73,'IE31A ')
      ENDIF
      IF (FR.LE.RES*EPS.AND.ITE.GE.NIT0) GOTO 200
C
      IF (BNOCON) THEN
       SIGMA1=FR*FR
      ELSE
       CALL LCP1(DR,DG,NEQ)
       CALL ID11A(DA,KCOL,KLD,KOP,DG,NEQ,OMEGA)
       CALL LSP1(DR,DG,NEQ,SIGMA1)
      ENDIF
C
      GAMMA=SIGMA1/SIGMA0
      SIGMA0=SIGMA1
      CALL LLC1(DG,DD,NEQ,-1D0,GAMMA)
100   CONTINUE
C
      WRITE (CPARAM,'(I15,2D25.16)') NIT,FR,RES
      CALL OMSG(71,'IE31A ')
      CALL OMSG(72,'IE31A ')
C
      IF (RES.GE.1D-70) THEN
       CAPPA=(FR/RES)**(1D0/NIT)
      ELSE
       CAPPA=0D0
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE31A ')
C
      IER=1
      GOTO 99999
C
200   IER=0
      IF (RES.GE.1D-70) RES=FR/RES
      WRITE (CPARAM,'(I15,2D25.16)') ITE,FR,RES
      CALL OMSG(72,'IE31A ')
C
      IF (ITE.EQ.0) THEN
       CAPPA=0D0
      ELSE
       CAPPA=RES**(1D0/ITE)
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE31A ')
C
99999 END
