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
* IE32n                                                                *
*                                                                      *
* Purpose  Solution of a linear system  A*X = B  using                 *
*          a preconditioned conjugate gradient method                  *
*          Single/single precision version                             *
*                                                                      *
* Subroutines/functions called  LSP2 , LLC2 , LL22 , LLI2 , LCL2 ,     *
*                               LAX2n , ID12n                          *
*                                                                      *
* Version from  11/12/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* VA       R*4    Matrix in storage technique n                        *
* KCOL     I*4                                                         *
* KLD      I*4                                                         *
* KOP      I*4                                                         *
* VX       R*4    Starting vector                                      *
* VB       R*4    Right hand side                                      *
* NEQ      I*4    Number of equations                                  *
* ITE      I*4    Minimum number of iterations                         *
* NIT      I*4    Maximum number of iterations                         *
* EPS      R*8    Desired precision, stop if !!RES!!/!!RES0!! < EPS    *
* OMEGA    R*8    0 <  OMEGA     No Preconditioning                    *
*                 0 >= OMEGA  -  SSOR Preconditioning                  *
* VR,VD    R*4    Workspace vectors of length NEQ                      *
* VD1,VG   R*4    For OMEGA < 0 , VG must be replaced by VR            *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* VX       R*4    Solution vector                                      *
* ITE      I*4    Number of iterations                                 *
* IER      I*4    Error indicator                                      *
*                 +1  Precision EPS not achieved after NIT iterations  *
*                                                                      *
************************************************************************
C
      SUBROUTINE IE323(VA,KDIA,KDIAS,NDIA,VX,VB,NEQ,NIT,ITE,EPS,OMEGA,
     *                 VR,VD,VD1,VG)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VA(*),KDIA(*),KDIAS(*)
      DIMENSION VX(*),VB(*),VR(*),VG(*),VD(*),VD1(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
      DATA FR/0D0/
C
      SUB='IE323'
      IF (ICHECK.GE.997) CALL OTRC('IE323 ','02/01/91')
C
      BMSG2=M.GE.2.OR.MT.GE.2
      NIT0=MAX(ITE,0)
C
      BNOCON=OMEGA.LT.0D0
C
      IF (ICHECK.GT.0) THEN
       CALL LLI2(VB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL2(VX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IE323 ')
        GOTO 99999
       ENDIF
      ENDIF
C
C *** Initialization
      CALL LAX23(VA,KDIA,KDIAS,NDIA,NEQ,VX,VR,1D0,0D0)
      CALL LLC2(VB,VR,NEQ,-1D0,1D0)
      CALL LL22(VR,NEQ,RES)
      IF (RES.LE.EPS) THEN
       ITE=0
       FR=RES
       GOTO 200
      ENDIF
C
      IF (BNOCON) THEN
       SIGMA0=RES*RES
      ELSE
       CALL LCP2(VR,VG,NEQ)
       CALL IA123(VA,VG,NEQ)
       CALL LSP2(VR,VG,NEQ,SIGMA0)
      ENDIF
C
      CALL LLC2(VG,VD,NEQ,-1D0,0D0)
C
C *** Iterative correction
      DO 100 ITE=1,NIT
C
      CALL LAX23(VA,KDIA,KDIAS,NDIA,NEQ,VD,VD1,1D0,0D0)
      CALL LSP2(VD,VD1,NEQ,ALPHA)
      ALPHA=SIGMA0/ALPHA
      CALL LLC2(VD,VX,NEQ,ALPHA,1D0)
      CALL LLC2(VD1,VR,NEQ,ALPHA,1D0)
C
      CALL LL22(VR,NEQ,FR)
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,FR
       CALL OMSG(73,'IE323 ')
      ENDIF
      IF (FR.LE.RES*EPS.AND.ITE.GE.NIT0) GOTO 200
C
      IF (BNOCON) THEN
       SIGMA1=FR*FR
      ELSE
       CALL LCP2(VR,VG,NEQ)
       CALL IA123(VA,VG,NEQ)
       CALL LSP2(VR,VG,NEQ,SIGMA1)
      ENDIF
C
      GAMMA=SIGMA1/SIGMA0
      SIGMA0=SIGMA1
      CALL LLC2(VG,VD,NEQ,-1D0,GAMMA)
100   CONTINUE
C
      WRITE (CPARAM,'(I15,2D25.16)') NIT,FR,RES
      CALL OMSG(71,'IE323 ')
      CALL OMSG(72,'IE323 ')
C
      IF (RES.GE.1D-70) THEN
       CAPPA=(FR/RES)**(1D0/NIT)
      ELSE
       CAPPA=0D0
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE323 ')
C
      IER=1
      GOTO 99999
C
200   IER=0
      IF (RES.GE.1D-70) RES=FR/RES
      WRITE (CPARAM,'(I15,2D25.16)') ITE,FR,RES
      CALL OMSG(72,'IE323 ')
C
      IF (ITE.EQ.0) THEN
       CAPPA=0D0
      ELSE
       CAPPA=RES**(1D0/ITE)
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE323 ')
C
99999 END
C
C
C
      SUBROUTINE IE324(VA,KDIA,KDIAS,NDIA,VX,VB,NEQ,NIT,ITE,EPS,OMEGA,
     *                 VR,VD,VD1,VG)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VA(*),KDIA(*),KDIAS(*)
      DIMENSION VX(*),VB(*),VR(*),VG(*),VD(*),VD1(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
      DATA FR/0D0/
C
      SUB='IE324'
      IF (ICHECK.GE.997) CALL OTRC('IE324 ','02/01/91')
C
      BMSG2=M.GE.2.OR.MT.GE.2
      NIT0=MAX(ITE,0)
C
      BNOCON=OMEGA.LT.0D0
C
      IF (ICHECK.GT.0) THEN
       CALL LLI2(VB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL2(VX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IE324 ')
        GOTO 99999
       ENDIF
      ENDIF
C
C *** Initialization
      CALL LAX24(VA,KDIA,KDIAS,NDIA,NEQ,VX,VR,1D0,0D0)
      CALL LLC2(VB,VR,NEQ,-1D0,1D0)
      CALL LL22(VR,NEQ,RES)
      IF (RES.LE.EPS) THEN
       ITE=0
       FR=RES
       GOTO 200
      ENDIF
C
      IF (BNOCON) THEN
       SIGMA0=RES*RES
      ELSE
       CALL LCP2(VR,VG,NEQ)
       CALL IA123(VA,VG,NEQ)
       CALL LSP2(VR,VG,NEQ,SIGMA0)
      ENDIF
C
      CALL LLC2(VG,VD,NEQ,-1D0,0D0)
C
C *** Iterative correction
      DO 100 ITE=1,NIT
C
      CALL LAX24(VA,KDIA,KDIAS,NDIA,NEQ,VD,VD1,1D0,0D0)
      CALL LSP2(VD,VD1,NEQ,ALPHA)
      ALPHA=SIGMA0/ALPHA
      CALL LLC2(VD,VX,NEQ,ALPHA,1D0)
      CALL LLC2(VD1,VR,NEQ,ALPHA,1D0)
C
      CALL LL22(VR,NEQ,FR)
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,FR
       CALL OMSG(73,'IE324 ')
      ENDIF
      IF (FR.LE.RES*EPS.AND.ITE.GE.NIT0) GOTO 200
C
      IF (BNOCON) THEN
       SIGMA1=FR*FR
      ELSE
       CALL LCP2(VR,VG,NEQ)
       CALL IA123(VA,VG,NEQ)
       CALL LSP2(VR,VG,NEQ,SIGMA1)
      ENDIF
C
      GAMMA=SIGMA1/SIGMA0
      SIGMA0=SIGMA1
      CALL LLC2(VG,VD,NEQ,-1D0,GAMMA)
100   CONTINUE
C
      WRITE (CPARAM,'(I15,2D25.16)') NIT,FR,RES
      CALL OMSG(71,'IE324 ')
      CALL OMSG(72,'IE324 ')
C
      IF (RES.GE.1D-70) THEN
       CAPPA=(FR/RES)**(1D0/NIT)
      ELSE
       CAPPA=0D0
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE324 ')
C
      IER=1
      GOTO 99999
C
200   IER=0
      IF (RES.GE.1D-70) RES=FR/RES
      WRITE (CPARAM,'(I15,2D25.16)') ITE,FR,RES
      CALL OMSG(72,'IE324 ')
C
      IF (ITE.EQ.0) THEN
       CAPPA=0D0
      ELSE
       CAPPA=RES**(1D0/ITE)
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE324 ')
C
99999 END
C
C
C
      SUBROUTINE IE327(VA,KCOL,KLD,VX,VB,NEQ,NIT,ITE,EPS,OMEGA,
     *                 VR,VD,VD1,VG)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VA(*),KCOL(*),KLD(*)
      DIMENSION VX(*),VB(*),VR(*),VG(*),VD(*),VD1(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
      DATA FR/0D0/
C
      SUB='IE327'
      IF (ICHECK.GE.997) CALL OTRC('IE327 ','11/12/89')
C
      BMSG2=M.GE.2.OR.MT.GE.2
      NIT0=MAX(ITE,0)
C
      BNOCON=OMEGA.LT.0D0
C
      IF (ICHECK.GT.0) THEN
       CALL LLI2(VB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL2(VX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IE327 ')
        GOTO 99999
       ENDIF
      ENDIF
C
C *** Initialization
      CALL LAX27(VA,KCOL,KLD,NEQ,VX,VR,1D0,0D0)
      CALL LLC2(VB,VR,NEQ,-1D0,1D0)
      CALL LL22(VR,NEQ,RES)
C
      IF (BNOCON) THEN
       SIGMA0=RES*RES
      ELSE
       CALL LCP2(VR,VG,NEQ)
       CALL ID127(VA,KCOL,KLD,VG,NEQ,OMEGA)
       CALL LSP2(VR,VG,NEQ,SIGMA0)
      ENDIF
C
      CALL LLC2(VG,VD,NEQ,-1D0,0D0)
C
C *** Iterative correction
      DO 100 ITE=1,NIT
C
      CALL LAX27(VA,KCOL,KLD,NEQ,VD,VD1,1D0,0D0)
      CALL LSP2(VD,VD1,NEQ,ALPHA)
      ALPHA=SIGMA0/ALPHA
      CALL LLC2(VD,VX,NEQ,ALPHA,1D0)
      CALL LLC2(VD1,VR,NEQ,ALPHA,1D0)
C
      CALL LL22(VR,NEQ,FR)
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,FR
       CALL OMSG(73,'IE327 ')
      ENDIF
      IF (FR.LE.RES*EPS.AND.ITE.GE.NIT0) GOTO 200
C
      IF (BNOCON) THEN
       SIGMA1=FR*FR
      ELSE
       CALL LCP2(VR,VG,NEQ)
       CALL ID127(VA,KCOL,KLD,VG,NEQ,OMEGA)
       CALL LSP2(VR,VG,NEQ,SIGMA1)
      ENDIF
C
      GAMMA=SIGMA1/SIGMA0
      SIGMA0=SIGMA1
      CALL LLC2(VG,VD,NEQ,-1D0,GAMMA)
100   CONTINUE
C
      WRITE (CPARAM,'(I15,2D25.16)') NIT,FR,RES
      CALL OMSG(71,'IE327 ')
      CALL OMSG(72,'IE327 ')
C
      IF (RES.GE.1D-70) THEN
       CAPPA=(FR/RES)**(1D0/NIT)
      ELSE
       CAPPA=0D0
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE327 ')
C
      IER=1
      GOTO 99999
C
200   IER=0
      IF (RES.GE.1D-70) RES=FR/RES
      WRITE (CPARAM,'(I15,2D25.16)') ITE,FR,RES
      CALL OMSG(72,'IE327 ')
C
      IF (ITE.EQ.0) THEN
       CAPPA=0D0
      ELSE
       CAPPA=RES**(1D0/ITE)
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE327 ')
C
99999 END
C
C
C
      SUBROUTINE IE328(VA,KCOL,KLD,VX,VB,NEQ,NIT,ITE,EPS,OMEGA,
     *                 VR,VD,VD1,VG)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VA(*),KCOL(*),KLD(*)
      DIMENSION VX(*),VB(*),VR(*),VG(*),VD(*),VD1(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
      DATA FR/0D0/
C
      SUB='IE328'
      IF (ICHECK.GE.997) CALL OTRC('IE328 ','11/12/89')
C
      BMSG2=M.GE.2.OR.MT.GE.2
      NIT0=MAX(ITE,0)
C
      BNOCON=OMEGA.LT.0D0
C
      IF (ICHECK.GT.0) THEN
       CALL LLI2(VB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL2(VX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IE328 ')
        GOTO 99999
       ENDIF
      ENDIF
C
C *** Initialization
      CALL LAX28(VA,KCOL,KLD,NEQ,VX,VR,1D0,0D0)
      CALL LLC2(VB,VR,NEQ,-1D0,1D0)
      CALL LL22(VR,NEQ,RES)
C
      IF (BNOCON) THEN
       SIGMA0=RES*RES
      ELSE
       CALL LCP2(VR,VG,NEQ)
       CALL ID128(VA,KCOL,KLD,VG,NEQ,OMEGA)
       CALL LSP2(VR,VG,NEQ,SIGMA0)
      ENDIF
C
      CALL LLC2(VG,VD,NEQ,-1D0,0D0)
C
C *** Iterative correction
      DO 100 ITE=1,NIT
C
      CALL LAX28(VA,KCOL,KLD,NEQ,VD,VD1,1D0,0D0)
      CALL LSP2(VD,VD1,NEQ,ALPHA)
      ALPHA=SIGMA0/ALPHA
      CALL LLC2(VD,VX,NEQ,ALPHA,1D0)
      CALL LLC2(VD1,VR,NEQ,ALPHA,1D0)
C
      CALL LL22(VR,NEQ,FR)
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,FR
       CALL OMSG(73,'IE328 ')
      ENDIF
      IF (FR.LE.RES*EPS.AND.ITE.GE.NIT0) GOTO 200
C
      IF (BNOCON) THEN
       SIGMA1=FR*FR
      ELSE
       CALL LCP2(VR,VG,NEQ)
       CALL ID128(VA,KCOL,KLD,VG,NEQ,OMEGA)
       CALL LSP2(VR,VG,NEQ,SIGMA1)
      ENDIF
C
      GAMMA=SIGMA1/SIGMA0
      SIGMA0=SIGMA1
      CALL LLC2(VG,VD,NEQ,-1D0,GAMMA)
100   CONTINUE
C
      WRITE (CPARAM,'(I15,2D25.16)') NIT,FR,RES
      CALL OMSG(71,'IE328 ')
      CALL OMSG(72,'IE328 ')
C
      IF (RES.GE.1D-70) THEN
       CAPPA=(FR/RES)**(1D0/NIT)
      ELSE
       CAPPA=0D0
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE328 ')
C
      IER=1
      GOTO 99999
C
200   IER=0
      IF (RES.GE.1D-70) RES=FR/RES
      WRITE (CPARAM,'(I15,2D25.16)') ITE,FR,RES
      CALL OMSG(72,'IE328 ')
C
      IF (ITE.EQ.0) THEN
       CAPPA=0D0
      ELSE
       CAPPA=RES**(1D0/ITE)
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE328 ')
C
99999 END
C
C
C
      SUBROUTINE IE32A(VA,KCOL,KLD,KOP,VX,VB,NEQ,NIT,ITE,EPS,OMEGA,
     *                 VR,VD,VD1,VG)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VA(*),KCOL(*),KLD(*),KOP(*)
      DIMENSION VX(*),VB(*),VR(*),VG(*),VD(*),VD1(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
      DATA FR/0D0/
C
      SUB='IE32A'
      IF (ICHECK.GE.997) CALL OTRC('IE32A ','11/12/89')
C
      BMSG2=M.GE.2.OR.MT.GE.2
      NIT0=MAX(ITE,0)
C
      BNOCON=OMEGA.LT.0D0
C
      IF (ICHECK.GT.0) THEN
       CALL LLI2(VB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL2(VX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IE32A ')
        GOTO 99999
       ENDIF
      ENDIF
C
C *** Initialization
      CALL LAX2A(VA,KCOL,KLD,KOP,NEQ,VX,VR,1D0,0D0)
      CALL LLC2(VB,VR,NEQ,-1D0,1D0)
      CALL LL22(VR,NEQ,RES)
C
      IF (BNOCON) THEN
       SIGMA0=RES*RES
      ELSE
       CALL LCP2(VR,VG,NEQ)
       CALL ID12A(VA,KCOL,KLD,KOP,VG,NEQ,OMEGA)
       CALL LSP2(VR,VG,NEQ,SIGMA0)
      ENDIF
C
      CALL LLC2(VG,VD,NEQ,-1D0,0D0)
C
C *** Iterative correction
      DO 100 ITE=1,NIT
C
      CALL LAX2A(VA,KCOL,KLD,KOP,NEQ,VD,VD1,1D0,0D0)
      CALL LSP2(VD,VD1,NEQ,ALPHA)
      ALPHA=SIGMA0/ALPHA
      CALL LLC2(VD,VX,NEQ,ALPHA,1D0)
      CALL LLC2(VD1,VR,NEQ,ALPHA,1D0)
C
      CALL LL22(VR,NEQ,FR)
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,FR
       CALL OMSG(73,'IE32A ')
      ENDIF
      IF (FR.LE.RES*EPS.AND.ITE.GE.NIT0) GOTO 200
C
      IF (BNOCON) THEN
       SIGMA1=FR*FR
      ELSE
       CALL LCP2(VR,VG,NEQ)
       CALL ID12A(VA,KCOL,KLD,KOP,VG,NEQ,OMEGA)
       CALL LSP2(VR,VG,NEQ,SIGMA1)
      ENDIF
C
      GAMMA=SIGMA1/SIGMA0
      SIGMA0=SIGMA1
      CALL LLC2(VG,VD,NEQ,-1D0,GAMMA)
100   CONTINUE
C
      WRITE (CPARAM,'(I15,2D25.16)') NIT,FR,RES
      CALL OMSG(71,'IE32A ')
      CALL OMSG(72,'IE32A ')
C
      IF (RES.GE.1D-70) THEN
       CAPPA=(FR/RES)**(1D0/NIT)
      ELSE
       CAPPA=0D0
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE32A ')
C
      IER=1
      GOTO 99999
C
200   IER=0
      IF (RES.GE.1D-70) RES=FR/RES
      WRITE (CPARAM,'(I15,2D25.16)') ITE,FR,RES
      CALL OMSG(72,'IE32A ')
C
      IF (ITE.EQ.0) THEN
       CAPPA=0D0
      ELSE
       CAPPA=RES**(1D0/ITE)
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE32A ')
C
99999 END
