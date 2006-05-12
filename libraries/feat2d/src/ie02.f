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
* IE02n                                                                *
*                                                                      *
* Purpose  Solution of a linear system  A*X = B  using                 *
*          a preconditioned conjugate gradient method                  *
*          Single/single precision version                             *
*                                                                      *
* Subroutines/functions called  LSP2 , LLC2 , LL22 , LLI2 , LCL2 ,     *
*                               LAX2n , ID12n                          *
*                                                                      *
* Version from  10/27/89                                               *
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
* NIT      I*4    Maximum number of iterations                         *
* EPS      R*8    Desired precision, stop if !!RES!! < EPS             *
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
      SUBROUTINE IE023(VA,KDIA,KDIAS,NDIA,VX,VB,NEQ,NIT,ITE,EPS,OMEGA,
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
      SUB='IE023'
      IF (ICHECK.GE.997) CALL OTRC('IE023 ','02/01/91')
C
      BMSG2=M.GE.2.OR.MT.GE.2
C
      BNOCON=OMEGA.LT.0D0
C
      IF (ICHECK.GT.0) THEN
       CALL LLI2(VB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL2(VX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IE023 ')
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
       CALL OMSG(73,'IE023 ')
      ENDIF
      IF (FR.LE.EPS) GOTO 200
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
      CALL OMSG(71,'IE023 ')
      CALL OMSG(72,'IE023 ')
C
      IF (RES.GE.1D-70) THEN
       CAPPA=(FR/RES)**(1D0/NIT)
      ELSE
       CAPPA=0D0
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE023 ')
C
      IER=1
      GOTO 99999
C
200   IER=0
      IF (RES.GE.1D-70) RES=FR/RES
      WRITE (CPARAM,'(I15,2D25.16)') ITE,FR,RES
      CALL OMSG(72,'IE023 ')
C
      IF (ITE.EQ.0) THEN
       CAPPA=0D0
      ELSE
       CAPPA=RES**(1D0/ITE)
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE023 ')
C
99999 END
C
C
C
      SUBROUTINE IE024(VA,KDIA,KDIAS,NDIA,VX,VB,NEQ,NIT,ITE,EPS,OMEGA,
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
      SUB='IE024'
      IF (ICHECK.GE.997) CALL OTRC('IE024 ','02/01/91')
C
      BMSG2=M.GE.2.OR.MT.GE.2
C
      BNOCON=OMEGA.LT.0D0
C
      IF (ICHECK.GT.0) THEN
       CALL LLI2(VB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL2(VX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IE024 ')
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
       CALL OMSG(73,'IE024 ')
      ENDIF
      IF (FR.LE.EPS) GOTO 200
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
      CALL OMSG(71,'IE024 ')
      CALL OMSG(72,'IE024 ')
C
      IF (RES.GE.1D-70) THEN
       CAPPA=(FR/RES)**(1D0/NIT)
      ELSE
       CAPPA=0D0
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE024 ')
C
      IER=1
      GOTO 99999
C
200   IER=0
      IF (RES.GE.1D-70) RES=FR/RES
      WRITE (CPARAM,'(I15,2D25.16)') ITE,FR,RES
      CALL OMSG(72,'IE024 ')
C
      IF (ITE.EQ.0) THEN
       CAPPA=0D0
      ELSE
       CAPPA=RES**(1D0/ITE)
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE024 ')
C
99999 END
C
C
C
      SUBROUTINE IE027(VA,KCOL,KLD,VX,VB,NEQ,NIT,ITE,EPS,OMEGA,
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
      SUB='IE027'
      IF (ICHECK.GE.997) CALL OTRC('IE027 ','10/27/89')
C
      BMSG2=M.GE.2.OR.MT.GE.2
C
      BNOCON=OMEGA.LT.0D0
C
      IF (ICHECK.GT.0) THEN
       CALL LLI2(VB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL2(VX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IE027 ')
        GOTO 99999
       ENDIF
      ENDIF
C
C *** Initialization
      CALL LAX27(VA,KCOL,KLD,NEQ,VX,VR,1D0,0D0)
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
       CALL OMSG(73,'IE027 ')
      ENDIF
      IF (FR.LE.EPS) GOTO 200
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
      CALL OMSG(71,'IE027 ')
      CALL OMSG(72,'IE027 ')
C
      IF (RES.GE.1D-70) THEN
       CAPPA=(FR/RES)**(1D0/NIT)
      ELSE
       CAPPA=0D0
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE027 ')
C
      IER=1
      GOTO 99999
C
200   IER=0
      IF (RES.GE.1D-70) RES=FR/RES
      WRITE (CPARAM,'(I15,2D25.16)') ITE,FR,RES
      CALL OMSG(72,'IE027 ')
C
      IF (ITE.EQ.0) THEN
       CAPPA=0D0
      ELSE
       CAPPA=RES**(1D0/ITE)
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE027 ')
C
99999 END
C
C
C
      SUBROUTINE IE028(VA,KCOL,KLD,VX,VB,NEQ,NIT,ITE,EPS,OMEGA,
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
      SUB='IE028'
      IF (ICHECK.GE.997) CALL OTRC('IE028 ','10/27/89')
C
      BMSG2=M.GE.2.OR.MT.GE.2
C
      BNOCON=OMEGA.LT.0D0
C
      IF (ICHECK.GT.0) THEN
       CALL LLI2(VB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL2(VX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IE028 ')
        GOTO 99999
       ENDIF
      ENDIF
C
C *** Initialization
      CALL LAX28(VA,KCOL,KLD,NEQ,VX,VR,1D0,0D0)
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
       CALL OMSG(73,'IE028 ')
      ENDIF
      IF (FR.LE.EPS) GOTO 200
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
      CALL OMSG(71,'IE028 ')
      CALL OMSG(72,'IE028 ')
C
      IF (RES.GE.1D-70) THEN
       CAPPA=(FR/RES)**(1D0/NIT)
      ELSE
       CAPPA=0D0
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE028 ')
C
      IER=1
      GOTO 99999
C
200   IER=0
      IF (RES.GE.1D-70) RES=FR/RES
      WRITE (CPARAM,'(I15,2D25.16)') ITE,FR,RES
      CALL OMSG(72,'IE028 ')
C
      IF (ITE.EQ.0) THEN
       CAPPA=0D0
      ELSE
       CAPPA=RES**(1D0/ITE)
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE028 ')
C
99999 END
C
C
C
      SUBROUTINE IE02A(VA,KCOL,KLD,KOP,VX,VB,NEQ,NIT,ITE,EPS,OMEGA,
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
      SUB='IE02A'
      IF (ICHECK.GE.997) CALL OTRC('IE02A ','10/27/89')
C
      BMSG2=M.GE.2.OR.MT.GE.2
C
      BNOCON=OMEGA.LT.0D0
C
      IF (ICHECK.GT.0) THEN
       CALL LLI2(VB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL2(VX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IE02A ')
        GOTO 99999
       ENDIF
      ENDIF
C
C *** Initialization
      CALL LAX2A(VA,KCOL,KLD,KOP,NEQ,VX,VR,1D0,0D0)
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
       CALL OMSG(73,'IE02A ')
      ENDIF
      IF (FR.LE.EPS) GOTO 200
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
      CALL OMSG(71,'IE02A ')
      CALL OMSG(72,'IE02A ')
C
      IF (RES.GE.1D-70) THEN
       CAPPA=(FR/RES)**(1D0/NIT)
      ELSE
       CAPPA=0D0
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE02A ')
C
      IER=1
      GOTO 99999
C
200   IER=0
      IF (RES.GE.1D-70) RES=FR/RES
      WRITE (CPARAM,'(I15,2D25.16)') ITE,FR,RES
      CALL OMSG(72,'IE02A ')
C
      IF (ITE.EQ.0) THEN
       CAPPA=0D0
      ELSE
       CAPPA=RES**(1D0/ITE)
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE02A ')
C
99999 END
