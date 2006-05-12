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
* IE22n                                                                *
*                                                                      *
* Purpose  Smoothing using a preconditioned conjugate gradient method  *
*          Single/single precision version                             *
*                                                                      *
* Subroutines/functions called  LSP2 , LLC2 , LL22 , LAX2n , ID12n     *
*                                                                      *
* Version from  10/30/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* VA       R*4    Matrix in storage technique n                        *
* KCOL     I*4                                                         *
* KLD      I*4                                                         *
* KOP      I*4                                                         *
* VX       R*4    Starting vector                                      *
* VB       R*4    Right hand side                                      *
* NEQ      I*4    Number of smoothing steps                            *
* OMEGA    R*8    0 <  OMEGA     No Preconditioning                    *
*                 0 >= OMEGA  -  SSOR-Preconditioning                  *
* VR,VD    R*4    Workspace vectors of length NEQ                      *
* VD1,VG   R*4    For OMEGA < 0 , VG must be replaced by VR            *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* VX       R*4    Smoothed vector                                      *
*                                                                      *
************************************************************************
C
      SUBROUTINE IE227(VA,KCOL,KLD,VX,VB,NEQ,NIT,OMEGA,VR,VD,VD1,VG)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VA(*),KCOL(*),KLD(*)
      DIMENSION VX(*),VB(*),VR(*),VG(*),VD(*),VD1(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /ERRCTL/,/CHAR/
C
      SUB='IE227'
      IF (ICHECK.GE.998) CALL OTRC('IE227 ','10/30/89')
C
      BNOCON=OMEGA.LT.0D0
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
      DO 1 ITE=1,NIT
C
      CALL LAX27(VA,KCOL,KLD,NEQ,VD,VD1,1D0,0D0)
      CALL LSP2(VD,VD1,NEQ,ALPHA)
      ALPHA=SIGMA0/ALPHA
      CALL LLC2(VD,VX,NEQ,ALPHA,1D0)
      CALL LLC2(VD1,VR,NEQ,ALPHA,1D0)
C
      CALL LL22(VR,NEQ,FR)
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
1     CONTINUE
C
      END
C
C
C
      SUBROUTINE IE228(VA,KCOL,KLD,VX,VB,NEQ,NIT,OMEGA,VR,VD,VD1,VG)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VA(*),KCOL(*),KLD(*)
      DIMENSION VX(*),VB(*),VR(*),VG(*),VD(*),VD1(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /ERRCTL/,/CHAR/
C
      SUB='IE228'
      IF (ICHECK.GE.998) CALL OTRC('IE228 ','10/30/89')
C
      BNOCON=OMEGA.LT.0D0
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
      DO 1 ITE=1,NIT
C
      CALL LAX28(VA,KCOL,KLD,NEQ,VD,VD1,1D0,0D0)
      CALL LSP2(VD,VD1,NEQ,ALPHA)
      ALPHA=SIGMA0/ALPHA
      CALL LLC2(VD,VX,NEQ,ALPHA,1D0)
      CALL LLC2(VD1,VR,NEQ,ALPHA,1D0)
C
      CALL LL22(VR,NEQ,FR)
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
1     CONTINUE
C
      END
C
C
C
      SUBROUTINE IE22A(VA,KCOL,KLD,KOP,VX,VB,NEQ,NIT,OMEGA,VR,VD,VD1,VG)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VA(*),KCOL(*),KLD(*),KOP(*)
      DIMENSION VX(*),VB(*),VR(*),VG(*),VD(*),VD1(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /ERRCTL/,/CHAR/
C
      SUB='IE22A'
      IF (ICHECK.GE.998) CALL OTRC('IE22A ','10/30/89')
C
      BNOCON=OMEGA.LT.0D0
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
      DO 1 ITE=1,NIT
C
      CALL LAX2A(VA,KCOL,KLD,KOP,NEQ,VD,VD1,1D0,0D0)
      CALL LSP2(VD,VD1,NEQ,ALPHA)
      ALPHA=SIGMA0/ALPHA
      CALL LLC2(VD,VX,NEQ,ALPHA,1D0)
      CALL LLC2(VD1,VR,NEQ,ALPHA,1D0)
C
      CALL LL22(VR,NEQ,FR)
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
1     CONTINUE
C
      END
