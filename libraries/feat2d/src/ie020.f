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
* IE020                                                                *
*                                                                      *
* Purpose  Solution of a linear system  A*X = B  using                 *
*          a preconditioned conjugate gradient method                  *
*          Single precision version                                    *
*                                                                      *
* Subroutines/functions called  , LLI2, LSP2 , LLC2 , LL22 , LCL2      *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* VX       R*4    Starting vector                                      *
* VB       R*4    Right hand side                                      *
* NEQ      I*4    Number of equations                                  *
* ITE      I*4    Minimum number of iterations (used for IREL=1 only)  *
* NIT      I*4    Maximum number of iterations                         *
* EPS      R*8    Desired precision                                    *
*                 IREL=0: Stop if !!RES!! < EPS                        *
*                 IREL=1: Stop if !!RES!!/!!RES0!! < EPS               *
*                         and a minimum of ITE iterations is performed *
* VAX0     SUBR   EXTERNAL Subroutine VAX0(VX,VAX,NEQ,A1,A2)           *
*                 Results  VAX := A1 * A * VX + A2 * VAX               *
* CG0C     SUBR   EXTERNAL Subroutine CG0C(VG,NEQ)                     *
*                 Results  VG := C**(-1) * VG  for the precondioning   *
*                 matrix C                                             *
* BNOCON   L*4    .TRUE.  No Preconditioning                           *
* VR,VD    R*4    Workspace vectors of length NEQ                      *
* VD1,VG   R*4    For BNOCON , VG must be replaced by VR               *
* IREL     I*4    See above                                            *
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
      SUBROUTINE IE020(VX,VB,NEQ,NIT,ITE,EPS,VAX0,CG0C,BNOCON,
     *                 VR,VD,VD1,VG,IREL)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VX(*),VB(*),VR(*),VG(*),VD(*),VD1(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
      DATA FR/0D0/
C
      SUB='IE020'
      IF (ICHECK.GE.997) CALL OTRC('IE020 ','01/02/89')
C
      BMSG2=M.GE.2.OR.MT.GE.2
      NIT0=MAX(ITE,0)
      BREL=IREL.EQ.1
C
      IF (ICHECK.GT.0) THEN
       CALL LLI2(VB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL2(VX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IE020 ')
        GOTO 99999
       ENDIF
      ENDIF
C
C *** Initialization
      CALL LL22(VB,NEQ,RBNORM)
      CALL VAX0(VX,VR,NEQ,1D0,0D0)
      CALL LLC2(VB,VR,NEQ,-1D0,1D0)
      CALL LL22(VR,NEQ,RES)
      IF (RES.LE.EPS.AND..NOT.BREL) THEN
       ITE=0
       FR=RES
       GOTO 200
      ENDIF
C
      IF (BNOCON) THEN
       SIGMA0=RES*RES
      ELSE
       CALL LCP2(VR,VG,NEQ)
       CALL CG0C(VG,NEQ)
       CALL LSP2(VR,VG,NEQ,SIGMA0)
      ENDIF
C
      CALL LLC2(VG,VD,NEQ,-1D0,0D0)
C
C *** Iterative correction
      DO 100 ITE=1,NIT
C
      CALL VAX0(VD,VD1,NEQ,1D0,0D0)
      CALL LSP2(VD,VD1,NEQ,ALPHA)
      ALPHA=SIGMA0/ALPHA
      CALL LLC2(VD,VX,NEQ,ALPHA,1D0)
      CALL LLC2(VD1,VR,NEQ,ALPHA,1D0)
C
      CALL LL22(VR,NEQ,FR)
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,FR
       CALL OMSG(73,'IE020 ')
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
       CALL LCP2(VR,VG,NEQ)
       CALL CG0C(VG,NEQ)
       CALL LSP2(VR,VG,NEQ,SIGMA1)
      ENDIF
C
      GAMMA=SIGMA1/SIGMA0
      SIGMA0=SIGMA1
      CALL LLC2(VG,VD,NEQ,-1D0,GAMMA)
100   CONTINUE
C
      WRITE (CPARAM,'(I15,2D25.16)') NIT,FR,RES
      CALL OMSG(71,'IE020 ')
      CALL OMSG(72,'IE020 ')
C
      IF (RES.GE.1D-70) THEN
       CAPPA=(FR/RES)**(1D0/NIT)
      ELSE
       CAPPA=0D0
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE020 ')
C
      IER=1
      GOTO 99999
C
200   IER=0
      IF (RES.GE.1D-70) RES=FR/RES
      WRITE (CPARAM,'(I15,2D25.16)') ITE,FR,RES
      CALL OMSG(72,'IE020 ')
C
      IF (ITE.EQ.0) THEN
       CAPPA=0D0
      ELSE
       CAPPA=RES**(1D0/ITE)
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE020 ')
C
99999 END
