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
* IE220                                                                *
*                                                                      *
* Purpose  Smoothing using a preconditioned conjugate gradient method  *
*          Single precision version                                    *
*                                                                      *
* Subroutines/functions called  LSP2 , LLC2 , LL22                     *
*                                                                      *
* Version from  10/30/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* VX       R*4    Starting vector                                      *
* VB       R*4    Right hand side                                      *
* NEQ      I*4    Number of equations                                  *
* NIT      I*4    Number of smoothing steps                            *
* VAX0     SUBR   EXTERNAL Subroutine VAX0(VX,VAX,NEQ,A1,A2)           *
*                 Results  VAX := A1 * A * VX + A2 * VAX               *
* CG0C     SUBR   EXTERNAL Subroutine CG0C(VG,NEQ)                     *
*                 Results  VG := C**(-1) * VG  for the precondioning   *
*                 matrix C                                             *
* BNOCON   L*4    .TRUE.   No Preconditioning                          *
* VR,VD    R*4    Workspace vectors of length NEQ                      *
* VD1,VG   R*4    For BNOCON , VG must be replaced by VR               *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* VX       R*4    Smoothed vector                                      *
*                                                                      *
************************************************************************
C
      SUBROUTINE IE220(VX,VB,NEQ,NIT,VAX0,CG0C,BNOCON,VR,VD,VD1,VG)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VX(*),VB(*),VR(*),VG(*),VD(*),VD1(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /ERRCTL/,/CHAR/
C
      SUB='IE220'
      IF (ICHECK.GE.998) CALL OTRC('IE220 ','10/30/89')
C
C *** Initialization
      CALL VAX0(VX,VR,NEQ,1D0,0D0)
      CALL LLC2(VB,VR,NEQ,-1D0,1D0)
      CALL LL22(VR,NEQ,RES)
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
      DO 1 ITE=1,NIT
C
      CALL VAX0(VD,VD1,NEQ,1D0,0D0)
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
       CALL CG0C(VG,NEQ)
       CALL LSP2(VR,VG,NEQ,SIGMA1)
      ENDIF
C
      GAMMA=SIGMA1/SIGMA0
      SIGMA0=SIGMA1
      CALL LLC2(VG,VD,NEQ,-1D0,GAMMA)
C
1     CONTINUE
C
      END
