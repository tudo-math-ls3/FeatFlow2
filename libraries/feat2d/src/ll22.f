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
* XLL22                                                                *
*                                                                      *
* Purpose  l2-norm of a single precision vector                        *
*                                                                      *
* Subroutines/functions called  LL22                                   *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LX       I*4    Number of vector                                     *
* NX       I*4    length                                               *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
*                                                                      *
************************************************************************
*                                                                      *
* LL22                                                                 *
*                                                                      *
* Purpose  l2-norm of a single precision vector                        *
*                                                                      *
*                                                                      *
* Subroutines/functions called   LSP2                                  *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* VX       R*4    Vector                                               *
* NX       I*4    Length                                               *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* XNORM    R*8    l2-norm                                              *
*                                                                      *
************************************************************************
C
      SUBROUTINE XLL22(LX,NX,XNORM)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EQUIVALENCE (DWORK(1),VWORK(1))
      SAVE /ERRCTL/,/CHAR/
      IF (ICHECK.GE.998) CALL OTRC('XLL22 ','01/02/89')
      IER=0
C
      IF (ICHECK.GE.3) THEN
       CALL ZTYPE(LX,ITYPE)
       CALL ZLEN(LX,ILEN)
       IF (NX.GT.ILEN) THEN
        CALL WERR(-121,'XLL22 ')
        GOTO 99999
       ELSE IF (ITYPE.NE.2) THEN
        CALL WERR(-170,'XLL22 ')
        GOTO 99999
       ENDIF
      ENDIF
C
      CALL LL22(VWORK(L(LX)),NX,XNORM)
C
99999 END
C
C
C
      SUBROUTINE LL22(VX,NX,XNORM)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('LL22  ','01/02/89')
C
      CALL LSP2(VX,VX,NX,XNORM)
      XNORM=SQRT(XNORM)
      END
