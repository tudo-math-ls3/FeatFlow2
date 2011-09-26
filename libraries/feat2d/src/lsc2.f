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
* XLSC2                                                                *
*                                                                      *
* Purpose  Scaling of a single precision vector                        *
*                                                                      *
* Subroutines/functions called  LSC2                                   *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LX       I*4    Number of vector                                     *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
*                                                                      *
************************************************************************
*                                                                      *
* LSC1                                                                 *
*                                                                      *
* Purpose  Scaling of a single precision vector                        *
*                                                                      *
*                                                                      *
* Subroutines/functions called   LCL2                                  *
* BLAS                           SSCAL                                 *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* VX       R*4    Vector                                               *
* NX       I*4    Length                                               *
* A        R*8    Scaling factor                                       *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* VX       R*4    Scaled vector                                        *
*                                                                      *
************************************************************************
C
      SUBROUTINE XLSC2(LX,NX,A)
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
C
      IF (ICHECK.GE.998) CALL OTRC('XLSC2 ','01/02/89')
      IER=0
C
      IF (ICHECK.GE.3) THEN
       CALL ZTYPE(LX,ITYPE)
       CALL ZLEN(LX,ILEN)
       IF (NX.GT.ILEN) THEN
        CALL WERR(-121,'XLSC2 ')
        GOTO 99999
       ELSE IF (ITYPE.NE.2) THEN
        CALL WERR(-170,'XLSC2 ')
        GOTO 99999
       ENDIF
      ENDIF
C
      CALL LSC2(VWORK(L(LX)),NX,A)
C
99999 END
C
C
C
C
      SUBROUTINE LSC2(VX,NX,A)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LSC2  ','01/05/89')
C
      IF (A.NE.1D0) THEN
       IF (A.EQ.0D0) THEN
        CALL LCL2(VX,NX)
       ELSE IF (A.EQ.-1D0) THEN
        DO 10 IX=1,NX
10      VX(IX)=-VX(IX)
       ELSE
        CALL SSCAL(NX,REAL(A),VX,1)
       ENDIF
      ENDIF
C
      END
