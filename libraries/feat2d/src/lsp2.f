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
* XLSP2                                                                *
*                                                                      *
* Purpose  Scalar product of two single precision vectors              *
*          Call  LSP2                                                  *
*                                                                      *
* Subroutines/functions called    ZTYPE, ZLEN, LSP2                    *
*                                                                      *
* Version from  06/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LX,LY    I*4    Numbers of vectors                                   *
* NX       I*4    Number of elements to be copied                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* IER      I*4    Error indicator                                      *
*                                                                      *
************************************************************************
*                                                                      *
* LSP2                                                                 *
*                                                                      *
* Purpose  Scalar product                                              *
*          Single precision version                                    *
*                                                                      *
* Subroutines/functions called   None                                  *
* BLAS                           DSDOT                                 *
*                                                                      *
* Version from  06/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* VX,VY    R*4    Input vectors                                        *
* NX       I*4    Number of elements                                   *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* SP       R*8    Scalar product                                       *
*                                                                      *
************************************************************************
C
      SUBROUTINE XLSP2(LX,LY,NX,SP)
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ERRCTL/ IER,ICHECK
      EQUIVALENCE (DWORK(1),VWORK(1))
      SAVE /OUTPUT/,/CHAR/,/ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('XLSP2 ','06/02/89')
      IER=0
C
      IF (ICHECK.GE.3) THEN
       CALL ZTYPE(LX,ITYPE1)
       CALL ZTYPE(LY,ITYPE2)
       CALL ZLEN(LX,ILEN1)
       CALL ZLEN(LY,ILEN2)
       IF (NX.GT.ILEN1.OR.NX.GT.ILEN2) THEN
        CALL WERR(-121,'XLSP2 ')
        GOTO 99999
       ELSE IF (ITYPE1.NE.2.OR.ITYPE2.NE.2) THEN
        CALL WERR(-170,'XLSP2 ')
        GOTO 99999
       ENDIF
      ENDIF
C
      CALL LSP2(VWORK(L(LX)),VWORK(L(LY)),NX,SP)
99999 END
C
C
C
      SUBROUTINE LSP2(VX,VY,NXY,SP)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VX(*),VY(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('LSP2  ','06/02/89')
C
C *** DSDOT uses double precision accumulation
      SP=DSDOT(NXY,VX,1,VY,1)
      END
