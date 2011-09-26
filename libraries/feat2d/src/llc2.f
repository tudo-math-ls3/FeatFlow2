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
* XLLC2                                                                *
*                                                                      *
* Purpose  Linear combination of two single precision vectors          *
*          Call  LLC2                                                  *
*                                                                      *
* Subroutines/functions called   ZLEN, ZTYPE, LLC2                     *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LX,LY    I*4    Number of input vectors                              *
* NX       I*4    Length of vectors                                    *
* A1,A2    R*8    VY := A1*VX + A2*VY                                  *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* IER      I*4    Error indicator                                      *
*                 -121 At least one of the vectors too small           *
*                      or invalid number                               *
*                 -170 Data type not single or double precision        *
*                                                                      *
************************************************************************
*                                                                      *
* LLC2                                                                 *
*                                                                      *
* Purpose  Linear combination of two double precision vectors          *
*                                                                      *
* Subroutines/functions called   None                                  *
* BLAS                           SAXPY,SCOPY,SSCAL                     *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* VX,VY    R*4    Input vectors                                        *
* NX       I*4    Length of vectors                                    *
* A1,A2    R*8    VY := A1*VX + A2*VY                                  *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* VY       R*4    Resulting vector                                     *
*                                                                      *
************************************************************************
C
      SUBROUTINE XLLC2(LX,LY,NX,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1))
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='XLLC2'
      IF (ICHECK.GE.998) CALL OTRC('XLLC2 ','01/02/89')
C
      IF (ICHECK.GE.3) THEN
       CALL ZLEN(LX,ILEN1)
       CALL ZLEN(LY,ILEN2)
       CALL ZTYPE(LX,ITYPE1)
       CALL ZTYPE(LY,ITYPE2)
       IF (NX.GT.ILEN1.OR.NX.GT.ILEN2) THEN
        CALL WERR(-121,'XLLC2 ')
        GOTO 99999
       ELSE IF (ITYPE1.NE.2.OR.ITYPE2.NE.2) THEN
        CALL WERR (-170,'XLLC2 ')
        GOTO 99999
       ENDIF
      ENDIF
C
      CALL LLC2(VWORK(L(LX)),VWORK(L(LY)),NX,A1,A2)
C
99999 END
C
C
C
      SUBROUTINE LLC2(VX,VY,NX,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VX(*),VY(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LLC2  ','01/02/89')
C
      IF (A2.EQ.0D0) THEN
       CALL SCOPY(NX,VX,1,VY,1)
       IF (A1.NE.1D0) CALL SSCAL(NX,REAL(A1),VY,1)
      ELSE IF (A2.EQ.1D0) THEN
       CALL SAXPY(NX,REAL(A1),VX,1,VY,1)
      ELSE
       A3=A1/A2
       CALL SAXPY(NX,REAL(A3),VX,1,VY,1)
       CALL SSCAL(NX,REAL(A2),VY,1)
      ENDIF
C
      END
