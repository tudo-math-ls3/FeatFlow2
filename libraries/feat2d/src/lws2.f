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
* XLWS2n                                                               *
*                                                                      *
* Purpose  Call  LWS2n                                                 *
*          Matrix stored in technique  n  (see Reference Manual)       *
*          Single/single precision version                             *
*                                                                      *
* Subroutines/functions called   ZTYPE, ZLEN, LWS2n                    *
*                                                                      *
* Version from  11/27/89                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LA       I*4    Number of matrix stored in technique  n              *
* NEQ      I*4    Number of equations (length of VX, VY)               *
* LX,LY    I*4    Numbers of vectors                                   *
* LLD      I*4    Numbers of pointer vectors corresponding to the      *
* LCOL     I*4    storage technique                                    *
* LOP      I*4                                                         *
*                                                                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* IER      I*4    Error indicator                                      *
*                                                                      *
************************************************************************
*                                                                      *
* LWS2n                                                                *
*                                                                      *
* Purpose  Weighted scalar product  (VA*VX,VY)                         *
*          Matrix stored in technique  n  (see Reference Manual)       *
*          Single/single precision version                             *
*                                                                      *
* Subroutines/functions called   None                                  *
*                                                                      *
* Version from  11/27/89                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* VA       R*4    Matrix stored in technique  n                        *
* NEQ      I*4    Number of equations (length of VX, VY)               *
* VX,VY    R*4    Input vectors                                        *
* KLD      I*4    Pointer vectors corresponding to the                 *
* KCOL     I*4    storage technique                                    *
* KOP      I*4                                                         *
*                                                                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* SP       R*8    Scalar product                                       *
*                                                                      *
************************************************************************
C
      SUBROUTINE XLWS23(LA,LDIA,LDIAS,NDIA,NEQ,LX,LY,SP)
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /OUTPUT/,/CHAR/,/ERRCTL/
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      IF (ICHECK.GE.998) CALL OTRC('XLWS23','11/19/90')
C
      IF (ICHECK.GE.3) THEN
       CALL ZTYPE(LX,ITYPE1)
       CALL ZTYPE(LY,ITYPE2)
       CALL ZLEN(LX,ILEN1)
       CALL ZLEN(LY,ILEN2)
       IF (NEQ.GT.ILEN1.OR.NEQ.GT.ILEN2) THEN
        CALL WERR(-121,'XLWS23')
        GOTO 99999
       ELSE IF (ITYPE1.NE.2.OR.ITYPE2.NE.2) THEN
        CALL WERR(-170,'XLWS23')
        GOTO 99999
       ENDIF
      ENDIF
C
      CALL LWS23(VWORK(L(LA)),KWORK(L(LDIA)),KWORK(L(LDIAS)),NDIA,NEQ,
     *           VWORK(L(LX)),VWORK(L(LY)),SP)
C
99999 END
C
C
C
      SUBROUTINE LWS23(VA,KDIA,KDIAS,NDIA,NEQ,VX,VY,SP)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KDIA(*),KDIAS(*),VX(*),VY(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LWS23 ','11/19/90')
C
      SP=0D0
      DO 1 IEQ=1,NEQ
1     SP=SP+DBLE(VX(IEQ))*DBLE(VA(IEQ))*DBLE(VY(IEQ))
      DO 2 IDIA=2,NDIA
      J1=KDIA(IDIA)
      IF (J1.GT.0) THEN
       I1=1
       I2=NEQ-J1
      ELSE
       I1=1-J1
       I2=NEQ
      ENDIF
      J0=KDIAS(IDIA)-I1
      DO 3 IROW=I1,I2
3     SP=SP+DBLE(VX(IROW))*DBLE(VA(IROW+J0))*DBLE(VY(IROW+J1))
2     CONTINUE
C
      END
C
C
C
      SUBROUTINE XLWS24(LA,LDIA,LDIAS,NDIA,NEQ,LX,LY,SP)
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /OUTPUT/,/CHAR/,/ERRCTL/
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      IF (ICHECK.GE.998) CALL OTRC('XLWS24','01/10/91')
C
      IF (ICHECK.GE.3) THEN
       CALL ZTYPE(LX,ITYPE1)
       CALL ZTYPE(LY,ITYPE2)
       CALL ZLEN(LX,ILEN1)
       CALL ZLEN(LY,ILEN2)
       IF (NEQ.GT.ILEN1.OR.NEQ.GT.ILEN2) THEN
        CALL WERR(-121,'XLWS24')
        GOTO 99999
       ELSE IF (ITYPE1.NE.2.OR.ITYPE2.NE.2) THEN
        CALL WERR(-170,'XLWS24')
        GOTO 99999
       ENDIF
      ENDIF
C
      CALL LWS24(VWORK(L(LA)),KWORK(L(LDIA)),KWORK(L(LDIAS)),NDIA,NEQ,
     *           VWORK(L(LX)),VWORK(L(LY)),SP)
C
99999 END
C
C
C
      SUBROUTINE LWS24(VA,KDIA,KDIAS,NDIA,NEQ,VX,VY,SP)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KDIA(*),KDIAS(*),VX(*),VY(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LWS24 ','01/10/91')
C
      SP=0D0
      DO 1 IEQ=1,NEQ
1     SP=SP+DBLE(VX(IEQ))*DBLE(VA(IEQ))*DBLE(VY(IEQ))
      DO 2 IDIA=2,NDIA
      J0=KDIAS(IDIA)-1
      J1= KDIA(IDIA)
      DO 3 IROW=1,NEQ-J1
3     SP=SP+DBLE(VX(IROW))*DBLE(VA(IROW+J0))*DBLE(VY(IROW+J1))
      DO 4 IROW=1+J1,NEQ
4     SP=SP+DBLE(VX(IROW))*DBLE(VA(IROW+J0-J1))*DBLE(VY(IROW-J1))
2     CONTINUE
C
      END
C
C
C
      SUBROUTINE XLWS27(LA,LCOL,LLD,NEQ,LX,LY,SP)
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /OUTPUT/,/CHAR/,/ERRCTL/
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      IF (ICHECK.GE.998) CALL OTRC('XLWS27','11/27/89')
C
      IF (ICHECK.GE.3) THEN
       CALL ZTYPE(LX,ITYPE1)
       CALL ZTYPE(LY,ITYPE2)
       CALL ZLEN(LX,ILEN1)
       CALL ZLEN(LY,ILEN2)
       IF (NEQ.GT.ILEN1.OR.NEQ.GT.ILEN2) THEN
        CALL WERR(-121,'XLWS27')
        GOTO 99999
       ELSE IF (ITYPE1.NE.2.OR.ITYPE2.NE.2) THEN
        CALL WERR(-170,'XLWS27')
        GOTO 99999
       ENDIF
      ENDIF
C
      CALL LWS27(VWORK(L(LA)),KWORK(L(LCOL)),KWORK(L(LLD)),NEQ,
     *           VWORK(L(LX)),VWORK(L(LY)),SP)
C
99999 END
C
C
C
      SUBROUTINE LWS27(VA,KCOL,KLD,NEQ,VX,VY,SP)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),VX(*),VY(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LWS27 ','11/27/89')
C
      SP=0D0
      DO 1 IEQ=1,NEQ
      AUX=0D0
      DO 2 ICOL=KLD(IEQ),KLD(IEQ+1)-1
2     AUX=AUX+DBLE(VA(ICOL))*DBLE(VX(KCOL(ICOL)))
      SP=SP+AUX*DBLE(VY(IEQ))
1     CONTINUE
C
      END
C
C
C
      SUBROUTINE XLWS28(LA,LCOL,LLD,NEQ,LX,LY,SP)
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /OUTPUT/,/CHAR/,/ERRCTL/
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      IF (ICHECK.GE.998) CALL OTRC('XLWS28','01/10/91')
C
      IF (ICHECK.GE.3) THEN
       CALL ZTYPE(LX,ITYPE1)
       CALL ZTYPE(LY,ITYPE2)
       CALL ZLEN(LX,ILEN1)
       CALL ZLEN(LY,ILEN2)
       IF (NEQ.GT.ILEN1.OR.NEQ.GT.ILEN2) THEN
        CALL WERR(-121,'XLWS28')
        GOTO 99999
       ELSE IF (ITYPE1.NE.2.OR.ITYPE2.NE.2) THEN
        CALL WERR(-170,'XLWS28')
        GOTO 99999
       ENDIF
      ENDIF
C
      CALL LWS28(VWORK(L(LA)),KWORK(L(LCOL)),KWORK(L(LLD)),NEQ,
     *           VWORK(L(LX)),VWORK(L(LY)),SP)
C
99999 END
C
C
C
      SUBROUTINE LWS28(VA,KCOL,KLD,NEQ,VX,VY,SP)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),VX(*),VY(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LWS28 ','01/10/91')
C
      SP=0D0
      DO 1 IEQ=1,NEQ
1     SP=SP+DBLE(VX(IEQ))*DBLE(VA(KLD(IEQ)))*DBLE(VY(IEQ))
      DO 2 IEQ=1,NEQ-1
      AUX=0D0
      DO 3 ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
      JCOL=KCOL(ICOL)
      AUX=AUX+DBLE(VA(ICOL))*DBLE(VX(JCOL))
      SP=SP+DBLE(VY(JCOL))*DBLE(VA(ICOL))*DBLE(VX(IEQ))
3     CONTINUE
      SP=SP+AUX*DBLE(VY(IEQ))
2     CONTINUE
C
      END
C
C
C
      SUBROUTINE XLWS29(LA,LCOL,LLD,NEQ,LX,LY,SP)
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /OUTPUT/,/CHAR/,/ERRCTL/
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      IF (ICHECK.GE.998) CALL OTRC('XLWS29','01/10/91')
C
      IF (ICHECK.GE.3) THEN
       CALL ZTYPE(LX,ITYPE1)
       CALL ZTYPE(LY,ITYPE2)
       CALL ZLEN(LX,ILEN1)
       CALL ZLEN(LY,ILEN2)
       IF (NEQ.GT.ILEN1.OR.NEQ.GT.ILEN2) THEN
        CALL WERR(-121,'XLWS29')
        GOTO 99999
       ELSE IF (ITYPE1.NE.2.OR.ITYPE2.NE.2) THEN
        CALL WERR(-170,'XLWS29')
        GOTO 99999
       ENDIF
      ENDIF
C
      CALL LWS29(VWORK(L(LA)),KWORK(L(LCOL)),KWORK(L(LLD)),NEQ,
     *           VWORK(L(LX)),VWORK(L(LY)),SP)
C
99999 END
C
C
C
      SUBROUTINE LWS29(VA,KCOL,KLD,NEQ,VX,VY,SP)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),VX(*),VY(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LWS29 ','01/10/91')
C
      SP=0D0
      DO 1 IEQ=1,NEQ
      AUX=0D0
      DO 2 ICOL=KLD(IEQ),KLD(IEQ+1)-1
2     AUX=AUX+DBLE(VA(ICOL))*DBLE(VX(KCOL(ICOL)))
      SP=SP+AUX*DBLE(VY(IEQ))
1     CONTINUE
C
      END
C
C
C
      SUBROUTINE XLWS2A(LA,LCOL,LLD,LOP,NEQ,LX,LY,SP)
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /OUTPUT/,/CHAR/,/ERRCTL/
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      IF (ICHECK.GE.998) CALL OTRC('XLWS2A','11/27/89')
C
      IF (ICHECK.GE.3) THEN
       CALL ZTYPE(LX,ITYPE1)
       CALL ZTYPE(LY,ITYPE2)
       CALL ZLEN(LX,ILEN1)
       CALL ZLEN(LY,ILEN2)
       IF (NEQ.GT.ILEN1.OR.NEQ.GT.ILEN2) THEN
        CALL WERR(-121,'XLWS2A')
        GOTO 99999
       ELSE IF (ITYPE1.NE.2.OR.ITYPE2.NE.2) THEN
        CALL WERR(-170,'XLWS2A')
        GOTO 99999
       ENDIF
      ENDIF
C
      CALL LWS2A(VWORK(L(LA)),KWORK(L(LCOL)),KWORK(L(LLD)),
     *           KWORK(L(LOP)),NEQ,VWORK(L(LX)),VWORK(L(LY)),SP)
C
99999 END
C
C
C
      SUBROUTINE LWS2A(VA,KCOL,KLD,KOP,NEQ,VX,VY,SP)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),KOP(*),VX(*),VY(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LWS2A ','11/27/89')
C
      SP=0D0
      DO 1 IEQ=1,NEQ
      JOP=KOP(IEQ)
      AUX=0D0
      DO 2 ICOL=KLD(JOP),KLD(JOP+1)-1
2     AUX=AUX+DBLE(VA(ICOL))*DBLE(VX(KCOL(ICOL)+IEQ))
      SP=SP+AUX*DBLE(VY(IEQ))
1     CONTINUE
C
      END
	
