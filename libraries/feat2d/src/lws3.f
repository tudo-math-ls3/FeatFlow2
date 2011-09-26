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
* XLWS3n                                                               *
*                                                                      *
* Purpose  Call LWS3n                                                  *
*          Matrix stored in technique  n  (see Reference Manual)       *
*          Single/double precision version                             *
*                                                                      *
* Subroutines/functions called   ZTYPE, ZLEN, LWS3n                    *
*                                                                      *
* Version from  11/27/89                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LA       I*4    Number of matrix stored in technique  n              *
* NEQ      I*4    Number of equations (length of DX, DY)               *
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
* LWS3n                                                                *
*                                                                      *
* Purpose  Weighted scalar product  (VA*DX,DY)                         *
*          Matrix stored in technique  n  (see Reference Manual)       *
*          Single/double precision version                             *
*                                                                      *
* Subroutines/functions called   None                                  *
*                                                                      *
* Version from  11/27/89                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* VA       R*4    Matrix stored in technique  n                        *
* NEQ      I*4    Number of equations (length of DX, DY)               *
* DX,DY    R*8    Input vectors                                        *
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
      SUBROUTINE XLWS33(LA,LDIA,LDIAS,NDIA,NEQ,LX,LY,SP)
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
      IF (ICHECK.GE.998) CALL OTRC('XLWS33','11/19/90')
C
      IF (ICHECK.GE.3) THEN
       CALL ZTYPE(LX,ITYPE1)
       CALL ZTYPE(LY,ITYPE2)
       CALL ZLEN(LX,ILEN1)
       CALL ZLEN(LY,ILEN2)
       IF (NEQ.GT.ILEN1.OR.NEQ.GT.ILEN2) THEN
        CALL WERR(-121,'XLWS33')
        GOTO 99999
       ELSE IF (ITYPE1.NE.1.OR.ITYPE2.NE.1) THEN
        CALL WERR(-170,'XLWS33')
        GOTO 99999
       ENDIF
      ENDIF
C
      CALL LWS33(VWORK(L(LA)),KWORK(L(LDIA)),KWORK(L(LDIAS)),NDIA,NEQ,
     *           DWORK(L(LX)),DWORK(L(LY)),SP)
C
99999 END
C
C
C
      SUBROUTINE LWS33(VA,KDIA,KDIAS,NDIA,NEQ,DX,DY,SP)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KDIA(*),KDIAS(*),DX(*),DY(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LWS33 ','11/19/90')
C
      SP=0D0
      DO 1 IEQ=1,NEQ
1     SP=SP+DX(IEQ)*DBLE(VA(IEQ))*DY(IEQ)
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
3     SP=SP+DX(IROW)*DBLE(VA(IROW+J0))*DY(IROW+J1)
2     CONTINUE
C
      END
C
C
C
      SUBROUTINE XLWS34(LA,LDIA,LDIAS,NDIA,NEQ,LX,LY,SP)
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
      IF (ICHECK.GE.998) CALL OTRC('XLWS34','01/10/91')
C
      IF (ICHECK.GE.3) THEN
       CALL ZTYPE(LX,ITYPE1)
       CALL ZTYPE(LY,ITYPE2)
       CALL ZLEN(LX,ILEN1)
       CALL ZLEN(LY,ILEN2)
       IF (NEQ.GT.ILEN1.OR.NEQ.GT.ILEN2) THEN
        CALL WERR(-121,'XLWS34')
        GOTO 99999
       ELSE IF (ITYPE1.NE.1.OR.ITYPE2.NE.1) THEN
        CALL WERR(-170,'XLWS34')
        GOTO 99999
       ENDIF
      ENDIF
C
      CALL LWS34(VWORK(L(LA)),KWORK(L(LDIA)),KWORK(L(LDIAS)),NDIA,NEQ,
     *           DWORK(L(LX)),DWORK(L(LY)),SP)
C
99999 END
C
C
C
      SUBROUTINE LWS34(VA,KDIA,KDIAS,NDIA,NEQ,DX,DY,SP)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KDIA(*),KDIAS(*),DX(*),DY(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LWS34 ','01/10/91')
C
      SP=0D0
      DO 1 IEQ=1,NEQ
1     SP=SP+DX(IEQ)*DBLE(VA(IEQ))*DY(IEQ)
      DO 2 IDIA=2,NDIA
      J0=KDIAS(IDIA)-1
      J1= KDIA(IDIA)
      DO 3 IROW=1,NEQ-J1
3     SP=SP+DX(IROW)*DBLE(VA(IROW+J0))*DY(IROW+J1)
      DO 4 IROW=1+J1,NEQ
4     SP=SP+DX(IROW)*DBLE(VA(IROW+J0-J1))*DY(IROW-J1)
2     CONTINUE
C
      END
C
C
C
      SUBROUTINE XLWS37(LA,LCOL,LLD,NEQ,LX,LY,SP)
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
      IF (ICHECK.GE.998) CALL OTRC('XLWS37','11/27/89')
C
      IF (ICHECK.GE.3) THEN
       CALL ZTYPE(LX,ITYPE1)
       CALL ZTYPE(LY,ITYPE2)
       CALL ZLEN(LX,ILEN1)
       CALL ZLEN(LY,ILEN2)
       IF (NEQ.GT.ILEN1.OR.NEQ.GT.ILEN2) THEN
        CALL WERR(-121,'XLWS37')
        GOTO 99999
       ELSE IF (ITYPE1.NE.1.OR.ITYPE2.NE.1) THEN
        CALL WERR(-170,'XLWS37')
        GOTO 99999
       ENDIF
      ENDIF
C
      CALL LWS37(VWORK(L(LA)),KWORK(L(LCOL)),KWORK(L(LLD)),NEQ,
     *           DWORK(L(LX)),DWORK(L(LY)),SP)
C
99999 END
C
C
C
      SUBROUTINE LWS37(VA,KCOL,KLD,NEQ,DX,DY,SP)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),DX(*),DY(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LWS37 ','11/27/89')
C
      SP=0D0
      DO 1 IEQ=1,NEQ
      AUX=0D0
      DO 2 ICOL=KLD(IEQ),KLD(IEQ+1)-1
2     AUX=AUX+DBLE(VA(ICOL))*DX(KCOL(ICOL))
      SP=SP+AUX*DY(IEQ)
1     CONTINUE
C
      END
C
C
C
      SUBROUTINE XLWS38(LA,LCOL,LLD,NEQ,LX,LY,SP)
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
      IF (ICHECK.GE.998) CALL OTRC('XLWS38','01/10/91')
C
      IF (ICHECK.GE.3) THEN
       CALL ZTYPE(LX,ITYPE1)
       CALL ZTYPE(LY,ITYPE2)
       CALL ZLEN(LX,ILEN1)
       CALL ZLEN(LY,ILEN2)
       IF (NEQ.GT.ILEN1.OR.NEQ.GT.ILEN2) THEN
        CALL WERR(-121,'XLWS38')
        GOTO 99999
       ELSE IF (ITYPE1.NE.1.OR.ITYPE2.NE.1) THEN
        CALL WERR(-170,'XLWS38')
        GOTO 99999
       ENDIF
      ENDIF
C
      CALL LWS38(VWORK(L(LA)),KWORK(L(LCOL)),KWORK(L(LLD)),NEQ,
     *           DWORK(L(LX)),DWORK(L(LY)),SP)
C
99999 END
C
C
C
      SUBROUTINE LWS38(VA,KCOL,KLD,NEQ,DX,DY,SP)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),DX(*),DY(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LWS38 ','01/10/91')
C
      SP=0D0
      DO 1 IEQ=1,NEQ
1     SP=SP+DX(IEQ)*DBLE(VA(KLD(IEQ)))*DY(IEQ)
      DO 2 IEQ=1,NEQ-1
      AUX=0D0
      DO 3 ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
      JCOL=KCOL(ICOL)
      AUX=AUX+DBLE(VA(ICOL))*DX(JCOL)
      SP=SP+DY(JCOL)*DBLE(VA(ICOL))*DX(IEQ)
3     CONTINUE
      SP=SP+AUX*DY(IEQ)
2     CONTINUE
C
      END
C
C
C
      SUBROUTINE XLWS39(LA,LCOL,LLD,NEQ,LX,LY,SP)
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
      IF (ICHECK.GE.998) CALL OTRC('XLWS39','01/10/91')
C
      IF (ICHECK.GE.3) THEN
       CALL ZTYPE(LX,ITYPE1)
       CALL ZTYPE(LY,ITYPE2)
       CALL ZLEN(LX,ILEN1)
       CALL ZLEN(LY,ILEN2)
       IF (NEQ.GT.ILEN1.OR.NEQ.GT.ILEN2) THEN
        CALL WERR(-121,'XLWS39')
        GOTO 99999
       ELSE IF (ITYPE1.NE.1.OR.ITYPE2.NE.1) THEN
        CALL WERR(-170,'XLWS39')
        GOTO 99999
       ENDIF
      ENDIF
C
      CALL LWS39(VWORK(L(LA)),KWORK(L(LCOL)),KWORK(L(LLD)),NEQ,
     *           DWORK(L(LX)),DWORK(L(LY)),SP)
C
99999 END
C
C
C
      SUBROUTINE LWS39(VA,KCOL,KLD,NEQ,DX,DY,SP)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),DX(*),DY(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LWS39 ','01/10/91')
C
      SP=0D0
      DO 1 IEQ=1,NEQ
      AUX=0D0
      DO 2 ICOL=KLD(IEQ),KLD(IEQ+1)-1
2     AUX=AUX+DBLE(VA(ICOL))*DX(KCOL(ICOL))
      SP=SP+AUX*DY(IEQ)
1     CONTINUE
C
      END
C
C
C
      SUBROUTINE XLWS3A(LA,LCOL,LLD,LOP,NEQ,LX,LY,SP)
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
      IF (ICHECK.GE.998) CALL OTRC('XLWS3A','11/27/89')
C
      IF (ICHECK.GE.3) THEN
       CALL ZTYPE(LX,ITYPE1)
       CALL ZTYPE(LY,ITYPE2)
       CALL ZLEN(LX,ILEN1)
       CALL ZLEN(LY,ILEN2)
       IF (NEQ.GT.ILEN1.OR.NEQ.GT.ILEN2) THEN
        CALL WERR(-121,'XLWS3A')
        GOTO 99999
       ELSE IF (ITYPE1.NE.1.OR.ITYPE2.NE.1) THEN
        CALL WERR(-170,'XLWS3A')
        GOTO 99999
       ENDIF
      ENDIF
C
      CALL LWS3A(VWORK(L(LA)),KWORK(L(LCOL)),KWORK(L(LLD)),
     *           KWORK(L(LOP)),NEQ,DWORK(L(LX)),DWORK(L(LY)),SP)
C
99999 END
C
C
C
      SUBROUTINE LWS3A(VA,KCOL,KLD,KOP,NEQ,DX,DY,SP)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),KOP(*),DX(*),DY(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LWS3A ','11/27/89')
C
      SP=0D0
      DO 1 IEQ=1,NEQ
      JOP=KOP(IEQ)
      AUX=0D0
      DO 2 ICOL=KLD(JOP),KLD(JOP+1)-1
2     AUX=AUX+DBLE(VA(ICOL))*DX(KCOL(ICOL)+IEQ)
      SP=SP+AUX*DY(IEQ)
1     CONTINUE
C
      END
