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
* XM010                                                                *
*                                                                      *
* Purpose  Allocate Workspace vector on DWORK                          *
*          Call M010                                                   *
*                                                                      *
* Subroutines/functions called   ZNEW, ZDISP, ZCPY, ZLEN, M010         *
*                                                                      *
* Version from  04/12/91                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LX       I*4    Number of starting vector, length >= KNEQ(NLMAX)     *
* LB       I*4    Number of right hand side, length >= KNEQ(NLMAX)     *
* YDAX     SUBR                                                        *
* YPROL    SUBR                                                        *
* YREST    SUBR   External subroutines                                 *
* YPRSM    SUBR                                                        *  
* YPOSM    SUBR                                                        *
* YEX      SUBR                                                        *
* YBC      SUBR                                                        *
* IDISP    I*4    =1: On return, vectors assigned to numbers           * 
*                     LX  and  LB are truncated to length  KNEQ(NLMAX) *
* Hint:   Use  IDISP=0  if  M010  is called several times              *
*         with the same parameters and  NNWORK  is large enough        *
*                                                                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* ITE      I*4    Number of iterations                                 *
* IER      I*4    Error indicator                                      *
*                                                                      *
************************************************************************
C
      SUBROUTINE XM010(LX,LB,KNEQ,NIT,ITE,EPS,
     *                 YDAX,YPROL,YREST,YPRSM,YPOSM,YEX,YBC,IDISP,IREL)
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B) 
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNLEV=9,NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION KNEQ(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL YDAX,YPROL,YREST,YPRSM,YPOSM,YEX,YBC
      SAVE /ERRCTL/,/CHAR/,/MGPAR/
C
      SUB='XM010 '
      IF (ICHECK.GE.997) CALL OTRC('XM010 ','04/12/91')
      IER=0
C
      CALL ZLEN(LX,ILENX)
      CALL ZLEN(LB,ILENB)
      IF (ICHECK.GT.0) THEN
       IF (ICYCLE.LT.0) THEN
        CALL WERR(-181,'XM010 ')
        GOTO 99999
       ENDIF
       IF (NLMAX.LT.NLMIN.OR.NLMIN.LE.0.OR.NLMAX.GT.NNLEV) THEN
        CALL WERR(-182,'XM010 ')
        GOTO 99999
       ENDIF
       IF (ILENX.LT.KNEQ(NLMAX).OR.ILENB.LT.KNEQ(NLMAX)) THEN
        CALL WERR(-121,'XM010 ')
        GOTO 99999
       ENDIF
       CALL ZTYPE(LX,ITYPE1)
       CALL ZTYPE(LB,ITYPE2)
       IF (ITYPE1.NE.1.OR.ITYPE2.NE.1) THEN
        CALL WERR(-170,'XM010 ')
        GOTO 99999
       ENDIF
      ENDIF
C
      IREQ=0
      DO 1 ILEV=NLMIN,NLMAX
1     IREQ=IREQ+KNEQ(ILEV)
      IF (ILENX.LT.IREQ) THEN
       CALL ZNEW(IREQ,-1,LX1,'DX    ')
       IF (IER.NE.0) GOTO 99999
       CALL ZCPY(LX,'DXOLD ',LX1,'DX    ')
       CALL ZDISP(0,LX,'DXOLD ')
       IF (IER.NE.0) GOTO 99999
       LX=LX1
      ENDIF
      IF (ILENB.LT.IREQ) THEN
       CALL ZNEW(IREQ,-1,LB1,'DB    ')
       IF (IER.NE.0) GOTO 99999
       CALL ZCPY(LB,'DBOLD ',LB1,'DB    ')
       CALL ZDISP(0,LB,'DBOLD ')
       IF (IER.NE.0) GOTO 99999
       LB=LB1
      ENDIF
      CALL ZNEW(IREQ,-1,LD,'DD    ')
      IF (IER.NE.0) GOTO 99999
C
      CALL ZNEW(NLMAX,-3,LOFFX,'KOFFX ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NLMAX,-3,LOFFB,'KOFFB ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NLMAX,-3,LOFFD,'KOFFD ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NLMAX,-3,LIT,'KIT   ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NLMAX,-3,LIT0,'KIT0  ')
      IF (IER.NE.0) GOTO 99999
C
      LLX=L(LX)-1
      LLB=L(LB)-1
      LLD=L(LD)-1
      DO 2 ILEV=NLMAX,NLMIN,-1
      KWORK(L(LOFFX)+ILEV-1)=LLX
      KWORK(L(LOFFB)+ILEV-1)=LLB
      KWORK(L(LOFFD)+ILEV-1)=LLD
      LLX=LLX+KNEQ(ILEV)
      LLB=LLB+KNEQ(ILEV)
2     LLD=LLD+KNEQ(ILEV)
C
      CALL M010(DWORK(1),DWORK(1),DWORK(1),
     *          KWORK(L(LOFFX)),KWORK(L(LOFFB)),KWORK(L(LOFFD)),
     *          KNEQ,NIT,ITE,EPS,
     *          YDAX,YPROL,YREST,YPRSM,YPOSM,YEX,YBC,
     *          KWORK(L(LIT0)),KWORK(L(LIT)),IREL)
C
      IER1=IER
      CALL ZDISP(0,LIT0,'KIT0  ')
      IF (IER.NE.0) GOTO 99999
      CALL ZDISP(0,LIT,'KIT   ')
      IF (IER.NE.0) GOTO 99999
      CALL ZDISP(0,LOFFD,'KOFFD ')
      IF (IER.NE.0) GOTO 99999
      CALL ZDISP(0,LOFFB,'KOFFB ')
      IF (IER.NE.0) GOTO 99999
      CALL ZDISP(0,LOFFX,'KOFFX ')
      IF (IER.NE.0) GOTO 99999
      CALL ZDISP(0,LD,'DD    ')
      IF (IER.NE.0) GOTO 99999
      IF (IDISP.EQ.1) THEN
       CALL ZDISP(KNEQ(NLMAX),LB,'DB    ')
       CALL ZDISP(KNEQ(NLMAX),LX,'DX    ')
      ENDIF
      IER=IER1
C
99999 END
