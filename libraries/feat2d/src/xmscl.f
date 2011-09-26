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
* XMSCL                                                                *
*                                                                      *
* Purpose  Make Clean multiple triangulations (multigrid version)      *
*                                                                      *
* Subroutines/functions called   ZDISP                                 *
*                                                                      *
* Version from  08/25/90                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* IER      I*4    Error indicator                                      *
*                 Set by ZDISP                                         *
*                                                                      *
************************************************************************
C
      SUBROUTINE XMSCL
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNLEV=9)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /MGTRD/,/MGTRA/,/ERRCTL/,/CHAR/
C
      SUB='XMSCL'
      IF (ICHECK.GE.997) CALL OTRC('XMSCL ','08/25/90')
      IER=0
C
      DO 1 ILEV=NNLEV,1,-1
C
      KNEL(ILEV) =0
      KNVT(ILEV) =0
      KNMT(ILEV) =0
      KNVEL(ILEV)=0
      KNVBD(ILEV)=0
C
      IF (KLMBDP(ILEV).NE.0) THEN
       CALL ZDISP(0,KLMBDP(ILEV),'DMBDP0')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLVBDP(ILEV).NE.0) THEN
       CALL ZDISP(0,KLVBDP(ILEV),'DVBDP0')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLBCT(ILEV).NE.0) THEN
       CALL ZDISP(0,KLBCT(ILEV),'KBCT0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLEBD(ILEV).NE.0) THEN
       CALL ZDISP(0,KLEBD(ILEV),'KEBD0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLVBD(ILEV).NE.0) THEN
       CALL ZDISP(0,KLVBD(ILEV),'KVBD0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLMM(ILEV).NE.0) THEN
       CALL ZDISP(0,KLMM(ILEV),'KMM0  ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLNPR(ILEV).NE.0) THEN
       CALL ZDISP(0,KLNPR(ILEV),'KNPR0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLMEL(ILEV).NE.0) THEN
       CALL ZDISP(0,KLMEL(ILEV),'KMEL0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLVEL(ILEV).NE.0) THEN
       CALL ZDISP(0,KLVEL(ILEV),'KVEL0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLADJ(ILEV).NE.0) THEN
       CALL ZDISP(0,KLADJ(ILEV),'KADJ0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLMID(ILEV).NE.0) THEN
       CALL ZDISP(0,KLMID(ILEV),'KMID0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLVERT(ILEV).NE.0) THEN
       CALL ZDISP(0,KLVERT(ILEV),'KVERT0')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLCMG(ILEV).NE.0) THEN
       CALL ZDISP(0,KLCMG(ILEV),'DCMG0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLCVG(ILEV).NE.0) THEN
       CALL ZDISP(0,KLCVG(ILEV),'DCVG0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
1     CONTINUE
C     
99999 END
      
      
