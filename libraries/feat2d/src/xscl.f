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
* XSCL                                                                 *
*                                                                      *
* Purpose  Make Clean subdivision                                      *
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
      SUBROUTINE XSCL
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNLEV=9)
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /TRIAD/,/TRIAA/,/ERRCTL/,/CHAR/
C
      SUB='XSCL'
      IF (ICHECK.GE.997) CALL OTRC('XSCL  ','08/25/90')
      IER=0
C
      NEL =0
      NVT =0
      NMT =0
      NVEL=0
      NVBD=0
C
      IF (LMBDP.NE.0) THEN
       CALL ZDISP(0,LMBDP,'DMBDP ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LVBDP.NE.0) THEN
       CALL ZDISP(0,LVBDP,'DVBDP ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LBCT.NE.0) THEN
       CALL ZDISP(0,LBCT,'KBCT  ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LEBD.NE.0) THEN
       CALL ZDISP(0,LEBD,'KEBD  ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LVBD.NE.0) THEN
       CALL ZDISP(0,LVBD,'KVBD  ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LMM.NE.0) THEN
       CALL ZDISP(0,LMM,'KMM   ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LNPR.NE.0) THEN
       CALL ZDISP(0,LNPR,'KNPR  ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LMEL.NE.0) THEN
       CALL ZDISP(0,LMEL,'KMEL  ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LVEL.NE.0) THEN
       CALL ZDISP(0,LVEL,'KVEL  ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LADJ.NE.0) THEN
       CALL ZDISP(0,LADJ,'KADJ  ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LMID.NE.0) THEN
       CALL ZDISP(0,LMID,'KMID  ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LVERT.NE.0) THEN
       CALL ZDISP(0,LVERT,'KVERT ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LCORMG.NE.0) THEN
       CALL ZDISP(0,LCORMG,'DCORMG')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LCORVG.NE.0) THEN
       CALL ZDISP(0,LCORVG,'DCORVG')
       IF (IER.NE.0) GOTO 99999
      ENDIF
C     
99999 END
      
      
