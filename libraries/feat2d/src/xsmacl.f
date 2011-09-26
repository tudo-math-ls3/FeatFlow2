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
* XSMACL                                                               *
*                                                                      *
* Purpose  Make clean information on macro-elements                    *
*                                                                      *
* Subroutines/functions called   ZDISP                                 *
*                                                                      *
* Version from  04/12/91                                               *
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
      SUBROUTINE XSMACL
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      COMMON /MACROD/ NMAEL,NMAVT,NMAEDG,NMAVE,NMAVEL,NMABCT,NMAVBD
      COMMON /MACROA/ LMACVG,LMACMG,LMAVT,LMAMID,LMAADJ,LMAVEL,LMAMEL,
     *                LMANPR,LMAMM,LMAVBD,LMAEBD,LMABCT,LMAVBP,LMAMBP,
     *                LMAVE
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /MACROD/,/MACROA/,/ERRCTL/,/CHAR/
C
      SUB='XSMACL'
      IF (ICHECK.GE.997) CALL OTRC('XSMACL','04/12/91')
      IER=0
C
      NMAEL =0
      NMAVT =0
      NMAEDG =0
      NMAVE =0
      NMAVEL =0
      NMABCT=0
      NMAVBD=0
C
      IF (LMAVE.NE.0) THEN
       CALL ZDISP(0,LMAVE,'KMAVE ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LMAMBP.NE.0) THEN
       CALL ZDISP(0,LMAMBP,'DMAMBP')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LMAVBP.NE.0) THEN
       CALL ZDISP(0,LMAVBP,'DMAVBP')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LMABCT.NE.0) THEN
       CALL ZDISP(0,LMABCT,'KMABCT')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LMAEBD.NE.0) THEN
       CALL ZDISP(0,LMAEBD,'KMAEBD')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LMAVBD.NE.0) THEN
       CALL ZDISP(0,LMAVBD,'KMAVBD')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LMAMM.NE.0) THEN
       CALL ZDISP(0,LMAMM,'KMAMM ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LMANPR.NE.0) THEN
       CALL ZDISP(0,LMANPR,'KMANPR ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LMAMEL.NE.0) THEN
       CALL ZDISP(0,LMAMEL,'KMAMEL')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LMAVEL.NE.0) THEN
       CALL ZDISP(0,LMAVEL,'KMAVEL')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LMAVE.NE.0) THEN
       CALL ZDISP(0,LMAVE,'KMAVE ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LMAADJ.NE.0) THEN
       CALL ZDISP(0,LMAADJ,'KMAADJ ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LMAMID.NE.0) THEN
       CALL ZDISP(0,LMAMID,'KMAMID  ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LMAVT.NE.0) THEN
       CALL ZDISP(0,LMAVT,'KMAVT ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LMACMG.NE.0) THEN
       CALL ZDISP(0,LMACMG,'DMACMG')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LMACVG.NE.0) THEN
       CALL ZDISP(0,LMACVG,'DMACVG')
       IF (IER.NE.0) GOTO 99999
      ENDIF
C     
99999 END
      
      
