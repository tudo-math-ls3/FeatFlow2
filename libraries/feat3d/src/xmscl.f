************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.3)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek                     *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
* based on the 2-D routine, modificated by P.Schreiber                 *
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
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
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
      KNET(ILEV) =0 
      KNAT(ILEV) =0 
      KNVE(ILEV) =0 
      KNEE(ILEV) =0 
      KNAE(ILEV) =0 
      KNVEL(ILEV) =0 
      KNEEL(ILEV) =0 
      KNVED(ILEV) =0 
      KNVAR(ILEV) =0
      KNEAR(ILEV) =0
      KNBCT(ILEV) =0
      KNVBD(ILEV) =0
      KNEBD(ILEV) =0 
      KNABD(ILEV) =0 
C
      IF (KLCVG(ILEV).NE.0) THEN
       CALL ZDISP(0,KLCVG(ILEV),'DCVG0')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLCMG(ILEV).NE.0) THEN
       CALL ZDISP(0,KLCMG(ILEV),'DCMG0')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLCAG(ILEV).NE.0) THEN
       CALL ZDISP(0,KLCAG(ILEV),'DCAG0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLVERT(ILEV).NE.0) THEN
       CALL ZDISP(0,KLVERT(ILEV),'KVERT0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLEDGE(ILEV).NE.0) THEN
       CALL ZDISP(0,KLEDGE(ILEV),'KEDG0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLAREA(ILEV).NE.0) THEN
       CALL ZDISP(0,KLAREA(ILEV),'KARE0  ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLADJ(ILEV).NE.0) THEN
       CALL ZDISP(0,KLADJ(ILEV),'KADJ0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLVEL(ILEV).NE.0) THEN
       CALL ZDISP(0,KLVEL(ILEV),'KVEL0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLEEL(ILEV).NE.0) THEN
       CALL ZDISP(0,KLEEL(ILEV),'KEEL0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLAEL(ILEV).NE.0) THEN
       CALL ZDISP(0,KLAEL(ILEV),'KAEL0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLVED(ILEV).NE.0) THEN
       CALL ZDISP(0,KLVED(ILEV),'KVED0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLAED(ILEV).NE.0) THEN
       CALL ZDISP(0,KLAED(ILEV),'KAED0')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLVAR(ILEV).NE.0) THEN
       CALL ZDISP(0,KLVAR(ILEV),'KVAR0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLEAR(ILEV).NE.0) THEN
       CALL ZDISP(0,KLEAR(ILEV),'KEAR0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLEVE(ILEV).NE.0) THEN
       CALL ZDISP(0,KLEVE(ILEV),'KEVE0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLAVE(ILEV).NE.0) THEN
       CALL ZDISP(0,KLAVE(ILEV),'KAVE0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLNPR(ILEV).NE.0) THEN
       CALL ZDISP(0,KLNPR(ILEV),'KNPR0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLBCT(ILEV).NE.0) THEN
       CALL ZDISP(0,KLBCT(ILEV),'KBCT0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLVBD(ILEV).NE.0) THEN
       CALL ZDISP(0,KLVBD(ILEV),'KVBD0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLEBD(ILEV).NE.0) THEN
       CALL ZDISP(0,KLEBD(ILEV),'KEBD0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (KLABD(ILEV).NE.0) THEN
       CALL ZDISP(0,KLABD(ILEV),'KABD0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
1     CONTINUE
C     
99999 END
      
      
