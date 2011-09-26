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
* SIVB                                                                 *
*                                                                      *
* Purpose  Determination of the minimum and maximum index on  KVBD     *
*          preceeding CALL of  SVEB  with  IPAR=1  and                 *
*          existence of DVBDP assumed                                  *
*                                                                      *
* Subroutines/functions called  None                                   *
*                                                                      *
* Version from  08/16/90                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* IBCT     I*4    Boundary component                                   *
* KBCT     I*4    Vertices on boundary component IBCT start at         *
*                 KBCT(IBCT)                                           *
* TMAX     R*8    EXTERNAL FUNCTION                                    *
* DVBDP    R*8    Parameters of boundary vertices                      *
* DPAR1    R*8    Minimum parameter value                              *
* DPAR2    R*8    Maximum parameter value                              *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* IVBD1    I*4    Minimum index - first node on first edge             *
* IVBD2    I*4    Maximum index - first node on last edge              *
* IER      I*4    Error indicator                                      *
*                 IER=-140  DPAR1 >= DPAR2                             *
*                 IER=-141  DPAR1 or DPAR2 exceed valid range          *
*                                                                      *
************************************************************************
C
      SUBROUTINE SIVB(IBCT,KBCT,TMAX,DVBDP,DPAR1,DPAR2,
     *                IVBD1,IVBD2)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (EPS=1D-8)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION KBCT(*),DVBDP(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='SIVB'
      IF (ICHECK.GE.997) CALL OTRC('SIVB  ','08/16/90')
C
      IER=0
      DTMAX=TMAX(IBCT)
      IF (ICHECK.GT.0) THEN
       IF (DPAR1.GE.DPAR2) CALL WERR(-140,'SIVB  ')
       IF ((DPAR1.LT.0D0).OR.(DPAR2.GT.DTMAX)) CALL WERR(-141,'SIVB  ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
      IF (DPAR1.EQ.0D0) THEN
       IVBD1=KBCT(IBCT)
      ELSE
       DO 10 IVBD1=KBCT(IBCT),KBCT(IBCT+1)-1
       DPAR=DVBDP(IVBD1)
       IF ((DPAR-DPAR1).GT.-EPS) GOTO 15
10     CONTINUE
      ENDIF
C
15    IF (DPAR2.EQ.DTMAX) THEN
       IVBD2=KBCT(IBCT+1)-1
      ELSE
       DO 20 IVBD2=KBCT(IBCT+1)-1,KBCT(IBCT),-1
       DPAR=DVBDP(IVBD2)
       IF ((DPAR-DPAR2).LT.EPS) GOTO 25
20     CONTINUE
25     IVBD2=IVBD2-1
      ENDIF
C
99999 END
