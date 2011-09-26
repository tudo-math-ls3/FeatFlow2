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
* XSMA2V                                                               *
*                                                                      *
* Purpose  Allocation of array KMAVEL                                  *
*          Call of S2MAV                                               *
*                                                                      *
* Subroutines/functions called  S2MAV, ZNEW, ZDISP, ZLEN               *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
*                                                                      *
************************************************************************
*                                                                      *
* S2V                                                                  *
*                                                                      *
* Purpose  Determine numbers of elements meeting at each vertex        *
*                                                                      *
* Subroutines/functions called   None                                  *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* KVERT    I*4    Numbers of vertices in each element, counterclockwise*
* KADJ     I*4    Number of neighbouring element, 0 for boundary edges *
* ICHK     I*4    See below                                            *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* KVEL     I*8    ICHK=1:                                              *
*                  KVEL(IVT,1) = #Elements meeting at vertex IVT       *
*                 ICHK=0:                                              *
*                  KVEL(IVT,.) = Numbers of elements meeting at IVT    *
*                                sorted in counterclockwise sense      *
*                                                                      *
************************************************************************
C
      SUBROUTINE XSMA2V
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /MACROD/ NMAEL,NMAVT,NMAEDG,NMAVE,NMAVEL,NMABCT,NMAVBD
      COMMON /MACROA/ LMACVG,LMACMG,LMAVT,LMAMID,LMAADJ,LMAVEL,LMAMEL,
     *                LMANPR,LMAMM,LMAVBD,LMAEBD,LMABCT,LMAVBP,LMAMBP,
     *                LMAVE
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /ERRCTL/,/CHAR/,/TRIAD/,/MACROD/,/MACROA/
C
      SUB='XSMA2V'
      IF (ICHECK.GE.997) CALL OTRC('XS2MAV','01/02/89')
C
      CALL ZLEN(LMAVEL,ILEN)
      IF (ILEN.LT.NMAVT) THEN
       IF (ILEN.GT.0) CALL ZDISP(0,LMAVEL,'KMAVEL')
       CALL ZNEW(NMAVT,3,LMAVEL,'KMAVEL')
       IF (IER.NE.0) GOTO 99999
      ELSE
       CALL ZCLEAR(LMAVEL,'KMAVEL')
      ENDIF
C
      NMAVEL=1
      CALL S2V(KWORK(L(LMAVT)),KWORK(L(LMAADJ)),KWORK(L(LMAVEL)),
     *         NVE,NMAEL,NMAVT,NMAVEL,1)
C
      CALL ZLEN(LMAVEL,ILEN)
      IF (ILEN.LT.NMAVEL*NMAVT) THEN
       CALL ZDISP(0,LMAVEL,'KMAVEL')
       CALL ZNEW (NMAVEL*NMAVT,3,LMAVEL,'KMAVEL')
       IF (IER.NE.0) GOTO 99999
      ELSE
       CALL ZCLEAR(LMAVEL,'KMAVEL')
      ENDIF
C
      CALL S2V(KWORK(L(LMAVT)),KWORK(L(LMAADJ)),KWORK(L(LMAVEL)),
     *         NVE,NMAEL,NMAVT,NMAVEL,0)
C
99999 END
