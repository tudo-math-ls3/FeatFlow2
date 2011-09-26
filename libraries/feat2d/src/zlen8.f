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
* ZLEN8                                                                *
*                                                                      *
* Purpose  Determine double word length of an array allocated on DWORK *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LNR      I*4    Number of array                                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* LENGTH   I*4    >0  Number of elements available                     *
*                  0  Array not yet allocated                          *
*                 -1  LNR invalid                                      *
*                                                                      *
************************************************************************
C
      SUBROUTINE ZLEN8(LNR,LENGTH)
C
      PARAMETER (NNARR=299)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /TABLE/  KTYPE(NNARR),KLEN(NNARR),KLEN8(NNARR),IFLAG
      SAVE /TABLE/
C
      IF (ICHECK.EQ.999) CALL OTRC('ZLEN  ','01/02/89')
C
      IF (LNR.GT.0.AND.LNR.LE.NNARR) THEN
       LENGTH=KLEN8(LNR)
      ELSE IF (LNR.EQ.0) THEN
       LENGTH=0
      ELSE
       LENGTH=-1
      ENDIF
C
      END
