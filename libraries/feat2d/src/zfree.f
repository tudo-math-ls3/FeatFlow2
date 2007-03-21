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
* ZFREE                                                                *
*                                                                      *
* Purpose  Calculate free space on workspace DWORK                     *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  08/14/90                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* ITYPE    I*4    Data type                                            *
*                  1  REAL*8     (Double precision)                    *
*                  2  REAL*4     (Single precision)                    *
*                  3  INTEGER*4                                        *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* IFREE    I*4    Number of free entries of type ITYPE                 *
* IER      I*4    Error indicator                                      *
*                 -106  Wrong value of ITYPE                           *
*                                                                      *
************************************************************************
C
      SUBROUTINE ZFREE(ITYPE,IFREE)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNARR=299)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /CHAR/,/ERRCTL/
C
      INTEGER ZILEND,ZVLEND
      EXTERNAL ZILEND,ZVLEND

      IF (ICHECK.GE.998) CALL OTRC('ZFREE ','08/14/90')
      IER=0
C
      IF (ITYPE.LT.1.OR.ITYPE.GT.3) THEN
C ***  Error *** Wrong value of ITYPE ***
       WRITE (CPARAM,'(I15)') ITYPE
       CALL WERR(-106,'ZFREE ')
       GOTO 99999
      ENDIF
C
C *** Calculate number of free elements on DWORK ***
      IFREE=NWORK-IWORK
      IF (ITYPE.EQ.2) IFREE=IFREE*ZVLEND()
      IF (ITYPE.EQ.3) IFREE=IFREE*ZILEND()
C
99999 END
