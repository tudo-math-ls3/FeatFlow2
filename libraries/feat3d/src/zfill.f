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
* ZFILL                                                                *
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
      SUBROUTINE ZFILL(IFILL,DFILL)
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
      IF (ICHECK.GE.998) CALL OTRC('ZFILL ','08/14/93')
      IER=0
      IFILL=0
C
C
C *** Calculate number of free elements on DWORK ***
      DO 10 IARR=1,NNARR
      CALL ZLEN8(IARR,ILEN)
      IFILL=IFILL+ILEN
      write(6,*) iarr,ilen,ifill
10    CONTINUE
C
      DFILL=DBLE(IFILL)/DBLE(IWORK)
C
99999 END
