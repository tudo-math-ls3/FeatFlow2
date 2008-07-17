************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.0)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller                               *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* ZTIME                                                                *
*                                                                      *
* Purpose  Return cpu time  (machine dependent)                        *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* New version from  06/22/94 (Monika Best) for IBM AIX                 *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* T        R*8    -1D0: Initialization                                 *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* T        R*8    Absolute cpu time in seconds                         *
*                                                                      *
************************************************************************
C
      SUBROUTINE ZTIME(T)
C
      DOUBLE PRECISION T
      INTEGER ITIME
C
      T=MCLOCK(ITIME) / 100.
C
      END
