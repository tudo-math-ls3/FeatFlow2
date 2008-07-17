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
* Version from  02/12/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* T        R*8    -1D0: Initialization                                 *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* T        R*8    Absolute cpu time                                    *
*                                                                      *
************************************************************************
C
      SUBROUTINE ZTIME(T)
C
      DOUBLE PRECISION T
      REAL ETIME,VTA
      DIMENSION VTA(2)
C
      T=0
C
      END
