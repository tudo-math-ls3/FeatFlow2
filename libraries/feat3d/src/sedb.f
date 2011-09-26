************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT3D  (Release 1.1)               *
*                                                                      *
* Authors: J. Harig, S. Turek                                          *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* SEDB                                                                *
*                                                                      *
* Purpose   Calculates the parameter value of the midpoint on the      *
*           boundary segment joining two vertices with given parameter *
*           values                                                     *
*                                                                      *
* Subroutines/functions called   TMAX                                  *
*                                                                      *
* Version from  01/02/93                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* Xi,Yi,Zi R*8    Parameter values                                     *
* IBCT            Boundary component                                   *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* XT,YT,ZT R*8    Parameter value of midpoint                          *
*                                                                      *
************************************************************************
C
      SUBROUTINE SEDB(X1,Y1,Z1,X2,Y2,Z2,XT,YT,ZT,IBCT)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('SEDB ','01/02/93')
C
      XT=0.5D0*(X1+X2)
      YT=0.5D0*(Y1+Y2)
      ZT=0.5D0*(Z1+Z2)
C
      END
