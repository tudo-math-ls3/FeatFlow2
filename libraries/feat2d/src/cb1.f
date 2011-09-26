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
* CB1                                                                  *
*                                                                      *
* Purpose  Set coordinates and weights for numerical integration       *
*          over the interval [-1,1]                                    *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* ICUB     I*4    Number of cubature formula                           *
*                 ICUB    NCUB    DEGREE                               *
*                   1       1        2      Gaussian formula           *
*                   2       2        2      Trapezoidal rule           *
*                   3       2        4      Gaussian formula           *
*                   4       3        4      Simpson rule               *
*                   5       3        6      Gaussian formula           *
*                   6       4        8      Gaussian formula           *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DXI      R*8    Cartesian coordinates of cubature points             *
* DOMEGA   R*8    Weights                                              *
* NCUBP    I*4    Number of cubature points                            *
* IER      I*4    Error indicator                                      *
*                 -119 Wrong value of ICUB                             *
*                                                                      *
************************************************************************
C
      SUBROUTINE CB1(ICUB)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNCUBP=36)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/CUB/
C
      SUB='CB1'
      IF (ICHECK.GE.997) CALL OTRC('CB1   ','01/02/89')
C
      IER=0
C
      GOTO (10,20,30,40,50,60) ICUB
C
      WRITE (CPARAM,'(I15)') ICUB
      CALL WERR(-119,'CB1   ')
      GOTO 99999
C
10    DXI(1,1)= 0D0
      DOMEGA(1)= 2D0
      NCUBP=1
      GOTO 99999
C
20    DXI(1,1)=-1D0
      DXI(2,1)= 1D0
      DOMEGA(1)= 1D0
      DOMEGA(2)= 1D0
      NCUBP=2
      GOTO 99999
C
30    DXI(1,1)=-0.577350269189626D0
      DXI(2,1)= 0.577350269189626D0
      DOMEGA(1)= 1D0
      DOMEGA(2)= 1D0
      NCUBP=2
      GOTO 99999
C
40    DXI(1,1)=-1D0
      DXI(2,1)= 0D0
      DXI(3,1)= 1D0
      DOMEGA(1)= 0.333333333333333D0
      DOMEGA(2)= 1.333333333333333D0
      DOMEGA(3)= 0.333333333333333D0
      NCUBP=3
      GOTO 99999
C
50    DXI(1,1)=-0.774596669241483D0
      DXI(2,1)= 0D0
      DXI(3,1)= 0.774596669241483D0
      DOMEGA(1)= 0.555555555555556D0
      DOMEGA(2)= 0.888888888888889D0
      DOMEGA(3)= 0.555555555555556D0
      NCUBP=3
      GOTO 99999
C
60    DXI(1,1)=-0.861136311594053D0
      DXI(2,1)=-0.339981043584856D0
      DXI(3,1)= 0.339981043584856D0
      DXI(4,1)= 0.861136311594053D0
      DOMEGA(1)= 0.347854845137454D0
      DOMEGA(2)= 0.652145154862546D0
      DOMEGA(3)= 0.652145154862546D0
      DOMEGA(4)= 0.347854845137454D0
      NCUBP=4
C
99999 END
